/* Title: Quantum Computer Simulation PHYM004
Author: Ewan-James Thomas
File_name: quantcomp.c
License: Public Domain
*/

/****************** PHYM004 PR2: Quantum Computer Simulation ******************/
//
// Purpose
// -------
//
// Usage
// -----
//  All gates apart from measure_register return type wavefunction and so the wavefunction of the entire
// system must be updates with each gate call .
// 
//  Compiled with make file by using command "make quant" on command line. Then execute the binary file with 
// "./quant" to yield output
// Example
// -------
//
/*      UPDATES
 Date         Version  Comments  #P.x.y means pending 
 ----         -------  --------
 13/01/21       0.0.1  Create file with intention to work on part 1: Building a quantum register
 14/01/21       0.0.2  Quantum register of 3 qubits created. Now working on quantum gates for register
 24/01/21       0.0.3  Measuring all Qbits in register now works with a different seed for rng each time based on TOD
 24/01/21       0.0.4  Started working on Hadamard Gate
 26/01/21       0.0.5  Hadamard gate for 3 qubits working, however not flexible wrt N
 28/01/21       0.0.6  Hadamard gate working for general N
 28/01/21       0.0.7  Same goes for phase gate
 28/01/21       0.0.8  Implemented Grovers algorithm, not sure if working...
 1/02/21        0.0.9  Grovers Algorithm corectly implemented.
 1/02/21        0.1.0  Need to Generalise for any number of qubits, N as i kepp getting memory problems
 1/02/21        0.1.0  Fixed for general N, although some refinement needed (line 96)
 16/02/21       P.1.1  Want to implement sparse matrices
 17/02/21       P.1.2  Implementing CNOT Gate
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <sys/time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>

#define BASIS 2
#define N 3 // Number of qubits defined
#define STATES_MAX 1024 //max of 10 qubits 

struct timeval tv;

void qubit_error(){
    printf("Please operate the gate on a valid qubit\n");
    exit(1);
}

typedef struct element{
    int a;
    int b;
} element;
//   Creates a wavefunction where all spins are "down" i.e. |00..0>.
//  A 1 in the ket represents a spin down. e.g. |010> means all qubits are spin down apart from qbit 2.  Note how 010 is 2 in binary.
//  Since 000 is 0 in binary, we set probability amplitude of state 0 to 1.
gsl_vector_complex* initWavefunctionSpinDown(int states){ 
    gsl_vector_complex* wavefunction = NULL;
    wavefunction = gsl_vector_complex_alloc(states);
    // Probability amplitude stored in "wavefuction" vector
    gsl_vector_complex_set(wavefunction, 0, GSL_COMPLEX_ONE);

    return wavefunction;
}

// Set the probability amplitudes of all states to 1/sqrt(basis^N) so equal probability of wavefunction collapsing into any state.
gsl_vector_complex* initWavefunctionEqualProb(int states){ //Initialising wf to "Equal Probability" of all states
    
    gsl_vector_complex* wavefunction = NULL;
    wavefunction = gsl_vector_complex_alloc(states);
    gsl_vector_complex_set_all(wavefunction, gsl_complex_rect(1/sqrt(8),0));

    return wavefunction;
}
 // Takes an integer and return binary representation in string format
char* intToBinary(int a){
    int bin = 0;
    int remainder, temp = 1;

    while(a != 0){
        remainder = a % 2;
        a /= 2;
        bin += remainder*temp;
        temp *= 10;
    }
    char *bin_str = (char *)malloc(N*sizeof(char));
    sprintf(bin_str, "%03d", bin);    
    return bin_str;
}

// To measure the quantum state we must use the discrete probability distribution provided by the 
// wavefuntion. The measurement function finds the probabilities of each state and observes the wf
// according to those probabilities using a random number generator. After measurement the wavefunction is "collapsed" and gives the 
// measurement over and over again.
void measureRegisterGate(gsl_vector_complex* wavefunction){
    int states = (int) wavefunction->size;

    double probabilities[states];
    for(int i = 0; i < states; i++){ //creating a list of the probabilities (non-normalised)
        probabilities[i] = GSL_REAL(gsl_vector_complex_get(wavefunction,i))*GSL_REAL(gsl_vector_complex_get(wavefunction,i)) + GSL_IMAG(gsl_vector_complex_get(wavefunction,i))*GSL_IMAG(gsl_vector_complex_get(wavefunction,i));
        //printf("%lg\n", probabilities[i]);
    }
    gsl_ran_discrete_t* lookup = gsl_ran_discrete_preproc(states, probabilities); // Preproc normalises the probabilities
    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937); // Mersene Twister Algorithm is used for high quality random numbers
    gettimeofday(&tv,NULL); // Get time of day in usec so that the seed changes giving a different stream of #'s each time
    unsigned long int seed = tv.tv_usec;    
    gsl_rng_set(r, seed); 
    size_t t = gsl_ran_discrete(r, lookup); // Choosing from the discrete probability distribution defined by the wavefunction 
    printf("Wavefunction collapsed into the state:\n|%s>\n", intToBinary(t));
    // Wavefunction collapsed so will only find system in the state from now on
    gsl_vector_complex_set_all(wavefunction, GSL_COMPLEX_ZERO);
    gsl_vector_complex_set(wavefunction, t, GSL_COMPLEX_ONE); // Set measured state to probability one so that if we measure again we get the same outcome
    // Free memory to avoid bloats
    gsl_ran_discrete_free(lookup);
    gsl_rng_free(r);
    return;
}


// This function will find the element of the tensor product for a given gate for one qubit
double findElementHad(char* inta, char* intb, int qubit){
    // Hadamard gate for single qubit used to calculate tensor product
    gsl_matrix_complex *hadamard_single = gsl_matrix_complex_alloc(BASIS, BASIS);
    gsl_matrix_complex_set_all(hadamard_single, gsl_complex_rect(1/sqrt(BASIS),0));
    gsl_matrix_complex_set(hadamard_single,1,1, gsl_complex_rect(-1/sqrt(BASIS),0));

    double value = 1.0;
    for(int i = 0; i < N; i++){
        if(inta[i] != intb[i] && i != qubit - 1){
            return 0.0; // Invokes Kronecker delta
        }
        value =  GSL_REAL(gsl_matrix_complex_get(hadamard_single, inta[qubit-1] - '0', intb[qubit -1] - '0'));

    }
    gsl_matrix_complex_free(hadamard_single);
    return value;
}

double findElementPhase(char* inta, char* intb, int qubit, double phi){
    // Hadamard gate for single qubit used to calculate tensor product
    gsl_matrix_complex *phase_single = gsl_matrix_complex_alloc(BASIS, BASIS);
    gsl_matrix_complex_set_identity(phase_single);
    gsl_matrix_complex_set(phase_single,1,1, gsl_complex_polar(1,phi));

    double value = 1.0;
    for(int i = 0; i < N; i++){
        if(inta[i] != intb[i] && i != qubit - 1){
            return 0.0; // Invokes Kronecker delta
        }
        value =  GSL_REAL(gsl_matrix_complex_get(phase_single, inta[qubit-1] - '0', intb[qubit -1] - '0'));

    }
    gsl_matrix_complex_free(phase_single);
    return value;
}

double findElementCnot(char* row, char* col, int target_qubit, int control_qubit, int num_qbits){
    // Defining the two qubit cnot gate that we will draw values from
    gsl_matrix *cnot = gsl_matrix_alloc(BASIS*BASIS, BASIS*BASIS);
    gsl_matrix_set(cnot, 0, 0, 1);
    gsl_matrix_set(cnot, 1, 1, 1);
    gsl_matrix_set(cnot, 2, 3, 1);
    gsl_matrix_set(cnot, 3, 2, 1);

    for(int i = 0; i < num_qbits; i++){
        // Employ the deltas first
        char char1 = row[i];
        char char2 = col[i];
        if(i != target_qubit-1 && i != control_qubit-1 && strcmp(&char1, &char2)!= 0){
            // If an element of row and col strings not a match on any element other that targ or contr then we know 
            // a delta will set the whole element to zero
            printf("%c\n", char1);
            printf("%c\n", char2);
            printf("%d\n", strcmp(&char1, &char2));
            return 0.0;
        }
    }
    // use strcat to put the control index first in the string that will be converted to 
    char str1[BASIS];
    char str2[BASIS];
    strcat(str1, &row[control_qubit-1]);
    strcat(str1, &row[target_qubit-1]);
    strcat(str2, &col[control_qubit-1]);
    strcat(str2, &col[target_qubit-1]);
    printf("i am a cat2");


    // use strtol to make a number that will give us a value from cnot  above
    long row_index = strtol(str1, NULL, 10);
    long col_index = strtol(str2, NULL, 10);
    double value = gsl_matrix_get(cnot, row_index, col_index);
    gsl_matrix_free(cnot);
    return value;
}
// The Hadamard gate sets qubits into a superposition of their basis states. This function does this 
// by allowing the user to specify which qubit we will be setting to a superposition. This will have 
// effects when it comes to measuring the whole register as sometimes that qubit will be spin up, 
// sometimes it will be spin down, which will change the overall state of the register.
gsl_vector_complex* hadamardGate(gsl_vector_complex* wavefunction, int qubit){
    if(qubit > N){
        printf("Please operate the gate on a valid qubit\n");
        exit(0);
    }
    // Will beome the NxN matrix for operation on whole register
    gsl_matrix_complex *hadamard = gsl_matrix_complex_alloc(wavefunction->size, wavefunction->size);
    gsl_matrix_complex_set_all(hadamard, GSL_COMPLEX_ZERO);

    for(int i = 0; i < wavefunction->size; i++){
        for(int j = 0; j < wavefunction->size; j++){
            double val = findElementHad(intToBinary(i), intToBinary(j), qubit);
            gsl_matrix_complex_set(hadamard, i , j, gsl_complex_rect(val,0));
        }
    }
    gsl_vector_complex* h_psi = gsl_vector_complex_alloc(wavefunction->size);
    gsl_vector_complex_set_zero(h_psi);
    // When implementing the sparse matrix version, the gate it real and the wf is complex so can split into real and im
    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, hadamard, wavefunction, GSL_COMPLEX_ZERO, h_psi);
    return h_psi;
}

gsl_vector_complex* phaseShiftGate(gsl_vector_complex *wavefunction, int qubit, float phase){
    if(qubit > N){
        printf("Please operate the gate on a valid qubit\n");
        exit(0);
    }
    // Will beome the NxN matrix for operation on whole register
    gsl_matrix_complex *phaseGate = gsl_matrix_complex_alloc(wavefunction->size, wavefunction->size);
    gsl_matrix_complex_set_all(phaseGate, GSL_COMPLEX_ZERO);

    for(int i = 0; i < wavefunction->size; i++){
        for(int j = 0; j < wavefunction->size; j++){
            double val = findElementPhase(intToBinary(i), intToBinary(j), qubit, phase); //This is causing some errors
            gsl_matrix_complex_set(phaseGate, i , j, gsl_complex_rect(val,0));
        }
    }
    gsl_vector_complex* r_psi = gsl_vector_complex_alloc(wavefunction->size);
    gsl_vector_complex_set_zero(r_psi);
    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, phaseGate, wavefunction, GSL_COMPLEX_ZERO, r_psi);
    return r_psi;
}

// Oracle gate used in grovers quantum search algorithm. Argument answer is the "Correct question" mentioned
// in paper
gsl_vector_complex* oracleGate(gsl_vector_complex* wavefunction, int answer){
    gsl_matrix_complex* oracleGate = gsl_matrix_complex_alloc(wavefunction->size, wavefunction->size);
    gsl_matrix_complex_set_identity(oracleGate);
    gsl_matrix_complex_set(oracleGate, answer-1, answer-1, gsl_complex_rect(-1,0)); //Minus one as index from 0
    
    gsl_vector_complex* o_psi = gsl_vector_complex_alloc(wavefunction->size);
    gsl_vector_complex_set_zero(o_psi);
    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, oracleGate, wavefunction, GSL_COMPLEX_ZERO, o_psi);

    return o_psi;
}
// Unity matrix of size 2^N*2^N with -1 in the 0,0th element
gsl_vector_complex* jGate(gsl_vector_complex* wavefunction){
    gsl_matrix_complex* jGate = gsl_matrix_complex_alloc(wavefunction->size, wavefunction->size);
    gsl_matrix_complex_set_identity(jGate);
    gsl_matrix_complex_set(jGate, 0, 0, gsl_complex_rect(-1,0));
    
    gsl_vector_complex* j_psi = gsl_vector_complex_alloc(wavefunction->size);
    gsl_vector_complex_set_zero(j_psi);
    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, jGate, wavefunction, GSL_COMPLEX_ZERO, j_psi);

    return j_psi;
}

gsl_vector_complex* groversBlock(gsl_vector_complex* wavefunction, int answer){
    gsl_vector_complex* b_psi = gsl_vector_complex_alloc(wavefunction->size);
    gsl_vector_complex_set_zero(b_psi);
    //First operation in the block is to apply a quantum oracle gate
    b_psi = oracleGate(wavefunction, answer); 
    //Then we apply a hadamard gate to each gate
    for(int i = 1; i < N+1; i++){
        b_psi = hadamardGate(b_psi, i); // +1 bc had takes in canonical qubit number.       
    }
    // Apply the j gate
    b_psi = jGate(b_psi);
    // Finally a hadamard gate on each qubit again
    for(int j = 1; j < N+1; j++){
        b_psi = hadamardGate(b_psi, j);        
    }
    return b_psi;
}

// gsl_vector_complex* cnotGate(gsl_vector_complex* wavefunction, int target_qubit, int control_qubit);
// //     There is a knronecker delta for every qubit that the gate doesn't operate on. So similar to 
// //  findElementHad and findElementPhase implementation of kronecker deltas. we have pqr and p'q'r' 
// //  that represent the binary number of the element. The control qubit (in example 2) has its bit placed
// //  at the beinging

void print_wf(gsl_vector_complex* wavefunction){
    for (int i = 0; i < wavefunction->size; i++){
        printf("%lg\n", GSL_REAL(gsl_vector_complex_get(wavefunction, i)));
    }
}
int main(){
    int states = (int)pow(BASIS, N);
    gsl_vector_complex* wavefunction = initWavefunctionSpinDown(states);
    //Putting system into equal super position of superposition all 2^N basis'
    wavefunction = hadamardGate(wavefunction, 1);
    wavefunction = hadamardGate(wavefunction, 2); 
    wavefunction = hadamardGate(wavefunction, 3);

    for(int i = 0; i < floor(M_PI_4*sqrt(pow(2,N))); i++){ // Needs to be called "floor(pi/4*sqrt(2^N))"" times for optimum output roughly 2 in our case
        wavefunction = groversBlock(wavefunction, 3); //Second argument is the basis state you want to be "right" in this case its |110>
    }
    measureRegisterGate(wavefunction);
    double val = findElementCnot(intToBinary(6),intToBinary(7), 3, 2, N);
    printf("%lg\n", val);

    return 0;
}
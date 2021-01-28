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
 Date         Version  Comments
 ----         -------  --------
 13/01/21       0.0.1  Create file with intention to work on part 1: Building a quantum register
 14/01/21       0.0.2  Quantum register of 3 qubits created. Now working on quantum gates for register
 24/01/21       0.0.3  Measuring all Qbits in register now works with a different seed for rng each time based on TOD
 24/01/21       0.0.4  Started working on Hadamard Gate
 26/01/21       0.0.5  Hadamard gate for 3 qubits working, however not flexible wrt N
 28/01/21       0.0.6  Hadamard gate working for general N
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

#define BASIS 2
#define N 3 // Number of qubits defined
#define STATES_MAX 1024 //max of 10 qubits 

// Bit map for 3 qubits and their corresponding states
const char *bit_rep[16] = {
    [ 0] = "000", [ 1] = "001", [ 2] = "010", [ 3] = "011",
    [ 4] = "100", [ 5] = "101", [ 6] = "110", [ 7] = "111",
}; 

struct timeval tv;

void qubit_error(){
    printf("Please operate the gate on a valid qubit\n");
    exit(1);
}

typedef struct element{
    int a;
    int b;
} element;

gsl_vector_complex* init_wavefunction_sd(int states){ //Initialising wf to all "Spin Down" 
    //Initialising the wf to |000> 
    gsl_vector_complex* wavefunction = NULL;
    wavefunction = gsl_vector_complex_alloc(states);
    gsl_vector_complex_set(wavefunction, 0, gsl_complex_rect(1,0));

    return wavefunction;
}

gsl_vector_complex* init_wavefunction_ep(int states){ //Initialising wf to "Equal Probability" of all states
    
    gsl_vector_complex* wavefunction = NULL;
    wavefunction = gsl_vector_complex_alloc(states);
    gsl_vector_complex_set_all(wavefunction, gsl_complex_rect(1/sqrt(8),0));

    return wavefunction;
}

// To measure the quantum state we must use the discrete probability distribution provided by the 
// wavefuntion. The measurement function finds the probabilities of each state and observes the wf
// according to those probabilities. After measurement the wavefunction is "collapsed" and gives the 
// measurement over and over again.
void measure_register_gate(gsl_vector_complex* wavefunction){
    int states = (int) wavefunction->size;

    double probabilities[states];
    for(int i = 0; i < states; i++){ //creating a list of the probabilities (non-normalised)
        probabilities[i] = GSL_REAL(gsl_vector_complex_get(wavefunction,i))*GSL_REAL(gsl_vector_complex_get(wavefunction,i)) + GSL_IMAG(gsl_vector_complex_get(wavefunction,i))*GSL_IMAG(gsl_vector_complex_get(wavefunction,i));
    }
    gsl_ran_discrete_t* lookup = gsl_ran_discrete_preproc(states, probabilities); // Preproc normalises the probabilities
    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937); // Mersene Twister Algorithm is used for high quality random numbers
    gettimeofday(&tv,NULL); // Get time of day in usec so that the seed changes giving a different stream of #'s each time
    unsigned long int seed = tv.tv_usec;    
    gsl_rng_set(r, seed); 
    size_t t = gsl_ran_discrete(r, lookup); // Choosing from the discrete probability distribution defined by the wavefunction 
    printf("Wavefunction collapsed into the state:\n|%s>\n", bit_rep[t]);
    // Wavefunction collapsed so will only find system in the state from now on
    gsl_vector_complex_set_all(wavefunction, GSL_COMPLEX_ZERO);
    gsl_vector_complex_set(wavefunction, t, GSL_COMPLEX_ONE); // Set measured state to probability one so that if we measure again we get the same outcome
    // Free memory to avoid bloats
    gsl_ran_discrete_free(lookup);
    gsl_rng_free(r);
    return;
}
// The Hadamard gate sets qubits into a superposition of their basis states. This function does this 
// by allowing the user to specify which qubit we will be setting to a superposition. This will have 
// effects when it comes to measuring the whole register as sometimes that qubit will be spin up, 
// sometimes it will be spin down, which will change the overall state of the register.


// // This function will find the element of the tensor product for a given gate for one qubit
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
    return value;
}

char* intToBinary(int a){ // Now works in regards to printing leading zeros
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

gsl_vector_complex* hadamard_gate(gsl_vector_complex* wavefunction, int qubit){
    if(qubit > N){
        printf("Please operate the gate on a valid qubit\n");
        exit(0);
    }
    // Will beome the NxN matrix for operation on whole register
    gsl_matrix_complex *hadamard = gsl_matrix_complex_alloc(wavefunction->size, wavefunction->size);
    gsl_matrix_complex_set_all(hadamard, GSL_COMPLEX_ZERO);

    for(int i = 0; i < wavefunction->size; i++){
        for(int j = 0; j < wavefunction->size; j++){
            double val = findElementHad(intToBinary(i), intToBinary(j), qubit); //This is causing some errors
            gsl_matrix_complex_set(hadamard, i , j, gsl_complex_rect(val,0));
        }
    }
    gsl_vector_complex* h_psi = gsl_vector_complex_alloc(wavefunction->size);
    gsl_vector_complex_set_zero(h_psi);
    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, hadamard, wavefunction, GSL_COMPLEX_ZERO, h_psi);
    return h_psi;
}

gsl_vector_complex* phase_shift_gate(gsl_vector_complex *wavefunction, int qubit, float phase){
    if(qubit > N){
        printf("Please operate the gate on a valid qubit\n");
        exit(0);
    }
    // Will beome the NxN matrix for operation on whole register
    gsl_matrix_complex *phase_gate = gsl_matrix_complex_alloc(wavefunction->size, wavefunction->size);
    gsl_matrix_complex_set_all(phase_gate, GSL_COMPLEX_ZERO);

    for(int i = 0; i < wavefunction->size; i++){
        for(int j = 0; j < wavefunction->size; j++){
            double val = findElementPhase(intToBinary(i), intToBinary(j), qubit, phase); //This is causing some errors
            gsl_matrix_complex_set(phase_gate, i , j, gsl_complex_rect(val,0));
        }
    }
    gsl_vector_complex* r_psi = gsl_vector_complex_alloc(wavefunction->size);
    gsl_vector_complex_set_zero(r_psi);
    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, phase_gate, wavefunction, GSL_COMPLEX_ZERO, r_psi);
    return r_psi;
}

// Oracle gate used in grovers quantum search algorithm. Argument answer is the "Correct question" mentioned
// in paper
gsl_vector_complex* oracle_gate(gsl_vector_complex* wavefunction, int answer){
    gsl_matrix_complex* oracle_gate = gsl_matrix_complex_alloc(wavefunction->size, wavefunction->size);
    gsl_matrix_complex_set_identity(oracle_gate);
    gsl_matrix_complex_set(oracle_gate, answer-1, answer-1, gsl_complex_rect(-1,0)); //Minus one as index from 0
    
    gsl_vector_complex* o_psi = gsl_vector_complex_alloc(wavefunction->size);
    gsl_vector_complex_set_zero(o_psi);
    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, oracle_gate, wavefunction, GSL_COMPLEX_ZERO, o_psi);

    return o_psi;
}

gsl_vector_complex* diffusion_gate(gsl_vector_complex* wavefunction){
    gsl_matrix_complex* diffusion_gate = gsl_matrix_complex_alloc(wavefunction->size, wavefunction->size);
    gsl_matrix_complex_set_identity(diffusion_gate);
    gsl_matrix_complex_set(diffusion_gate, 0, 0, gsl_complex_rect(-1,0));
    
    gsl_vector_complex* j_psi = gsl_vector_complex_alloc(wavefunction->size);
    gsl_vector_complex_set_zero(j_psi);
    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, diffusion_gate, wavefunction, GSL_COMPLEX_ZERO, j_psi);

    return j_psi;
}

gsl_vector_complex* groversBlock(gsl_vector_complex* wavefunction, int answer){
    gsl_vector_complex* b_psi = gsl_vector_complex_alloc(wavefunction->size);
    gsl_vector_complex_set_zero(b_psi);
    //First operation in the block is to apply a quantum oracle gate
    b_psi = oracle_gate(wavefunction, answer); 
    //Then we apply a hadamard gate to each gate
    for(int i = 0; i < N; i++){
        b_psi = hadamard_gate(wavefunction, i);        
    }
    // Apply the diffusion gate
    b_psi = diffusion_gate(wavefunction);
    // Finally a hadamard gate on each qubit again
    for(int i = 0; i < N; i++){
        b_psi = hadamard_gate(wavefunction, i);        
    }
    return b_psi;
}

void print_wf(gsl_vector_complex* wavefunction){
    for (int i = 0; i < wavefunction->size; i++){
        printf("%lg\n", GSL_REAL(gsl_vector_complex_get(wavefunction, i)));
    }
}
int main(){
    int states = (int)pow(BASIS, N);
    gsl_vector_complex* wavefunction = init_wavefunction_sd(states);
    //Putting system into equal super position of superposition all 2^N basis'
    // wavefunction = hadamard_gate(wavefunction, 1);
    // wavefunction = hadamard_gate(wavefunction, 2); 
    // wavefunction = hadamard_gate(wavefunction, 3); 
    // wavefunction = phase_shift_gate(wavefunction, 1,  3.14159);
    // wavefunction = hadamard_gate(wavefunction, 1); 
    for(int i = 0; i < 4; i++){
        wavefunction = groversBlock(wavefunction, 7);
    }
    measure_register_gate(wavefunction);
    return 0;
}
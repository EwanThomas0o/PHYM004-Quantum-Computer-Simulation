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
// system must be updates with each gate call 
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
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
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
double findElement(char* inta, char* intb, int qubit){
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
 
        // if(i == (qubit - 1)){
        //     int row = inta[i] - '0'; //takes from ASCII to the actual int rep of the char
        //     int col = intb[i] - '0';
        //     value = GSL_REAL(gsl_matrix_complex_get(hadamard_single, row, col));
        //     return value;
        // }
    }
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
    printf("%s\n", bin_str);
    
    return bin_str;
}

gsl_vector_complex* hadamard_gate(gsl_vector_complex* wavefunction, int qubit){
    // Will beome the NxN matrix for operation on whole register
    gsl_matrix_complex *hadamard = gsl_matrix_complex_alloc(wavefunction->size, wavefunction->size);
    gsl_matrix_complex_set_all(hadamard, GSL_COMPLEX_ZERO);
    gsl_matrix_complex_set_identity(hadamard);

    // Hadamard gate for single qubit used to calculate tensor product
    gsl_matrix_complex *hadamard_single = gsl_matrix_complex_alloc(BASIS, BASIS);
    gsl_matrix_complex_set_all(hadamard_single, gsl_complex_rect(1/sqrt(BASIS),0));
    gsl_matrix_complex_set(hadamard_single,1,1, gsl_complex_rect(-1/sqrt(BASIS),0));

    if (qubit == 1){
        for(int i = 0; i < wavefunction->size; i++){
            for(int j = 0; j < wavefunction->size; j++){
                double val = findElement(intToBinary(i), intToBinary(j), qubit); //This is causing some errors
                gsl_matrix_complex_set(hadamard, i , j, gsl_complex_rect(val,0));
            }
        }
        // gsl_matrix_complex_set(hadamard, 4,4, gsl_complex_rect(-1.0,0));
        // gsl_matrix_complex_set(hadamard, 5,5, gsl_complex_rect(-1.0,0));
        // gsl_matrix_complex_set(hadamard, 6,6, gsl_complex_rect(-1.0,0));
        // gsl_matrix_complex_set(hadamard, 7,7, gsl_complex_rect(-1.0,0));
        // gsl_matrix_complex_set(hadamard, 4,0, gsl_complex_rect(1.0,0));
        // gsl_matrix_complex_set(hadamard, 5,1, gsl_complex_rect(1.0,0));
        // gsl_matrix_complex_set(hadamard, 6,2, gsl_complex_rect(1.0,0));
        // gsl_matrix_complex_set(hadamard, 7,3, gsl_complex_rect(1.0,0));
        // gsl_matrix_complex_set(hadamard, 0,4, gsl_complex_rect(1.0,0));
        // gsl_matrix_complex_set(hadamard, 1,5, gsl_complex_rect(1.0,0));
        // gsl_matrix_complex_set(hadamard, 2,6, gsl_complex_rect(1.0,0));
        // gsl_matrix_complex_set(hadamard, 3,7, gsl_complex_rect(1.0,0));
        // gsl_matrix_complex_scale(hadamard, gsl_complex_rect(1/sqrt(BASIS),0));
    }
    if (qubit == 2){
        gsl_matrix_complex_set(hadamard, 2,2, gsl_complex_rect(-1.0,0));
        gsl_matrix_complex_set(hadamard, 3,3, gsl_complex_rect(-1.0,0));
        gsl_matrix_complex_set(hadamard, 6,6, gsl_complex_rect(-1.0,0));
        gsl_matrix_complex_set(hadamard, 7,7, gsl_complex_rect(-1.0,0));
        gsl_matrix_complex_set(hadamard, 0,2, gsl_complex_rect(1.0,0));
        gsl_matrix_complex_set(hadamard, 1,3, gsl_complex_rect(1.0,0));
        gsl_matrix_complex_set(hadamard, 4,6, gsl_complex_rect(1.0,0));
        gsl_matrix_complex_set(hadamard, 5,7, gsl_complex_rect(1.0,0));
        gsl_matrix_complex_set(hadamard, 2,0, gsl_complex_rect(1.0,0));
        gsl_matrix_complex_set(hadamard, 3,1, gsl_complex_rect(1.0,0));
        gsl_matrix_complex_set(hadamard, 6,4, gsl_complex_rect(1.0,0));
        gsl_matrix_complex_set(hadamard, 7,5, gsl_complex_rect(1.0,0));
        gsl_matrix_complex_scale(hadamard, gsl_complex_rect(1/sqrt(BASIS),0));
    }
    if (qubit == 3){
        gsl_matrix_complex_set(hadamard, 1,1, gsl_complex_rect(-1.0,0));
        gsl_matrix_complex_set(hadamard, 3,3, gsl_complex_rect(1.0,0));
        gsl_matrix_complex_set(hadamard, 5,5, gsl_complex_rect(-1.0,0));
        gsl_matrix_complex_set(hadamard, 7,7, gsl_complex_rect(-1.0,0));
        gsl_matrix_complex_set(hadamard, 1,0, gsl_complex_rect(1.0,0));
        gsl_matrix_complex_set(hadamard, 0,1, gsl_complex_rect(1.0,0));
        gsl_matrix_complex_set(hadamard, 2,3, gsl_complex_rect(1.0,0));
        gsl_matrix_complex_set(hadamard, 3,2, gsl_complex_rect(1.0,0));
        gsl_matrix_complex_set(hadamard, 4,5, gsl_complex_rect(1.0,0));
        gsl_matrix_complex_set(hadamard, 5,4, gsl_complex_rect(1.0,0));
        gsl_matrix_complex_set(hadamard, 6,7, gsl_complex_rect(1.0,0));
        gsl_matrix_complex_set(hadamard, 7,6, gsl_complex_rect(1.0,0));
        gsl_matrix_complex_scale(hadamard, gsl_complex_rect(1/sqrt(BASIS),0));
    }
    else{
        //qubit_error();
    }
    gsl_vector_complex* h_psi = gsl_vector_complex_alloc(wavefunction->size);
    gsl_vector_complex_set_zero(h_psi);
    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, hadamard, wavefunction, GSL_COMPLEX_ZERO, h_psi);
    return h_psi;
}

void phase_shift_gate(gsl_vector_complex *wavefunction, int qubit, float phase){
    if(qubit == 1){
    }
    if(qubit == 2){
    }
    if(qubit == 3){
    }
    else{
        //qubit_error();
    }
    return;
}

void print_wf(gsl_vector_complex* wavefunction){
    for (int i = 0; i < wavefunction->size; i++){
        printf("%lg\n", GSL_REAL(gsl_vector_complex_get(wavefunction, i)));
    }
}
int main(){
    int states = (int)pow(BASIS, N);
    gsl_vector_complex* wavefunction = init_wavefunction_sd(states);

    measure_register_gate(wavefunction);
    wavefunction  = hadamard_gate(wavefunction, 1); 
    measure_register_gate(wavefunction);
    return 0;
}
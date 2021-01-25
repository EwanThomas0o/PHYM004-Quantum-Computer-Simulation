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
//
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
void hadamard_gate(gsl_vector_complex* wavefunction, int qubit){
    gsl_matrix_complex *hadamard = gsl_matrix_alloc(wavefunction->size, wavefunction->size);
    gsl_matrix__set_identity(hadamard);
    if (qubit == 1){
        gsl_matrix_set(hadamard, 4,4, -1);
        gsl_matrix_set(hadamard, 5,5, -1);
        gsl_matrix_set(hadamard, 6,6, -1);
        gsl_matrix_set(hadamard, 7,7, -1);
        gsl_matrix_set(hadamard, 0,4, 1);
        gsl_matrix_set(hadamard, 1,5, 1);
        gsl_matrix_set(hadamard, 2,6, 1);
        gsl_matrix_set(hadamard, 3,7, 1);
        gsl_matrix_set(hadamard, 3,7, 1);
        gsl_matrix_set(hadamard, 4,0, 1);
        gsl_matrix_set(hadamard, 5,1, 1);
        gsl_matrix_set(hadamard, 6,2, 1);
        gsl_matrix_set(hadamard, 7,3, 1);
        gsl_matrix_scale(hadamard, 1/sqrt(BASIS));
    }
    if (qubit == 2){
        gsl_matrix_set(hadamard, 2,2, -1);
        gsl_matrix_set(hadamard, 3,3, -1);
        gsl_matrix_set(hadamard, 6,6, -1);
        gsl_matrix_set(hadamard, 7,7, -1);
        gsl_matrix_set(hadamard, 0,2, 1);
        gsl_matrix_set(hadamard, 1,3, 1);
        gsl_matrix_set(hadamard, 2,6, 1);
        gsl_matrix_set(hadamard, 4,6, 1);
        gsl_matrix_set(hadamard, 5,7, 1);
        gsl_matrix_set(hadamard, 2,0, 1);
        gsl_matrix_set(hadamard, 3,1, 1);
        gsl_matrix_set(hadamard, 6,4, 1);
        gsl_matrix_set(hadamard, 7,5, 1);
        gsl_matrix_scale(hadamard, 1/sqrt(BASIS));
    }
    if (qubit == 3){
        gsl_matrix_set(hadamard, 1,1, -1);
        gsl_matrix_set(hadamard, 3,3, -1);
        gsl_matrix_set(hadamard, 5,5, -1);
        gsl_matrix_set(hadamard, 7,7, -1);
        gsl_matrix_set(hadamard, 0,4, 1);
        gsl_matrix_set(hadamard, 1,5, 1);
        gsl_matrix_set(hadamard, 2,6, 1);
        gsl_matrix_set(hadamard, 3,7, 1);
        gsl_matrix_set(hadamard, 3,7, 1);
        gsl_matrix_set(hadamard, 4,0, 1);
        gsl_matrix_set(hadamard, 5,1, 1);
        gsl_matrix_set(hadamard, 6,2, 1);
        gsl_matrix_set(hadamard, 7,3, 1);
        gsl_matrix_scale(hadamard, 1/sqrt(BASIS));
    }
    return;
}

int main(){
    int states = (int)pow(BASIS, N);
    gsl_vector_complex* wavefunction = init_wavefunction_ep(states);
    measure_register_gate(wavefunction);

    return 0;
}
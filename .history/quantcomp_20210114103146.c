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

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>

#define BASIS 2
#define N 3 // Number of qubits defined
#define STATES_MAX 1024 //max of 10 qubits 


gsl_vector_complex* init_wavfunction_sd(int states){ //Initialising wf to all sping down "sd"
    //Initialising the wf to |000> 
    gsl_vector_complex* wavefunction = NULL;
    wavefunction = gsl_vector_complex_alloc(states);
    gsl_vector_complex_set(wavefunction, 1, gsl_complex_rect(1,0));

    return wavefunction;
}

void measure_register_gate(gsl_vector_complex* wavefunction){
    int states = (int) wavefunction->size;

    double probabilities[states];
    for(int i = 0; i < states; i++){
        probabilities[i] = GSL_REAL(gsl_vector_complex_get(wavefunction,i))*GSL_REAL(gsl_vector_complex_get(wavefunction,i)) + GSL_IMAG(gsl_vector_complex_get(wavefunction,i))*GSL_IMAG(gsl_vector_complex_get(wavefunction,i));
    }
    gsl_ran_discrete_t* lookup = gsl_ran_discrete_preproc(states, probabilities);
    gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
    for(int l = 0; l < 20; l++){
        size_t t = gsl_ran_discrete(r, lookup);
        printf("%zu\n", t);
    }
    gsl_ran_discrete_free(lookup);
    
    
    
    return;
}


int main(){
    int states = (int)pow(BASIS, N);
    gsl_vector_complex* wavefunction = init_wavfunction_sd(states);
    measure_register_gate(wavefunction);

    return 0;
}
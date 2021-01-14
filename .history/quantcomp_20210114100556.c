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


gsl_matrix_complex* init_wavfunction(){

    int states = (int)pow(BASIS, N);
    gsl_matrix_complex* wavefunction = NULL;
    wavefunction = gsl_matrix_alloc(1,states);
    gsl_matrix_complex_set(wavefunction, 1,1, gsl_complex_rect(1,0)); //Initialises the wf to |000> 

    return wavefunction;
}


int main(){
    int states = (int)pow(BASIS, N);
    gsl_matrix_complex* wavefunction = init_wavfunction();
    
    double probabilities[states];
    for (int j = 0; j < states; j ++){
        probabilities[j] = wavefunction[j].real*wavefunction[j].real + wavefunction[j].imag*wavefunction[j].imag;
    }

    gsl_ran_discrete_t* lookup = gsl_ran_discrete_preproc(states, probabilities);
    gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
    for(int l = 0; l < 20; l++){
        size_t t = gsl_ran_discrete(r, lookup);
        printf("%zu\n", t);
    }
    double pdf =  gsl_ran_discrete_pdf(7, lookup); //gives normalised probability
    //printf("%lg\n", pdf);
    gsl_ran_discrete_free(lookup);
;    return 0;
}
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

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_randist.h>

#define BASIS 2
#define N 3 // Number of qubits defined
#define STATES_MAX 1024 //max of 10 qubits 

typedef struct{
    double real;  //TODO: Change to gsl complex type
    double imag;
}complex; //creating homemade complex dtype

complex* init_wavfunction(){
    int states = (int)pow(BASIS,N);

    complex *wavefunction;
    wavefunction = malloc(states*sizeof(complex));
    
    for(int i = 0; i < states; i++){
        wavefunction[i].real = 1/sqrt(states); //setting equal probability of each state
        wavefunction[i].imag = 0;
        //printf("%lg+%lgi\n", wavefunction[i].real, wavefunction[i].imag);
    }

    return wavefunction;
}


int main(){
    int states = (int)pow(BASIS, N);
    complex* wavefunction = init_wavfunction();
    
    double probabilities[states];
    for (int j = 0; j < states; j ++){
        probabilities[j] = wavefunction[j].real*wavefunction[j].real + wavefunction[j].imag*wavefunction[j].imag;
    }

    gsl_ran_discrete_t* lookup = gsl_ran_discrete_preproc(states, probabilities);
    gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
    size_t t = gsl_ran_discrete(r, lookup);
    double pdf =  gsl_ran_discrete_pdf(2, lookup);
    printf("%lg\n", pdf);
;    return 0;
}
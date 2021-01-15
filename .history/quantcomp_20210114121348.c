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

// Bit map for 3 qubits and their corresponding states
const char *bit_rep[16] = {
    [ 0] = "000", [ 1] = "001", [ 2] = "010", [ 3] = "011",
    [ 4] = "100", [ 5] = "101", [ 6] = "110", [ 7] = "111",
}; 


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
    for(int i = 0; i < states; i++){
        probabilities[i] = GSL_REAL(gsl_vector_complex_get(wavefunction,i))*GSL_REAL(gsl_vector_complex_get(wavefunction,i)) + GSL_IMAG(gsl_vector_complex_get(wavefunction,i))*GSL_IMAG(gsl_vector_complex_get(wavefunction,i));
    }
    gsl_ran_discrete_t* lookup = gsl_ran_discrete_preproc(states, probabilities);
    gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
    for(int l = 0; l < 5; l++){
        size_t t = gsl_ran_discrete(r, lookup);
        printf("Wavefunction collapsed into the state:\n|%s>\n", bit_rep[t]);
    }
    //gsl_vector_complex_set_all(wavefunction, GSL_COMPLEX_ZERO);
    //gsl_vector_complex_set(wavefunction, t, GSL_COMPLEX_ONE); // Set that state to probability one so that if we measure again we get the same outcome
    //double pdf =  gsl_ran_discrete_pdf(1, lookup); //gives normalised probability
    //printf("%lg\n", pdf);
    gsl_ran_discrete_free(lookup);
    
    
    return;
}
// The Hadamard gate sets qubits into a superposition of their basis states. This function does this 
// by allowing the user to specify which qubit we will be setting to a superposition. This will have 
// effects when it comes to measuring the whole register as sometimes that qubit will be spin up, 
// sometimes it will be spin down, which will change the overall state of the register. 
void hadamard_gate(gsl_vector_complex* wavefunction, int qubit){

}

int main(){
    int states = (int)pow(BASIS, N);
    gsl_vector_complex* wavefunction = init_wavefunction_ep(states);
    measure_register_gate(wavefunction);

    gsl_vector_complex* wavefunction2 = init_wavefunction_ep(states);
    measure_register_gate(wavefunction2);



    return 0;
}
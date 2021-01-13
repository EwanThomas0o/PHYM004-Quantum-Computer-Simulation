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

#define N 3 // Number of qubits defined
#define STATES_MAX 1024 //max of 10 qubits 

typedef struct{
    double real;
    double imag;
}complex;

int main(){
    int states = (int)pow(2,N);

    complex *wavefunction;
    wavefunction = malloc(states*sizeof(complex));
    
    for(int i = 0; i < states; i++){
        wavefunction[i].real = 1/sqrt(states); //setting equal probability of each state
        wavefunction[i].imag = 0;
        printf("%lg+%lgi\n", wavefunction[i].real, wavefunction[i].imag);
    }
    return 0;
}
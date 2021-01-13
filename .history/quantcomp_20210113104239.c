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
#include <complex.h>

#define N 3
#define STATES_MAX 1024

typedef struct{
    bool normed;
    double complex wf[];
} wavefunction;

int main(){
    int states = (int)pow(2,N);
    double complex wavefunction[states];
    //printf("%d\n", states);
    for(int i = 0; i< states; i++){
        //normalising the wavefunction with eqaul probability
        wavefunction[i] = 1/sqrt(states);
        printf("%lg\n", wavefunction[i]);
    }
    return 0;
}
/* Title: Quantum Computer Simulation PHYM004
Author: Ewan-James Thomas
File_name: quantcomp.c
License: Public Domain & GNU licensing 
*/

/****************** PHYM004 PR2: Quantum Computer Simulation ******************/
//
// Purpose
// -------
//
// Usage
// -----
//  All gates apart from measure_register return type wavefunction and so the wavefunction of the entire
//  system must be updates with each gate call .
// 
//  Compiled with make file by using command "make quant -q 5" on command line. -q flag defines number of qubits.
//  NOTE: Shor's Alg only works with 7 qubits
//  Then execute the binary file with "make do"
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
 28/01/21       0.0.7  Same goes for phase gate
 28/01/21       0.0.8  Implemented Grovers algorithm, not sure if working...
 1/02/21        0.0.9  Grovers Algorithm corectly implemented.
 1/02/21        0.1.0  Need to Generalise for any number of qubits, N as i keep getting memory problems
 1/02/21        0.1.0  Fixed for general N, although some refinement needed
 17/02/21       0.1.1  Implementing CNOT Gate. Error found. Will squash tomorrow.
 18/02/21       0.1.2  CNOT implemented. NOT 100% correct.
 18/02/21       0.1.3  CNOT implemented successfully
 22/02/21       0.1.4  Trying to implement sparse matrices to speed up
 22/02/21       0.1.4  using vector views etc to improve usability and speed (no need for memcpy)
 22/02/21       0.1.5  Sparse matrices implemented, Need to create complex sparse multiplier
 22/02/21        ""    Some progress making complex sp matrix multiplier, not working 100% yet tho
 23/02/21        ""    Change matrices to ints (rather than doubles to save memory)
 1/03/21         ""    Handbuild complex sparse matrix multiplier seems to be working in test program but sets wf to zero here...
 2/03/21        0.1.6  Phase gate now works! Was setting nz elements in spmatrix to 0!
 8/03/21        0.2.0  CPhase working. Now on to implementing the non-quantum section of Shor's Alg
 9/03/21         ""    Working on step 4 of shors algorithm still. Some bugs need fixing.
 9/03/21        0.3.0  Output of x~ is now correct for a=7, C=15, N=7
 10/03/21       0.4.0  Now working on non-quantum steps, getting denominator from x-tilde and trying to find real period
 11/03/21       0.4.1  Added a few functions to assist with finding the period, no success yet
 12/03/21       0.5.0  Can use quantum method to factor :)
 22/03/21       1.0.0  Implemented xmalloc (malloc wrapper)
 22/03/21       1.0.0  Comments on all fucntions informing user of usage.
 22/03/21       1.0.0  Reason for getting 1 and 15 is bc testP must be wrong see step 5 in paper
 23/03/21       1.1.0  Pass wavefunction and gates by reference to swapsies so constantly updating wf
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
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>
#include "quantcomp.h"

#define BASIS 2
#define N 9 // Number of qubits defined
#define STATES_MAX 1024 //max of 10 qubits 
#define L 3
#define M 4


struct timeval tv;

typedef struct primeFactors{
    int a;
    int b;
}primeFactors;

void* xmalloc(int bits){
    void* ptr = malloc(bits);
    
    if (ptr == NULL)
    {
        printf("Memory exhausted");
        exit(0);
    }
    
    return ptr;
}
void qubit_error(){
    printf("Please operate the gate on a valid qubit\n");
    exit(1);
}

gsl_vector_complex* createRegister(int qubits){
    if(qubits < 1)
    {
        fprintf(stderr, "Error: Qubit number must be greater than 0\n");
        exit(0);
    }
    if(qubits > 14){
        fprintf(stderr, "Warning: Computations with more than 14 qubits are memory and time intensive");
    }
    gsl_vector_complex* quantumRegister = gsl_vector_complex_alloc(pow(2, qubits));
    gsl_vector_complex_set_basis(quantumRegister, 0);

    return quantumRegister;
}
// Simply prints the probability amplitudes of each state.
//
// Arguments
// ---------
// [1] Wavefunction -> Contains probability amplitudes to be printed.
// 
// Returns
// -------
// void
void print_wf(gsl_vector_complex* wavefunction){
    for (int i = 0; i < wavefunction->size; i++){
        printf("%lg + %lgi\n", GSL_REAL(gsl_vector_complex_get(wavefunction, i)), GSL_IMAG(gsl_vector_complex_get(wavefunction, i)));
    }printf("\n");
}
// Hand built function based on cblas functions to multiply non sparse complex matrix and vecotrs and function used to
// multiply sparse matrices in COO format. Possible extensions in the future could be to allow this func to manage CCS and CRS 
// formated matrices
//
// Argumnents
// ----------
// [1] CBLAS_TRANSPOSE_t -> is the matrix to be transposed or not, legacy feature.
// [2] A -> sparse complex matrix to be multiplied
// [3] x -> complex vector to be multiplied
// [4] y -> Where the procuct of A and x is stored.
//
// Returns
// -------
// [1] y -> The product of the matrix A and vector x
gsl_vector_complex* myMulFunc(const CBLAS_TRANSPOSE_t TransA,
                 gsl_spmatrix_complex *A, gsl_vector_complex *x,
                 gsl_vector_complex *y)
{
  const size_t P = A->size1;
  const size_t K = A->size2;

  if ((TransA == CblasNoTrans && K != x->size) ||
      (TransA == CblasTrans && P != x->size))
    {
        return NULL;
      //GSL_ERROR("invalid length of x vector", GSL_EBADLEN);
    }
  else if ((TransA == CblasNoTrans && P != y->size) ||
           (TransA == CblasTrans && K != y->size))
    {
        return NULL; 
        //GSL_ERROR("invalid length of y vector", GSL_EBADLEN);
    }
  else
    {
      size_t incX;
      size_t lenX, lenY;
      double *X, *Y;
      double *Ad;
      int *Ai, *Aj;
      int p;

      if (TransA == CblasNoTrans)
        {
          lenX = K;
          lenY = P;
        }
      else
        {
          lenX = P;
          lenY = K;
        }

      /* form y := op(A)*x */
      Ad = A->data;
      X = x->data;
      incX = x->stride;
      Y = y->data;

      if (GSL_SPMATRIX_ISTRIPLET(A))
        {
          if (TransA == CblasNoTrans)
            {
              Ai = A->i;
              Aj = A->p;
            }
          else
            {
              Ai = A->p;
              Aj = A->i;
            }

          for (p = 0; p < (int) A->nz; ++p)
            {
            double ar = Ad[2*Ai[p]];
            double ai = Ad[2*Ai[p]+1];

            double xr = X[2*p*incX];
            double xi = X[2*p*incX+1];

            // printf("Multiplying %.3lg with %.3lg and %.3lg with %.3lg for real part of element\n", ar, xr, ai, xi);
            // printf("Multiplying %.3lg with %.3lg and %.3lg with %lg for imag part of element\n", ar, xi, ai, xr);

            gsl_vector_complex_set(y, Aj[p], gsl_complex_rect((ar * xr - ai * xi), (ar * xi + ai * xr)) );

            }
        }
      else
        {
          return NULL;
          //GSL_ERROR("unsupported matrix type", GSL_EINVAL);
        }

      return y;
    }
}

//   Creates a wavefunction where all spins are "down" i.e. |00..0>.
//  A 1 in the ket represents a spin down. e.g. |010> means all qubits are spin down apart from qbit 2.  Note how 010 is 2 in binary.
//  Since 000 is 0 in binary, we set probability amplitude of state 0 to 1.
//  Arguments
//  ---------
// [1] states -> The number of total states for the system. This depends on the number of quibits defined and the basis vector.
//
//  Returns
//  ---------
// [1] wavefunction of gsl_vector_complex type initialised to probability of 1 for all qubits in spin down state |000>
void initWavefunctionSpinDown(gsl_vector_complex* wavefunction){ 
    // Probability amplitude stored in "wavefuction" vector
    gsl_vector_complex_set_basis(wavefunction, 0);

}

// Set the probability amplitudes of all states to 1/sqrt(basis^N) so equal probability of wavefunction collapsing into any state.
//  Arguments
//  ---------
// [1] states -> The number of total states for the system. This depends on the number of quibits defined and the basis vector.]
//
//  Returns
//  ---------
//  [1] wavefunction of gsl_vector_complex type initialised so that there is equal chance to find the system in every state.
gsl_vector_complex* initWavefunctionEqualProb(int states){ //Initialising wf to "Equal Probability" of all states
    
    gsl_vector_complex* wavefunction = NULL;

    wavefunction = gsl_vector_complex_alloc(states);

    gsl_vector_complex_set_all(wavefunction , gsl_complex_rect(1/sqrt(pow(BASIS, N)), 0));

    return wavefunction;
}

// Initialises the wavefunction of the system to |00...01> as required by Shor's Algorithm
// Arguments
// ---------
// [1] States -> Number of possible quantum states 
//
// Returns
// -------
// [1] wavefunction -> wf initialised to |00...01> 

gsl_vector_complex* initWavefunctionShors(int states){ //Initialising wf to "Equal Probability" of all states
    
    gsl_vector_complex* wavefunction = NULL;

    wavefunction = gsl_vector_complex_alloc(states);

    gsl_vector_complex_set(wavefunction, 1, GSL_COMPLEX_ONE);
    return wavefunction;
}

// Takes an integer and return binary representation in string format
//  Arguments
//  ---------
// [1] a -> The integer number that is to be turned into a binary string
//
//  Returns
//  ---------
//  [1] a string that is the binary representation of the integer a
char* intToBinary(int a){

    int bin = 0;

    int remainder, temp = 1;

    while(a != 0){

        remainder = a % 2;
        a /= 2;
        bin += remainder*temp;
        temp *= 10;
    }
    char *bin_str = (char *) xmalloc(N*sizeof(char));

    sprintf(bin_str, "%09d", bin);   

    return bin_str;

}

// Uses pointers to update the wavefunction, rather than using memcopy which has complexity 0(n^2)
//
// Arguments
// ---------
// [1] new_state -> The state that we wavefunction will be updated to
// [2] wavefunction -> The wavefuction that describes the system that is to be updated

// Returns
// -------
// [1] An updated wavefunction of type gsl_vector_complex
void swapsies(gsl_vector_complex* new_state, gsl_vector_complex* wavefunction){
    //update the old wf with the new 

    gsl_vector_complex_swap(new_state, wavefunction);
    gsl_vector_complex_free(new_state);
}

// A function that alllows one to multiply a complex vector and a real sparse matrix, involves splitting vector into two "views"
//
// Arguments
// ---------
// [1] vector -> The vector that is to be multipled
// [2] matrix -> Matrix to be multiplied
//
// Returns
// -------
// [1] mv -> matrix*vector result of type gsl_vector_complex
gsl_vector_complex* complexVecMultRealMat(gsl_vector_complex* vector, gsl_spmatrix* matrix){
    
    // Splitting vector up into real and imag
    gsl_vector_view real = gsl_vector_complex_real(vector); // Temp object on the stack so not mem intensive

    gsl_vector_view imag = gsl_vector_complex_imag(vector);


    gsl_vector_complex* mv = gsl_vector_complex_alloc(vector->size);

    gsl_vector_view new_real = gsl_vector_complex_real(mv);

    gsl_vector_view new_imag = gsl_vector_complex_imag(mv);

    gsl_spblas_dgemv(CblasNoTrans, 1.0, matrix, &real.vector, 0.0, &new_real.vector);

    gsl_spblas_dgemv(CblasNoTrans, 1.0, matrix, &imag.vector, 0.0, &new_imag.vector);

    return mv;
}
// Get unitary matrix in sparse format COO
//
// Arguments
// ---------
// [1] matrix -> Matrix to become unit amtrix
//
// Returns
// -------
// [1] matrix -> unit matrix in sparse form
gsl_spmatrix* spmatrixIdentity(gsl_spmatrix* matrix){
    
    int i = 0;
    int j = 0;

    while(i < matrix->size1 && j < matrix->size2)
    {    
        gsl_spmatrix_set(matrix, i ,j, 1);
        i++;
        j++;
    
    }

    return matrix;
}
//  Arguments
//  ---------
// [1] matrix -> The matrix that is to be printed
void print_matrix(gsl_spmatrix_complex* matrix){

    for(int i = 0; i < matrix->size1;i++){
    
        for(int j = 0; j < matrix->size2; j++){
    
            printf("%.3lg + %.3lgi\t", GSL_REAL(gsl_spmatrix_complex_get(matrix, i, j)), GSL_IMAG(gsl_spmatrix_complex_get(matrix, i, j)));
    
        }printf("\n");
    
    }
}


// To measure the quantum state we must use the discrete probability distribution provided by the 
// wavefuntion. The measurement function finds the probabilities of each state and observes the wf
// according to those probabilities using a random number generator. After measurement the wavefunction is "collapsed" and gives the 
// measurement over and over again.

//  Arguments
//  ---------
// [1] Wavefunction -> Defines the state of the system with probability amplitudes of possible configurations
char* measureRegisterGate(gsl_vector_complex* wavefunction){
    int states = (int) wavefunction->size;
    double probabilities[states];
    
    for(int i = 0; i < states; i++){ //creating a list of the probabilities (non-normalised)
        
        probabilities[i] = GSL_REAL(gsl_vector_complex_get(wavefunction,i))*GSL_REAL(gsl_vector_complex_get(wavefunction,i)) + GSL_IMAG(gsl_vector_complex_get(wavefunction,i))*GSL_IMAG(gsl_vector_complex_get(wavefunction,i));
    
    }

    gsl_ran_discrete_t* lookup = gsl_ran_discrete_preproc(states, probabilities); // Preproc normalises the probabilities given by wf
    
    gettimeofday(&tv,NULL); // Get time of day in usec so that the seed changes giving a different stream of #'s each time
    
    unsigned long int seed = tv.tv_usec; // Seed depends on time of day so repeates every 24 hours.  
    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937); // Mersene Twister Algorithm is used for high quality random numbers
    gsl_rng_set(r, seed); 
   
    size_t t = gsl_ran_discrete(r, lookup); // Choosing from the discrete probability distribution defined by the wavefunction probability amplitudes
   
    // Wavefunction collapsed so will only find system in the state from now on
    printf("Wavefunction collapsed into the state:\n|%s>\n", intToBinary(t));
    
    gsl_vector_complex_set_all(wavefunction, GSL_COMPLEX_ZERO);
    gsl_vector_complex_set(wavefunction, t, GSL_COMPLEX_ONE); // Set measured state to probability one so that if we measure again we get the same outcome
    
    // Free memory to avoid bloats
    gsl_ran_discrete_free(lookup);
    gsl_rng_free(r);
    
    return intToBinary(t);

}

// This function will find the element of the tensor product for a given gate for one qubit

//  Arguments
//  ---------
// [1] inta -> the row of the desired element in base 2 
// [2] intb -> the col of the desired element in base 2
// [3] qubit -> specifies which qubit the Hadamard gate will be acting upon.

//  Returns
//  ---------
//  [1] The value of the abth element of the corrresponding hadamard gate
double findElementHad(char* inta, char* intb, int qubit){
    
    // Hadamard gate for single qubit used to calculate tensor product
    
    gsl_matrix *hadamard_single = gsl_matrix_alloc(BASIS, BASIS);
    gsl_matrix_set_all(hadamard_single, 1/sqrt(BASIS));
    gsl_matrix_set(hadamard_single,1,1, -1/sqrt(BASIS));

    double value = 1.0;

    for(int i = 0; i < N; i++){

        if(inta[i] != intb[i] && i != qubit - 1){
            
            // Invokes Kronecker delta

            return 0.0; 

        }

        value =  gsl_matrix_get(hadamard_single, inta[qubit-1] - '0', intb[qubit -1] - '0');

    }

    gsl_matrix_free(hadamard_single);

    return value;

}

//  This function calculates values of the phase matrix for a system of arbitrary size in an element-wise method.

//  Arguments
//  ---------
// [1] inta -> the row of the desired element in base 2 
// [2] intb -> the col of the desired element in base 2
// [3] qubit -> specifies which qubit the Phase gate will be acting upon.
// [4] phi -> The angle of rotation

//  Returns
//  ---------
//  [1] The value for the abth element of the corresponding phase matrix
gsl_complex findElementPhase(char* inta, char* intb, int qubit, double phi){
    
    // Phase gate for single qubit used to calculate tensor product
    gsl_matrix_complex *phase_single = gsl_matrix_complex_alloc(BASIS, BASIS);
    gsl_matrix_complex_set_identity(phase_single);
    gsl_matrix_complex_set(phase_single,1,1, gsl_complex_polar(1,phi));

    gsl_complex value = gsl_complex_rect(1,0);

    for(int i = 0; i < N; i++){

        if(inta[i] != intb[i] && i != qubit - 1){

            return GSL_COMPLEX_ZERO; // Invokes Kronecker delta
        
        }
        
        value =  gsl_matrix_complex_get(phase_single, inta[qubit-1] - '0', intb[qubit -1] - '0');

    }
    
    gsl_matrix_complex_free(phase_single);
    
    return value;
}

//  This function calculates values of the cnot matrix for a system of arbitrary size in an element-wise method.

//  Arguments
//  ---------
// [1] row -> the row of the desired element in base 2 
// [2] col -> the col of the desired element in base 2
// [3] control_qubit -> specifies which qubit acts as the control
// [4] target_qubit -> specifies target qubit
// [5] num_qubits -> total number of qubits

//  Returns
//  ---------
//  [1] The value of the row-columnth element of the corresponding CNOT gate
double findElementCnot(char* row, char* col, int control_qubit, int target_qubit, int num_qbits){
    // Defining the two qubit cnot gate that we will draw values from
    gsl_matrix *cnot = gsl_matrix_alloc(BASIS*BASIS, BASIS*BASIS);
    
    gsl_matrix_set_zero(cnot);
    gsl_matrix_set(cnot, 0, 0, 1);
    gsl_matrix_set(cnot, 1, 1, 1);
    gsl_matrix_set(cnot, 2, 3, 1);
    gsl_matrix_set(cnot, 3, 2, 1);

    for(int i = 0; i < num_qbits; i++){

        // Employ the deltas first

        char char1 = row[i];
        char char2 = col[i];

        if(i != target_qubit-1 && i != control_qubit-1 && char1!=char2){

            // If an element of row and col strings not a match on any element other that targ or contr then we know 
            // a delta will set the whole element to zero

            return 0.0;
        }
   
    }
   
    // If deltas dont cause the value to be zero then we find corresponding element from cnot matrix
   
    char str1[] = "00";
    char str2[] = "00";
   
    // Pushing the control qubit to the front of row and col string
   
    str1[0] = row[control_qubit-1];
    str1[1] = row[target_qubit-1];
    str2[0] = col[control_qubit-1];
    str2[1] = col[target_qubit-1];
   
    // Need to go from binary to dec to find the element
   
    long row_index = strtol(str1, NULL, 2);
   
    long col_index = strtol(str2, NULL, 2);
   
    double value = gsl_matrix_get(cnot, row_index, col_index);
   
    gsl_matrix_free(cnot);
   
    return value;
}

// Gate used to calculate f(x) on the l register
// Arguments
// ---------
// [1] control -> which qubit is the control qubit 
// [2] a -> random number between 1 and C
// [2] C -> Composite number to be factorised
// [3] wavefunction -> The wavefunction that describes the system that is to be updated.
//
// Returns
// -------
// [1] -> updated wavefunction after multiplying by f(x)
void amodcGate(int control, int a, int C, gsl_vector_complex* wavefunction){ // only need to cycle through rows as permutation matrix
    gsl_vector_complex* newState = gsl_vector_complex_alloc(wavefunction->size);

    double A;
    double A0 = (double) (a % C);
    double A1 = (double) ((int)pow((double)a, 2.0) % C);
    double A2 = (double) ((int)pow((double)a, 4.0) % C);

    if (control == 3)
    {
        A = A0;
    }
    if (control == 2)
    {
        A = A1;
    }
    if (control == 1)
    {
        A = A2;
    }
    
    gsl_spmatrix* amodx = gsl_spmatrix_alloc(128, 128);

    for(int k = 0; k < amodx->size2; k++){
        
        char* binK = intToBinary(k);
        
        if(binK[control-1] - '0' == 0){

            gsl_spmatrix_set(amodx, k, k, 1.0);

        }

        if (binK[control-1] - '0' != 0)
        {
            
            char* f = xmalloc(4);
            strlcpy(f, binK+3, 4+1); // Taking the substring m3m2m1m0 from l2l1l0m3m2m1m0
                        
            if(  (int)strtol(f, NULL, 2)  >= C){
               
                gsl_spmatrix_set(amodx, k, k, 1.0);
            
            }
            else if(binK[control-1] - '0' != 0 && (int)strtol(f, NULL, 2) < C)
            {
                int fprime = (int)A*strtol(f, NULL, 2) % C; // f' = f*A*mod(C) 

                char* binFprime = xmalloc(N); // a char is one byte

                binFprime = intToBinary(fprime);

                char *l = xmalloc(N);
                
                strlcpy(l, binK, 3+1);
                
                strncat(l, binFprime+3, 4);
                
                int j = (int)strtol(l, NULL, 2);

                gsl_spmatrix_set(amodx, j, k, 1.0);
                
                free(l);
                free(binFprime);
                free(f);
            }
            free(binK);

        }
    }
    newState = complexVecMultRealMat(wavefunction, amodx);
    swapsies(newState, wavefunction);

}

//  This function calculates values of the cphase matrix for a system of arbitrary size in an element-wise method.

//  Arguments
//  ---------
// [1] row -> the row of the desired element in base 2 
// [2] col -> the col of the desired element in base 2
// [3] control_qubit -> specifies which qubit acts as the control
// [4] target_qubit -> specifies target qubit
// [5] num_qubits -> total number of qubits

//  Returns
//  ---------
//  [1] The value of the row-columnth element of the corresponding cphase gate
gsl_complex findElementCphase(char* row, char* col, int control_qubit, int target_qubit, int num_qbits, double phase){
    // Defining the two qubit cphase gate that we will draw values from
    gsl_spmatrix_complex *cphase = gsl_spmatrix_complex_alloc(BASIS*BASIS, BASIS*BASIS);
    
    gsl_spmatrix_complex_set_zero(cphase);
    gsl_spmatrix_complex_set(cphase, 0, 0, GSL_COMPLEX_ONE);
    gsl_spmatrix_complex_set(cphase, 1, 1, GSL_COMPLEX_ONE);
    gsl_spmatrix_complex_set(cphase, 2, 2, GSL_COMPLEX_ONE);
    gsl_spmatrix_complex_set(cphase, 3, 3, gsl_complex_polar(1, phase));

    for(int i = 0; i < num_qbits; i++){

        // Employ the deltas first

        char char1 = row[i];
        char char2 = col[i];

        if(i != target_qubit-1 && i != control_qubit-1 && char1!=char2){

            // If an element of row and col strings not a match on any element other that targ or contr then we know 
            // a delta will set the whole element to zero

            return GSL_COMPLEX_ZERO; // may not need to return this as using spmatrix and these are all initialised to zero anyway!
        }
   
    }
   
    // If deltas dont cause the value to be zero then we find corresponding element from cphase matrix
   
    char str1[] = "00";
    char str2[] = "00";
   
    // Pushing the control qubit to the front of row and col string
   
    str1[0] = row[control_qubit-1];
    str1[1] = row[target_qubit-1];
    str2[0] = col[control_qubit-1];
    str2[1] = col[target_qubit-1];
   
    // Need to go from binary to dec to find the element location
   
    long row_index = strtol(str1, NULL, 2);
   
    long col_index = strtol(str2, NULL, 2);
   
    gsl_complex value = gsl_spmatrix_complex_get(cphase, row_index, col_index);

   
    gsl_spmatrix_complex_free(cphase);
   
    return value;
}

// Controlled phase gate. Architecture based of Cnot gate
// Arguments
// ---------
// [1] wavefunction -> The wavefunction that describes the system that is to be updated.
// [2] control -> defines control qubit
// [3] target -> defines target qubit
// [4] phase -> defined the angle \theta used in e^{i\theta}
//
// Returns
// -------
// [1] -> wavefunction -> updated wavefunction
void CphaseGate(gsl_vector_complex* wavefunction, int control, int target, double phase){

    /* Will need to only put val if val is non-zero due to sparse nature, need to use myMulFunc in order to multiply
    complex matrix by complex vector */

    gsl_spmatrix_complex* cphasegate = gsl_spmatrix_complex_alloc(wavefunction->size, wavefunction->size);
    for(int i = 0; i < wavefunction->size; i++){
        for(int j = 0; j < wavefunction->size; j++){
           
            gsl_complex val = findElementCphase(intToBinary(i), intToBinary(j), control, target, N, phase);
           
            if(GSL_REAL(val) == 0 && GSL_IMAG(val) == 0){
                
                /* Do Nothing */
                
            }
            else{
                gsl_spmatrix_complex_set(cphasegate, i , j, val);
            }
        }
    }
    gsl_vector_complex* r_psi = gsl_vector_complex_alloc(wavefunction->size);

    myMulFunc(CblasNoTrans, cphasegate, wavefunction, r_psi); // complex sparse multiplier, works with COO only

    /* Maybe use the same block of memory for the matrix and vectors? */
    swapsies(r_psi, wavefunction);

}

//  The Hadamard gate operates on the register by setting qubits into a superposition of their basis states. 
//  This function does this by allowing the user to specify which qubit we will be setting to a superposition. 
//  This will have effects when it comes to measuring the whole register as sometimes that qubit will be spin up, 
//  sometimes it will be spin down, which will change the overall state of the register.

//  Arguments
//  ---------
// [1] Wavefunction -> Defines the state of the system with probability amplitudes of possible configurations
// [2] qubit -> specifies which qubit the Hadamard gate will be acting upon.
//  Returns
//  ---------
//  [1] Wavefunction after entire register has been acted on by desired hadamard gate.
void hadamardGate(gsl_vector_complex* wavefunction, int qubit){
    
    if(qubit > N){
    
        printf("Please operate the gate on a valid qubit\n");
    
        exit(0);
    
    }
    
    // Will beome the NxN matrix for operation on whole register
    gsl_spmatrix *hadamard = gsl_spmatrix_alloc(wavefunction->size, wavefunction->size);

    for(int i = 0; i < wavefunction->size; i++){

        for(int j = 0; j < wavefunction->size; j++){

            double val = findElementHad(intToBinary(i), intToBinary(j), qubit);
            gsl_spmatrix_set(hadamard, i , j, val);
       
        }

    }

    gsl_vector_complex* h_psi = gsl_vector_complex_alloc(wavefunction->size);

    h_psi = complexVecMultRealMat(wavefunction, hadamard);
    
    swapsies(h_psi, wavefunction); ; //updates wf to new state using pointers rather than memcopy
}

//  A phase gate does not alter the probabilites of finding the the system in a given state, 
//  however it does alter the relative phase between the states. Quantum phase cannot be measured
//  directly, however it can have effects if more gates are subsequently applied before register measuremnt.

//  Arguments
//  ---------
// [1] Wavefunction -> Defines the state of the system with probability amplitudes of possible configurations
// [2] qubit -> specifies which qubit the phase shift gate is operating on.
// [3] phase -> the angle of rotation 

//  Returns
//  ---------
//  [1] Wavefunction after entire register has been acted upon by desired phase shift gate
void phaseShiftGate(gsl_vector_complex *wavefunction, int qubit, double phase){
    if(qubit > N){
        printf("Please operate the gate on a valid qubit\n");
        exit(0);
    }
    // Will beome the NxN matrix for operation on whole register
    gsl_spmatrix_complex *phaseGate = gsl_spmatrix_complex_alloc(wavefunction->size, wavefunction->size);
    

    for(int i = 0; i < wavefunction->size; i++){
        for(int j = 0; j < wavefunction->size; j++){
            gsl_complex val = findElementPhase(intToBinary(i), intToBinary(j), qubit, phase);
            if(GSL_REAL(val) == 0 && GSL_IMAG(val) == 0){
                /* Do Nothing */
            }
            else{
                gsl_spmatrix_complex_set(phaseGate, i , j, val);
            }
        }
    }
    gsl_vector_complex* r_psi = gsl_vector_complex_alloc(wavefunction->size);

    print_matrix(phaseGate);
    myMulFunc(CblasNoTrans, phaseGate, wavefunction, r_psi); // complex sparse multiplier, works with COO only
    
    swapsies(r_psi, wavefunction);
}

// Quantum oracle gate used in grovers quantum search algorithm. Argument answer is the "Correct question" mentioned
// in paper
// Arguments
// ---------
// [1] wavefunction -> 
// [2] answer -> which state is the "correct" states
//
// Returns
// -------
// [1] wavefunction -> updated wf 
void oracleGate(gsl_vector_complex* wavefunction, int answer){
    
    gsl_spmatrix* oracleGate = gsl_spmatrix_alloc(wavefunction->size, wavefunction->size);
    
    spmatrixIdentity(oracleGate);
    gsl_spmatrix_set(oracleGate, answer-1, answer-1, -1); //Minus one as index from 0
    
    gsl_vector_complex* o_psi = gsl_vector_complex_alloc(wavefunction->size);
    o_psi = complexVecMultRealMat(wavefunction, oracleGate);

    swapsies(o_psi, wavefunction);

}

// Unity matrix of size 2^N*2^N with -1 in the 0,0th element
//
// Arguments
// ---------
// [1] wavefunction -> Defines the state of the system with probability amplitudes of possible configurations
//
// Returns
// -------
// [1] wavefunction -> wavefunction after muliplication with j gate
void jGate(gsl_vector_complex* wavefunction){
    gsl_spmatrix* jGate = gsl_spmatrix_alloc(wavefunction->size, wavefunction->size);
    spmatrixIdentity(jGate);
    gsl_spmatrix_set(jGate, 0, 0, -1);
    
    gsl_vector_complex* j_psi = gsl_vector_complex_alloc(wavefunction->size);
    j_psi = complexVecMultRealMat(wavefunction, jGate);

    swapsies(j_psi, wavefunction);
}

// An ensemble of gates that are used within grovers block

// Arguments
// ---------
// [1] wavefunction -> Defines the state of the system with probability amplitudes of possible configurations
// [2] Answer -> The right choice in grovers alg

// Returns
// -------
// [1] wavefunction -> wavefunction after muliplication with j gate
void groversBlock(gsl_vector_complex* wavefunction, int answer){
    gsl_vector_complex* b_psi = gsl_vector_complex_alloc(wavefunction->size);
    gsl_vector_complex_set_zero(b_psi);
    //First operation in the block is to apply a quantum oracle gate
    oracleGate(wavefunction, answer); 
    //Then we apply a hadamard gate to each gate
    for(int i = 1; i < N+1; i++){
        hadamardGate(b_psi, i); // +1 bc had takes in canonical qubit number.       
    }
    // Apply the j gate
    jGate(b_psi);
    // Finally a hadamard gate on each qubit again
    for(int j = 1; j < N+1; j++){
        hadamardGate(b_psi, j);        
    }
    swapsies(b_psi, wavefunction);
}

// The controlled-not gate operates on two qubits. The target and the control. If the control is |1>, the target quibit is flipped
// irrespective of its value. If the control qubit is |0> the target remains unchanged. 
//
// Arguments
// ---------
// [1] Wavefunction -> Defines the state of the system with probability amplitudes of possible configurations
// [2] Control -> Specifies the control qubit (Min: 0 ; Max: 2^N-1 )
// [3] target -> Specifies the target qubit (Min: 0 ; Max: 2^N-1 )
//
// Returns
// -------
// [1] Wavefunction -> Wavefunction after manipulation by the cnot gate
void cnotGate(gsl_vector_complex* wavefunction, int control, int target){
    gsl_vector_complex* c_psi = gsl_vector_complex_alloc(wavefunction->size);
    gsl_vector_complex_set_zero(c_psi);

    gsl_spmatrix *cnot = gsl_spmatrix_alloc(wavefunction->size, wavefunction->size);
    
    for(int i = 0; i < wavefunction->size; i++){
        for(int j = 0; j < wavefunction->size; j++){
            gsl_spmatrix_set(cnot, i, j, findElementCnot(intToBinary(i), intToBinary(j), control, target, N));
        }
    }
    c_psi = complexVecMultRealMat(wavefunction, cnot);
    swapsies(c_psi, wavefunction);
}

// Given a composite number, i.e. the product of two primes, Shors alg is able to factorise this composite number faster than any classic algorithm. Use a high degree of encapsulation here.
//
// Arguments
// ---------
// [1] composite_number -> number that is to be factored
//
// Returns
// -------
// [1] Factor
int isPower(int number){
    
    if (number % 2 == 0) // Checking if number is even or power of two
    {
        return 2;
    }
    
    for(int i = 3; i < (int)(number / 3); i++){ //Starting with seeing if 3^n == number for a whole number n, then to 4^n etc
        
        double n = log(number) / log(i);
        
        if(ceil(n) == n){ //Checking if power of a smaller number
            return pow(i,n);
        }
    }
    return 1;
}

// Utilising Euclids algorithm (recursive) to find GCD
//
// Arguments
// ---------
// [1] n1 -> number 1
// [2] n2 -> number 2
//
// Returns
// -------
// [1] greatest common divisor of n1 and n2
int greatestCommonDivisor(int n1, int n2){

    if(n2 != 0){
        return greatestCommonDivisor(n2, n1 % n2); //Call until n2 is set to zero
    }
    else
    {
        return n1;
    }
}

// Given a string, this function will reverse it. Used in finding xtilde
// Arguments
// ---------
// [1] string -> to be reversed
//
// Returns
// -------
// [1] rev -> reversed string

char* reverseString(const char* string){
    
    char* rev;
    int begin, end, count;
    
    count = strlen(string);

    rev = xmalloc(count);

    end = count - 1;

    for(begin = 0; begin < count; begin++)
    {
        rev[begin] = string[end];
        end--;
    }

    // Not forgetting the terminating charecter
    rev[begin] = '\0'; 

    return rev; 
}
// Obtains the values of qubits in the x register by measuring the wf, then turns result into the IQFT.
// Arguments
// ---------
// [1] wavefunction -> Describes state of the system.
//
// Returns
// -------
// [1] xtilde -> The IQFT of the x register
int readsXReg(gsl_vector_complex* wavefunction){
    //Wavefunction needs to be collapsed in order for this to work properly! Add a measure gate just to be sure it is collapsed

    char* state;
    char* binXTransformed;
    // After wf is collapsed we measure which state it's in, turn this number into binary, then read the first 3 digits in reverse
    // and convert back into base 10 to find the fourier trans of x
    state = measureRegisterGate(wavefunction);

    
    binXTransformed = xmalloc(3); // use xmalloc??
    if(binXTransformed == NULL){
        printf("Malloc failed");
    }

    strlcpy(binXTransformed, state, 4); // 4 to include terminating charecter for security as apposed to strncpy()
    
    // Need to reverse binXTransformed
    binXTransformed = reverseString(binXTransformed);

    // printf("%s\n", binXTransformed);
    
    int xTilde = (int)strtol(binXTransformed, NULL, 2);
    
    printf("x tilde = %d\n", xTilde);
    
    free(binXTransformed);
    
    return xTilde;

}

// Gives numerator and demoninator of a fraction
// Arguments
// ---------
// [1] number -> decimal number
// [2] ->
//
// Returns
// -------
// [1] ints -> returns a type containing a pair of ints. In this case they represent the numerator and denominator 
primeFactors decimalToFraction(double decimalNum){

    primeFactors ints;

    double intVal = floor(decimalNum);
    double decVal = decimalNum - intVal;
    double pVal = 1000000000;
    int gcd = greatestCommonDivisor(round(decVal * pVal), pVal);

   double num = round(decVal * pVal) / gcd;

   double deno = pVal / gcd;

   ints.a = intVal * deno + num; //top
   ints.b = deno; //bottom

   return ints;
}

//  Checks to see whether the period generated is suitable according to rules set out by 
// D. Candela
// Arguments
// ---------
// [1] a -> Random number chosen between 1 and C
// [2] p -> period attempt
// [3] C -> Composite number
//
// Return
// ------
// [1] 0 if p is not suitable, in this case go back and pick a new a!
// [2] 1 if p is all good for usage
int testP(int a, int p, int C){

    double congruence = (double)(((int)pow((double)a, (double)p/2)) % C);
    printf("congruence = %lg\n", congruence);
    double minusOneCongruence = 14;
    // printf("%lg", congruence);
    
    if((p % 2 == 0 /*is even*/ && congruence == minusOneCongruence) || (p % 2 != 0) /* ensuring a^p/2 != -1 mod(C)*/ ){
        // If this is the case we need to pick a new "a" or rand, as i've called it in shors().
        return 1;
    }
    else
    {
        return 0;
    }
    
}

// Chooses a random int between 1 and C.
// Arguments
// ---------
// [1] compositeNumber -> known as C in rest of code. The number that will be factored.
//
// Returns
// -------
// [1] rand -> Random int between 1 and compositeNumber 
int getRand(int compositeNumber){
    gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); // High quality random number generator

    gettimeofday(&tv,NULL);
    unsigned long int seed = tv.tv_usec; // Seed depends on time of day so repeates every 24 hours.  
    
    gsl_rng_set(r, seed);
    
    int rand = 0;
    while (rand < 2)
    {
        rand = gsl_rng_uniform_int(r, compositeNumber); // Goes until rand is in range 1 < rand < C
    }
    // printf("rand %d\n", rand);
    gsl_rng_free(r);
    
    return rand;
}

// All gates needed to set up Shor's algorithm
// Arguments
// ---------
// [1] wavefunction -> Description of system to be updated
// [2] rand -> the "a" for a given attempt of shor's alg
// [3] compositeNumber -> The number to be factored.

// Returns
// -------
// [1] wavefunction -> updated wavefunction 
void shorsBlock(gsl_vector_complex* wavefunction, int rand, int compositeNumber){
    // Applying hadamard gates to the x register
    wavefunction = initWavefunctionShors(wavefunction->size);
    for(int xQubitHad = 1; xQubitHad < M; xQubitHad++){
     
        hadamardGate(wavefunction, xQubitHad);
    
    }
    // Multiplying f register by f(x)
    for(int xQubitCGate = L; xQubitCGate > 0; xQubitCGate--){
    
        amodcGate(xQubitCGate, rand, compositeNumber, wavefunction);
    
    }
    // IQFT block--------------------------------------------
    hadamardGate(wavefunction, 1); 
    CphaseGate(wavefunction, 1, 2, M_PI_2);
    CphaseGate(wavefunction, 1, 3, M_PI_4);
    hadamardGate(wavefunction, 2);
    CphaseGate(wavefunction, 2, 3, M_PI_2);
    hadamardGate(wavefunction, 3);


}

// Stitches together multiple functions in order to execute full shor's algorithm
// Arguments
// ---------
// [1] wavefunction ->
// [2] compositeNumber ->
//
// Returns
// -------
// [1] Factors -> non-trivial factors of the composite number 
primeFactors shors(gsl_vector_complex* wavefunction, int compositeNumber){
    
    wavefunction = initWavefunctionShors(wavefunction->size);
    primeFactors factors;
    factors.a = 0;
    factors.b = 0;

    // Step 1 in SA
    int isp = isPower(compositeNumber);
    if(isp != 1){
        factors.a = isp;
        factors.b = (int) compositeNumber / isp;
        printf("factor a = %d \n factor b = %d", factors.a, factors.b);
        return factors;
    
    }
    
    // Step 2 in SA
    // if not a factor of two or a power of another number, pick a random number between 1 and compositeNumber
    int rand = getRand(compositeNumber);
    // Now we have our "a",  we can now carry out the quantum part of shors algorithm
    int gcf = greatestCommonDivisor(rand, compositeNumber);
    
    if( gcf > 1){
        factors.a = gcf;
        factors.b = (int)compositeNumber/gcf;
        printf("factor a = %d \nfactor b = %d\n", factors.a, factors.b);

        return factors;
    }
    // Finding the period p
    shorsBlock(wavefunction, rand, compositeNumber);
    // Measure the wavefunction to collapse is and observe the IQFT of x  
    // print_wf(wavefunction);  
    double omega = (double) readsXReg(wavefunction) / pow(2,L);
    // need an if omega = 0 to reset as cannot obtain a period guess from this as not expressable as a fraction.
    printf("omega try 1 = %lg\n", omega);
    //Now we need to get p from omega! omega = s/p, and then try some factors of s/p i.e. if omega = 0.5 then omega = 1/2 = 2/4 = 4/8, so we try 2,4,8 to fit aP con 1 mod c
    // Remember we need the smallest value of p that satisfies this rule!
    
    // Hopefully this while will ensure we never measure a zero on the x reg. The 0 contains no infomation about the period
    int i = 2;
    while(omega == 0){
        shorsBlock(wavefunction, rand, compositeNumber);
        omega = (double) readsXReg(wavefunction) / pow(2,L);
        printf("omega try %d = %lg\n",i, omega);
        i++;
    }

    // omega = 0.75 = 3/4 = 6/8 
    // omega = 0.25 = 1/4 = 2/8
    // Now we have an omega in the form s/p, we try the values of p proposed.
    int p = decimalToFraction(omega).b; //Extract the denominator from the fraction
    printf("Denominator = %d\n", p);
    // Checking for true period
    for (size_t i = 1; i < N; i++)
    {
        double periodTrial = i*p;
        printf("rand = %d\n", rand);
        double k = (double)((int)pow((double)rand, periodTrial) % compositeNumber);
        double oneCongruence = 1 % compositeNumber;

        printf("Check to see if 1 mod(15) = %d^%lg mod(15)\nk  = %lg and 1 mod(15) = %lg\n",rand, periodTrial, k, oneCongruence);
        
        if(k == oneCongruence) //How we check for congruence, ensuring we have the right period.
        {
            int test = testP(rand, periodTrial, compositeNumber);
            if(test == 0){
                factors.a = greatestCommonDivisor( pow(rand, (periodTrial)/2) + 1, compositeNumber);
                factors.b = greatestCommonDivisor( pow(rand, (periodTrial)/2) - 1, compositeNumber);
                printf("a = %d b = %d\n", factors.a, factors.b);
                return factors;
            }
            else
            {   // Going back and choosing another a
                factors = shors(wavefunction, compositeNumber);
            }
        }        
    }
    return factors;
}


// int main(){
//     int states = (int)pow(BASIS, N);
//     gsl_vector_complex* wavefunction = initWavefunctionShors(states);
//     //Putting system into equal super position of superposition all 2^N basis'
//     // wavefunction = hadamardGate(wavefunction, 1);
//     // wavefunction = hadamardGate(wavefunction, 2); 
//     // wavefunction = hadamardGate(wavefunction, 3);
// //     // wavefunction = cnotGate(wavefunction, 1,2 );
// //     // wavefunction = CphaseGate(wavefunction, 2, 1, M_PI_4);

// //     // for(int i = 0; i < floor(M_PI_4*sqrt(pow(2,N))); i++){ // Needs to be called "floor(pi/4*sqrt(2^N))"" times for optimum output roughly 2 in our case
// //     //     wavefunction = groversBlock(wavefunction, 7); //Second argument is the basis state you want to be "right" in this case its |110>
// //     // }
//     // wavefunction = amodcGate(3, 7, 15, wavefunction);
//     // wavefunction = amodcGate(2, 7, 15, wavefunction);
//     // wavefunction = amodcGate(1, 7, 15, wavefunction);

// // // // IQFT block
//     // wavefunction = hadamardGate(wavefunction, 1);
//     // wavefunction = CphaseGate(wavefunction, 1, 2, M_PI_2);
//     // wavefunction = CphaseGate(wavefunction, 1, 3, M_PI_4);
//     // wavefunction = hadamardGate(wavefunction, 2);
//     // wavefunction = CphaseGate(wavefunction, 2, 3, M_PI_2);
//     // wavefunction = hadamardGate(wavefunction, 3);

// //     // wavefunction = phaseShiftGate(wavefunction, 3,  3.14159);
//     // shors(wavefunction, 15);
//     // print_wf(wavefunction);
//     // readsXReg(wavefunction);
//     // measureRegisterGate(wavefunction);
//     return 0;
// }
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
 Date         Version  Comments  #P.x.y means pending 
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
 1/02/21        0.1.0  Need to Generalise for any number of qubits, N as i kepp getting memory problems
 1/02/21        0.1.0  Fixed for general N, although some refinement needed (line 96)
 17/02/21       0.1.1  Implementing CNOT Gate. Error found. Will squash tomorrow.
 18/02/21       0.1.2  CNOT implemented. NOT 100% correct.
 18/02/21       0.1.3  CNOT implemented successfully
 22/02/21       0.1.4  Trying to implement sparse matrices to speed up
 22/02/21       0.1.4  using vector views etc to improve usability and speed (no need for memcpy)
 22/02/21       0.1.5  Sparse matrices implemented, Need to create complex sparse multiplier
 22/02/21        ""    Some progress making complex sp matrix multiplier, not working 100% yet tho
 23/02/21        ""    Change matrices to ints (rather than doubles to save memory) ALMOST THERE FGS
 1/03/21         ""    Handbuild complex sparse matrix multiplier seems to be working in test program but sets wf to zero here...
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

#define BASIS 2
#define N 3 // Number of qubits defined
#define STATES_MAX 1024 //max of 10 qubits 

struct timeval tv;

void qubit_error(){
    printf("Please operate the gate on a valid qubit\n");
    exit(1);
}

gsl_vector_complex* myMulFunc(const CBLAS_TRANSPOSE_t TransA,
                 gsl_spmatrix_complex *A, gsl_vector_complex *x,
                 gsl_vector_complex *y)
{
  const size_t M = A->size1;
  const size_t L = A->size2;

  if ((TransA == CblasNoTrans && L != x->size) ||
      (TransA == CblasTrans && M != x->size))
    {
        return NULL;
      //GSL_ERROR("invalid length of x vector", GSL_EBADLEN);
    }
  else if ((TransA == CblasNoTrans && M != y->size) ||
           (TransA == CblasTrans && L != y->size))
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
          lenX = L;
          lenY = M;
        }
      else
        {
          lenX = M;
          lenY = L;
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
gsl_vector_complex* initWavefunctionSpinDown(int states){ 
    
    gsl_vector_complex* wavefunction = NULL;
    
    wavefunction = gsl_vector_complex_alloc(states);
    // Probability amplitude stored in "wavefuction" vector
    gsl_vector_complex_set(wavefunction, 0, GSL_COMPLEX_ONE);

    return wavefunction;
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

    gsl_vector_complex_set_all(wavefunction , gsl_complex_rect(1/sqrt(8), 0));

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
    char *bin_str = (char *) malloc(N*sizeof(char));

    sprintf(bin_str, "%03d", bin);   

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
gsl_vector_complex* swapsies(gsl_vector_complex* new_state, gsl_vector_complex* wavefunction){
    //update the old wf with the new 
    gsl_vector_complex** ptr1 = &new_state;

    gsl_vector_complex** ptr2 = &wavefunction;

    gsl_vector_complex* temp;

    temp = *ptr1;
    *ptr1 = *ptr2;
    *ptr2 = temp;

    gsl_vector_complex_free(new_state);

    return wavefunction;
}

// A function that alllows one to multiply a complex vector and a real sparse matrix, involves splitting vector into two "views"
//
// Arguments
// ---------
// [1] vector -> The vector that is to be multipled
// [2] matrix -> MAtrix to be multiplied
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
// get unitary matrix in sparse format COO
//
// Arguments
// ---------
// [1] matrix -> Matrix to become unit amtrix
//
// Returns
// -------
// [1] matrix -> unit matrix
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
void measureRegisterGate(gsl_vector_complex* wavefunction){
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
    
    return;

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
gsl_vector_complex* hadamardGate(gsl_vector_complex* wavefunction, int qubit){
    
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
    
    return swapsies(h_psi, wavefunction); ; //updates wf to new state using pointers rather than memcopy
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
gsl_vector_complex* phaseShiftGate(gsl_vector_complex *wavefunction, int qubit, float phase){ // NEEDS COMPLEX SPARSE MULTIPLIER!
    if(qubit > N){
        printf("Please operate the gate on a valid qubit\n");
        exit(0);
    }
    // Will beome the NxN matrix for operation on whole register
    gsl_spmatrix_complex *phaseGate = gsl_spmatrix_complex_alloc(wavefunction->size, wavefunction->size);
    

    for(int i = 0; i < wavefunction->size; i++){
        for(int j = 0; j < wavefunction->size; j++){
            gsl_complex val = findElementPhase(intToBinary(i), intToBinary(j), qubit, phase); //This is causing some errors
            gsl_spmatrix_complex_set(phaseGate, i , j, val);
        }
    }
    gsl_vector_complex* r_psi = gsl_vector_complex_alloc(wavefunction->size);

    print_matrix(phaseGate);
    myMulFunc(CblasNoTrans, phaseGate, wavefunction, r_psi); //no sparse equivalent so must build one
    
    return swapsies(r_psi, wavefunction);
}

// Oracle gate used in grovers quantum search algorithm. Argument answer is the "Correct question" mentioned
// in paper
gsl_vector_complex* oracleGate(gsl_vector_complex* wavefunction, int answer){
    
    gsl_spmatrix* oracleGate = gsl_spmatrix_alloc(wavefunction->size, wavefunction->size);
    
    spmatrixIdentity(oracleGate);
    gsl_spmatrix_set(oracleGate, answer-1, answer-1, -1); //Minus one as index from 0
    
    gsl_vector_complex* o_psi = gsl_vector_complex_alloc(wavefunction->size);
    o_psi = complexVecMultRealMat(wavefunction, oracleGate);

    return swapsies(o_psi, wavefunction);

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
gsl_vector_complex* jGate(gsl_vector_complex* wavefunction){
    gsl_spmatrix* jGate = gsl_spmatrix_alloc(wavefunction->size, wavefunction->size);
    spmatrixIdentity(jGate);
    gsl_spmatrix_set(jGate, 0, 0, -1);
    
    gsl_vector_complex* j_psi = gsl_vector_complex_alloc(wavefunction->size);
    j_psi = complexVecMultRealMat(wavefunction, jGate);

    return swapsies(j_psi, wavefunction);
}

// An ensemble of gates that are used within grovers block

// Arguments
// ---------
// [1] wavefunction -> Defines the state of the system with probability amplitudes of possible configurations
// [2] Answer -> The right choice in grovers alg

// Returns
// -------
// [1] wavefunction -> wavefunction after muliplication with j gate
gsl_vector_complex* groversBlock(gsl_vector_complex* wavefunction, int answer){
    gsl_vector_complex* b_psi = gsl_vector_complex_alloc(wavefunction->size);
    gsl_vector_complex_set_zero(b_psi);
    //First operation in the block is to apply a quantum oracle gate
    b_psi = oracleGate(wavefunction, answer); 
    //Then we apply a hadamard gate to each gate
    for(int i = 1; i < N+1; i++){
        b_psi = hadamardGate(b_psi, i); // +1 bc had takes in canonical qubit number.       
    }
    // Apply the j gate
    b_psi = jGate(b_psi);
    // Finally a hadamard gate on each qubit again
    for(int j = 1; j < N+1; j++){
        b_psi = hadamardGate(b_psi, j);        
    }
    return b_psi;
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
gsl_vector_complex* cnotGate(gsl_vector_complex* wavefunction, int control, int target){
    gsl_vector_complex* c_psi = gsl_vector_complex_alloc(wavefunction->size);
    gsl_vector_complex_set_zero(c_psi);

    gsl_spmatrix *cnot = gsl_spmatrix_alloc(wavefunction->size, wavefunction->size);
    
    for(int i = 0; i < wavefunction->size; i++){
        for(int j = 0; j < wavefunction->size; j++){
            gsl_spmatrix_set(cnot, i, j, findElementCnot(intToBinary(i), intToBinary(j), control, target, N));
        }
    }
    c_psi = complexVecMultRealMat(wavefunction, cnot);
    return swapsies(c_psi, wavefunction);
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
        printf("%lg + %lgi\n", gsl_vector_complex_get(wavefunction, i).dat[0], gsl_vector_complex_get(wavefunction, i).dat[1]);
    }
}

int main(){
    int states = (int)pow(BASIS, N);
    gsl_vector_complex* wavefunction = initWavefunctionSpinDown(states);
    //Putting system into equal super position of superposition all 2^N basis'
    wavefunction = hadamardGate(wavefunction, 1);
    wavefunction = hadamardGate(wavefunction, 2); 
    wavefunction = hadamardGate(wavefunction, 3);
    // wavefunction = cnotGate(wavefunction, 2, 3);
    // wavefunction = cnotGate(wavefunction, 2, 1);
    // Putting into cat state.

    // for(int i = 0; i < floor(M_PI_4*sqrt(pow(2,N))); i++){ // Needs to be called "floor(pi/4*sqrt(2^N))"" times for optimum output roughly 2 in our case
    //     wavefunction = groversBlock(wavefunction, 7); //Second argument is the basis state you want to be "right" in this case its |110>
    // }
    // print_wf(wavefunction);

    wavefunction = phaseShiftGate(wavefunction, 3,  3.14159);

    print_wf(wavefunction);


    return 0;
}
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
#include <stdatomic.h>

#define ATOMIC sig_atomic_t
#define FUNCTION 
#define TYPE

#define N 5
#define BASIS 2

/*
gsl_spblas_dgemv()
  Multiply a sparse matrix and a vector
Inputs: alpha - scalar factor
        A     - sparse matrix
        x     - dense vector
        beta  - scalar factor
        y     - (input/output) dense vector
Return: y = alpha*op(A)*x + beta*y
*/

int myMulFunc(const CBLAS_TRANSPOSE_t TransA,
                 gsl_spmatrix_complex *A, gsl_vector_complex *x,
                 gsl_vector_complex *y)
{
  const size_t M = A->size1;
  const size_t L = A->size2;

  if ((TransA == CblasNoTrans && L != x->size) ||
      (TransA == CblasTrans && M != x->size))
    {
      GSL_ERROR("invalid length of x vector", GSL_EBADLEN);
    }
  else if ((TransA == CblasNoTrans && M != y->size) ||
           (TransA == CblasTrans && L != y->size))
    {
      GSL_ERROR("invalid length of y vector", GSL_EBADLEN);
    }
  else
    {
      size_t j;
      size_t incX, incY;
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
      incY = y->stride;

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
            // printf("Multiplying %.3lg with %.3lg and %.3lg with %.3lg for imag part of element\n", ar, xi, ai, xr);
            gsl_vector_complex_set(y, Aj[p], gsl_complex_rect((ar * xr - ai * xi), (ar * xi + ai * xr)) );

            }
        }
      else
        {
          GSL_ERROR("unsupported matrix type", GSL_EINVAL);
        }

      return 0;
    }
} /* gsl_spblas_dgemv() */
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

gsl_vector_complex* initWavefunctionEqualProb(int states){ //Initialising wf to "Equal Probability" of all states
    
    gsl_vector_complex* wavefunction = NULL;

    wavefunction = gsl_vector_complex_alloc(states);

    gsl_vector_complex_set_all(wavefunction , gsl_complex_rect(1/sqrt(8), 0));

    return wavefunction;
}

void print_matrix(gsl_spmatrix_complex* matrix){

    for(int i = 0; i < matrix->size1;i++){
    
        for(int j = 0; j < matrix->size2; j++){
    
            printf("%.3lg + %.3lgi\t", GSL_REAL(gsl_spmatrix_complex_get(matrix, i, j)), GSL_IMAG(gsl_spmatrix_complex_get(matrix, i, j)));
    
        }printf("\n");
    
    }
}
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
gsl_vector_complex* phaseShiftGate(gsl_vector_complex *wavefunction, int qubit, float phase){ // NEEDS COMPLEX SPARSE MULTIPLIER!
    if(qubit > N){
        printf("Please operate the gate on a valid qubit\n");
        exit(0);
    }
    // Will beome the NxN matrix for operation on whole register
    gsl_spmatrix_complex *phaseGate = gsl_spmatrix_complex_alloc(wavefunction->size, wavefunction->size);
    gsl_vector_complex* r_psi = gsl_vector_complex_alloc(wavefunction->size);
    

    for(int i = 0; i < wavefunction->size; i++){
        for(int j = 0; j < wavefunction->size; j++){
            gsl_complex val = findElementPhase(intToBinary(i), intToBinary(j), qubit, phase); //This is causing some errors
            gsl_spmatrix_complex_set(phaseGate, i, j, val);
            gsl_spmatrix_complex_get(phaseGate, i, j);
        }
    }

    print_matrix(phaseGate);
    myMulFunc(CblasNoTrans, phaseGate, wavefunction, r_psi); //no sparse equivalent so must build one
    
    return r_psi;
}

int main(){
    
    int states = 4;

    gsl_vector_complex* vector = gsl_vector_complex_alloc(states);
    gsl_spmatrix_complex* matrix = gsl_spmatrix_complex_alloc(states, states);
    
    int i,j = 0;
    while( i < states){
      while( j < states){
        if (j % 2 == 1 && i % 2 ==1){
          gsl_spmatrix_complex_set(matrix, i, j, gsl_complex_rect(-1,0));
        }
        else
        {
          gsl_spmatrix_complex_set(matrix, i, j, gsl_complex_rect(1,0));
        }
        i++;
        j++;
      }
    }
    gsl_vector_complex* vector2 = gsl_vector_complex_alloc(states);

    gsl_vector_complex_set(vector, 0, gsl_complex_rect(2,4));
    gsl_vector_complex_set(vector, 1, gsl_complex_rect(1,5));
    gsl_vector_complex_set(vector, 3, gsl_complex_rect(2,6));

    vector2 = phaseShiftGate(vector, 3, 3.14159); // There's something wrong with my phase shift gate
    //myMulFunc(CblasNoTrans, matrix, vector, vector2);


    for(int j = 0; j < vector2->size; j++){
        printf("%lg + %lgi\n", GSL_REAL(gsl_vector_complex_get(vector2, j)), GSL_IMAG(gsl_vector_complex_get(vector2, j)));
    }

    for(int k = 0; k < 2*matrix->nz; k++){
      printf("%lg\n", matrix->data[k]);
    }

    print_matrix(matrix);

    return 0;
    //complex numbers in sp matricies are not stored as complex data types! but they are split up into an array of doubles!


}

// Y[Ai[p] * incY] += alpha * Ad[p] * X[j * incX];
// Need this in terms of imag # stored in sparse format
// One complex # = gsl_complex_rect(y[0], y[1])

// step one, remove alpha as we wont be using it, it will always be one.
// Will only be using COO storage as it's good enough



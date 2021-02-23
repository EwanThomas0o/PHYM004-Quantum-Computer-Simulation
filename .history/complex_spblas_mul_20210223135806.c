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

int
myMulFunc(const CBLAS_TRANSPOSE_t TransA,
                 gsl_spmatrix_complex *A, gsl_vector_complex *x,
                 gsl_vector_complex *y)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if ((TransA == CblasNoTrans && N != x->size) ||
      (TransA == CblasTrans && M != x->size))
    {
      GSL_ERROR("invalid length of x vector", GSL_EBADLEN);
    }
  else if ((TransA == CblasNoTrans && M != y->size) ||
           (TransA == CblasTrans && N != y->size))
    {
      GSL_ERROR("invalid length of y vector", GSL_EBADLEN);
    }
  else
    {
      size_t j;
      size_t incX, incY;
      size_t lenX, lenY;
      gsl_complex *X, *Y;
      double *Ad;
      int *Ap, *Ai, *Aj;
      int p;

      if (TransA == CblasNoTrans)
        {
          lenX = N;
          lenY = M;
        }
      else
        {
          lenX = M;
          lenY = N;
        }

      /* form y := op(A)*x */
      Ap = A->p;
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
            //   Y[Ai[p] * incY].dat[0] = Ad[2*p] * X[Aj[p] * incX];
            //   Y[Ai[p] * incY].dat[1] = Ad[2*p+1] * 
            Y[Ai[p] * incY] = gsl_complex_mul(gsl_complex_rect(Ad[2*p], Ad[2*p+1]), X[Aj[p]* incX]);
            }
        }
      else
        {
          GSL_ERROR("unsupported matrix type", GSL_EINVAL);
        }

      return GSL_SUCCESS;
    }
} /* gsl_spblas_dgemv() */

int main(){
    
    int states = 2;

    gsl_vector_complex* vector = gsl_vector_complex_alloc(states);
    gsl_spmatrix_complex* matrix = gsl_spmatrix_complex_alloc(states, states);
    gsl_vector_complex* vector2 = gsl_vector_complex_alloc(states);


    gsl_spmatrix_complex_set(matrix, 0,0, gsl_complex_rect(1,3));
    gsl_spmatrix_complex_set(matrix, 1,1, gsl_complex_rect(1,3));

    gsl_vector_complex_set(vector, 0, gsl_complex_rect(2,4));
    gsl_vector_complex_set(vector, 1, gsl_complex_rect(1,5));


    myMulFunc(CblasNoTrans, matrix, vector, vector2);

    for(int i = 0; i < matrix->nz; i++){
        printf("%lg + %lgi\n", matrix->data[2*i], matrix->data[2*i+1]);
    }

    gsl_complex z1 = gsl_complex_rect(1,2);
    z1 = gsl_complex_add(z1, gsl_complex_rect(1,2));
    printf("%lg + %lgi\n", GSL_REAL(z1), GSL_IMAG(z1));
    printf("%zu\n", vector->stride);

    for(int j = 0; j < vector2->size; j++){
        printf("%lg + %lgi\n", GSL_REAL(gsl_vector_complex_get(vector2, j)), GSL_IMAG(gsl_vector_complex_get(vector2, j)));
    }


    return 0;
    //complex numbers in sp matricies are not stored as complex data types! but they are split up into an array of doubles!


}

// Y[Ai[p] * incY] += alpha * Ad[p] * X[j * incX];
// Need this in terms of imag # stored in sparse format
// One complex # = gsl_complex_rect(y[0], y[1])

// step one, remove alpha as we wont be using it, it will always be one.
// Will only be using COO storage as it's good enough



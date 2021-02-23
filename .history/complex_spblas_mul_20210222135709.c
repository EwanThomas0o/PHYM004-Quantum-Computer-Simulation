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

int myMulFunc(const CBLAS_TRANSPOSE_t TransA, const double alpha,
                 const gsl_spmatrix_complex *A, const gsl_vector_complex *x,
                 const double beta, gsl_vector_complex *y)
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
      double *X, *Y;
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

      /* form y := beta*y */

      Y = y->data;
      incY = y->stride;

      if (beta == 0.0)
        {
          size_t jy = 0;
          for (j = 0; j < lenY; ++j)
            {
              Y[jy] = 0.0;
              jy += incY;
            }
        }
      else if (beta != 1.0)
        {
          size_t jy = 0;
          for (j = 0; j < lenY; ++j)
            {
              Y[jy] *= beta;
              jy += incY;
            }
        }

      if (alpha == 0.0)
        return GSL_SUCCESS;

      /* form y := alpha*op(A)*x + y */
      Ap = A->p;
      Ad = A->data;
      X = x->data;
      incX = x->stride;

      if ((GSL_SPMATRIX_ISCCS(A) && (TransA == CblasNoTrans)) ||
          (GSL_SPMATRIX_ISCRS(A) && (TransA == CblasTrans)))
        {
          Ai = A->i;

          for (j = 0; j < lenX; ++j)
            {
              for (p = Ap[j]; p < Ap[j + 1]; ++p)
                {
                  Y[Ai[p] * incY] += alpha * Ad[p] * X[j * incX];
                }
            }
        }
      else if ((GSL_SPMATRIX_ISCCS(A) && (TransA == CblasTrans)) ||
               (GSL_SPMATRIX_ISCRS(A) && (TransA == CblasNoTrans)))
        {
          Ai = A->i;

          for (j = 0; j < lenY; ++j)
            {
              for (p = Ap[j]; p < Ap[j + 1]; ++p)
                {
                  Y[j * incY] += alpha * Ad[p] * X[Ai[p] * incX];
                }
            }
        }
      else if (GSL_SPMATRIX_ISTRIPLET(A))
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
              Y[Ai[p] * incY] += alpha * Ad[p] * X[Aj[p] * incX];
            }
        }
      else
        {
          GSL_ERROR("unsupported matrix type", GSL_EINVAL);
        }

      return GSL_SUCCESS;
    }
} /* gsl_spblas_dgemv() */

// int spComplexMul(gsl_vector_complex * a, gsl_spmatrix_complex* b)
// {
//   const size_t N = a->size;

//   if (b->size1 != N)
//     {
//       GSL_ERROR ("lenght of vector must have same length as matrix columns", GSL_EBADLEN);
//     }
//   else 
//     {
//       const size_t stride_a = a->stride;
//       const size_t stride_b = b->stride1;
//       const size_t stride_c = b->stride2

//       size_t i;

//       for (i = 0; i < N; i++)
//         {
//           ATOMIC ar = a->data[2 * i * stride_a];
//           ATOMIC ai = a->data[2 * i * stride_a + 1];
          
//           ATOMIC br = b->data[2 * i * stride_b];
//           ATOMIC bi = b->data[2 * i * stride_b + 1];

//           a->data[2 * i * stride_a] = ar * br - ai * bi;
//           a->data[2 * i * stride_a + 1] = ar * bi + ai * br;
//         }
      
//       return GSL_SUCCESS;
//     }
// }

int main(){
    
    int states = 8;

    gsl_vector_complex* vector = gsl_vector_complex_alloc(states);
    gsl_spmatrix_complex* matrix = gsl_spmatrix_complex_alloc(states, states);
    gsl_vector_complex* vector2 = gsl_vector_complex_alloc(states);

    myMulFunc(CblasNoTrans, 1.0, matrix, vector, 0.0, vector2);

    return 0;


}

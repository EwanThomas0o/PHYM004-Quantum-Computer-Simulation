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
#include <gsl/gsl_math.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_spmatrix.h>

int main(){

    // Creating the complex sparse matrices
    gsl_spmatrix_complex* matrix_one = NULL;
    matrix_one = gsl_spmatrix_complex_alloc(5,5);
    gsl_sp_matrix_complex_set(matrix_one, 2, 3, gsl_complex_rect(1.0, 2.0) ); 

    gsl_spmatrix* matrix_two = NULL;
    matrix_two = gsl_spmatrix_complex_alloc(5,5);
    gsl_sp_matrix_complex_set (matrix_two, 4, 2, gsl_complex_rect(1.0, 4.0) );

    // Creating 4 matrices, real and imaginary sparse matrices for matrix_one and matrix_two. So that we can use CBLAS to multiply.
    gsl_spmatrix* real_matrix_one = NULL;
    gsl_spmatrix* imag_matrix_one = NULL;
    gsl_spmatrix* real_matrix_two = NULL;
    gsl_spmatrix* imag_matrix_two = NULL;

    real_matrix_one = gsl_spmatrix_complex_alloc(matrix_one->size1, matrix_one->size2);
    imag_matrix_one = gsl_spmatrix_complex_alloc(matrix_one->size1, matrix_one->size2);
    real_matrix_two = gsl_spmatrix_complex_alloc(matrix_one->size1, matrix_one->size2);
    imag_matrix_one = gsl_spmatrix_complex_alloc(matrix_one->size1, matrix_one->size2);

    //Looping through non-zero elements of original complex spmatrices and assigning values to the corresponding real and imaginary spmatrices
    for(int j = 0 ; j < matrix_one->nz; j++){
        gsl_spmatrix_set(real_matrix_one, matrix_one->i[j], matrix_one->p[j], GSL_REAL(gsl_spmatrix_complex_get(matrix_one, matrix_one->i[j], matrix_one->p[j])) );
        gsl_spmatrix_set(imag_matrix_one, matrix_one->i[j], matrix_one->p[j], GSL_IMAG(gsl_spmatrix_complex_get(matrix_one, matrix_one->i[j], matrix_one->p[j])) );
    }

    for(int k = 0; k < matrix_two->nz; k++){
        gsl_spmatrix_set(real_matrix_two, matrix_two->i[k], matrix_two->p[k], GSL_REAL(gsl_spmatrix_complex_get(matrix_two, matrix_two->i[k], matrix_two->p[k])) );
        gsl_spmatrix_set(real_matrix_two, matrix_two->i[k], matrix_two->p[k], GSL_IMAG(gsl_spmatrix_complex_get(matrix_two, matrix_two->i[k], matrix_two->p[k])) );
    }

    // Now we should be able to use CBLAS in order to multiply elements
    // Need 4 more matrices to store answers as CBLAS  doesn't allow in place multiplication
    
    // Real(A) * Real(B)  -|
    //                     |-- Both form the real part of the answer
    // Imag(A) * Imag(B)  -|

    // Real(A) * Imag(B)  -|
    //                     |-- These form the imaginary part of the answer
    // Imag(A) * Real(B)  -|

    gsl_spmatrix* real_real = gsl_spmatrix_alloc(matrix_one->size1, matrix_two->size2);
    gsl_spmatrix* imag_imag = gsl_spmatrix_alloc(matrix_one->size1, matrix_two->size2);
    gsl_spmatrix* real_imag = gsl_spmatrix_alloc(matrix_one->size1, matrix_two->size2);
    gsl_spmatrix* imag_real = gsl_spmatrix_alloc(matrix_one->size1, matrix_two->size2);

    gsl_spblas_dgemm(1, real_matrix_one, real_matrix_two, real_real);
    gsl_spblas_dgemm(1, real_matrix_one, imag_matrix_two, real_imag);
    gsl_spblas_dgemm(1, imag_matrix_one, real_matrix_two, imag_real);
    gsl_spblas_dgemm(1, imag_matrix_one, imag_matrix_two, imag_imag);




    return 0;
}
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

    gsl_spmatrix_complex* matrix_one = NULL;
    matrix_one = gsl_spmatrix_complex_alloc(5,5);
    gsl_sp_matrix_complex_set(matrix_one, 2, 3, gsl_complex_rect(1.0, 2.0) );

    gsl_spmatrix* matrix_two = NULL;
    matrix_two = gsl_spmatrix_complex_alloc(5,5);
    gsl_sp_matrix_complex_set (matrix_two, 4, 2, gsl_complex_rect(1.0, 4.0) );

    gsl_matrix* real_matrix_one = NULL;
    gsl_matrix* imag_matrix_one = NULL;
    gsl_matrix* real_matrix_two = NULL;
    gsl_matrix* imag_matrix_two = NULL;

    real_matrix_one = gsl_matrix_complex_alloc(matrix_one->size1, matrix_one->size2);
    imag_matrix_one = gsl_matrix_complex_alloc(matrix_one->size1, matrix_one->size2);
    real_matrix_two = gsl_matrix_complex_alloc(matrix_one->size1, matrix_one->size2);
    imag_matrix_one = gsl_matrix_complex_alloc(matrix_one->size1, matrix_one->size2);

    gsl_matrix_set_all(real_matrix_one, 0.0);
    gsl_matrix_set_all(imag_matrix_one, 0.0);
    gsl_matrix_set_all(real_matrix_two, 0.0);
    gsl_matrix_set_all(imag_matrix_one, 0.0);


    for(int j = 0 ; j < matrix_one->nz; j++){
        gsl_spmatrix_set(real_matrix_one, matrix_one->i[j], matrix_one->p[j], GSL_REAL(gsl_spmatrix_complex_get(matrix_one, matrix_one->i[j], matrix_one->p[j])) );
        gsl_spmatrix_set(imag_matrix_one, matrix_one->i[j], matrix_one->p[j], GSL_IMAG(gsl_spmatrix_complex_get(matrix_one, matrix_one->i[j], matrix_one->p[j])) );
    }

    for(int k = 0; k < matrix_two->nz; k++){
        gsl_spmatrix_set(real_matrix_two, matrix_two->i[j], matrix_two->p[j], GSL_REAL(gsl_spmatrix_complex_get(matrix_two, matrix_two->i[j], matrix_two->p[j])) );
        gsl_spmatrix_set(real_matrix_two, matrix_two->i[j], matrix_two->p[j], GSL_IMAG(gsl_spmatrix_complex_get(matrix_two, matrix_two->i[j], matrix_two->p[j])) );

    }

    return 0;
}
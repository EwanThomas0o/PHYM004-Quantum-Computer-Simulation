#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_complex.h>

int main(){
    const size_t size = 8;
    gsl_spmatrix_complex *A = gsl_spmatrix_complex_alloc(size,size);

    gsl_spmatrix_complex_set(A, 2, 1, GSL_COMPLEX_ONE);
}


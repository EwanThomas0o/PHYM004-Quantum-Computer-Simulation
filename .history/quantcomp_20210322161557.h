#ifndef __QUANTCOMP__
#define __QUANTCOMP__
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>

gsl_vector_complex* initWavefunctionSpinDown(int states);
gsl_vector_complex* hadamardGate(gsl_vector_complex* wavefunction, int qubit);
void print_wf(gsl_vector_complex* register);


#endif 
#ifndef __QUANTCOMP__
#define __QUANTCOMP__
//May include extern C to protect against name mangling of c++
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>

void initWavefunctionSpinDown(gsl_vector_complex* wavefunction);
gsl_vector_complex* createRegister(int qubits);
void hadamardGate(gsl_vector_complex* wavefunction, int qubit);
void print_wf(gsl_vector_complex* wavefunction);
char* measureRegisterGate(gsl_vector_complex* wavefunction);
void swapsies(gsl_vector_complex* new_state, gsl_vector_complex* wavefunction);
void groversBlock(gsl_vector_complex* wavefunction, int answer);
void CphaseGate(gsl_vector_complex* wavefunction, int control, int target, double phase);
void phaseShiftGate(gsl_vector_complex *wavefunction, int qubit, double phase);
void oracleGate(gsl_vector_complex* wavefunction, int answer);

#endif 
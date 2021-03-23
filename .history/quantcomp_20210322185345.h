#ifndef __QUANTCOMP__
#define __QUANTCOMP__


gsl_vector_complex* initWavefunctionSpinDown(int states);
gsl_vector_complex* hadamardGate(gsl_vector_complex* wavefunction, int qubit);
void print_wf(gsl_vector_complex* wavefunction);
char* measureRegisterGate(gsl_vector_complex* wavefunction);


#endif 
#include "quantcomp.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>

int main(int argc, char **argv){
    
    gsl_vector_complex* quantum_register = initWavefunctionSpinDown(4); // 2 Qubits
    hadamardGate(quantum_register, 1);
    print_wf(quantum_register);   

    return 0;
}
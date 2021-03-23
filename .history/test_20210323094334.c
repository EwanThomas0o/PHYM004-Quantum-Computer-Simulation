#include "quantcomp.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>


int main(){
    
    gsl_vector_complex* quantum_register = initWavefunctionSpinDown(4); // 2 Qubits
    
    print_wf(quantum_register);
    
    hadamardGate(quantum_register, 1);

    print_wf(quantum_register);

    measureRegister(quantum_register);


    return 0;
}
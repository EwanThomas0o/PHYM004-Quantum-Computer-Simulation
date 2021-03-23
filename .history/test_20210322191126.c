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
    
    gsl_vector_complex* quantum_register1 = initWavefunctionSpinDown(4); // 2 Qubits
    gsl_vector_complex* quantum_register2 = initWavefunctionSpinDown(4); // 2 Qubits
    gsl_vector_complex_set(quantum_register2, 1, GSL_COMPLEX_ONE);

    print_wf(quantum_register1);
    print_wf(quantum_register2);

    gsl_vector_complex** ptr1 = &quantum_register1;

    gsl_vector_complex** ptr2 = &quantum_register2;

    gsl_vector_complex* temp;

    temp = *ptr1;
    *ptr1 = *ptr2;
    *ptr2 = temp;

    print_wf(quantum_register1);
    print_wf(quantum_register2);


    return 0;
}
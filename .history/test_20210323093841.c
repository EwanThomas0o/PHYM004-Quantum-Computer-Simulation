#include "quantcomp.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>

void updateReg(gsl_vector_complex* new_state, gsl_vector_complex* wavefunction){
    //update the old wf with the new 
    gsl_vector_complex_swap(new_state, wavefunction);
    gsl_vector_complex_free(newstate);


}

int main(){
    
    gsl_vector_complex* quantum_register1 = initWavefunctionSpinDown(4); // 2 Qubits
    gsl_vector_complex* quantum_register2 = initWavefunctionSpinDown(4); // 2 Qubits
    gsl_vector_complex_set(quantum_register2, 1, gsl_complex_rect(1, 3.14159));

    print_wf(quantum_register1);
    print_wf(quantum_register2);

    swaper(quantum_register1, quantum_register2);


    print_wf(quantum_register2);


    return 0;
}
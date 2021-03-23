#include "quantcomp.h"

int main(int argc, char **argc){
    
    gsl_vector_complex* quantum_register = initWavefunctionSpinDown(4); // 2 Qubits
    hadamardGate(quantum_register, 1);
    print_wf(quantum_register);   
     
    return 0;
}
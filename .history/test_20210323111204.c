#include "quantcomp.h"

int main(){
    
     // 2 Qubits all spin down |00>
    gsl_vector_complex* quantum_register = createRegister(9);

    print_wf(quantum_register);
    
    hadamardGate(quantum_register, 1);

    print_wf(quantum_register);

    measureRegisterGate(quantum_register);

    return 0;
}
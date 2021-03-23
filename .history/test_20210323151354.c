#include "quantcomp.h"

int main(){
    
     // 2 Qubits all spin down |00>
    gsl_vector_complex* quantum_register = createRegister(3);

    print_wf(quantum_register);

    groversAlgorithm(quantum_register, 7);

    print_wf(quantum_register);

    printf("%s\n", intToBinary(10));

    measureRegisterGate(quantum_register);

    return 0;
}
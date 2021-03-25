#include "quantcomp.h"

int main(){
    
     // 3 Qubits all spin down |00>
    gsl_vector_complex* quantum_register = createRegister(3);

    // print_wf(quantum_register);

    groversAlgorithm(quantum_register, 8); //TODO: not allows stated out of bounds of register

    print_wf(quantum_register);


    measureRegisterGate(quantum_register);

    return 0;
}
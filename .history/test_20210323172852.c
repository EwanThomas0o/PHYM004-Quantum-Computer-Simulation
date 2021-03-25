#include "quantcomp.h"

int main(){
    
     // 2 Qubits all spin down |00>
    gsl_vector_complex* quantum_register = createRegister(3);

    // print_wf(quantum_register);

    groversAlgorithm(quantum_register, 7); //TODO: not allows stated out of bounds of register

    // print_wf(quantum_register);


    measureRegisterGate(quantum_register);

    return 0;
}
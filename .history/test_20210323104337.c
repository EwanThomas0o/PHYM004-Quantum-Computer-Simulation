#include "quantcomp.h"

int main(){
    
    gsl_vector_complex* quantum_register = createRegister(2); // 2 Qubits

    print_wf(quantum_register);
    
    hadamardGate(quantum_register, 1);



    measureRegisterGate(quantum_register);


    return 0;
}
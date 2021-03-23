#include "quantcomp.h"

void groversAlgorithm(gsl_vector_complex* quantum_register){
    int qubits = log2(quantum_register->size);

    for(int i = 1; i <= qubits; i++)
    {
        hadamardGate(quantum_register, i);
    }

    for(int j = 0; j < floor(M_PI_4*sqrt(pow(2,qubits))); j++){ // Needs to be called "floor(pi/4*sqrt(2^N))"" times for optimum output roughly 2 in our case
        groversBlock(quantum_register, 7); //Second argument is the basis state you want to be "right" in this case its |110>
    }
}

int main(){
    
     // 2 Qubits all spin down |00>
    gsl_vector_complex* quantum_register = createRegister(3);

    groversAlgorithm(quantum_register);

    print_wf(quantum_register);

    return 0;
}
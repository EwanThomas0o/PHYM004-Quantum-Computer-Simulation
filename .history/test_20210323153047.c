#include "quantcomp.h"

char* intToBinaryTest(int a, gsl_vector_complex* quantum_register){

    int qubits = log2(quantum_register->size);

    int bin[qubits];

    int remainder, temp = 1;
    
    int index = 0;

    while(a != 0){

        remainder = a % 2;
        a /= 2;
        bin += remainder*temp;
        temp *= 10;
        printf("%d\n", bin[i]);
    }
    char *bin_str = (char *) malloc(qubits*sizeof(char));


    return bin_str;

}

int main(){
    
     // 2 Qubits all spin down |00>
    gsl_vector_complex* quantum_register = createRegister(3);

    print_wf(quantum_register);

    groversAlgorithm(quantum_register, 7);

    print_wf(quantum_register);

    intToBinaryTest(3, quantum_register);

    measureRegisterGate(quantum_register);

    return 0;
}
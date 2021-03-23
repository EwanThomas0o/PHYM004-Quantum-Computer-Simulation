#include "quantcomp.h"

char* intToBinaryTest(int a, gsl_vector_complex* quantum_register){

    int qubits = log2(quantum_register->size);

    int* bin = calloc(qubits, sizeof(int));

    int remainder, temp = 1;
    
    int index = 0;

    while(a != 0){

        remainder = a % 2;
        a /= 2;
        bin[qubits-index+1] = remainder*temp;
        // temp *= 10;
        index++;
    }
    char *bin_str = (char *) malloc(qubits*sizeof(char));
    for(int l = 0; l < qubits; l++){
        printf("%d", bin[l]);
    }printf("\n");
    free(bin);
    return bin_str;

}

int main(){
    
     // 2 Qubits all spin down |00>
    gsl_vector_complex* quantum_register = createRegister(3);

    // print_wf(quantum_register);

    // groversAlgorithm(quantum_register, 7);

    // print_wf(quantum_register);

    intToBinaryTest(3, quantum_register);

    // measureRegisterGate(quantum_register);

    return 0;
}
#include "quantcomp.h"

char* intToBinaryTest(int a, gsl_vector_complex* quantum_register){

    int qubits = log2(wavefunction->size);

    int bin = 0;

    int remainder, temp = 1;

    while(a != 0){

        remainder = a % 2;
        a /= 2;
        bin += remainder*temp;
        temp *= 10;
        print("%d", bin);
    }
    char *bin_str = (char *) xmalloc(qubits*sizeof(char));


    return bin_str;

}

int main(){
    
     // 2 Qubits all spin down |00>
    gsl_vector_complex* quantum_register = createRegister(3);

    print_wf(quantum_register);

    groversAlgorithm(quantum_register, 7);

    print_wf(quantum_register);

    intToBinaryTest(3 quantum_register);

    measureRegisterGate(quantum_register);

    return 0;
}
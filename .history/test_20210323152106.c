#include "quantcomp.h"

char* intToBinaryTest(int a){

    int bin = 0;

    int remainder, temp = 1;

    while(a != 0){

        remainder = a % 2;
        a /= 2;
        bin += remainder*temp;
        temp *= 10;
    }
    char *bin_str = (char *) xmalloc(N*sizeof(char));

    sprintf(bin_str, "%03d", bin);

    return bin_str;

}

int main(){
    
     // 2 Qubits all spin down |00>
    gsl_vector_complex* quantum_register = createRegister(3);

    print_wf(quantum_register);

    groversAlgorithm(quantum_register, 7);

    print_wf(quantum_register);

    intToBinaryTest();

    measureRegisterGate(quantum_register);

    return 0;
}
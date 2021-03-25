#include "quantcomp.h"
#define GNUPLOT_EXE    "gnuplot"
#define GNUPLOT_SCRIPT "myprog_gnuplot.script"
#define GNUPLOT_DATA   "myprog_gnuplot.dat"

int main(){
    
    for(int i = 0; i < 3; i++){
        // 3 Qubits all spin down |00>
        gsl_vector_complex* quantum_register = createRegister(3);

        groversAlgorithm(quantum_register, 8);

        measureRegisterGate(quantum_register);

        gsl_vector_complex_free(quantum_register);
    }

    snprintf(command, sizeof(command), "%s %s", GNUPLOT_EXE, GNUPLOT_SCRIPT );
	system( command );

    return 0;
}
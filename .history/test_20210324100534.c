#include "quantcomp.h"
#define GNUPLOT_EXE    "gnuplot"
#define GNUPLOT_SCRIPT "myprog_gnuplot.script"
#define GNUPLOT_DATA   "myprog_gnuplot.dat"
static struct timeval stop, start;

int main(){

    long avgs[9];
    long unsigned times[3];
    
    for(int i = 0; i < 3; i++){

        // 3 Qubits all spin down |00>
        gsl_vector_complex* quantum_register = createRegister(3);

        gettimeofday(&start, NULL);
        groversAlgorithm(quantum_register, 8);
        gettimeofday(&stop, NULL);

        measureRegisterGate(quantum_register);

        gsl_vector_complex_free(quantum_register);
    }

    snprintf(command, sizeof(command), "%s %s", GNUPLOT_EXE, GNUPLOT_SCRIPT );
	system( command );

    return 0;
}
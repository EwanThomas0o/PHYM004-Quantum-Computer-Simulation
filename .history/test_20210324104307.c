#include "quantcomp.h"
#define GNUPLOT_EXE    "gnuplot"
#define GNUPLOT_SCRIPT "myprog_gnuplot.script"
#define GNUPLOT_DATA   "myprog_gnuplot.dat"
#define USEC_IN_SEC 1000000
static struct timeval stop, start;

int main(){

    long avgs[9];
    long unsigned times[3];
    
    for(int j = 1; j <= 9; j++){
        
        for(int i = 0; i < 3; i++){

            // 3 Qubits all spin down |00>
            gsl_vector_complex* quantum_register = createRegister(i);

            gettimeofday(&start, NULL);
            groversAlgorithm(quantum_register, 1);
            gettimeofday(&stop, NULL);

            times[i] = (stop.tv_sec - start.tv_sec) * USEC_IN_SEC + stop.tv_usec - start.tv_usec

            measureRegisterGate(quantum_register);

            gsl_vector_complex_free(quantum_register);
        }
    }

    snprintf(command, sizeof(command), "%s %s", GNUPLOT_EXE, GNUPLOT_SCRIPT );
	system( command );

    return 0;
}
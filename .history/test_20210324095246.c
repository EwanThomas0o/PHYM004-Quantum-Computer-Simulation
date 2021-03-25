#include "quantcomp.h"

int main(){
    
     // 3 Qubits all spin down |00>
    gsl_vector_complex* quantum_register = createRegister(3);

    // print_wf(quantum_register);

    groversAlgorithm(quantum_register, 8); //TODO: not allows stated out of bounds of register

    print_wf(quantum_register);


    measureRegisterGate(quantum_register);

    FILE *pipe_gp = popen("gnuplot -p", "w");
    fputs("set terminal png \n",pipe_gp);
    fputs("set output 'grovers_scaling.png' \n",pipe_gp);
    fputs("set xlabel 'f' \n",pipe_gp);
    fputs("set xrange [0:100] \n",pipe_gp);
    fputs("set yrange [0:100] \n",pipe_gp);
    fputs("plot 'data.temp' u 1:2 w circles lc rgb 'pink' notitle \n",pipe_gp);
    pclose(pipe_gp);

    return 0;
}
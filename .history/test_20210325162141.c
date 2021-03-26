/****************** PHYM004 PR2: Quantum Computer Simulation -> Test environment ******************/
// USAGE
// -----
// First use "make" to compile and link the header file
// Then use "make do" to execute the binary
//
// Input 1 to run an analysis on Grover's algorithm
// Input 2 to run shors Algorithm with 7 qubits
// Input 3 to run code in the sandbox environment

// Gate available in sandbox mode
// ------------------------------
// Hadamard Gate
// Phase Gate
// Oracle Gate
// Grovers Block
// Controlled not gate
// Controlled phase gate
// J gate

// Utility functions available
// ---------------------------
// print_wf
// measureRegisterGate

// Getting Started
// ---------------
// declare your register as type "gsl_vector_complex*" and initialise with createRegister



#include "quantcomp.h"
#define GNUPLOT_EXE    "gnuplot"
#define GNUPLOT_SCRIPT "grovers.script"
#define GNUPLOT_DATA   "grovers_time.dat"
#define USEC_IN_SEC 1000000
#include <sys/time.h>

static struct timeval stop, start;

int main(){

    int mode;
    printf("Choose which mode you'd like to enter\n");
    printf("Mode 1 = Analysis of Grovers Algorithm\n");
    printf("Mode 2 = How many qubits?\n");
    printf("Mode 3 = Shor's Algorithm using 7 Qubits for the number 15\n");
    printf("Mode 4 = Sandbox mode which compiles code in the sandbox area\n");
    scanf("%d", &mode);
    
    if(mode == 1)
    {
        char command[PATH_MAX];
        long avgs[9];
        FILE* fp = fopen("grovers_time.dat", "w+");
        if(fp == NULL){
            printf("Error: failed to open file");
            exit(0);
        }
        
        for(int j = 1; j <= 8; j++){
            
            long unsigned times[3];
            
            for(int i = 0; i < 3; i++){
                
                // j Qubits all spin down |00..0>
                gsl_vector_complex* quantum_register = createRegister(j);
                
                //time shors alg
                gettimeofday(&start, NULL);
                groversAlgorithm(quantum_register, 1);
                gettimeofday(&stop, NULL);
                
                //put in times
                times[i] = (stop.tv_sec - start.tv_sec) * USEC_IN_SEC + stop.tv_usec - start.tv_usec;

                gsl_vector_complex_free(quantum_register);
            }

            //Find avg
            double avg = 0;
            double sum = 0;
            for(int num = 0; num < 3; num++){
                sum += times[num];
            }
            avg = sum/3;

            //Put avg in list of avgerages
            avgs[j] = avg;


            //For every time we input a new value into avgs, we need to write to a .dat file
            fprintf(fp, "%d\t%ld\n", j, avgs[j]);

        }
        fclose(fp);

        snprintf(command, sizeof(command), "%s %s", GNUPLOT_EXE, GNUPLOT_SCRIPT );
        system( command );
    }

    if(mode == 2)
    {
        char command[PATH_MAX];
        long avgs[9];
        FILE* fp = fopen("sim_time.dat", "w+");
        if(fp == NULL){
            printf("Error: failed to open file");
            exit(0);
        }
        
        for(int j = 1; j <= 8; j++){
            
            long unsigned times[3];
            
            for(int i = 0; i < 3; i++){
                
                // j Qubits all spin down |00..0>
                gettimeofday(&start, NULL);
                gsl_vector_complex* quantum_register = createRegister(j);
                //time shors creation of reg and simple operation
                hadamardGate(quantum_register, 1);
                gettimeofday(&stop, NULL);
                
                //put in times
                times[i] = (stop.tv_sec - start.tv_sec) * USEC_IN_SEC + stop.tv_usec - start.tv_usec;

                gsl_vector_complex_free(quantum_register);
            }

            //Find avg
            double avg = 0;
            double sum = 0;
            for(int num = 0; num < 3; num++){
                sum += times[num];
            }
            avg = sum/3;

            //Put avg in list of avgerages
            avgs[j] = avg;


            //For every time we input a new value into avgs, we need to write to a .dat file
            fprintf(fp, "%d\t%ld\n", j, avgs[j]);

        }
        fclose(fp);

        snprintf(command, sizeof(command), "%s %s", GNUPLOT_EXE, GNUPLOT_SCRIPT );
        system( command );

    }

    if(mode == 3)
    {
        gsl_vector_complex* reg = createRegister(7);
        shors(reg,15);
    }

    if(mode == 4)
    {   
        // Example of grovers algorithm and printing probability amplitudes.
        gsl_vector_complex* reg = createRegister(3);
        groversAlgorithm(reg, 2);
        print_wf(reg);
        measureRegisterGate(reg);
        
    }

    return 0;
}
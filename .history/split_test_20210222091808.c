#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <sys/time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>

#define BASIS 2
#define N 3
// Program to test multiplying real matrix by complex vector and then saving result as complex vector



gsl_vector_complex* initWavefunctionEqualProb(int states){ //Initialising wf to "Equal Probability" of all states
    
    gsl_vector_complex* wavefunction = NULL;
    wavefunction = gsl_vector_complex_alloc(states);
    gsl_vector_complex_set_all(wavefunction, gsl_complex_rect(1/sqrt(8),0));

    return wavefunction;
}

double findElementHad(char* inta, char* intb, int qubit){
    // Hadamard gate for single qubit used to calculate tensor product
    gsl_matrix *hadamard_single = gsl_matrix_alloc(BASIS, BASIS);
    gsl_matrix_set_all(hadamard_single, gsl_rect(1/sqrt(BASIS),0));
    gsl_matrix_set(hadamard_single,1,1, gsl_rect(-1/sqrt(BASIS),0));

    double value = 1.0;
    for(int i = 0; i < N; i++){
        if(inta[i] != intb[i] && i != qubit - 1){
            return 0.0; // Invokes Kronecker delta
        }
        value =  gsl_matrix_get(hadamard_single, inta[qubit-1] - '0', intb[qubit -1] - '0');

    }
    gsl_matrix_free(hadamard_single);
    return value;
}

gsl_vector_complex* hadamardGate(gsl_vector_complex* wavefunction, int qubit){
    if(qubit > N){
        printf("Please operate the gate on a valid qubit\n");
        exit(0);
    }
    // Will beome the NxN matrix for operation on whole register
    gsl_matrix *hadamard = gsl_matrix_alloc(wavefunction->size, wavefunction->size);
    gsl_matrix_set_all(hadamard, GSL_POSZERO);

    for(int i = 0; i < wavefunction->size; i++){
        for(int j = 0; j < wavefunction->size; j++){
            double val = findElementHad(intToBinary(i), intToBinary(j), qubit);
            gsl_matrix_set(hadamard, i , j, val); // Hadamard matrix created
        }
    }
    // Can use gsl_vector_views to split wf into real and imaginary
    gsl_vector_view real = gsl_vector_complex_real(wavefunction);
    gsl_vector_view imag = gsl_vector_complex_imag(wavefunction);

    gsl_vector_complex* h_psi = gsl_vector_complex_alloc(wavefunction->size);
    gsl_vector_view new_real = gsl_vector_complex_real(h_psi);
    gsl_vector_view new_imag = gsl_vector_complex_imag(h_psi);

    gsl_blas_dgemv(CblasNoTrans, 1.0, hadamard, &real.vector, 0.0, &new_real.vector);
    gsl_blas_dgemv(CblasNoTrans, 1.0, hadamard, &imag.vector, 0.0, &new_imag.vector);
    return h_psi;

}


int main(){




    return 0;
}
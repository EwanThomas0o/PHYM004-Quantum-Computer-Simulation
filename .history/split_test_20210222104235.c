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
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_spmatrix.h>


#define BASIS 2
#define N 3

struct timeval tv;
// Program to test multiplying real matrix by complex vector and then saving result as complex vector



gsl_vector_complex* initWavefunctionEqualProb(int states){ //Initialising wf to "Equal Probability" of all states
    
    gsl_vector_complex* wavefunction = NULL;
    wavefunction = gsl_vector_complex_alloc(states);
    gsl_vector_complex_set_all(wavefunction, gsl_complex_rect(1/sqrt(8),0));

    return wavefunction;
}
gsl_vector_complex* initWavefunctionSpinDown(int states){ 
    gsl_vector_complex* wavefunction = NULL;
    wavefunction = gsl_vector_complex_alloc(states);
    // Probability amplitude stored in "wavefuction" vector
    gsl_vector_complex_set(wavefunction, 0, GSL_COMPLEX_ONE);

    return wavefunction;
}

char* intToBinary(int a){
    int bin = 0;
    int remainder, temp = 1;

    while(a != 0){
        remainder = a % 2;
        a /= 2;
        bin += remainder*temp;
        temp *= 10;
    }
    char *bin_str = (char *)malloc(N*sizeof(char));
    sprintf(bin_str, "%03d", bin);    
    return bin_str;
}

gsl_vector_complex* swapsies(gsl_vector_complex* new_state, gsl_vector_complex* wavefunction){
    //update the old wf with the new 
    gsl_vector_complex** ptr1 = &new_state;

    gsl_vector_complex** ptr2 = &wavefunction;

    gsl_vector_complex* temp;

    temp = *ptr1;
    *ptr1 = *ptr2;
    *ptr2 = temp;

    gsl_vector_complex_free(new_state);

    return wavefunction;
}

double findElementHad(char* inta, char* intb, int qubit){
    // Hadamard gate for single qubit used to calculate tensor product
    gsl_matrix *hadamard_single = gsl_matrix_alloc(BASIS, BASIS);
    gsl_matrix_set_all(hadamard_single, 1/sqrt(BASIS));
    gsl_matrix_set(hadamard_single,1,1, -1/sqrt(BASIS));

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
    gsl_spmatrix *hadamard = gsl_spmatrix_alloc(wavefunction->size, wavefunction->size);

    for(int i = 0; i < wavefunction->size; i++){
        for(int j = 0; j < wavefunction->size; j++){
            double val = findElementHad(intToBinary(i), intToBinary(j), qubit);
            gsl_spmatrix_set(hadamard, i , j, val); // Hadamard matrix created
        }
    }
    // Can use gsl_vector_views to split wf into real and imaginary
    gsl_vector_view real = gsl_vector_complex_real(wavefunction);
    gsl_vector_view imag = gsl_vector_complex_imag(wavefunction);

    gsl_vector_complex* h_psi = gsl_vector_complex_alloc(wavefunction->size);
    gsl_vector_view new_real = gsl_vector_complex_real(h_psi); // On the stack so not mem intensive
    gsl_vector_view new_imag = gsl_vector_complex_imag(h_psi);

    gsl_spblas_dgemv(CblasNoTrans, 1.0, hadamard, &real.vector, 0.0, &new_real.vector);
    gsl_spblas_dgemv(CblasNoTrans, 1.0, hadamard, &imag.vector, 0.0, &new_imag.vector);
    
    gsl_vector_complex** ptr1 = &h_psi;
    gsl_vector_complex**ptr2 = &wavefunction;
    gsl_vector_complex* temp;

    temp = *ptr1;
    *ptr1 = *ptr2;
    *ptr2 = temp;
    // gsl_vector_complex_memcpy(wavefunction, h_psi);
    gsl_vector_complex_free(h_psi);
    
    return swapsies(h_psi, wavefunction);

}


void measureRegisterGate(gsl_vector_complex* wavefunction){
    int states = (int) wavefunction->size;

    double probabilities[states];
    for(int i = 0; i < states; i++){ //creating a list of the probabilities (non-normalised)
        probabilities[i] = GSL_REAL(gsl_vector_complex_get(wavefunction,i))*GSL_REAL(gsl_vector_complex_get(wavefunction,i)) + GSL_IMAG(gsl_vector_complex_get(wavefunction,i))*GSL_IMAG(gsl_vector_complex_get(wavefunction,i));
        //printf("%lg\n", probabilities[i]);
    }
    gsl_ran_discrete_t* lookup = gsl_ran_discrete_preproc(states, probabilities); // Preproc normalises the probabilities
    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937); // Mersene Twister Algorithm is used for high quality random numbers
    gettimeofday(&tv,NULL); // Get time of day in usec so that the seed changes giving a different stream of #'s each time
    unsigned long int seed = tv.tv_usec;    
    gsl_rng_set(r, seed); 
    size_t t = gsl_ran_discrete(r, lookup); // Choosing from the discrete probability distribution defined by the wavefunction 
    printf("Wavefunction collapsed into the state:\n|%s>\n", intToBinary(t));
    // Wavefunction collapsed so will only find system in the state from now on
    gsl_vector_complex_set_all(wavefunction, GSL_COMPLEX_ZERO);
    gsl_vector_complex_set(wavefunction, t, GSL_COMPLEX_ONE); // Set measured state to probability one so that if we measure again we get the same outcome
    // Free memory to avoid bloats
    gsl_ran_discrete_free(lookup);
    gsl_rng_free(r);
    return;
}

void printWavefunction(gsl_vector_complex* wavefunction){
    
    for(int i = 0; i < wavefunction->size; i++){
        printf("%lg\n", gsl_vector_complex_get(wavefunction, i).dat[0]);
    }
}

int main(){
    int states = (int)pow(BASIS, N);

    gsl_vector_complex* wavefunction = initWavefunctionSpinDown(states);

    wavefunction = hadamardGate(wavefunction, 1);

    printWavefunction(wavefunction);


    return 0;
}
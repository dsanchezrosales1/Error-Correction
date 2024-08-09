#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "LDPC_functions.c"

#define m 21600
#define n 43200
#define k 172800

#define Na 555
#define Ma 3893
#define Nf 556
#define Mf 3896


int main() {
    // Read cs_list from file
    cs_list_t* cs_list = read_cs_list("../TXT_files/cs_list.txt", k);

    // Read cs_index from file
    int* cs_index = read_int_list("../TXT_files/cs_index.txt", m);
    cs_index[m] = k+1;

    // Read LL_index from file
    int* LL_index = read_int_list("../TXT_files/LL_index.txt", n);
    LL_index[n] = k+1;

    // Read H Matrix/array from file
    int* H = read_int_list("../TXT_files/H_array.txt", k);

    /////////////////////////////////////////////////////////////////

    int n_runs = 10;
    float QBER[] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};

    int c[m];
    int b[m];
    int d[m];

    int LL_reg_i[k];
    int y1[n];

    //Generate Log Lookup Tables
    int a2f[Ma];
    int f2a[Mf];
    generate_log_tables(Ma, Na, Mf, Nf, a2f, f2a);

    // How many loops are there for different QBER?
    int num_QBER = sizeof(QBER) / sizeof(QBER[0]);

    //Generate time tracking strcutures
    struct timeval start_time, end_time;

    //Start first loop ieterating over all QBER values (1-9%)
    for (int q = 0; q < num_QBER; q++) {
        // Initialize the Log Likelihood Register with initial value, a_init
        int a_init = round(Na * log(1 - 2 * QBER[q]));
        int f_init = round(Nf * log((1 - QBER[q]) / QBER[q]));

        int LL_reg_i[k];

        for (int i = 0; i < k; i++) {
            LL_reg_i[i] = a_init;
        }

        printf("QBER = %f \n", QBER[q]);

        //Keep track of all the successful runs for a single QBER value
        int wins = 0;

        //Start second loop iterating over n_runs for each QBER value
        for (int i = 0; i < n_runs; i++) {
            //For each QBER, LL_reg must be initialized to the same value (a_init) for that QBER
            int LL_reg[k];
            memcpy(LL_reg, LL_reg_i, sizeof(LL_reg_i));

            //Generate file names to read
            char name_A[100], name_B[100];
            sprintf(name_A, "../TXT_files/Test_files/QBER_%d_A_%d.txt", (int)(QBER[q] * 100), i + 1);
            sprintf(name_B, "../TXT_files/Test_files/QBER_%d_B_%d.txt", (int)(QBER[q] * 100), i + 1);

            // Read Alice & Bob's data from files
            int* Alice = read_int_list(name_A, n);
            int* Bob = read_int_list(name_B, n);

            //Generate Alice's and Bob's syndromes by multiplyting the matrix  by their raw data
            matrix_multiply(Alice, H, k, m, c);
            matrix_multiply(Bob, H, k, m, b);

            //Generate d, an array describing the parity of c and b (if theyr equal, d is all 0's)
            for (int i = 0; i < m; i++) {
                d[i] = c[i]^b[i];
            }

            

        ////////// Run algorithm
            //Get current time before running algorithm
            gettimeofday(&start_time, NULL);
            //Initialize the number of iterations to zero
            int itr = 0;
            //Start loop to error corrcect Bob's data (maximum of 31 loops)
            for (int l = 0; l < 31; l++) {
                cs_msgs_2_bits(cs_list, cs_index, LL_reg, m, d, Ma);
                bit_msgs_2_cs(LL_reg, LL_index, n, a2f, f2a, Bob, y1, f_init, Mf);
                int success = converged(m, cs_list, cs_index, y1, c);
                itr++;
                if (success) {
                    //Get current time to see how long it took to converge one data set and print it
                    gettimeofday(&end_time, NULL);
                    double time_spent = (double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_usec - start_time.tv_usec) / 1000000.0;
                    printf("Iterations: %d     Time: %f seconds\n", itr, time_spent);
                    break;
                }
            }
            
            // Check if Alice equals y1( Bob's corrected data)
            int match = 1;
            for (int j = 0; j < n; j++) {
                if (Alice[j] != y1[j]) {
                    match = 0;
                    break;
                }
            }
            //Add 1 to wins to see how many data sets converged
            if (match) {
                wins += 1;
            }

            free(Alice);
            free(Bob);
        }

        printf("The fraction of successful runs is: %f \n\n", (float)wins / n_runs);
    }
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}

//LDPC belief propagation code functions for error correction
//Daniel Sanchez Rosales 
//Aliro Technologies
//The code presented here is an implementation of the paper shown here: https://drive.google.com/file/d/128vyeSn9-xBCvmXX7zK9qA6vryDDNCJO/view?usp=sharing
//I have documented the LDPC encoding and decoding protocols in detail at: https://docs.google.com/document/d/1VQ7SOWKeOjMaxbkn2DDh5PuEIR-P28c1O_bzaorVuBA/edit?usp=sharing
////This document will be known as the documentation in the comments of this code.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//The following functions are used to read in TXT data file which were generated in MATLAB.
//Because I used MATLAB to generate the files, you will often see an index offset of -1.
////This is because MATLAB indexing starts at 1 and in C it starts at 0.
//Addtionally, the mapping required for this code to work was done in MATLAB and it is not shown here. 
////Please refer to Pg. 5-6 in the documetation and to the paper above.

// Define the cs_list structure
typedef struct {
    int LL_ptr;
    int bit;
} cs_list_t;

// Function to read cs_list from a file
cs_list_t* read_cs_list(const char* filename, int count) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening cs_list file");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for cs_list
    cs_list_t* cs_list = (cs_list_t*)malloc(count * sizeof(cs_list_t));
    if (cs_list == NULL) {
        perror("Memory allocation error for cs_list");
        exit(EXIT_FAILURE);
    }

    // Skip the header line
    char buffer[32];
    fgets(buffer, sizeof(buffer), file);

    // Read the actual data
    int i = 0;
    while (fgets(buffer, sizeof(buffer), file) != NULL && i < count) {
        sscanf(buffer, "%d,%d", &cs_list[i].bit, &cs_list[i].LL_ptr);
        i++;
    }

    fclose(file);
    return cs_list;
}

// Function to read an integer list from a file
int* read_int_list(const char* filename, int count) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening int list file");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the integer list
    int* list = (int*)malloc(count * sizeof(int));
    if (list == NULL) {
        perror("Memory allocation error for int list");
        exit(EXIT_FAILURE);
    }


    // Read the actual data
    int i = 0;
    char buffer[16];
    while (fgets(buffer, sizeof(buffer), file) != NULL && i < count) {
        sscanf(buffer, "%d", &list[i]);
        i++;
    }

    fclose(file);
    return list;
}

//H is a sparse matrix, and as such it contains mostly zeros. 
//To reduce the amount fo memory required to store H, I save only the location of every one on each row.
////H conains 8 ones on every row, thus with 21600 rows we have 172800 values.
////Since each row contains 43200 columns, we can represent every value in H with 16 bits, for a total of 345600 bytes in memory.
////See documentation Pg.4-6
//To do matrix multiplication (H*Alice or H*bob) to generate CS or DS, I multiply only the indeces where H == 1.
////Additionally, the multiplication sums are done bitwie using XOR
void matrix_multiply(int *info_bits, int *H, int k, int m, int *result) {
    int acc_sum = 0;
    int idx = 0;

    for (int i = 0; i < k; i++) {
        acc_sum += info_bits[H[i]-1]; //The -1 is because H was generated in MATLAB, which indexes from 1 instead of 0.

        if (i % 8 == 7) {
            result[idx] = acc_sum % 2; //Use bitwise AND instead ***ASK MATT
            idx++; 
            acc_sum = 0;
        }
    }
}

//These are the lookup tables needed to avoid calculating the log values during execution
//See Pg. 4 in the documentation. These are done using predetermined values of Na, Ma, Nf, and Mf.
//These values were optimized in the paper documented above and are Na=555,Ma=3893, Nf=556, Mf=3896.
void generate_log_tables(int Ma, int Na, int Mf, int Nf, int* a2f, int* f2a) {
    // Generate a2f table
    for (int i = 1; i <= Ma; i++) {
        double a = exp(-((double)i) / Na);
        double f = (1 + a) / (1 - a);
        double z = log(f);
        int j = round(Nf * z + 0.5);
        if (j <= 0) {
            j = 1;
        } else if (j > Mf) {
            j = Mf;
        }
        a2f[i-1] = j;
    }

    // Generate f2a table
    for (int i = 1; i <= Mf; i++) {
        double f = exp((double)i / Nf);
        double a = (f - 1) / (f + 1);
        double z = -log(a);
        int j = round(Na * z + 0.5);
        if (j <= 0) {
            j = 1;
        } else if (j > Ma) {
            j = Ma;
        }
        f2a[i-1] = j;
    }
}

//This function implements message paasing from checksum nodes to bit nodes (step 3 in message passing algorith)
////See Pg. 3 in documentation.
//Agian, because I used MATLAB to generate cs_list, cs_index, LL_index, and H, you will see an index offset of -1.
void cs_msgs_2_bits(cs_list_t* cs_list, int* cs_index, int* LL_reg, int m, int* d, int Ma) {
    for (int i = 0; i < m; i++) {
        int sign = d[i];
        int big_alpha = 0;
        int j1 = cs_index[i]-1;
        int j2 = cs_index[i+1]-1;

        for (int j = j1; j < j2; j++) {
            int a1 = LL_reg[cs_list[j].LL_ptr -1];
            if (a1 < 0) {
                sign = 1 - sign;
                big_alpha = big_alpha - a1;
            } else {
                big_alpha = big_alpha + a1;
            }
        }

        for (int j = j1; j < j2; j++) {
            int a1 = LL_reg[cs_list[j].LL_ptr -1];
            int p_sign, p_alpha;
            if (a1 < 0) {
                p_sign = 1 - sign;
                p_alpha = big_alpha + a1;
            } else {
                p_sign = sign;
                p_alpha = big_alpha - a1;
            }
            if (p_alpha <= 0) {
                p_alpha = 1;
            }
            if (p_alpha > Ma) {
                p_alpha = Ma;
            }
            if (p_sign == 0) {
                LL_reg[cs_list[j].LL_ptr -1] = p_alpha;
            } else {
                LL_reg[cs_list[j].LL_ptr -1] = -p_alpha;
            }
        }
    }
}


//This function implements message paasing bit nodes to checksum nodes (step 2 in message passing algorith)
////This function also implements steps 4 and 5 in the message passage algorithm (see Pg. 2-4 of the documentation)
//Agian, because I used MATLAB to generate cs_list, cs_index, LL_index, and H, you will see an index offset of -1.
void bit_msgs_2_cs(int* LL_reg, int* LL_index, int n, int* a2f, int* f2a, int* y, int* y1, int f_init, int Mf) {
    for (int i = 0; i < n; i++) {
        int j1 = LL_index[i] - 1;
        int j2 = LL_index[i+1] - 1;
        int f_tot = f_init;

        for (int j = j1; j < j2; j++) {
            int u = LL_reg[j];
            if (u > 0) {
                u = a2f[u];
            } else {
                u = -a2f[-u];
            }
            LL_reg[j] = u;
            f_tot += u;
        }

        for (int j = j1; j < j2; j++) {
            int k = f_tot - LL_reg[j];
            int p_sign;
            if (k < 0) {
                p_sign = 1;
                k = -k;
            } else {
                p_sign = 0;
            }
            if (k < 1) {
                k = 1;
            }
            if (k > Mf) {
                k = Mf-1;
            }
            if (p_sign == 1) {
                LL_reg[j] = -f2a[k];
            } else {
                LL_reg[j] = f2a[k];
            }
        }

        if (f_tot < 0) {
            y1[i] = 1 - y[i];
        } else {
            y1[i] = y[i];
        }
    }
}

//This function checks whether Bob's new updated info bits satisfy the CS = DS condition of the algorithm
////See protocol figure in Pg. 3 of documentation
//Agian, because I used MATLAB to generate cs_list, cs_index, LL_index, and H, you will see an index offset of -1.
int converged(int m, cs_list_t* cs_list, int* cs_index, int* y1, int* c) {
    int success = 1;

    for (int i = 0; i < m; i++) {
        int sum = 0;
        int j1 = cs_index[i] - 1;
        int j2 = cs_index[i+1] - 1;

        for (int j = j1; j < j2; j++) {
            sum ^= y1[cs_list[j].bit-1];
        }

        if (sum != c[i]) {
            success = 0;
            break;
        }
    }

    return success;
}
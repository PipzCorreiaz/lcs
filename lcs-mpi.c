#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define ALPHA_SIZE 4

#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(index,p,n) (((p)*((index)+1)-1)/(n))


void print_new_matrix(int* matrix, int size1, int size2) {
    int i, j;
    for (i = 0; i < size1 + 1; i++) {
        for (j = 0; j < size2 + 1; j++) {
            int c = matrix[i * (size2 + 1) + j];
            printf("%d", c);
        }

        printf("\n");
    }

    printf("\n");
}

int letter_index(char alpha[4], char letter) {
    int i;
    for (i = 0; i < ALPHA_SIZE; i++) {
        if (alpha[i] == letter) {
            return i;
        }
    }

    return -1;
}

short cost(int x) {
    int i, n_iter = 20;
    double dcost = 0;
    for (i = 0; i < n_iter; i++) {
        dcost += pow(sin((double) x), 2) + pow(cos((double) x), 2);
    }

    return (short) (dcost / n_iter + 0.1);
}

int block_width(int strlen) {
    if (strlen < 50) {
        return 3;
    } if (strlen < 500) {
        return 100;
    } if (strlen < 10000) {
        return 4000;
    } if (strlen < 20000) {
        return 9000;
    } else {
        return 15000;
    }
}

int blocks_per_proc(int strlen, int block_width) {
    float x = strlen / (block_width * 1.0);
    int res = x;

    if (x - res > 0.5) {
        res++;
    }
    return res;
}

int last_block_size(int strlen, int block_width) {
    float x = strlen / (block_width * 1.0);
    int res = x;

    if (x - res > 0.5) {
        return strlen - (block_width * res);
    } else {
        return block_width + (strlen - (block_width * res));
    }
}


int main(int argc, char *argv[]) {
    int me, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    MPI_Status status;
    FILE *file;
    int seq1_size = 0;
    char *seq1, *seq2;
    int seq2_size = 0;

    file = fopen(argv[1], "r");
    if (file != NULL) {
        fscanf(file, "%d %d", &seq1_size, &seq2_size);

        seq1 = (char *) calloc(seq1_size + 1, sizeof(char));
        seq2 = (char *) calloc(seq2_size + 1, sizeof(char));

        fscanf(file, "%s %s", seq1, seq2);

    } else {
        printf("%s\n", "Unable to open file");
        exit(-1);
    }

    fclose(file);

    char C[ALPHA_SIZE] = { 'A', 'C', 'G', 'T' };
    int *P;

    P = (int *) calloc(ALPHA_SIZE * (seq2_size + 1), sizeof(int));


    int i, j, k, b;
    int gap = seq2_size + 1;
    int block_low = BLOCK_LOW(me, nprocs, ALPHA_SIZE);
    int block_high = BLOCK_HIGH(me, nprocs, ALPHA_SIZE);
    int block_size = BLOCK_SIZE(me, nprocs, ALPHA_SIZE);
    int b_width = block_width(seq2_size + 1);
    int b_per_proc = blocks_per_proc(seq2_size + 1, b_width);
    int last_b_size = last_block_size(seq2_size + 1, b_width);


    MPI_Barrier(MPI_COMM_WORLD);

    /* Compute P table */
    #pragma omp parallel for private(j)
    for (i = 0; i < ALPHA_SIZE; i++) {
        int index = i * gap;
        char c = C[i];
        for (j = 1; j < seq2_size + 1; j++) {
            char b = seq2[j - 1];
            if (b == c) {
                P[index + j] = j;
            } else {
                P[index + j] = P[index + j - 1];
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    block_low = BLOCK_LOW(me, nprocs, seq1_size + 1);
    block_high = BLOCK_HIGH(me, nprocs, seq1_size + 1);
    block_size = BLOCK_SIZE(me, nprocs, seq1_size + 1);
    int primer = BLOCK_OWNER(0, nprocs, seq1_size + 1);
    int laster = BLOCK_OWNER(seq1_size, nprocs, seq1_size + 1);
    int array_size;

    if (me == primer) {
        array_size = block_size;
    } else {
        array_size = block_size + 1;
    }

    int *S = (int*) calloc(array_size * (seq2_size + 1), sizeof(int));

    if (S == NULL) {
        printf("Not enough memory. Learning how to program might help...\n");
        exit(-1);
    }

    /*printf("me: %d, nprocs: %d, bl: %d, bh: %d, bs: %d\n", me, nprocs, block_low, block_high, block_size);*/

    for(k = 0; k < b_per_proc; k++) {
        int columns = b_width;
        int current_column = k * columns;

        if (k + 1 == b_per_proc) {
            columns = last_b_size;
        }
        if (me != primer && block_size > 0) {
            MPI_Recv(S + current_column, columns, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        }

        b = block_low == 0? 1 : block_low;
        for (i = 1; b < block_high + 1 && i < seq1_size + 1; i++, b++) {
            int index = i * gap;

            #pragma omp parallel for
            for (j = 0; j < columns; j++) {
                if (k != 0 || j != 0) {
                    char xi = seq1[b - 1];
                    int column = current_column + j;
                    char yj = seq2[column - 1];
                    int l_value = letter_index(C, xi);
                    int p_value = P[l_value * gap + column];

                    int t = (0 - p_value) < 0? 1 : 0;
                    int s;

                    if (p_value == 0) {
                        s = (0 - (S[index - gap + column] - t * S[index - gap])) < 0? 1 : 0;
                    } else {
                        s = (0 - (S[index - gap + column] - t * S[index - gap + p_value - 1])) < 0? 1 : 0;
                    }

                    if (xi == yj) cost(i);
                    S[index + column] = S[index - gap + column] + t * (s ^ 1);
                }
            }
        }

        if (block_size > 0 && block_high < seq1_size) {
            int next = BLOCK_OWNER(block_high + 1, nprocs, seq1_size + 1);
            if (me == primer) {
                MPI_Send(&S[block_high * gap + current_column], columns, MPI_INT, next, 0, MPI_COMM_WORLD);
            } else {
                MPI_Send(&S[block_size * gap + current_column], columns, MPI_INT, next, 0, MPI_COMM_WORLD);
            }
        }
    }

    /*print_new_matrix(S, array_size - 1, seq2_size);*/
    i = array_size - 1;
    j = seq2_size;
    char xi, yj;
    int last_cell;
    int len;
    int *lcs;


    if (me != laster) {
        MPI_Recv(&len, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        lcs = (int*) calloc(len + 1, sizeof(int));
        MPI_Recv(lcs, len + 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        j = lcs[0];
        last_cell = S[i * gap + j];
    } else {
        last_cell = S[i * gap + j];
        len = last_cell;
        lcs = (int*) calloc(len + 1, sizeof(int));
    }

    while(i > 0 && last_cell > 0) {
        int index = i * gap + j;
        if (me == primer) {
            xi = seq1[(block_low + i) - 1];
        } else {
            xi = seq1[(block_low + i) - 2];
        }
        yj = seq2[j - 1];
        if (xi == yj) {
            lcs[last_cell] = xi;
            last_cell--;
            i--;
            j--;
        } else if (S[index - gap] > S[index - 1]) {
            i--;
        } else {
            j--;
        }
    }

    if (last_cell > 0) {
        lcs[0] = j;
        int previous = BLOCK_OWNER(block_low - 1, nprocs, seq1_size + 1);
        MPI_Send(&len, 1, MPI_INT, previous, 0, MPI_COMM_WORLD);
        MPI_Send(lcs, len + 1, MPI_INT, previous, 0, MPI_COMM_WORLD);
    } else {
        printf("%d\n", len);
        for (i = 1; i < len + 1; i++) {
            printf("%c", lcs[i]);
        }

        printf("\n");
    }

/*

        while(last_cell > 0) {
            int index = i * gap + j;
            xi = seq1[i - 1];
            yj = seq2[j - 1];
            if (xi == yj) {
                lcs[last_cell - 1] = xi;
                last_cell--;
                i--;
                j--;
            } else if (S[index - gap] > S[index - 1]) {
                i--;
            } else {
                j--;
            }
        }

        printf("%d\n%s\n", len, lcs);

    }
*/

    free(S);
    free(P);
    free(seq1);
    free(seq2);

    MPI_Finalize();

    return 0;
}

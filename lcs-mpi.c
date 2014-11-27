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
    return 500;
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

    /*int i, j;
    int array[10];*/
    MPI_Status status;
    /*int buffer[10];*/
    /*int **matrix = (int**) calloc(nprocs, sizeof(int*));*/
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
    int buffer_max_size = ALPHA_SIZE / nprocs + 1;
    int *block = (int*) calloc(buffer_max_size * (seq2_size + 1), sizeof(int));
    int *buffer = (int*) calloc(buffer_max_size * (seq2_size + 1), sizeof(int));
    int b_width = block_width(seq2_size + 1);
    int b_per_proc = blocks_per_proc(seq2_size + 1, b_width);
    int last_b_size = last_block_size(seq2_size + 1, b_width);


    MPI_Barrier(MPI_COMM_WORLD);

    /* Compute P table */
    if (block_size > 0) {
        for (i = block_low, k = 0; i < block_high + 1; i++, k++) {
            /*int index = i * gap;*/
            int block_index = k * gap;
            char c = C[i];
            for (j = 1; j < seq2_size + 1; j++) {
                char b = seq2[j - 1];
                if (b == c) {
                    /*P[index + j] = j;*/
                    block[block_index + j] = j;
                } else {
                    /*P[index + j] = P[index + j - 1];*/
                    block[block_index + j] = block[block_index + j - 1];
                }
            }
        }

        for (j = 0; j < nprocs; j++) {
            MPI_Send(block, buffer_max_size * (seq2_size + 1), MPI_INT, j, 0, MPI_COMM_WORLD);
        }
    }

    int min = nprocs > ALPHA_SIZE? ALPHA_SIZE : nprocs;

    for (i = 0; i < min; i++) {
        MPI_Recv(buffer, buffer_max_size * (seq2_size + 1), MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        int source = status.MPI_SOURCE;
        int block_l = BLOCK_LOW(source, nprocs, ALPHA_SIZE);
        int size_b = BLOCK_SIZE(source, nprocs, ALPHA_SIZE);

        for (j = block_l, k = 0; k < size_b * (seq2_size + 1); k++) {
            int p_index = j * gap + k;
            P[p_index] = buffer[k];
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    block_low = BLOCK_LOW(me, nprocs, seq1_size + 1);
    block_high = BLOCK_HIGH(me, nprocs, seq1_size + 1);
    block_size = BLOCK_SIZE(me, nprocs, seq1_size + 1);
    int master = BLOCK_OWNER(0, nprocs, seq1_size + 1);
    int array_size = (seq1_size + 1) / nprocs + 1;

    int *S;
    if (me == master) {
        S = (int*) calloc((seq1_size + 1) * (seq2_size + 1), sizeof(int));
    } else {
        S = (int*) calloc((array_size + 1) * (seq2_size + 1), sizeof(int));
    }

    if (S == NULL) {
        printf("Not enough memory. Learning how to program might help...\n");
        exit(-1);
    }



    for(k = 0; k < b_per_proc; k++) {
        int columns = b_width;
        int current_column = k * columns;

        if (k + 1 == b_per_proc) {
            columns = last_b_size;
        }
        if (me != master && block_size > 0) {
            MPI_Recv(S + current_column, columns, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        }
        
        /*printf("me: %d, bl: %d, bs: %d\n", me, block_low, block_size);*/
        b = block_low == 0? 1 : block_low;
        for (i = 1; b < block_high + 1 && i < seq1_size + 1; i++, b++) {
            int index = i * gap;
            
            #pragma omp parallel for
            for (j = 0; j < columns; j++) {
                if (k != 0 || j != 0) {
                    char xi = seq1[b - 1];
                    int column = current_column + j;
                    char yj = seq2[column - 1];
                    /*printf("me %d, (i, column): (%d, %d) index - gap + column: %d yj: %c\n", me, i, column, index - gap + column, yj);*/
                    int l_value = letter_index(C, xi);
                    int p_value = P[l_value * gap + column];

                    int t = (0 - p_value) < 0? 1 : 0;
                    int s;

                    /*printf("i: %d, b: %d, p_value: %d, lol: %d\n", i, b, p_value, index - gap + p_value - 1);*/
                    if (p_value == 0) {
                        s = (0 - (S[index - gap + column] - t * S[index - gap])) < 0? 1 : 0;
                    } else {
                        s = (0 - (S[index - gap + column] - t * S[index - gap + p_value - 1])) < 0? 1 : 0;
                    }

                    if (xi == yj) cost(i);
                    S[index + column] = S[index - gap + column] + t * (s ^ 1);
                }
            }
            /*print_new_matrix(S + index, 0, seq2_size);*/
        }

        /*if (me != master) print_new_matrix(S, array_size, seq2_size);*/
        if (block_size > 0 && block_high < seq1_size) {
            int next = BLOCK_OWNER(block_high + 1, nprocs, seq1_size + 1);
            if (me == master) {
                MPI_Send(&S[block_high * gap + current_column], columns, MPI_INT, next, 0, MPI_COMM_WORLD);
            } else {
                MPI_Send(&S[block_size * gap + current_column], columns, MPI_INT, next, 0, MPI_COMM_WORLD);
            }
            /*print_new_matrix(S +(block_size - 1) * gap, 0, seq2_size);*/
        }
    }

    if (me != master) {
        MPI_Send(&S[gap], array_size * (seq2_size + 1), MPI_INT, master, 0, MPI_COMM_WORLD);
    }

/*
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0) {
        printf("me: %d\n", me);
        print_new_matrix(S, array_size, seq2_size);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 1) {
        printf("me: %d\n", me);
        print_new_matrix(S, array_size, seq2_size);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 2) {
        printf("me: %d\n", me);
        print_new_matrix(S, array_size, seq2_size);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 3) {
        printf("me: %d\n", me);
        print_new_matrix(S, array_size, seq2_size);
    }
*/

    if (me == master) {
        for (i = 0; i < nprocs; i++) {
            if (i != master) {
                int *array = (int *) malloc(array_size * (seq2_size + 1) * sizeof(int));
                MPI_Recv(array, array_size * (seq2_size + 1), MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
                int source = status.MPI_SOURCE;
                int block_l = BLOCK_LOW(source, nprocs, seq1_size + 1);
                int size_b = BLOCK_SIZE(source, nprocs, seq1_size + 1);

                for (j = block_l, k = 0; k < size_b * (seq2_size + 1); k++) {
                    int index = j * gap + k;
                    S[index] = array[k];
                }

                free(array);
            }
        }

        /*print_new_matrix(S, seq1_size, seq2_size);*/

        i = seq1_size;
        j = seq2_size;
        char xi, yj;
        int last_cell = S[i * gap + j];
        int len = last_cell;
        char lcs[len + 1];
        lcs[len] = '\0';

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


    free(S);
    free(P);
    free(seq1);
    free(seq2);

    MPI_Finalize();

    return 0;
}

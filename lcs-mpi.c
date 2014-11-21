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

    int *matrix = (int*) calloc((seq1_size + 1) * (seq2_size + 1), sizeof(int));
    if (matrix == NULL) {
        printf("Not enough memory. Learning how to program in C might help...\n");
        exit(-1);
    }

    int i, j, k;
    int gap = seq2_size + 1;
    int block_low = BLOCK_LOW(me, nprocs, ALPHA_SIZE);
    int block_high = BLOCK_HIGH(me, nprocs, ALPHA_SIZE);
    int block_size = BLOCK_SIZE(me, nprocs, ALPHA_SIZE);
    int buffer_max_size = ALPHA_SIZE / nprocs + 1;
    int *block = (int*) calloc(buffer_max_size * (seq2_size + 1), sizeof(int));
    int *buffer = (int*) calloc(buffer_max_size * (seq2_size + 1), sizeof(int));


    MPI_Barrier(MPI_COMM_WORLD);

    /* Compute P table */
    if(block_size > 0) {
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


    free(matrix);
    free(P);
    free(seq1);
    free(seq2);

    MPI_Finalize();

    return 0;
}

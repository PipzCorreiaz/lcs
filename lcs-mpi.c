#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    int me, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    int i, j;
    int array[10];
    MPI_Status status;
    int buffer[10];
    int **matrix = (int**) calloc(nprocs, sizeof(int*));
    for (i = 0; i < nprocs; i++) {
        matrix[i] = (int*) calloc(10, sizeof(int));
    }

    for (i = 0; i < 10; i++) {
        array[i] = me;
    }

    for (i = 0; i < nprocs; i++) {
        if (i != me) {
            MPI_Send(&array, 10, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }

    for (i = 0; i < nprocs; i++) {
        if (me != i) {
            MPI_Recv(&buffer, 10, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            for (j = 0; j < 10; j++) {
                matrix[status.MPI_SOURCE][j] = buffer[j];
            }
        } else {
            for (j = 0; j < 10; j++) {
                matrix[me][j] = array[j];
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (!me) {
        printf("Hi from node %d of %d\n", me, nprocs);
        for (i = 0; i < nprocs; i++) {
            for (j = 0; j < 10; j++) {
                printf("%d ", matrix[i][j]);
            }
            printf("\n");
        }
    }

    MPI_Finalize();

    return 0;
}

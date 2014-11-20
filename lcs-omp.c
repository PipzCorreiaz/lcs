#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ALPHA_SIZE 4

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


int main(int argc, char const *argv[]) {
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

    char C[ALPHA_SIZE] = { 'A', 'T', 'C', 'G' };
    int *P;

    P = (int *) calloc(ALPHA_SIZE * (seq2_size + 1), sizeof(int));

	int *matrix = (int*) calloc((seq1_size + 1) * (seq2_size + 1), sizeof(int));
    if (matrix == NULL) {
    	printf("Not enough memory. Learning how to program in C might help...\n");
    	exit(-1);
    }

    int i, j;
    int gap = seq2_size + 1;


    /* Compute P table */
    #pragma omp parallel for private(j)
    for (i = 0; i < ALPHA_SIZE; i++) {
    	int index = i * gap;
    	char c = C[i];
    	for (j = 0; j < seq2_size + 1; j++) {
    		char b = (j - 2) < 0 ? '#' : seq2[j - 2];

    		if (j == 0) {
    			P[index] = 0;
    		} else if (b == c) {
    			P[index + j] = j;
    		} else {
				P[index + j] = P[index + j - 1];
    		}
    	}
    }

	for (i = 1; i < seq1_size + 1; i++) {
		int index = i * gap;
        #pragma omp parallel for
		for (j = 1; j < seq2_size + 1; j++) {
			char xi = seq1[i - 1];
            int l_value = letter_index(C, xi);
            int p_value = P[l_value * gap + j];

            int t = (0 - p_value) < 0? 1 : 0;
            int s = (0 - (matrix[index - gap + j] - t * matrix[index - gap + p_value - 1])) < 0? 1 : 0;

            matrix[index + j] = matrix[index - gap + j] + t * (s ^ cost(i));

		}
	}

    i = seq1_size;
    j = seq2_size;
    char xi, yj;
    int last_cell = matrix[i * gap + j];
    int len = last_cell;
    char lcs[len + 1];
    lcs[len] = '\0';

    while(i > 0 && j > 0) {
        int index = i * gap + j;
        xi = seq1[i - 1];
        yj = seq2[j - 1];

        if (xi == yj) {
            last_cell--;
            lcs[last_cell] = xi;
            i--;
            j--;
        } else {
            int current = matrix[index];
            int top = matrix[index - gap];
            int left = matrix[index - 1];
            int right_of_top = matrix[index - gap + 1];

            if (right_of_top > top) {
                top++;
            }
            if (current > left) {
                left++;
            }

            if (left >= top) {
                j--;
            } else {
                i--;
            }
        }
    }

    printf("%d\n%s\n", len, lcs);

    free(matrix);
    free(seq1);
    free(seq2);


	return 0;
}

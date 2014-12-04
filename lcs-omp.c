#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ALPHA_SIZE 26

char C(int index) {
    return 'A' + index;
}

int letter_index(char letter) {
    return letter - 'A';
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

    int *P = (int *) calloc(ALPHA_SIZE * (seq2_size + 1), sizeof(int));

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
        char c = C(i);
        for (j = 1; j < seq2_size + 1; j++) {
            char b = seq2[j - 1];
            if (b == c) {
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
            char yj = seq2[j - 1];
            int l_value = letter_index(xi);
            int p_value = P[l_value * gap + j];

            int t = (0 - p_value) < 0? 1 : 0;
            int s = (0 - (matrix[index - gap + j] - t * matrix[index - gap + p_value - 1])) < 0? 1 : 0;

            if (xi == yj) cost(i);
            matrix[index + j] = matrix[index - gap + j] + t * (s ^ 1);
		}
	}

    i = seq1_size;
    j = seq2_size;
    char xi, yj;
    int last_cell = matrix[i * gap + j];
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
        } else if (matrix[index - gap] > matrix[index - 1]) {
            i--;
        } else {
            j--;
        }
    }

    printf("%d\n%s\n", len, lcs);

    free(matrix);
    free(P);
    free(seq1);
    free(seq2);

	return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


typedef struct matrix_cell {
	int value;
	char letter;
	struct matrix_cell* parent;
} cell;

void print_matrix(cell** matrix, int size1, int size2) {
    int i, j;
	for (i = 0; i < size1 + 1; i++) {
		for (j = 0; j < size2 + 1; j++) {
			printf("%d", matrix[i][j].value);
		}

        printf("\n");
	}
}

void print_new_matrix(cell* matrix, int size1, int size2) {
    int i, j;
    for (i = 0; i < size1 + 1; i++) {
        for (j = 0; j < size2 + 1; j++) {
            cell* c = &matrix[i * (size1 + 1) + j];
            printf("%d", c->value);
            c = NULL;
        }

        printf("\n");
    }
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

	cell *matrix;

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

	matrix = (cell*) malloc((seq1_size + 1) * (seq2_size + 1) * sizeof(cell));

	cell *current_cell;
	cell *parent;
	cell *top_cell;
	cell *left_cell;

    int i, j;
    int gap = seq2_size + 1;

    if (matrix == NULL) {
    	printf("Not enough memory. Learning how to program in C might help...\n");
    	exit(-1);
    }

	for (i = 1; i < seq1_size + 1; i++) {
		for (j = 1; j < seq2_size + 1; j++) {
			char xi = seq1[i - 1];
			char yj = seq2[j - 1];
            current_cell = &matrix[i * gap + j];
			/* current_cell = &matrix[i][j]; */

			if (xi != yj) {
                top_cell = &matrix[(i - 1) * gap + j];
                left_cell = &matrix[i * gap + (j - 1)];
				/* top_cell = &matrix[i - 1][j];
				left_cell = &matrix[i][j - 1]; */

				if (top_cell->value > left_cell->value) {
					current_cell->parent = top_cell;
					current_cell->value = top_cell->value;
				} else {
					current_cell->parent = left_cell;
					current_cell->value = left_cell->value;
				}

				/*current_cell->letter = '\0';*/
			} else {
				/* parent = &matrix[i - 1][j - 1]; */
                parent = &matrix[(i - 1) * gap + (j - 1)];
				current_cell->value = parent->value + cost(i);
				current_cell->letter = xi;
				current_cell->parent = parent;
			}
		}
	}

	current_cell = NULL;
	parent = NULL;
	top_cell = NULL;
	left_cell = NULL;

	/* cell* last_cell = &matrix[seq1_size][seq2_size]; */
    cell *last_cell = &matrix[seq1_size * gap + seq2_size];
	int len = last_cell->value;
	char lcs[len + 1];
	char letter;
	int pos;

	lcs[len] = '\0';

	while (last_cell != NULL) {
		pos = last_cell->value - 1;
		letter = last_cell->letter;
		if (letter != '\0') {
			lcs[pos] = letter;
		}
		last_cell = last_cell->parent;
	}

    printf("%d\n%s\n", len, lcs);


	for (i = 0; i < seq1_size + 1; i++) {
	    for (j = 0; j < seq2_size + 1; j++) {
	      matrix[i * gap + j].parent = NULL;
	    }
	}

	free(matrix);
	free(seq1);
	free(seq2);


	return 0;
}

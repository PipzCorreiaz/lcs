#include <stdio.h>
#include <stdlib.h>
#include <math.h>


typedef struct matrix_cell {
	int value;
	char letter;
	struct matrix_cell* parent;
} cell;

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

	int *matrix;

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

	matrix = (int*) calloc((seq1_size + 1) * (seq2_size + 1), sizeof(int));
	if (matrix == NULL) {
		printf("Not enough memory. Learning how to program in C might help...\n");
		exit(-1);
	}

	int i, j;
	const int gap = seq2_size + 1;
	const int max_diagonals = abs(seq1_size - seq2_size) + 1;
	const int shortest_seq_size = seq1_size < seq2_size? seq1_size : seq2_size;


	int k, max_j = 0, diagonals_counter = max_diagonals;
	for (i = 1, k = 1; i < seq1_size + 1 && k < seq2_size + 1; ) {
		if (i < shortest_seq_size) {
			max_j = i;
		} else {
			if (diagonals_counter > 0) {
				max_j = shortest_seq_size;
				diagonals_counter--;
			} else {
				max_j--;
			}
		}

		for (j = k; j < max_j + k; j++) {
			int u = (i - j + k);
			int index = u * gap + j;
			char xi = seq1[u - 1];
			char yj = seq2[j - 1];

			if (xi != yj) {
				int top_cell = matrix[index - gap];
				int left_cell = matrix[index - 1];

				if (top_cell > left_cell) {
					matrix[index] = top_cell;
				} else {
					matrix[index] = left_cell;
				}

			} else {
				matrix[index] = matrix[index - gap - 1] + cost(i);
			}
		}

		if (i < seq1_size) {
			i++;
		} else {
			k++;
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
		xi = seq1[i-1];
		yj = seq2[j-1];
		if(xi == yj) {
			lcs[last_cell-1] = xi;
			last_cell--;
			i--;
			j--;
		}
		else {
			if(matrix[(i-1) * gap + j] > matrix[i * gap + (j - 1)]) {
				i--;
			}
			else {
				j--;
			}
		}
	}
	
	printf("%d\n%s\n", len, lcs);

	

	free(matrix);
	free(seq1);
	free(seq2);


	return 0;
}

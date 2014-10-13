#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>

typedef struct matrix_cell {
	int value;
	char letter;
	struct matrix_cell* parent;
} cell;

void print_matrix(cell** matrix, int size1, int size2) {
	for (int i = 0; i < size1 + 1; i++) {
		for (int j = 0; j < size2 + 1; j++) {
			std::cout << matrix[i][j].value;
		}

		std::cout << std::endl;
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
	std::ifstream file;

	int seq1_size = 0;
	std::string seq1;
	int seq2_size = 0;
	std::string seq2;

	cell **matrix;

	file.open(argv[1]);
	if (file.is_open()) {
		file >> seq1_size;
		file >> seq2_size;
		file >> seq1;
		file >> seq2;
	} else {
		std::cout << "Unable to open file"; 
	}

	file.close();

	// Init matrix

	matrix = (cell**)calloc(seq1_size + 1, sizeof(cell*));
	for (int i = 0; i < seq1_size + 1; i++) {
		matrix[i] = (cell*)calloc(seq2_size + 1, sizeof(cell));
	}


	for (int i = 1; i < seq1_size + 1; i++) {
		for (int j = 1; j < seq2_size + 1; j++) {
			char xi = seq1[i - 1];
			char yj = seq2[j - 1];
			cell* current_cell = &matrix[i][j];

			if (xi == yj) {
				cell* parent = &matrix[i - 1][j - 1];
				current_cell->value = parent->value + cost(i);
				current_cell->letter = xi;
				current_cell->parent = parent;
			} else {
				cell* top_cell = &matrix[i][j - 1];
				cell* left_cell = &matrix[i - 1][j];

				if (top_cell->value > left_cell->value) {
					current_cell->parent = top_cell;
					current_cell->value = top_cell->value;
				} else {
					current_cell->parent = left_cell;
					current_cell->value = left_cell->value;
				}

				current_cell->letter = '\0';
			}
		}
	}

	print_matrix(matrix, seq1_size, seq2_size);

	// cell* last_cell = matrix[seq1_size][seq2_size];
	// std::string lcs = "";

	// while (last_cell == NULL) {

	// }
	
	return 0;
}

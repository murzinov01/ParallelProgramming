#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define MATRIX_SIZE 16

int main(int argc, char **argv) {
  int rank, size;
  int i, j;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double matrix[MATRIX_SIZE][MATRIX_SIZE];
  double vector[MATRIX_SIZE];
  double result[MATRIX_SIZE];
  double local_matrix[MATRIX_SIZE / size][MATRIX_SIZE];
  double local_result[MATRIX_SIZE / size];
  MPI_Status status;

  // Initialize matrix and vector
  int n = 0;
  if (rank == 0) {
    for (i = 0; i < MATRIX_SIZE; i++) {
      for (j = 0; j < MATRIX_SIZE; j++) {
        matrix[i][j] = ++n;
      }
    }

    for (i = 0; i < MATRIX_SIZE; i++) {
      vector[i] = i;
    }
  }

  // Share vector
  MPI_Bcast(vector, MATRIX_SIZE, MPI_DOUBLE,
            0, MPI_COMM_WORLD);
  // Scatter matrix rows
  MPI_Scatter(matrix, MATRIX_SIZE * MATRIX_SIZE / size, MPI_DOUBLE,
              local_matrix, MATRIX_SIZE * MATRIX_SIZE / size, MPI_DOUBLE,
              0, MPI_COMM_WORLD);

  // Perform local matrix-vector multiplication
  for (i = 0; i < MATRIX_SIZE / size; ++i) {
    local_result[i] = 0.0;
    for (j = 0; j < MATRIX_SIZE; j++) {
      local_result[i] += local_matrix[i][j] * vector[j];
    }
  }

  // Gather local results
  MPI_Gather(local_result, MATRIX_SIZE / size, MPI_DOUBLE,
             result, MATRIX_SIZE / size, MPI_DOUBLE,
             0, MPI_COMM_WORLD);

  // Print result
  if (rank == 0) {
    printf("Result:\n");
    for (i = 0; i < MATRIX_SIZE; i++) {
      printf("%f\n", result[i]);
    }
  }

  MPI_Finalize();
  return 0;
}

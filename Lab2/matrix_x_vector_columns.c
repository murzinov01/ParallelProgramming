#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define MATRIX_SIZE 32

int main(int argc, char **argv) {
  int rank, size;
  int i, j;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double matrix[MATRIX_SIZE][MATRIX_SIZE];
  double vector[MATRIX_SIZE];
  double result_from_processes[size][MATRIX_SIZE];
  double result[MATRIX_SIZE];
  double local_matrix[MATRIX_SIZE / size][MATRIX_SIZE];
  double local_result[MATRIX_SIZE];
  double local_vector[MATRIX_SIZE / size];
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

  // Transpose matrix
  double temp;
  for (i = 0; i < MATRIX_SIZE; i++) {
    for (j = i + 1; j < MATRIX_SIZE; j++) {
      temp = matrix[i][j];
      matrix[i][j] = matrix[j][i];
      matrix[j][i] = temp;
    }
  }

  // Scatter vector
  MPI_Scatter(vector, MATRIX_SIZE / size, MPI_DOUBLE,
              local_vector, MATRIX_SIZE / size, MPI_DOUBLE,
              0, MPI_COMM_WORLD);
  // Scatter matrix columns
  MPI_Scatter(matrix, MATRIX_SIZE * MATRIX_SIZE / size, MPI_DOUBLE,
              local_matrix, MATRIX_SIZE * MATRIX_SIZE / size, MPI_DOUBLE,
              0, MPI_COMM_WORLD);

  // Perform local matrix-vector multiplication
  for (i = 0; i < MATRIX_SIZE; ++i) {
    local_result[i] = 0;
    for (j = 0; j < MATRIX_SIZE / size; j++) {
      local_result[i] += local_matrix[j][i] * local_vector[j];
    }
  }

  // Gather local results
  MPI_Gather(local_result, MATRIX_SIZE, MPI_DOUBLE,
             result_from_processes, MATRIX_SIZE, MPI_DOUBLE,
             0, MPI_COMM_WORLD);

  // Sum vectors from processes
  for (int j = 0; j < MATRIX_SIZE; ++j) {
    result[j] = 0.0;
    for (int i = 0; i < size; ++i) {
      result[j] += result_from_processes[i][j];
    }
  }

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

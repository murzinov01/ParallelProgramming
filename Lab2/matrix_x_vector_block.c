#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define MATRIX_SIZE 16
#define MATRIX_TAG 1
#define VECTOR_TAG 2
#define RESULT_TAG 3

int main(int argc, char **argv) {
  int rank, size, procceses;
  int i, j;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procceses);

  if (procceses < 4) {
    printf("procceses number is incompatible with matrix size\n");
    exit(-1);
  }
  size = (int)sqrt(procceses);
  if (MATRIX_SIZE % size) {
    printf("procceses number is incompatible with matrix size\n");
    exit(-1);
  }

  double matrix[MATRIX_SIZE][MATRIX_SIZE];
  double vector[MATRIX_SIZE];
  double result_from_processes[size][MATRIX_SIZE];
  double result[MATRIX_SIZE];

  double local_matrix[MATRIX_SIZE / size][MATRIX_SIZE / size];
  double local_vector[MATRIX_SIZE / size];
  double local_result[MATRIX_SIZE / size];

  double matrix_block[MATRIX_SIZE / size][MATRIX_SIZE / size];
  double vector_block[MATRIX_SIZE / size];
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

    // Send Matrix blocks
    // choose block
    int send_to = 0;
    for (i = 0; i < size; ++i) {
      for (j = 0; j < size; ++j) {
        // fill block
        for (int k = 0; k < MATRIX_SIZE / size; ++k) {
          for (int m = 0; m < MATRIX_SIZE / size; ++m) {
            int offset_k = i * MATRIX_SIZE / size;
            int offset_m = j * MATRIX_SIZE / size;
            matrix_block[k][m] = matrix[offset_k + k][offset_m + m];
          }
        }
        MPI_Send(matrix_block, pow(MATRIX_SIZE / size, 2), MPI_DOUBLE,
                send_to++, MATRIX_TAG, MPI_COMM_WORLD);
      }
    }
    // Send Vector blocks
    // choose block
    for (i = 0; i < size; ++i) {
      // fill block
      for (int k = 0; k < MATRIX_SIZE / size; ++k) {
        int offset = i * MATRIX_SIZE / size;
        vector_block[k] = offset + k;
      }
      for (j = 0; j < size; ++j) {
        MPI_Send(vector_block, MATRIX_SIZE / size, MPI_DOUBLE,
                i + size * j, VECTOR_TAG, MPI_COMM_WORLD);
      }
    }
  }

  // Receive Matrix block
  MPI_Recv(local_matrix, pow(MATRIX_SIZE / size, 2), MPI_DOUBLE,
           0, MATRIX_TAG, MPI_COMM_WORLD, &status);
  MPI_Recv(local_vector, MATRIX_SIZE / size, MPI_DOUBLE,
           0, VECTOR_TAG, MPI_COMM_WORLD, &status);

  // Perform local matrix-vector multiplication
  for (i = 0; i < MATRIX_SIZE / size; ++i) {
    local_result[i] = 0.0;
    for (j = 0; j < MATRIX_SIZE / size; j++) {
      local_result[i] += local_matrix[i][j] * local_vector[j];
    }
  }

  // Send result
  MPI_Send(local_result, MATRIX_SIZE / size, MPI_DOUBLE,
           0, RESULT_TAG, MPI_COMM_WORLD);

  // initialize result
  for (i = 0; i < MATRIX_SIZE; ++i) {
    result[i] = 0.0;
  }

  // Receive result
  if (rank == 0) {
    for (i = 0; i < procceses; ++i) {
      MPI_Recv(local_result, MATRIX_SIZE / size, MPI_DOUBLE,
               i, RESULT_TAG, MPI_COMM_WORLD, &status);
      int vector_result_offset = (i / size) * (MATRIX_SIZE / size);
      for (int k = 0; k < MATRIX_SIZE / size; ++k) {
        result[vector_result_offset + k] += local_result[k];
      }
    }
    printf("Result:\n");
    for (i = 0; i < MATRIX_SIZE; i++) {
      printf("%f\n", result[i]);
    }
  }

  MPI_Finalize();
  return 0;
}

#include <cstdio>
#include "utils.hpp"

/// Pretty-print a general matrix in a traditional format.
void print_mat(double **mat, int dim_rows, int dim_cols) {
  // first, handle the column labels
  printf("     ");
  for (int i = 0; i < dim_cols; i++)
    printf("%11d", i+1);
  printf("\n");
  // then, handle the row labels
  for (int i = 0; i < dim_rows; i++) {
    printf("%5d", i+1);
    // print the matrix data
    for (int j = 0; j < dim_cols; i++)
      printf("%11.7f", &mat[i][j]);
    printf("\n");
  }
}

#include <cstdio>
#include "utils.hpp"

/// Pretty-print a general matrix in a traditional format.
void print_mat(double **mat, int dim_rows, int dim_cols) {
  // first, handle the column labels
  printf("     ");
  for (int i = 0; i < dim_cols; i++)
    printf("%12d", i+1);
  printf("\n");
  // then, handle the row labels
  for (int i = 0; i < dim_rows; i++) {
    printf("%5d", i+1);
    // print the matrix data
    for (int j = 0; j < dim_cols; j++)
      printf("%12.7f", mat[i][j]);
    printf("\n");
  }
}

int compound_index_2(int i, int j) {
  if (i < j)
    return i*(i+1)/2 + j;
  else
    return j*(j+1)/2 + i;
}

int compound_index_4(int i, int j, int k, int l) {
  int ij = compound_index_2(i,j);
  int kl = compound_index_2(k,l);
  return compound_index_2(ij,kl);
}

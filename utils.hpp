#ifndef UTILS_HPP
#define UTILS_HPP

void tred2(int n, double **a, double *d, double *e, int matz);
void tqli(int n, double *d, double **z, double *e, int matz, double toler);
void eigsort(double *d, double **v, int n);
double *init_array(int size);
double **init_matrix(int n,int m);
void free_matrix(double **array, int size);
// void print_mat(double **a, int m, int n, FILE *out);
void diag(int nm, int n, double **array, double *e_vals, int matz, double **e_vecs, double toler);

void print_mat(double **mat, int dim_rows, int dim_cols);

#endif /* UTILS_HPP */

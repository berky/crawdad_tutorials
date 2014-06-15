#include <cstdio>
#include <cmath>
#include "../utils.hpp"

int main()
{
  int i, j, k, l;
  double val;
  int mu, nu, lam, sig;

  /** Step 1: Nuclear Repulsion Energy
   * Read the nuclear repulsion energy from the enuc.dat file.
   */

  FILE *enuc_file;
  enuc_file = fopen("h2o_sto3g_enuc.dat", "r");
  double Vnn;
  fscanf(enuc_file, "%lf", &Vnn);
  fclose(enuc_file);

  /** Step 2: One-Electron Integrals
   * Read the AO-basis overlap, kinetic energy, and nuclear attraction
   * and store them in appropriately constructed matrices. Then form the
   * "core Hamiltonian":
   *  H_{\mu\nu}^{core} = T_{\mu\nu} + V_{\mu\nu}.
   * Note that the one-electron integrals provided include only the
   * permutationally unique integrals, but you should store the full matrices
   * for convenience. Note also that the AO indices on the integrals in the
   * files start with 1 rather than 0.
   */

  int NBasis = 7;
  double **S_AO = new double* [NBasis];
  double **T_AO = new double* [NBasis];
  double **V_AO = new double* [NBasis];
  double **H_AO_Core = new double* [NBasis];
  for (int i = 0; i < NBasis; i++) {
    S_AO[i] = new double[NBasis];
    T_AO[i] = new double[NBasis];
    V_AO[i] = new double[NBasis];
    H_AO_Core[i] = new double[NBasis];
  }

  FILE *S_AO_file, *T_AO_file, *V_AO_file;
  S_AO_file = fopen("h2o_sto3g_s.dat", "r");
  T_AO_file = fopen("h2o_sto3g_t.dat", "r");
  V_AO_file = fopen("h2o_sto3g_v.dat", "r");

  while (fscanf(S_AO_file, "%d %d %lf", &i, &j, &val) != EOF)
    S_AO[i-1][j-1] = S_AO[j-1][i-1] = val;
  while (fscanf(T_AO_file, "%d %d %lf", &i, &j, &val) != EOF)
    T_AO[i-1][j-1] = T_AO[j-1][i-1] = val;
  while (fscanf(V_AO_file, "%d %d %lf", &i, &j, &val) != EOF)
    V_AO[i-1][j-1] = V_AO[j-1][i-1] = val;

  fclose(S_AO_file);
  fclose(T_AO_file);
  fclose(V_AO_file);

  printf("AO Overlap Integrals [S_AO]:\n");
  print_mat(S_AO, NBasis, NBasis);
  printf("AO Kinetic Energy Integrals [T_AO]:\n");
  print_mat(T_AO, NBasis, NBasis);
  printf("AO Nuclear Attraction Integrals [V_AO]:\n");
  print_mat(V_AO, NBasis, NBasis);

  for (int i = 0; i < NBasis; i++)
    for (int j = 0; j < NBasis; j++)
      H_AO_Core[i][j] = H_AO_Core[j][i] = (T_AO[i][j] + V_AO[i][j]);

  printf("AO Core Hamiltonian [H_AO_Core]:\n");
  print_mat(H_AO_Core, NBasis, NBasis);

  /** Step #3: Two-Electron Integrals
   * Read the two-electron repulsion integrals from the eri.dat file. The
   * integrals in this file are provided in Mulliken notation over real AO
   * basis functions:
   *  ...
   * Hence, the integrals obey the eight-fold permutational symmetry
   * relationships:
   *  ...
   * and only the permutationally unique integrals are provided in the file,
   * with the restriction that, for each integral, the following relationships
   * hold:
   *  ...
   * where
   *  ...
   * Note that the two-electron integrals may be store efficiently in a
   * one-dimensional array and the above relationship used to map between
   * given \mu, \nu, \lambda, and \sigma indices and a "compound index" defined
   * as:
   *  ...
   */

  /**
   * This block uses the full NBasis*NBasis matrix allocation.
   */
  double ****ERI_AO = new double*** [NBasis];
  for (int i = 0; i < NBasis; i++) {
    ERI_AO[i] = new double** [NBasis];
    for (int j = 0; j < NBasis; j++) {
      ERI_AO[i][j] = new double* [NBasis];
      for (int k = 0; k < NBasis; k++) {
	ERI_AO[i][j][k] = new double[NBasis];
      }
    }
  }

  FILE *ERI_AO_file;
  ERI_AO_file = fopen("h2o_sto3g_eri.dat", "r");

  while (fscanf(ERI_AO_file, "%d %d %d %d %lf", &i, &j, &k, &l, &val) != EOF) {
    mu = i-1; nu = j-1; lam = k-1; sig = l-1;
    ERI_AO[mu][nu][lam][sig] 
      = ERI_AO[mu][nu][sig][lam] 
      = ERI_AO[nu][mu][lam][sig] 
      = ERI_AO[nu][mu][sig][lam] 
      = ERI_AO[lam][sig][mu][nu] 
      = ERI_AO[lam][sig][nu][mu] 
      = ERI_AO[sig][lam][mu][nu] 
      = ERI_AO[sig][lam][nu][mu] = val;
  }

  fclose(ERI_AO_file);

  /** Step #4: Build the Orthogonalization Matrix
   * Diagonalize the overlap matrix:
   *  ...
   * where \mathbf{L}_{S} is the matrix of eigenvectors (columns) and 
   * \mathbf{\Lambda}_{S} is the diagonal matrix of corresponding eigenvalues.
   * Build the symmetric orthogonalization matrix using:
   */

  double* Lam_S_AO = new double[NBasis];
  double **Lam_S_AO_mat = new double* [NBasis];
  double **L_S_AO = new double* [NBasis];
  double **Lam_sqrt_inv_AO = new double* [NBasis];
  double **Symm_Orthog = new double* [NBasis];
  for (int i = 0; i < NBasis; i++) {
    Lam_S_AO_mat[i] = new double[NBasis];
    L_S_AO[i] = new double[NBasis];
    Lam_sqrt_inv_AO[i] = new double[NBasis];
    Symm_Orthog[i] = new double[NBasis];
  }
  
  diag(NBasis, NBasis, S_AO, Lam_S_AO, true, L_S_AO, 1.0e-13);

  for (int i = 0; i < NBasis; i++) {
    for (int j = 0; j < NBasis; j++) {
      Lam_S_AO_mat[i][j] = 0.0;
    }
    Lam_S_AO_mat[i][i] = Lam_S_AO[i];
  }

  printf("matrix of eigenvectors (columns) [L_S_AO]:\n");
  print_mat(L_S_AO, NBasis, NBasis);
  printf("diagonal matrix of corresponding eigenvalues [Lam_S_AO]:\n");
  print_mat(Lam_S_AO_mat, NBasis, NBasis);

  // take the square root of the inverse of Lam_S_AO_mat (element-wise)
  for (int i = 0; i < NBasis; i++) {
    for (int j = 0; j < NBasis; j++) {
      if (i == j)
	Lam_sqrt_inv_AO[i][j] = sqrt(1.0/Lam_S_AO_mat[i][j]);
      else
	Lam_sqrt_inv_AO[i][j] = 0.0;
    }
  }

  // build the symmetric orthogonalization matrix as L_S_AO * Lam_sqrt_inv_AO * L_S_AO.t()
  /* steps:
   * 1. allocate tmp
   * 2. tmp = A * B.t + tmp
   * 3. C = B * tmp + C
   * 4. free tmp
   */
  double** tmp = init_matrix(NBasis, NBasis);
  mmult(Lam_sqrt_inv_AO, 0, L_S_AO, 1, tmp, NBasis, NBasis, NBasis);
  mmult(L_S_AO, 0, tmp, 0, Symm_Orthog, NBasis, NBasis, NBasis);
  free_matrix(tmp, NBasis);

  printf("Symmetric Orthogonalization Matrix [S^-1/2]:\n");
  print_mat(Symm_Orthog, NBasis, NBasis);

  /// Clean up after ourselves...
  for (int i = 0; i < NBasis; i++) {
    delete[] S_AO[i]; delete[] T_AO[i]; delete[] V_AO[i]; delete[] H_AO_Core[i];
  }
  delete[] S_AO; delete[] T_AO; delete[] V_AO; delete[] H_AO_Core;

  for (int i = 0; i < NBasis; i++) {
    for (int j = 0; j < NBasis; j++) {
      for (int k = 0; k < NBasis; k++) {
	delete[] ERI_AO[i][j][k];
      }
      delete[] ERI_AO[i][j];
    }
    delete[] ERI_AO[i];
  }
  delete[] ERI_AO;
  for (int i = 0; i < NBasis; i++) {
    delete[] L_S_AO[i];
    delete[] Lam_S_AO_mat[i];
    delete[] Lam_sqrt_inv_AO[i];
    delete[] Symm_Orthog[i];
  }
  delete[] L_S_AO; delete[] Lam_S_AO; delete[] Lam_S_AO_mat;
  delete[] Lam_sqrt_inv_AO; delete[] Symm_Orthog;

  return 0;
}

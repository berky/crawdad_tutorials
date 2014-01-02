#include <cstdio>
#include <cmath>
#include "utils.hpp"

int main()
{
  int i, j, k, l;
  double val;

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

  /// Clean up after ourselves...
  for (int i = 0; i < NBasis; i++) {
    delete[] S_AO[i]; delete[] T_AO[i]; delete[] V_AO[i]; delete[] H_AO_Core[i];
  }
  delete[] S_AO; delete[] T_AO; delete[] V_AO; delete[] H_AO_Core;

  return 0;
}

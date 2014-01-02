#include <cstdio>
#include <cmath>
#include "utils.hpp"

int main()
{
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
  for (int i = 0; i < NBasis; i++) {
    S_AO[i] = new double[NBasis];
    T_AO[i] = new double[NBasis];
    V_AO[i] = new double[NBasis];
  }

  FILE *S_AO_file, *T_AO_file, *V_AO_file;
  S_AO_file = fopen("h2o_sto3g_s.dat", "r");
  T_AO_file = fopen("h2o_sto3g_t.dat", "r");
  V_AO_file = fopen("h2o_sto3g_v.dat", "r");

  /// Clean up after ourselves...
  for (int i = 0; i < NBasis; i++) {
    delete[] S_AO[i]; delete[] T_AO[i]; delete[] V_AO[i];
  }
  delete[] S_AO; delete[] T_AO; delete[] V_AO;

  return 0;
}

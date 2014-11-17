#include <cstdio>
#include <cmath>
#include <cstdlib>
#include "../utils.hpp"

double calc_elec_energy(double** P, double** H, double** F, int NBasis) {
  double energy = 0.0;
  for (int mu = 0; mu < NBasis; mu++)
    for (int nu = 0; nu < NBasis; nu++)
      energy += P[mu][nu] * (H[mu][nu] + F[mu][nu]);
  return energy;
}

void build_fock(double** F, double** P, double** H, double**** ERI, int NBasis) {
  for (int mu = 0; mu < NBasis; mu++)
    for (int nu = 0; nu < NBasis; nu++) {
      F[mu][nu] = H[mu][nu];
      for (int lm = 0; lm < NBasis; lm++)
        for (int sg = 0; sg < NBasis; sg++)
          F[mu][nu] += P[lm][sg] * (2*ERI[mu][nu][lm][sg] -
                                    ERI[mu][lm][nu][sg]);
    }
}

void build_fock_prime(double** F_prime, double** Symm_Orthog, double** F, int NBasis) {
  double** tmp = init_matrix(NBasis, NBasis);
  mmult(F, 0, Symm_Orthog, 0, tmp, NBasis, NBasis, NBasis);
  mmult(Symm_Orthog, 1, tmp, 0, F_prime, NBasis, NBasis, NBasis);
  free_matrix(tmp, NBasis);
}

double rmsd_density(double** P_new, double** P_old, int NBasis) {
  double rmsd = 0.0;
  for (int mu = 0; mu < NBasis; mu++)
    for (int nu = 0; nu < NBasis; nu++)
      rmsd += pow(P_new[mu][nu] - P_old[mu][nu], 2);
  return sqrt(rmsd);
}

int main()
{
  int i, j, k, l;
  double val;
  int mu, nu, lm, sg;

  /**
   * Step 1: Nuclear Repulsion Energy
   */

  FILE *enuc_file;
  enuc_file = fopen("h2o_sto3g_enuc.dat", "r");
  double vnn;
  fscanf(enuc_file, "%lf", &vnn);
  fclose(enuc_file);

  /**
   * Step 2: One-Electron Integrals
   */

  int NElec = 10;
  int NOcc = NElec / 2;
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

  /**
   * Step #3: Two-Electron Integrals
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
    mu = i-1; nu = j-1; lm = k-1; sg = l-1;
    ERI_AO[mu][nu][lm][sg]
      = ERI_AO[mu][nu][sg][lm]
      = ERI_AO[nu][mu][lm][sg]
      = ERI_AO[nu][mu][sg][lm]
      = ERI_AO[lm][sg][mu][nu]
      = ERI_AO[lm][sg][nu][mu]
      = ERI_AO[sg][lm][mu][nu]
      = ERI_AO[sg][lm][nu][mu] = val;
  }

  fclose(ERI_AO_file);

  double** F = init_matrix(NBasis, NBasis);

  /**
   * Step #4: Build the Orthogonalization Matrix
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

  /// take the square root of the inverse of Lam_S_AO_mat (element-wise)
  for (int i = 0; i < NBasis; i++) {
    for (int j = 0; j < NBasis; j++) {
      if (i == j)
        Lam_sqrt_inv_AO[i][j] = sqrt(1.0/Lam_S_AO_mat[i][j]);
      else
        Lam_sqrt_inv_AO[i][j] = 0.0;
    }
  }

  /** build the symmetric orthogonalization matrix as L_S_AO * Lam_sqrt_inv_AO * L_S_AO.t()
   * steps:
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

  /**
   * Step #5: Build the Initial (Guess) Density
   */

  /**
   * F_prime = Symm_Orthog.t() * H * Symm_Orthog;
   * steps:
   * 1. allocate tmp
   * 2. tmp = H * Symm_Orthog + tmp
   * 3. F_prime = Symm_Orthog.t() * tmp + F_prime
   * 4. free tmp
   */
  double** F_prime = init_matrix(NBasis, NBasis);
  //zero_matrix(F_prime, NBasis, NBasis);
  tmp = init_matrix(NBasis, NBasis);
  //zeros_matrix(tmp, NBasis, NBasis);
  mmult(H_AO_Core, 0, Symm_Orthog, 0, tmp, NBasis, NBasis, NBasis);
  mmult(Symm_Orthog, 1, tmp, 0, F_prime, NBasis, NBasis, NBasis);
  free_matrix(tmp, NBasis);

  printf("Initial (guess) Fock Matrix [F_prime_0_AO]:\n");
  print_mat(F_prime, NBasis, NBasis);

  /**
   * Diagonalize the Fock Matrix
   */

  double* eps_vec = init_array(NBasis);
  double** eps_mat = init_matrix(NBasis, NBasis);
  double** C_prime = init_matrix(NBasis, NBasis);
  diag(NBasis, NBasis, F_prime, eps_vec, true, C_prime, 1e-13);
  for (int i = 0; i < NBasis; i++) {
    for (int j = 0; j < NBasis; j++)
      eps_mat[i][j] = 0.0;
    eps_mat[i][i] = eps_vec[i];
  }

  printf("Initial MO Coefficients [C_prime_0_AO]:\n");
  print_mat(C_prime, NBasis, NBasis);
  printf("Initial Orbital Energies [eps_0_AO]:\n");
  print_mat(eps_mat, NBasis, NBasis);

  /**
   * Transform the eigenvectors into the original (non-orthogonal) AO basis
   */

  double** C = init_matrix(NBasis, NBasis);
  mmult(Symm_Orthog, 0, C_prime, 0, C, NBasis, NBasis, NBasis);
  printf("Initial MO Coefficients (non-orthogonal) [C_0_AO]:\n");
  print_mat(C, NBasis, NBasis);

  /**
   * Build the density matrix using the occupied MOs
   */

  double** D = init_matrix(NBasis, NBasis);
  double** D_old = init_matrix(NBasis, NBasis);
  mmult(C, 0, C, 1, D, NBasis, NBasis, NOcc);
  printf("Initial Density Matrix [D_0]:\n");
  print_mat(D, NBasis, NBasis);

  /**
   * Step 6: Compute the Initial SCF Energy
   */

  double thresh_E = 1.0e-15;
  double thresh_D = 1.0e-7;
  int iter = 0;
  int max_iter = 1024;
  double E_total, E_elec_old, E_elec_new, delta_E, rmsd_D;

  E_elec_new = calc_elec_energy(D, H_AO_Core, H_AO_Core, NBasis);
  E_total = E_elec_new + vnn;
  delta_E = E_total;
  printf("%4d %20.12f %20.12f %20.12f\n",
         iter, E_elec_new, E_total, delta_E);
  iter++;

  /**
   * Start the SCF Iterative Procedure
   */

  while (iter < max_iter) {
    build_fock(F, D, H_AO_Core, ERI_AO, NBasis);
    build_fock_prime(F_prime, Symm_Orthog, F, NBasis);
    diag(NBasis, NBasis, F_prime, eps_vec, true, C_prime, 1e-13);
    // transform MO coefficients into original non-orthogonal AO basis
    mmult(Symm_Orthog, 0, C_prime, 0, C, NBasis, NBasis, NBasis);
    D_old = D;
    // build the density matrix
    mmult(C, 0, C, 1, D, NBasis, NBasis, NOcc);
    E_elec_old = E_elec_new;
    E_elec_new = calc_elec_energy(D, H_AO_Core, F, NBasis);
    E_total = E_elec_new + vnn;
    if (iter == 1)
      printf("%4d %20.12f %20.12f %20.12f\n",
             iter, E_elec_new, E_total, delta_E);
    else
      printf("%4d %20.12f %20.12f %20.12f %20.12f\n",
             iter, E_elec_new, E_total, delta_E, rmsd_D);
    delta_E = E_elec_new - E_elec_old;
    rmsd_D = rmsd_density(D, D_old, NBasis);
    if (delta_E < thresh_E && rmsd_D < thresh_D) {
      printf("Convergence achieved.\n");
      break;
    }
    F = F_prime;
    iter++;
  };

  /// Clean up after ourselves...
  for (int i = 0; i < NBasis; i++) {
    for (int j = 0; j < NBasis; j++) {
      for (int k = 0; k < NBasis; k++) {
        delete[] ERI_AO[i][j][k];
      }
      delete[] ERI_AO[i][j];
    }
    delete[] S_AO[i]; delete[] T_AO[i]; delete[] V_AO[i]; delete[] H_AO_Core[i];
    delete[] L_S_AO[i];
    delete[] Lam_S_AO_mat[i];
    delete[] Lam_sqrt_inv_AO[i];
    delete[] Symm_Orthog[i];
    delete[] ERI_AO[i];
  }
  delete[] S_AO; delete[] T_AO; delete[] V_AO; delete[] H_AO_Core;
  delete[] ERI_AO;
  delete[] L_S_AO; delete[] Lam_S_AO; delete[] Lam_S_AO_mat;
  delete[] Lam_sqrt_inv_AO; delete[] Symm_Orthog;

  free_matrix(F_prime, NBasis);
  free(eps_vec);
  free_matrix(C_prime, NBasis);
  free_matrix(eps_mat, NBasis);
  free_matrix(C, NBasis);
  free_matrix(D, NBasis);
  //free_matrix(D_old, NBasis);
  free_matrix(F, NBasis);

  return 0;
}

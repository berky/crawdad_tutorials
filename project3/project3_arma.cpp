#include <cstdio>
#include <cmath>
#include <iostream>
#include "../utils.hpp"

using namespace std;

/*!
 * @brief Calculate the Hartree-Fock electronic energy.
 */
double calc_elec_energy(const arma::mat &P,
                        const arma::mat &H,
                        const arma::mat &F) {
  return arma::accu(P%(H+F));
}

/*!
 * @brief Form the density matrix from the MO coefficients.
 */
void make_density(arma::mat &P,
                  const arma::mat &C,
                  const int NOcc) {
  P = C.cols(0, NOcc-1) * C.cols(0, NOcc-1).t();
}

/*!
 * @brief Build the Fock matrix from the density, one-electron,
 * and two-electron integrals.
 */
void build_fock(arma::mat &F,
                const arma::mat &P,
                const arma::mat &H,
                const arma::vec &ERI) {
  for (int mu = 0; mu < H.n_rows; mu++) {
    for (int nu = 0; nu < H.n_cols; nu++) {
      F(mu,nu) = H(mu,nu);
      for (int lam = 0; lam < P.n_rows; lam++)
        for (int sig = 0; sig < P.n_cols; sig++)
          F(mu,nu) += P(lam,sig) * (2*ERI(idx4(mu,nu,lam,sig)) - ERI(idx4(mu,lam,nu,sig)));
    }
  }
}

/*!
 * @brief Calculate the RMS deviation between two density matrices.
 */
double rmsd_density(const arma::mat &P_new,
                    const arma::mat &P_old) {
  return sqrt(arma::accu(arma::pow((P_new - P_old), 2)));
}

/*!
 * @brief Perform Hartree "damping" by mixing a fraction of old density
 * with the new density to aid convergence.
 */
void mix_density(arma::mat &P_new,
                 const arma::mat &P_old,
                 const double alpha) {
  // alpha must be in the range [0,1)
  P_new = ((1.0-alpha)*P_new) + (alpha*P_old);
}

int main()
{
  int i, j, k, l;
  double val;
  int mu, nu, lam, sig;

  cout.width(12);
  cout.precision(7);
  cout.setf(ios::fixed);

  FILE *enuc_file;
  enuc_file = fopen("h2o_sto3g_enuc.dat", "r");
  double Vnn;
  fscanf(enuc_file, "%lf", &Vnn);
  fclose(enuc_file);

  printf("Nuclear Repulsion Energy: %12f\n", Vnn);

  int NElec = 10;
  int NOcc = NElec/2;
  int NBasis = 7;
  int M = idx4(NBasis, NBasis, NBasis, NBasis);

  arma::mat S = arma::mat(NBasis, NBasis);
  arma::mat T = arma::mat(NBasis, NBasis);
  arma::mat V = arma::mat(NBasis, NBasis);
  arma::mat H = arma::mat(NBasis, NBasis);
  arma::mat F = arma::mat(NBasis, NBasis, arma::fill::zeros);
  arma::mat F_prime = arma::mat(NBasis, NBasis, arma::fill::zeros);
  arma::mat D = arma::mat(NBasis, NBasis, arma::fill::zeros);
  arma::mat D_old = arma::mat(NBasis, NBasis, arma::fill::zeros);
  arma::mat C = arma::mat(NBasis, NBasis);

  arma::vec eps_vec = arma::vec(NBasis);
  arma::mat eps_mat = arma::mat(NBasis, NBasis);
  arma::mat C_prime = arma::mat(NBasis, NBasis);

  arma::vec Lam_S_vec = arma::vec(NBasis);
  arma::mat Lam_S_mat = arma::mat(NBasis, NBasis);
  arma::mat L_S = arma::mat(NBasis, NBasis);

  FILE *S_file, *T_file, *V_file;
  S_file = fopen("h2o_sto3g_s.dat", "r");
  T_file = fopen("h2o_sto3g_t.dat", "r");
  V_file = fopen("h2o_sto3g_v.dat", "r");

  while (fscanf(S_file, "%d %d %lf", &i, &j, &val) != EOF)
    S(i-1, j-1) = S(j-1, i-1) = val;
  while (fscanf(T_file, "%d %d %lf", &i, &j, &val) != EOF)
    T(i-1, j-1) = T(j-1, i-1) = val;
  while (fscanf(V_file, "%d %d %lf", &i, &j, &val) != EOF)
    V(i-1, j-1) = V(j-1, i-1) = val;

  fclose(S_file);
  fclose(T_file);
  fclose(V_file);

  cout << "AO Overlap Integrals [S_AO]:" << endl; print_arma_mat(S);
  cout << "AO Kinetic Energy Integrals [T_AO]:" << endl; print_arma_mat(T);
  cout << "AO Nuclear Attraction Integrals [V_AO]:" << endl; print_arma_mat(V);

  H = T + V;
  cout << "AO Core Hamiltonian [H_AO_Core]:" << endl; print_arma_mat(H);

  arma::vec ERI = arma::vec(M, arma::fill::zeros);

  FILE *ERI_file;
  ERI_file = fopen("h2o_sto3g_eri.dat", "r");

  while (fscanf(ERI_file, "%d %d %d %d %lf", &i, &j, &k, &l, &val) != EOF) {
    mu = i-1; nu = j-1; lam = k-1; sig = l-1;
    ERI(idx4(mu,nu,lam,sig)) = val;
  }

  fclose(ERI_file);

  double thresh_E = 1.0e-15;
  double thresh_D = 1.0e-7;
  int iter = 0;
  int max_iter = 1024;
  double E_total, E_elec_old, E_elec_new, delta_E, rmsd_D;

  /*!
   * Build the Orthogonalization Matrix
   */
  arma::eig_sym(Lam_S_vec, L_S, S);
  Lam_S_mat = arma::diagmat(Lam_S_vec);
  cout << "matrix of eigenvectors (columns) [L_S_AO]:" << endl; print_arma_mat(L_S);
  cout << "diagonal matrix of corresponding eigenvalues [Lam_S_AO]:" << endl; print_arma_mat(Lam_S_mat);

  arma::mat Lam_sqrt_inv = arma::sqrt(arma::inv(Lam_S_mat));
  arma::mat Symm_Orthog = L_S * Lam_sqrt_inv * L_S.t();
  cout << "Symmetric Orthogonalization Matrix [S^-1/2]:" << endl; print_arma_mat(Symm_Orthog);

  /*!
   * Build the Initial (Guess) Density
   */
  F_prime = Symm_Orthog.t() * H * Symm_Orthog;
  cout << "Initial (guess) Fock Matrix [F_prime_0_AO]:" << endl; print_arma_mat(F_prime);

  /*!
   * Diagonalize the Fock Matrix
   */
  arma::eig_sym(eps_vec, C_prime, F_prime);
  eps_mat = arma::diagmat(eps_vec);
  cout << "Initial MO Coefficients [C_prime_0_AO]:" << endl; print_arma_mat(C_prime);
  cout << "Initial Orbital Energies [eps_0_AO]:" << endl; print_arma_mat(eps_mat);

  /*!
   * Transform the eigenvectors into the original (non-orthogonal) AO basis
   */
  C = Symm_Orthog * C_prime;
  cout << "Initial MO Coefficients (non-orthogonal) [C_0_AO]:" << endl; print_arma_mat(C);

  /*!
   * Build the density matrix using the occupied MOs
   */
  make_density(D, C, NOcc);
  cout << "Initial Density Matrix [D_0]:" << endl; print_arma_mat(D);

  /*!
   * Compute the Initial SCF Energy
   */
  E_elec_new = calc_elec_energy(D, H, H);
  E_total = E_elec_new + Vnn;
  delta_E = E_total;
  // printf("%4c %20c %20c %20c %20c\n",
  //        "Iter",
  //        "E_elec",
  //        "E_tot",
  //        "delta_E",
  //        "rmsd_D");
  printf("%4d %20.12f %20.12f %20.12f\n",
         iter, E_elec_new, E_total, delta_E);
  iter++;

  /*!
   * Start the SCF iterative procedure
   */
  while (iter < max_iter) {
    build_fock(F, D, H, ERI);
    if (iter == 1) {
      cout << "First iteration Fock matrix:" << endl; print_arma_mat(F);
    }
    F_prime = Symm_Orthog.t() * F * Symm_Orthog;
    arma::eig_sym(eps_vec, C_prime, F_prime);
    C = Symm_Orthog * C_prime;
    D_old = D;
    make_density(D, C, NOcc);
    E_elec_old = E_elec_new;
    E_elec_new = calc_elec_energy(D, H, F);
    E_total = E_elec_new + Vnn;
    if (iter == 1)
      printf("%4d %20.12f %20.12f %20.12f\n",
             iter, E_elec_new, E_total, delta_E);
    else
      printf("%4d %20.12f %20.12f %20.12f %20.12f\n",
             iter, E_elec_new, E_total, delta_E, rmsd_D);
    delta_E = E_elec_new - E_elec_old;
    rmsd_D = rmsd_density(D, D_old);
    if (delta_E < thresh_E && rmsd_D < thresh_D) {
      printf("Convergence achieved.\n");
      break;
    }
    F = F_prime;
    iter++;
  };

  arma::mat F_MO = C.t() * F * C;
  printf("Fock matrix in MO basis:\n"); print_arma_mat(F_MO);

  // Save the TEIs and MO coefficients/energies to disk
  // for use in other routines.
  ERI.save("ERI.mat", arma::arma_ascii);
  C.save("C.mat", arma::arma_ascii);
  F_MO.save("F_MO.mat", arma::arma_ascii);

  return 0;
}

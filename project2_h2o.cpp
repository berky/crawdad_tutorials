#include <cstdio>
#include <cmath>
#include "molecule.hpp"
#include "utils.hpp"
#include "constants.hpp"

int main()
{
  /** Step 1: Read the Coordinate Data
   *
   */
  FILE *xyzfile;
  xyzfile = fopen("h2o_geom.txt", "r");
  int natom;
  fscanf(xyzfile, "%d", &natom);
  Molecule mol(natom, 0);
  for (int i = 0; i < natom; i++)
    fscanf(xyzfile, "%lf %lf %lf %lf", &mol.zvals[i], &mol.geom[i][0], &mol.geom[i][1], &mol.geom[i][2]);
  fclose(xyzfile);

  /** Step 2: Read the Cartesian Hessian Data
   *
   */
  FILE *hessfile;
  hessfile = fopen("h2o_hessian.txt", "r");
  int hessatom;
  fscanf(hessfile, "%d", &hessatom);
  if (fabs(natom-hessatom) > 0) {
    printf("The number of atoms doesn't match.\n");
    return -1;
  }
  double **H = new double* [3*natom];
  for (int i = 0; i < 3*natom; i++)
    H[i] = new double[3*natom];
  for (int i = 0; i < 3*natom; i++)
    for (int j = 0; j < natom; j++)
      fscanf(hessfile, "%lf %lf %lf", &H[i][3*j], &H[i][3*j+1], &H[i][3*j+2]);
  fclose(hessfile);

  printf("Hessian:\n");
  print_mat(H, 3*natom, 3*natom);
  printf("\n");

  /** Step 3: Mass-Weight the Hessian Matrix
   * Divide each element of the Hessian matrix by the product of square roots of the masses of the atoms associated with the given coordinates:
   *  \vect{F}_{M}^{ij} = \frac{F_{ij}}{\sqrt{m_{i}m_{j}}}
   */
  double **Hmw = new double* [3*natom];
  for (int i = 0; i < 3*natom; i++)
    Hmw[i] = new double[3*natom];
  double mi, mj, mimj;
  for (int i = 0; i < natom; i++) {
    for (int j = 0; j < natom; j++) {
      mi = masses[(int)mol.zvals[i]]; mj = masses[(int)mol.zvals[j]];
      mimj = sqrt(mi*mj);
      Hmw[i*natom+0][j*natom+0] = H[i*natom+0][j*natom+0]/mimj;
      Hmw[i*natom+0][j*natom+1] = H[i*natom+0][j*natom+1]/mimj;
      Hmw[i*natom+0][j*natom+2] = H[i*natom+0][j*natom+2]/mimj;
      Hmw[i*natom+1][j*natom+0] = H[i*natom+1][j*natom+0]/mimj;
      Hmw[i*natom+1][j*natom+1] = H[i*natom+1][j*natom+1]/mimj;
      Hmw[i*natom+1][j*natom+2] = H[i*natom+1][j*natom+2]/mimj;
      Hmw[i*natom+2][j*natom+0] = H[i*natom+2][j*natom+0]/mimj;
      Hmw[i*natom+2][j*natom+1] = H[i*natom+2][j*natom+1]/mimj;
      Hmw[i*natom+2][j*natom+2] = H[i*natom+2][j*natom+2]/mimj;
    }
  }

  printf("Mass-weighted Hessian:\n");
  print_mat(Hmw, 3*natom, 3*natom);
  printf("\n");

  /** Step 4: Diagonalize the Mass-Weighted Hessian Matrix
   * Compute the eigenvalues of the mass-weighted Hessian:
   *  \vect{F}^{M}\vect{L} = \vect{L}\vect{\Lambda}
   */
  double *evals = new double[3*natom];
  for (int i = 0; i < 3*natom; i++) evals[i] = 0.0;
  double **evecs = new double* [3*natom];
  for (int i = 0; i < 3*natom; i++) evecs[i] = new double[3*natom];
  diag(3*natom, 3*natom, Hmw, evals, false, evecs, 1e-19);
  for (int i = 0; i < 3*natom; i++) delete[] evecs[i];
  delete[] evecs;

  printf("Mass-weighted Hessian eigenvalues:\n");
  for (int i = 0; i < 3*natom; i++)
    printf("%12.10f\n", evals[i]);
  printf("\n");

  /** Step 5: Compute the Harmonic Vibrational Frequencies
   * The vibrational frequencies are proportional to the square root of the eigenvalues of the mass-weighted Hessian:
   *  \omega_{i} = \textrm{constant} \times \sqrt{\lambda_{i}}
   */

  printf("Harmonic vibrational frequences [cm]^-1:\n");
  for (int i = 0; i < 3*natom; i++)
    printf("%10.4f\n", sqrt(evals[i])*vib_constant);
  printf("\n");

  /// Clean up after ourselves...
  for (int i = 0; i < 3*natom; i++) {
    delete[] H[i];
    delete[] Hmw[i];
  }
  delete[] H;
  delete[] Hmw;
  delete[] evals;

  return 0;
}



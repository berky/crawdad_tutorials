#include <cstdio>
#include <cmath>
#include "../molecule.hpp"

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
  for (int i = 0; i < mol.natom; i++)
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
    return -1
  }
  double **H = new double* [3*natom];
  for (int i = 0; i < 3*natom; i++)
    new double[3*natom];
  for (int i = 0; i < 3*natom; i++)
    for (int j = 0; j < natom; j++)
      fscanf(hessfile, "%lf %lf %lf", &H[i][3*j], &H[i][3*j+1], &H[i][3*j+2]);
  fclose(hessfile);

  /** Step 3: Mass-Weight the Hessian Matrix
   * Divide each element of the Hessian matri by the product of square roots of the masses of the atoms associated with the given coordinates:
   *  \vect{F}_{M}^{ij} = \frac{F_{ij}}{\sqrt{m_{i}m_{j}}}
   */
  

  /** Step 4: Diagonalize the Mass-Weighted Hessian Matrix
   * Compute the eigenvalues of the mass-weighted Hessian:
   *  \vect{F}^{M}\vect{L} = \vect{L}\vect{\Lambda}
   */
  

  /** Step 5: Compute the Harmonic Vibrational Frequencies
   * The vibrational frequencies are proportional to the square root of the eigenvalues of the mass-weighted Hessian:
   *  \omega_{i} = \textrm{constant} \times \sqrt{\lambda_{i}}
   */

  return 0;
}



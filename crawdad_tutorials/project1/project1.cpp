#include <cstdio>
#include "molecule.hpp"

int main()
{
  /// Step 1: Read the Coordinate Data from Input
  FILE *xyzfile;
  xyzfile = fopen("acetaldehyde.txt", "r");
  int natom;
  fscanf(xyzfile, "%d", &natom);
  Molecule mol(natom, 0);
  for (int i = 0; i < mol.natom; i++)
    fscanf(xyzfile, "%lf %lf %lf %lf", &mol.zvals[i], &mol.geom[i][0], &mol.geom[i][1], &mol.geom[i][2]);
  fclose(xyzfile);

  printf("\nNumber of atoms: %d\n", mol.natom);
  printf("Input Cartesian coordinates:\n");
  mol.print_geom();

  /** Step 2: Bond Lengths
   * Calculate the interatomic distances using the expression:
   * R_{ij} = \sqrt{(x_{i}-x_{j})^2 + (y_{i}-y_{j})^2 + (z_{i}-z_{j})^2}
   * where x, y, and z are Cartesian coordinates and i and j denote atomic indices.
   */
  mol.print_bonds();

  /** Step 3: Bond Angles
   * \mathrm{cos }\phi_{ijk} = \vect{\tilde{e}_{ji}} \cdot \vect{\tilde{e}_{jk}}
   * where
   *  \vect{\tilde{e}_{ij}} = (e_{ij}^{x}, e_{ij}^{y}, e_{ij}^{z}
   *  e_{ij}^{x} = -(x_{i} - x_{j})/R_{ij}
   *  e_{ij}^{y} = -(y_{i} - y_{j})/R_{ij}
   *  e_{ij}^{z} = -(z_{i} - z_{j})/R_{ij}
   */
  mol.print_angles();

  /** Step 4: Out-of-Plane Angles
   * Calculate all possible out-of-plane angles. For example, the angle 
   * \theta_{ijkl} for atom \mathbf{i} out of the plane containing atoms
   * \mathbf{j-k-l} (with \mathbf{k} as the central atom, connected to
   * \mathbf{i}) is given by:
   *  \mathrm{sin }\theta_{ijkl} = \frac{\vect{\tilde{e}_{kj}}\times\vect{\tilde{e}_{kl}}}{\mathrm{sin }\phi_{jkl}} \cdot \vect{\tilde{e}_{ki}}
   */
  mol.print_oop_angles();

  /** Step 5: Torsion/Dihedral Angles
   * Calculate all possible torsional angles. For examples, the torsional angle
   * \tau_{ijkl} for the atom connectivity \mathbf{i-j-k-l} is given by:
   *  \mathrm{cos }\tau_{ijkl} = \frac{(\vect{\tilde{e}_{ij}}\times\vect{\tilde{e}_{jk}}) \cdot (\vect{\tilde{e}_{jk}}\times\vect{\tilde{e}_{kl}})}{\mathrm{sin }\phi_{ijk}\mathrm{sin }\phi_{jkl}}
   * Can you also determine the sign of the torsional angle?
   */
  mol.print_torsion_angles();

  /** Step 6: Center-of-Mass Translation
   * Find the center of mass of the molecule:
   *  X_{c.m.} = \frac{\sum_{i} m_{i} x_{i}}{\sum_{i} m_{i}} etc.
   * where m_{i} is the mass of atom i and the summation runs over all atoms in the molecule.
   * Translate the input coordinates of the molecule to the center-of-mass.
   */
  mol.print_com();

  /** Step 7: Principal Moments of Inertia
   * Calculate elements of the moment of inertia tensor.
   * Diagonal:
   *  I_{xx} = \sum_{i} m_{i} (y_{i}^{2} + z_{i}^{2})
   *  I_{yy} = \sum_{i} m_{i} (z_{i}^{2} + x_{i}^{2})
   *  I_{zz} = \sum_{i} m_{i} (x_{i}^{2} + y_{i}^{2})
   * Off-diagonal:
   *  I_{xy} = \sum_{i} m_{i}x_{i}y_{i} = I_{yx}
   *  I_{xz} = \sum_{i} m_{i}x_{i}z_{i} = I_{zx}
   *  I_{yz} = \sum_{i} m_{i}y_{i}z_{i} = I_{zy}
   * Diagonalize the inertia tensor to obtain the principal moments of inertia:
   *  I_{a} \leq I_{b} \leq I_{c}
   * Report the moments of inertial in amu bohr^{2}, amu \AA^{2}, and g cm^{2}.
   * Based on the relative values of the principal moments, determine the molecular rotor type: linear, oblate, prolate, asymmetric.
   */

  /** Step 8: Rotational Constants
   * Compute the rotational constants in cm^{-1} and MHz:
   *  A = \frac{h}{8\pi^{2}cI_{a}}
   *  B = \frac{h}{8\pi^{2}cI_{b}}
   *  C = \frac{h}{8\pi^{2}cI_{c}}
   *  A \geq B \geq C
   */

  return 0;
}

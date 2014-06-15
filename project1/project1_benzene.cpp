#include <cstdio>
#include "../molecule.hpp"

int main()
{
  FILE *xyzfile;
  xyzfile = fopen("benzene.txt", "r");
  int natom;
  fscanf(xyzfile, "%d", &natom);
  Molecule mol(natom, 0);
  for (int i = 0; i < mol.natom; i++)
    fscanf(xyzfile, "%lf %lf %lf %lf", &mol.zvals[i], &mol.geom[i][0], &mol.geom[i][1], &mol.geom[i][2]);
  fclose(xyzfile);
  
  printf("\nNumber of atoms: %d\n", mol.natom);
  printf("Input Cartesian coordinates:\n");
  mol.print_geom();

  mol.print_bonds();
  mol.print_angles();
  mol.print_oop_angles();
  mol.print_torsion_angles();
  mol.print_com();
  mol.print_moi();
  mol.print_rot_const();

  return 0;
}

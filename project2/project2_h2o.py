from __future__ import division

import math
import numpy as np
import scipy.linalg as spl

from crawdad_tutorials.utils import *

"""
Step 0: Read the Coordinate Data (piratechem)
"""
import piratechem as pc
h2o = pc.molecule.Molecule("water")
pc.read.xyz("h2o_geom.xyz", h2o)
natom = h2o.num_atoms()

"""
Step 1: Read the Coordinate Data
"""

"""
Step 2: Read the Cartesian Hessian Data
"""
hessfile = open("h2o_hessian.txt", "r")
hessatom = int(hessfile.readline())
H = np.zeros((3*natom, 3*natom))
for i in range(0,3*natom):
    for j in range(0,natom):
        line = hessfile.readline().split()
        H[i][3*j] = float(line[0])
        H[i][3*j+1] = float(line[1])
        H[i][3*j+2] = float(line[2])

print("Hessian:")
print_mat(H)
print("")

"""
Step 3: Mass-Weight the Hessian Matrix
"""
Hmw = np.zeros((3*natom, 3*natom))
for i in range(0,natom):
    for j in range(0,natom):
        mi = h2o.atoms[i].m
        mj = h2o.atoms[j].m
        mimj = math.sqrt(mi*mj)
        Hmw[i*natom+0][j*natom+0] = H[i*natom+0][j*natom+0]/mimj
        Hmw[i*natom+0][j*natom+1] = H[i*natom+0][j*natom+1]/mimj
        Hmw[i*natom+0][j*natom+2] = H[i*natom+0][j*natom+2]/mimj
        Hmw[i*natom+1][j*natom+0] = H[i*natom+1][j*natom+0]/mimj
        Hmw[i*natom+1][j*natom+1] = H[i*natom+1][j*natom+1]/mimj
        Hmw[i*natom+1][j*natom+2] = H[i*natom+1][j*natom+2]/mimj
        Hmw[i*natom+2][j*natom+0] = H[i*natom+2][j*natom+0]/mimj
        Hmw[i*natom+2][j*natom+1] = H[i*natom+2][j*natom+1]/mimj
        Hmw[i*natom+2][j*natom+2] = H[i*natom+2][j*natom+2]/mimj

print("Mass-weighted Hessian:")
print_mat(Hmw)
print("")

"""
Step 4: Diagonalize the Mass-Weighted Hessian Matrix
"""
omega = spl.eigvalsh(Hmw)

print("Mass-weighted Hessian eigenvalues:")
for i in omega:
    print("{0:12.10f}".format(i))
print("")

"""
Step 5: Compute the Harmonic Vibrational Frequencies
"""
freq = np.sqrt(omega) * vib_constant

print("Harmonic vibrational frequencies [cm]^-1:")
for i in freq:
    print("{0:10.4f}".format(i))
print("")

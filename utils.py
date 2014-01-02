#!/usr/bin/env python2

from __future__ import division
import math

amu2kg = 1.6605402e-27
bohr2m = 0.529177249e-10
hartree2joule = 4.35974434e-18
planck = 6.6260755e-34
pi = math.acos(-1.0)
planckbar = planck/(2*pi)
speed_of_light = 299792458
avogadro = 6.0221413e+23
rot_constant = planck/(8*pi*pi*speed_of_light)
vib_constant = math.sqrt((avogadro*hartree2joule*1000)/(bohr2m*bohr2m))/(2*pi*speed_of_light*100)

def print_mat(mat):
    """Pretty-print a general NumPy matrix in a traditional format."""
    dim_rows, dim_cols = mat.shape
    # first, handle the column labels
    print(" " * 5),
    for i in range(dim_cols):
        print("{0:11d}".format(i+1)),
    print("")
    # then, handle the row labels
    for i in range(dim_rows):
        print("{0:5d}".format(i+1)),
        # print the matrix data
        for j in range(dim_cols):
            print("{0:11.7f}".format(mat[i][j])),
        print("")

def compound_index_2(i,j):
    if i > j:
        return i*(i+1)/2 + j
    else:
        return j*(j+1)/2 + i

def compound_index_4(i,j,k,l):
    ij = compound_index_2(i,j)
    kl = compound_index_2(k,l)
    compound_index_2(ij,kl)

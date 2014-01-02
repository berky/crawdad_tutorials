#!/usr/bin/env python2

from __future__ import division

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

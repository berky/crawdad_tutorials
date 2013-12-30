#!/usr/bin/env python2

def print_mat(mat):
    """Pretty-print a general NumPy matrix in a traditional format."""
    dim_rows, dim_cols = mat.shape
    width = 10
    # first, handle the column labels
    print(" " * width),
    for i in range(dim_cols):
        print("{0:10d}".format(i+1)),
    print("")
    # then, handle the row labels
    for i in range(dim_rows):
        print("{0:10d}".format(i+1)),
        # print the matrix data
        for j in range(dim_cols):
            print("{0:10f}".format(mat[i][j])),
        print("")

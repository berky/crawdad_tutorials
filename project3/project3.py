from __future__ import division

import numpy as np
import scipy as sp
import scipy.linalg as spl

from crawdad_tutorials.utils import *

def read_vnn():
    vnn_handle = open(stub + "enuc.dat")
    vnn = float(vnn_handle.read())
    vnn_handle.close()
    return vnn

def read_s():
    s_handle = open(stub + "s.dat")
    s_file = s_handle.readlines()
    s_handle.close()
    matsize = int(s_file[-1].split()[0])
    s = np.zeros(shape=(matsize, matsize))
    for line in s_file:
        tmp1, tmp2, tmp3 = line.split()
        mu, nu, s_mu_nu = int(tmp1)-1, int(tmp2)-1, float(tmp3)
        s[mu][nu] = s[nu][mu] = s_mu_nu
    return s

def read_t():
    t_handle = open(stub + "t.dat")
    t_file = t_handle.readlines()
    t_handle.close()
    matsize = int(t_file[-1].split()[0])
    t = np.zeros(shape=(matsize, matsize))
    for line in t_file:
        tmp1, tmp2, tmp3 = line.split()
        mu, nu, t_mu_nu = int(tmp1)-1, int(tmp2)-1, float(tmp3)
        t[mu][nu] = t[nu][mu] = t_mu_nu
    return t

def read_v():
    v_handle = open(stub + "v.dat")
    v_file = v_handle.readlines()
    v_handle.close()
    matsize = int(v_file[-1].split()[0])
    v = np.zeros(shape=(matsize, matsize))
    for line in v_file:
        tmp1, tmp2, tmp3 = line.split()
        mu, nu, v_mu_nu = int(tmp1)-1, int(tmp2)-1, float(tmp3)
        v[mu][nu] = v[nu][mu] = v_mu_nu
    return v

def read_eri():
    eri_handle = open(stub + "eri.dat")
    eri_file = eri_handle.readlines()
    eri_handle.close()
    matsize = int(eri_file[-1].split()[0])
    # be very inefficient with how we store these for now -- use all
    # four indices
    eri = np.zeros(shape=(matsize, matsize, matsize, matsize))
    for line in eri_file:
        tmp1, tmp2, tmp3, tmp4, tmp5 = line.split()
        mu = int(tmp1)-1
        nu = int(tmp2)-1
        lm = int(tmp3)-1
        sg = int(tmp4)-1
        eri_val = float(tmp5)
        eri[mu][nu][lm][sg] \
            = eri[mu][nu][sg][lm] \
            = eri[nu][mu][lm][sg] \
            = eri[nu][mu][sg][lm] \
            = eri[lm][sg][mu][nu] \
            = eri[lm][sg][nu][mu] \
            = eri[sg][lm][mu][nu] \
            = eri[sg][lm][nu][mu] \
            = eri_val
    return eri

def calc_elec_energy(D, H, F):
    return np.sum(D*(H+F))

def build_density(D, C):
    for mu in range(NBasis):
        for nu in range(NBasis):
            for m in range(NOcc):
                D[mu][nu] += C[mu][m] * C[nu][m]    

def build_fock(F, D, H, ERI):
    NBasis = len(D)
    for mu in range(NBasis):
        for nu in range(NBasis):
            F[mu][nu] = H[mu][nu]
            for lm in range(NBasis):
                for sg in range(NBasis):
                    F[mu][nu] += D[mu][nu] * (2*ERI[mu, nu, lm, sg] -
                                              ERI[mu, lm, nu, sg])

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--stub", dest="stub", type=str, default="h2o_sto3g")
    args = parser.parse_args()
    stub = args.stub + "_"

    ##########################################################################
    # Step #1: Nuclear Repulsion Energy

    vnn = read_vnn()

    print("Nuclear Repulsion Energy: {0:12f}".format(vnn))

    ##########################################################################
    # Step #2: One-Electron Integrals
    
    S = read_s()
    T = read_t()
    V = read_v()

    print("AO Overlap Integrals [S]:")
    print_mat(S)
    print("AO Kinetic Energy Integrals [T]:")
    print_mat(T)
    print("AO Nuclear Attraction Integrals [V]:")
    print_mat(V)
    
    H = T + V
    print("AO Core Hamiltonian [H]:")
    print_mat(H)

    NElec = 10
    NOcc = NElec//2
    NBasis = S.shape[0]

    ##########################################################################
    # Step #3: Two-Electron Integrals

    ERI = read_eri()

    ##########################################################################
    # Step #4: Build the Orthogonalization Matrix

    Lam_S, L_S = spl.eigh(S)
    Lam_S = Lam_S * np.eye(len(Lam_S))
    
    print("matrix of eigenvectors (columns) [L_S]:")
    print_mat(L_S)
    print("diagonal matrix of corresponding eigenvalues [Lam_S]:")
    print_mat(Lam_S)
    
    Lam_sqrt_inv = np.sqrt(spl.inv(Lam_S))
    Symm_Orthog = np.dot(L_S, np.dot(Lam_sqrt_inv, L_S.T))
    
    print("Symmetric Orthogonalization Matrix [Symm_Orthog]:")
    print_mat(Symm_Orthog)
    
    ##########################################################################
    # Step #5: Build the (Inital) Guess Density

    F_prime = np.dot(Symm_Orthog.T, np.dot(H, Symm_Orthog))

    print("Initial (guess) Fock Matrix [F_prime]:")
    print_mat(F_prime)

    # Diagonalize the Fock matrix

    eps, C_prime = spl.eigh(F_prime)
    eps = eps * np.eye(len(eps))

    print("Initial MO Coefficients [C_prime]:")
    print_mat(C_prime)
    print("Initial Orbital Energies [eps]:")
    print_mat(eps)

    # Transform the eigenvectors into the original (non-orthogonal) AO basis

    C = np.dot(Symm_Orthog, C_prime)

    print("Initial MO Coefficients (non-orthogonal) [C]:")
    print_mat(C)

    # Build the density matrix using the occupied MOs

    D = np.zeros((NBasis, NBasis))
    build_density(D, C)

    print("Initial Density Matrix [D]:")
    print_mat(D)

    ##########################################################################
    # Step #6: Compute the Initial SCF Energy

    E_elec_new = calc_elec_energy(D, H, H)
    E_total = E_elec_new + vnn

    print("Initial Electronic Energy: {0:20.12f}".format(E_elec_new))
    print("     Initial Total Energy: {0:20.12f}".format(E_total))

    ##########################################################################
    # Step #7: Compute the New Fock Matrix

    F = np.empty((NBasis,NBasis))

    iter = 1
    max_iter = 100

    while iter < max_iter:
        build_fock(F, D, H, ERI)
        if iter == 1:
            print("First iteration Fock matrix [F]:")
            print_mat(F)
        iter += 1

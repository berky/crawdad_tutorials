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

def read_s_ao():
    s_ao_handle = open(stub + "s.dat")
    s_ao_file = s_ao_handle.readlines()
    s_ao_handle.close()
    matsize = int(s_ao_file[-1].split()[0])
    s_ao = np.zeros(shape=(matsize, matsize))
    for line in s_ao_file:
        tmp1, tmp2, tmp3 = line.split()
        mu, nu, s_mu_nu = int(tmp1)-1, int(tmp2)-1, float(tmp3)
        s_ao[mu][nu] = s_ao[nu][mu] = s_mu_nu
    return s_ao

def read_t_ao():
    t_ao_handle = open(stub + "t.dat")
    t_ao_file = t_ao_handle.readlines()
    t_ao_handle.close()
    matsize = int(t_ao_file[-1].split()[0])
    t_ao = np.zeros(shape=(matsize, matsize))
    for line in t_ao_file:
        tmp1, tmp2, tmp3 = line.split()
        mu, nu, t_mu_nu = int(tmp1)-1, int(tmp2)-1, float(tmp3)
        t_ao[mu][nu] = t_ao[nu][mu] = t_mu_nu
    return t_ao

def read_v_ao():
    v_ao_handle = open(stub + "v.dat")
    v_ao_file = v_ao_handle.readlines()
    v_ao_handle.close()
    matsize = int(v_ao_file[-1].split()[0])
    v_ao = np.zeros(shape=(matsize, matsize))
    for line in v_ao_file:
        tmp1, tmp2, tmp3 = line.split()
        mu, nu, v_mu_nu = int(tmp1)-1, int(tmp2)-1, float(tmp3)
        v_ao[mu][nu] = v_ao[nu][mu] = v_mu_nu
    return v_ao

def read_eri_ao():
    eri_ao_handle = open(stub + "eri.dat")
    eri_ao_file = eri_ao_handle.readlines()
    eri_ao_handle.close()
    matsize = int(eri_ao_file[-1].split()[0])
    # be very inefficient with how we store these for now -- use all
    # four indices
    eri_ao = np.zeros(shape=(matsize, matsize, matsize, matsize))
    for line in eri_ao_file:
        tmp1, tmp2, tmp3, tmp4, tmp5 = line.split()
        mu = int(tmp1)-1
        nu = int(tmp2)-1
        lm = int(tmp3)-1
        sg = int(tmp4)-1
        eri_val = float(tmp5)
        eri_ao[mu][nu][lm][sg] \
            = eri_ao[mu][nu][sg][lm] \
            = eri_ao[nu][mu][lm][sg] \
            = eri_ao[nu][mu][sg][lm] \
            = eri_ao[lm][sg][mu][nu] \
            = eri_ao[lm][sg][nu][mu] \
            = eri_ao[sg][lm][mu][nu] \
            = eri_ao[sg][lm][nu][mu] \
            = eri_val
    return eri_ao

def calc_elec_energy(P, H, F):
    return np.sum(P*(H+F))

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
    
    S_AO = read_s_ao()
    T_AO = read_t_ao()
    V_AO = read_v_ao()

    print("AO Overlap Integrals [S_AO]:")
    print_mat(S_AO)
    print("AO Kinetic Energy Integrals [T_AO]:")
    print_mat(T_AO)
    print("AO Nuclear Attraction Integrals [V_AO]:")
    print_mat(V_AO)
    
    H_AO_Core = T_AO + V_AO
    print("AO Core Hamiltonian [H_AO_Core]:")
    print_mat(H_AO_Core)

    NElec = 10
    NOcc = NElec//2
    NBasis = S_AO.shape[0]

    ##########################################################################
    # Step #3: Two-Electron Integrals

    ERI_AO = read_eri_ao()

    ##########################################################################
    # Step #4: Build the Orthogonalization Matrix

    Lam_S_AO, L_S_AO = spl.eigh(S_AO)
    Lam_S_AO = Lam_S_AO * np.eye(len(Lam_S_AO))
    
    print("matrix of eigenvectors (columns) [L_S_AO]:")
    print_mat(L_S_AO)
    print("diagonal matrix of corresponding eigenvalues [Lam_S_AO]:")
    print_mat(Lam_S_AO)
    
    Lam_sqrt_inv_AO = np.sqrt(spl.inv(Lam_S_AO))
    Symm_Orthog = np.dot(L_S_AO, np.dot(Lam_sqrt_inv_AO, L_S_AO.T))
    
    print("Symmetric Orthogonalization Matrix [S^-1/2]:")
    print_mat(Symm_Orthog)
    
    ##########################################################################
    # Step #5: Build the (Inital) Guess Density

    F_prime_0_AO = np.dot(Symm_Orthog.T, np.dot(H_AO_Core, Symm_Orthog))

    print("Initial (guess) Fock Matrix [F_prime_0_AO]:")
    print_mat(F_prime_0_AO)

    # Diagonalize the Fock matrix

    eps_0_AO, C_prime_0_AO = spl.eigh(F_prime_0_AO)
    eps_0_AO = eps_0_AO * np.eye(len(eps_0_AO))

    print("Initial MO Coefficients [C_prime_0_AO]:")
    print_mat(C_prime_0_AO)
    print("Initial Orbital Energies [eps_0_AO]:")
    print_mat(eps_0_AO)

    # Transform the eigenvectors into the original (non-orthogonal) AO basis

    C_0_AO = np.dot(Symm_Orthog, C_prime_0_AO)

    print("Initial MO Coefficients (non-orthogonal) [C_0_AO]:")
    print_mat(C_0_AO)

    # Build the density matrix using the occupied MOs

    D_0 = np.zeros((NBasis,NBasis))

    for mu in range(NBasis):
        for nu in range(NBasis):
            for m in range(NOcc):
                D_0[mu][nu] += C_0_AO[mu][m] * C_0_AO[nu][m]

    print("Initial Density Matrix [D_0]:")
    print_mat(D_0)

    ##########################################################################
    # Step #6: Compute the Initial SCF Energy

    E_elec_0 = calc_elec_energy(D_0, H_AO_Core, H_AO_Core)
    E_total_0 = E_elec_0 + vnn

    print("Initial Electronic Energy: {0:20.12f}".format(E_elec_0))
    print("     Initial Total Energy: {0:20.12f}".format(E_total_0))

    ##########################################################################
    # Step #7: Compute the New Fock Matrix

    F_AO = np.empty((NBasis,NBasis))
    build_fock(F_AO, D_0, H_AO_Core, ERI_AO)

    print("First iteration Fock matrix:")
    print_mat(F_AO)
    

#!/usr/bin/env python2

from __future__ import division
import numpy as np
import scipy as sp
import scipy.linalg as spl
from utils import *

def read_vnn():
    vnn_handle = open(stub + "enuc.dat", "r")
    vnn = float(vnn_handle.read())
    return vnn

def read_s_ao():
    s_ao_handle = open(stub + "s.dat", "r")
    s_ao_file = s_ao_handle.readlines()
    matsize = int(s_ao_file[-1].split()[0])
    s_ao = np.zeros(shape=(matsize, matsize))
    for line in s_ao_file:
        tmp1, tmp2, tmp3 = line.split()
        mu, nu, s_mu_nu = int(tmp1)-1, int(tmp2)-1, float(tmp3)
        s_ao[mu][nu] = s_ao[nu][mu] = s_mu_nu
    return s_ao

def read_t_ao():
    t_ao_handle = open(stub + "t.dat", "r")
    t_ao_file = t_ao_handle.readlines()
    matsize = int(t_ao_file[-1].split()[0])
    t_ao = np.zeros(shape=(matsize, matsize))
    for line in t_ao_file:
        tmp1, tmp2, tmp3 = line.split()
        mu, nu, t_mu_nu = int(tmp1)-1, int(tmp2)-1, float(tmp3)
        t_ao[mu][nu] = t_ao[nu][mu] = t_mu_nu
    return t_ao

def read_v_ao():
    v_ao_handle = open(stub + "v.dat", "r")
    v_ao_file = v_ao_handle.readlines()
    matsize = int(v_ao_file[-1].split()[0])
    v_ao = np.zeros(shape=(matsize, matsize))
    for line in v_ao_file:
        tmp1, tmp2, tmp3 = line.split()
        mu, nu, v_mu_nu = int(tmp1)-1, int(tmp2)-1, float(tmp3)
        v_ao[mu][nu] = v_ao[nu][mu] = v_mu_nu
    return v_ao

def calc_elec_energy(P, H, F):
    return np.sum(P*(H+F))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--stub", dest="stub", type=str, default="h2o_sto3g")
    args = parser.parse_args()
    stub = args.stub + "_"

    ##########################################################################
    # Step #1: Nuclear Repulsion Energy

    Vnn = read_vnn()
    print("Nuclear Repulsion Energy: {0:12f}".format(Vnn))

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

    def read_eri_ao():
        eri_ao_handle = open(stub + "eri.dat", "r")
        eri_ao_file = eri_ao_handle.readlines()
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
    # Form an initial (guess) Fock matrix in the orthonormal AO basis using the core Hamiltonian as a guess:
    #  {\mathbf F}'_0 \equiv {\mathbf {\tilde S}}^{-1/2} {\mathbf H}^{\rm core} {\mathbf S}^{-1/2}

    F_prime_0_AO = np.dot(Symm_Orthog.T, np.dot(H_AO_Core, Symm_Orthog))

    print("Initial (guess) Fock Matrix [F_prime_0_AO]:")
    print_mat(F_prime_0_AO)

    # Diagonalize the Fock matrix:
    #  {\mathbf F}'_0 {\mathbf C}'_0 = {\mathbf C}'_0 \epsilon_0.
    # Note that the \epsilon_{0} matrix contains the initial orbital energies.

    eps_0_AO, C_prime_0_AO = spl.eigh(F_prime_0_AO)
    eps_0_AO = eps_0_AO * np.eye(len(eps_0_AO))

    print("Initial MO Coefficients [C_prime_0_AO]:")
    print_mat(C_prime_0_AO)
    print("Initial Orbital Energies [eps_0_AO]:")
    print_mat(eps_0_AO)

    # Transform the eigenvectors into the original (non-orthogonal) AO basis:
    #  C_{0} = \mathbf{S}^{1/2}\mathbf{C}_{0}^{'}

    C_0_AO = np.dot(Symm_Orthog, C_prime_0_AO)

    print("Initial MO Coefficients (non-orthogonal) [C_0_AO]:")
    print_mat(C_0_AO)

    # Build the density matrix using the occupied MOs:
    #  D_{\mu\nu}^{0} = \sum_{m}^{occ} (\mathbf{C}_{0})_{\mu}^{m} (\mathbf{C}_{0})_{\nu}^{m}
    # where m indexes the columns of the coefficient matrices, and the summation includes only the occupied spatial MOs.

    D_0 = np.zeros((NBasis,NBasis))

    for mu in range(NBasis):
        for nu in range(NBasis):
            for m in range(NOcc):
                D_0[mu][nu] += C_0_AO[mu][m] * C_0_AO[nu][m]
    # D_0 = C_0_AO[:,:NOcc] * C_0_AO[:,:NOcc]

    print("Initial Density Matrix [D_0]:")
    print_mat(D_0)

    ##########################################################################

    # Step #6: Compute the Initial SCF Energy
    # The SCF electronic energy may be computed using the density matrix as:
    #  E_{elec}^{0} = \sum_{\mu\nu}^{AO} D_{\mu\nu}^{0} (H_{\mu\nu}^{core} + F_{\mu\nu})
    # The total energy is the sum of the electronic energy and the nuclear repulsion energy:
    #  E_{total}^{0} = E_{elec}^{0} + E_{nuc}
    # where 0 denotes the initial SCF iteration.

    E_elec_0 = calc_elec_energy(D_0, H_AO_Core, F_prime_0_AO)
    E_total_0 = E_elec_0 + Vnn

    print("Initial Electronic Energy: {0:20.12f}".format(E_elec_0))
    print("     Initial Total Energy: {0:20.12f}".format(E_total_0))

    ##########################################################################

    # Step #7: Compute the New Fock Matrix
    # Start the SCF iterative procedure by building a new Fock matrix using the 
    #  previous iteration's density as:
    #   F_{\mu\nu} = H_{\mu\nu}^{core} + \sum_{\lambda\sigma}^{AO}
    #    D_{\lambda\sigma}^{i-1} [2(\mu\nu|\lambda\sigma) - (\mu\lambda|\nu\sigma)]
    #  where the double summation runs over all AOs and i-1 denotes the density
    #   for the last iteration.

    F_AO = np.empty((NBasis,NBasis))
    for i in range(NBasis):
        for j in range(NBasis):
            F_AO[i][j] = H_AO_Core[i][j]
            for k in range(NBasis):
                for l in range(NBasis):
                    F_AO[i][j] += D_0[k][l] * (2*ERI_AO[i][j][k][l] -
                                                     ERI_AO[i][k][j][l])
    
    print("First iteration Fock matrix:")
    print_mat(F_AO)
    

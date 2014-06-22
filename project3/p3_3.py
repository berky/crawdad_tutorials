# -*- coding: utf-8 -*-
"""
Created on Tue Jun 17 17:54:55 2014

@author: kevin
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 14:04:36 2014

@author: kevin
"""

import numpy as np

K = 7
N = 10

enucfile = open('enuc.dat', 'r')
e_nuc = float(enucfile.read())
enucfile.close()

svtfile = [open('s.dat', 'r'), open('v.dat', 'r'), open('t.dat', 'r')]
svtraw = [svtfile[i].readlines() for i in range(3)]
for i in range(3):
    svtfile[i].close

erifile = open('eri.dat', 'r')
eriraw = erifile.readlines()
erifile.close

eri = np.zeros((K, K, K, K))

svth = np.zeros((4, K, K))
s, v, t, h = svth[0], svth[1], svth[2], svth[3]

for idx in range(3):
    for line in svtraw[idx]:
        l = line.split()
        i, j = int(l[0]) - 1, int(l[1]) - 1
        svth[idx][i][j] = svth[idx][j][i] = l[2]
    
for i in range(K):
    for j in range(i+1):
        h[i][j] = h[j][i] = t[i][j] + v[i][j]
        
for line in eriraw:
    ln = line.split()
    [i, j, k, l] = [int(ln[m]) - 1 for m in range(4)]
    eri[i][j][k][l] = eri[j][i][k][l] = eri[i][j][l][k] = eri[j][i][l][k] =\
    eri[l][k][j][i] = eri[l][k][i][j] = eri[k][l][j][i] = eri[k][l][i][j] =\
    ln[4]

#ls, Ls = np.linalg.eig(s)
ls, Ls = np.linalg.eigh(s)
ldiag_s = np.zeros((K,K))
ldiag_invsqrt_s = np.zeros((K,K))

for i in range(K):
    ldiag_s[i][i] = ls[i]
    ldiag_invsqrt_s[i][i] = np.sqrt(1 / ls[i])

s_invsqrt = np.dot(np.dot(Ls, ldiag_invsqrt_s), Ls.transpose())

def addz(x, dim = 2):
    if dim == 2:
        return np.append(x, np.zeros((1, K, K)), 0)
    if dim == 1:
        return np.append(x, np.zeros((1, K)), 0)
    if dim == 0:
        return np.append(x, np.zeros((1, 1)), 0)

f = np.zeros((1,K,K))
fprime = np.zeros((1,K,K))
j_coul = np.zeros((1,K,K))
k_exch = np.zeros((1,K,K))
#dens = np.identity(K).reshape(1,K,K)
dens = np.zeros((1,K,K))
#dens = np.ones((1,K,K))
cprime = np.zeros((1,K,K))
c = np.zeros((1,K,K))
e_i = np.zeros((1,K))
e_elec = np.zeros((1,1))


tol_log10 = -14
i_scf = 0
tol = 10 ** tol_log10
d_dens_rms = 1

goal = 5
str_out = ['{:>5} {:>25} {:>25} {:>25}\n'.format('iter', 'E_elec', 'E_tot', u'D_rms') ]
bench = ['{:>12} {:>5} {:>25}\n'.format(u'-log[D_rms]', 'iter', 'E_elec')]

while d_dens_rms > tol:
    
    if i_scf > 0:
        f = addz(f)
        fprime = addz(fprime)
        j_coul = addz(j_coul)
        k_exch = addz(k_exch)
        dens = addz(dens)
        cprime = addz(cprime)
        c = addz(c)
        e_i = addz(e_i, 1)
        e_elec = addz(e_elec, 0)
    
    j_coul[i_scf] = np.einsum('kl,ijkl', dens[i_scf-1], eri)
    k_exch[i_scf] = np.einsum('jl,ijkl', dens[i_scf-1], eri)    
    f[i_scf] = h + 2 * (j_coul[i_scf] - 0.5 * k_exch[i_scf]) 
    
    fprime[i_scf] = np.dot(np.dot(s_invsqrt.transpose(), f[i_scf]), s_invsqrt)

    e_i[i_scf], cprime[i_scf] = np.linalg.eigh(fprime[i_scf])
    
    c[i_scf] = np.dot(s_invsqrt, cprime[i_scf])
        
    for i in range(K):
        for j in range(i+1):
            dens[i_scf][i][j] = dens[i_scf][j][i] = \
            sum([c[i_scf][i][m] * c[i_scf][j][m] for m in range(int(N/2))])

    
    e_elec[i_scf] = np.einsum('ij,ij', dens[i_scf], h + f[i_scf])
    
    if i_scf > 0:
        d_dens = dens[i_scf] - dens[i_scf - 1]
    else:
        d_dens = dens[i_scf]
    
    d_dens_rms = np.sqrt(np.einsum('ij,ij', d_dens, d_dens))
    
    str_out.append('{:>5d} {:>25.14f} {:>25.14f} {:>25.14f} \n'.format(i_scf, e_elec[i_scf][0], e_elec[i_scf][0] + e_nuc, d_dens_rms))

    if int(abs(np.log10(d_dens_rms))) > goal:
        bench.append('{:>12d} {:>5d} {:>25.14f} \n'.format(goal+1, i_scf, e_elec[i_scf][0]))         
        goal = int(abs(np.log10(d_dens_rms)))
        
    i_scf += 1
    

output = file('SCF_out.txt', 'w')
output.seek(0)
output.writelines(bench)
output.writelines(str_out)
output.close()    

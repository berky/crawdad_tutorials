function print_mat(mat)
    nrow, ncol = size(mat)
    # first, handle the column labels
    @printf("     ")
    for i = 1:ncol
        @printf("%12d", i)
    end
    @printf("\n")
    # then, handle the row labels
    for i = 1:nrow
        @printf("%5d", i)
        # print the matrix data
        for j = 1:ncol
            @printf("%12.7f", mat[i, j])
        end
        @printf("\n")
    end
end

function print_vec(vec, ncol)
    # print ncol elements of vec before moving to the next row
    counter = 0
    veclen = length(vec)
    while counter < veclen
        # do stuff!
        println("Hey! Listen!")
    end
end

function idx2(i, j)
    if i > j
        idx = div(i*(i+1), 2) + j
    else
        idx = div(j*(j+1), 2) + i
    end
    return idx
end

function idx4(i, j, k, l)
    ij = idx2(i, j)
    kl = idx2(k, l)
    return idx2(ij, kl)
end

function calc_elec_energy(P, H, F)
    return trace(P * (H + F))
end

function make_density(P, C, NOcc)
    P = C[:, 1:NOcc] * transpose(C[:, 1:NOcc])
end

function build_fock(F, P, H, ERI)
    r, c = size(H)
    # F = zeros(r, c)
    for mu = 1:r
        for nu = 1:c
            F[mu,nu] = H[nu]
            for lm = 1:r
                for sg = 1:c
                    F[mu,nu] += P[lm,sg] * (2*ERI[] - ERI[])
                end
            end
        end
    end
    return F
end


function rmsd_density(P_new, P_old)
    return sqrt(sum((P_new - P_old)^2))
end

function mix_density(P_new, P_old, alpha)
    # return ((1.0-alpha)*P_new) + (alpha*P_old)
    P_new = ((1.0-alpha)*P_new) + (alpha*P_old)
end

function build_error_matrix(F, D, S)
    return (F*D*S) - (S*D*F)
end

function build_B_matrix(e)
    NErr = length(e)
    B = Array(Int64, (NErr + 1, NErr + 1))
    B[NErr + 1, NErr + 1] = 0.0
    for a = 1:NErr
        B[a, NErr + 1] = B[NErr + 1, a] = -1.0
        for b = 1:a
            B[a, b] = B[b, a] = dot(transpose(e[a]), e[b])
        end
    end
    return B
end

function build_extrap_fock(F_extrap, diis_coeffs, diis_fock_vec)
    
end

function build_diis_zero_vec(len)
    diis_zero_vec = zeros(len)
    diis_zero_vec[len] = -1.0
    return diis_zero_vec
end

function main()
    NElec = 10
    NOcc = NElec/2
    NBasis = 7
    M = idx4(NBasis, NBasis, NBasis, NBasis)

    Vnn = readdlm("h2o_sto3g_enuc.dat")[1]
    S = Array(Float64, (NBasis, NBasis))
    T = Array(Float64, (NBasis, NBasis))
    V = Array(Float64, (NBasis, NBasis))
    open("h2o_sto3g_s.dat") do file_s
        for line in eachline(file_s)
            sline = split(line)
            i, j = int(sline[1]), int(sline[2])
            S[i, j] = S[j, i] = float(sline[3])
        end
    end
    open("h2o_sto3g_t.dat") do file_t
        for line in eachline(file_t)
            sline = split(line)
            i, j = int(sline[1]), int(sline[2])
            T[i, j] = T[j, i] = float(sline[3])
        end
    end
    open("h2o_sto3g_v.dat") do file_v
        for line in eachline(file_v)
            sline = split(line)
            i, j = int(sline[1]), int(sline[2])
            V[i, j] = V[j, i] = float(sline[3])
        end
    end
    H = T + V

    println("AO Overlap Integrals:"); print_mat(S);
    println("AO Kinetic Energy Integrals:"); print_mat(T);
    println("AO Nuclear Attraction Integrals:"); print_mat(V);
    println("AO Core Hamiltonian:"); print_mat(H);

    ERI = zeros(M)

    open("h2o_sto3g_eri.dat") do file_eri
        for line in eachline(file_eri)
            sline = split(line)
            mu, nu, lm, sg = map(int, sline[1:4])
            ERI[idx4(mu, nu, lm, sg)] = float(sline[5])
        end
    end

    thresh_E = 1.0e-15
    thresh_D = 1.0e-7
    iter = 0
    max_iter = 1024
    # E_total, E_elec_old, E_elec_new, delta_E, rmsd_D

    S_eigendecomp = eigfact(S)
    Lam_S = diagm(S_eigendecomp.values)
    L_S = S_eigendecomp.vectors
    Lam_sqrt_inv = sqrt(inv(Lam_S))
    Symm_Orthog = L_S * Lam_sqrt_inv * transpose(L_S)
    F_prime = transpose(Symm_Orthog) * H * Symm_Orthog
    F_prime_eigendecomp = eigfact(F_prime)
    epsilon = F_prime_eigendecomp.values
    C_prime = F_prime_eigendecomp.vectors
    C = Symm_Orthog * C_prime
    println("matrix of eigenvectors (columns)"); print_mat(L_S)
    println("diagonal matrix of corresponding eigenvalues"); print_mat(Lam_S)
    println("Symmetric Orthogonalization Matrix"); print_mat(Symm_Orthog)
    println("Initial (guess) Fock Matrix"); print_mat(F_prime)
    println("Initial MO Coefficients"); print_mat(C_prime)
    println("Initial Orbital Energies"); print_mat(diagm(epsilon))
    println("Initial MO Coefficients (non-orthogonal)"); print_mat(C)
    D = Array(Float64, (NBasis, NBasis))
    D = make_density(D, C, NOcc)
    println("Initial Density Matrix"); print_mat(D)

    E_elec_new = calc_elec_energy(D, H, H)
    E_total = E_elec_new + Vnn
    delta_E = E_total
    @printf("%4d %20.12f %20.12f %20.12f\n",
            iter, E_elec_new, E_total, delta_E)

    F = Array(Float64, (NBasis, NBasis))

    # The main loop of the SCF iterative procedure.
    while iter < max_iter
        F = build_fock(F, D, H, ERI)
        if iter == 1
            @printf("First iteration Fock matrix:\n"); print_mat(F)
        end
        F_prime = transpose(Symm_Orthog) * F * Symm_Orthog
        F_prime_eigendecomp = eigfact(F_prime)
        epsilon = F_prime_eigendecomp.values
        C_prime = F_prime_eigendecomp.vectors
        C = Symm_Orthog * C_prime
        D_old = D
        D = make_density(D, C, NOcc)
        E_elec_old = E_elec_new
        E_elec_new = calc_elec_energy(D, H, F)
        E_total = E_elec_new + Vnn
        delta_E = E_elec_new - E_elec_old
        rmsd_D = rmsd_density(D, D_old)
        if iter == 1
            @printf("%4d %20.12f %20.12f %20.12f\n",
                    iter, E_elec_new, E_total, delta_E)
        else
            @printf("%4d %20.12f %20.12f %20.12f %20.12f\n",
                    iter, E_elec_new, E_total, delta_E, rmsd_D)
        end
        if (delta_E < thresh_E && rmsd_D < thresh_D)
            @printf("Convergence achieved.\n")
            break
        end
        F = F_prime
        iter += 1
    end # while iter < max_iter
end

main()

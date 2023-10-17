using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf
using ActiveSpaceSolvers

C = npzread("mo_coeffs.npy")
h0 = npzread("ints_h0.npy")
h1 = npzread("ints_h1.npy")
h2 = npzread("ints_h2.npy")
ints = InCoreInts(h0, h1, h2)

Pa = npzread("Pa.npy")
Pb = npzread("Pb.npy")
@printf(" Input energy:    %12.8f\n", compute_energy(ints, RDM1(Pa, Pb)))

init_fspace =  [(5, 0), (3, 3), (5, 0), (3, 3), (3, 3), (3, 3), (3, 3), (3, 3), (3, 3)]
clusters    =  [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [11, 12, 13, 14, 15, 16], [17, 18, 19, 20, 21, 22, 23, 24, 25, 26], [27, 28, 29, 30, 31, 32], [33, 34, 35, 36, 37, 38], [39, 40, 41, 42, 43, 44], [45, 46, 47, 48, 49, 50], [51, 52, 53, 54, 55, 56], [57, 58, 59, 60, 61, 62]]


clusters = [MOCluster(i, collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)

d1 = RDM1(n_orb(ints))


ansatze = [FCIAnsatz(10, 5, 0), 
           FCIAnsatz(6,  3, 3),
           FCIAnsatz(10, 0, 5),
           FCIAnsatz(6,  3, 3),
           FCIAnsatz(6,  3, 3),
           FCIAnsatz(6,  3, 3),
           FCIAnsatz(6,  3, 3),
           FCIAnsatz(6,  3, 3),
           FCIAnsatz(6,  3, 3)]

e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, ansatze, d1,
                           maxiter_oo   = 100, 
                           maxiter_ci   = 200, 
                           maxiter_d1   = 200,
                           orb_hessian  = true,
                           verbose      = 0, 
                           tol_oo       = 1e-7, 
                           tol_d1       = 1e-9, 
                           tol_ci       = 1e-11, 
                           sequential   = false, 
                           diis_start   = 1,
                           max_ss_size  = 24)





ints = orbital_rotation(ints, U)
C = C*U

npzwrite("Ccmf.npy", C)

@save "data_cmf.jld2" clusters init_fspace ints d1 e_cmf U 

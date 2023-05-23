using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf

C = npzread("/home/arnabbachhar/FermiCG-data/bimetallics/cr2_morokuma/26__3d4s4d_2ps_3d4s4d/mo_coeffs.npy")
h0 = npzread("/home/arnabbachhar/FermiCG-data/bimetallics/cr2_morokuma/26__3d4s4d_2ps_3d4s4d/ints_h0.npy")
h1 = npzread("/home/arnabbachhar/FermiCG-data/bimetallics/cr2_morokuma/26__3d4s4d_2ps_3d4s4d/ints_h1.npy")
h2 = npzread("/home/arnabbachhar/FermiCG-data/bimetallics/cr2_morokuma/26__3d4s4d_2ps_3d4s4d/ints_h2.npy")
ints = InCoreInts(h0, h1, h2)

Pa = npzread("/home/arnabbachhar/FermiCG-data/bimetallics/cr2_morokuma/26__3d4s4d_2ps_3d4s4d/Pa.npy")
Pb = npzread("/home/arnabbachhar/FermiCG-data/bimetallics/cr2_morokuma/26__3d4s4d_2ps_3d4s4d/Pb.npy")
@printf(" Input energy:    %12.8f\n", compute_energy(ints, RDM1(Pa, Pb)))


init_fspace=[(6, 3), (4, 4), (6, 3)]
 clusters   =  [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], [12, 13, 14, 15], [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]]



clusters = [MOCluster(i, collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)

d1 = RDM1(n_orb(ints))


# # Do CMF
e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace, RDM1(n_orb(ints)),
    verbose=0, sequential=false, max_iter_oo=20)

ints = orbital_rotation(ints, U)
C = C*U# Do CMF


e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, d1,
                           maxiter_oo   = 700, 
                           maxiter_ci   = 200, 
                           maxiter_d1   = 200, 
                           verbose      = 0, 
                           tol_oo       = 1e-6, 
                           tol_d1       = 1e-9, 
                           tol_ci       = 1e-11, 
                           sequential   = true, 
                           alpha        = .1,
                           diis_start   = 1,
                           max_ss_size  = 36)

ints = orbital_rotation(ints, U)
C = C*U

npzwrite("Ccmf_26_cr2.npy", C)

@save "data_cmf_26_cr2.jld2" clusters init_fspace ints d1 e_cmf U 

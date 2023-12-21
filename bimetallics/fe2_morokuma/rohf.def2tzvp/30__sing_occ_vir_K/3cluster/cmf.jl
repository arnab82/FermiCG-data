using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf

C = npzread("../mo_coeffs.npy")
h0 = npzread("../ints_h0.npy")
h1 = npzread("../ints_h1.npy")
h2 = npzread("../ints_h2.npy")
ints = InCoreInts(h0, h1, h2)

Pa = npzread("../Pa.npy")
Pb = npzread("../Pb.npy")
@printf(" Input energy:    %12.8f\n", compute_energy(ints, RDM1(Pa, Pb)))

#clusters =  [[1], [2], [3], [4], [5, 6, 7, 8, 9, 10, 11, 12, 13, 14], [15, 16, 17], [18, 19, 20, 21, 22, 23, 24, 25, 26, 27], [28], [29], [30]]
#init_fspace =  [(1, 1), (1, 1), (1, 1), (1, 1), (5, 0), (3, 3), (0, 5), (1, 1), (1, 1), (1, 1)]

#clusters    = [[2, 3, 4], [5, 6, 7, 8, 9, 10, 11, 12, 13, 14], [1, 15, 16, 17], [18, 19, 20, 21, 22, 23, 24, 25, 26, 27], [28, 29, 30]]
#init_fspace = [(3, 3), (5, 0), (4, 4), (0, 5), (3, 3)]

clusters    = [[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], [1, 15, 16, 17], [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]]
init_fspace = [(8, 3), (4, 4), (3, 8)]


clusters = [MOCluster(i, collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)

d1 = RDM1(n_orb(ints))

ansatze=[FCIAnsatz(length(ci),init_fspace[ci.idx][1],init_fspace[ci.idx][2]) for ci in clusters]

# # Do CMF
#e_cmf, U, d1 = ClusterMeanField.cmf_oo_newton(ints, clusters, init_fspace, ansatze, rdm1, maxiter_oo = 400,
#                                              tol_oo=1e-6, 
#                                              tol_d1=1e-9, 
#                                              tol_ci=1e-11,
#                                              verbose=4, 
#                                              zero_intra_rots = true,
#                                              sequential=true)
e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace, RDM1(n_orb(ints)),
                                       verbose=2, sequential=true, max_iter_oo=100)

ints = orbital_rotation(ints, U)
C = C*U

npzwrite("Ccmf.npy", C)

@save "data_cmf.jld2" clusters init_fspace ints d1 e_cmf U 

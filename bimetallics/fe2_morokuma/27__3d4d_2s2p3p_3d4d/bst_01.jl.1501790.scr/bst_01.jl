using QCBase
using ClusterMeanField
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2


@load "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/27__3d4d_2s2p3p_3d4d/data_cmf_27.jld2"

M = 300

init_fspace = FermiCG.FockConfig([(5,0), (4, 4), (0,5)])

cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3], init_fspace, max_roots=M, verbose=1);

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=6
# start by defining P/Q spaces
p_spaces = Vector{ClusterSubspace}()

for ci in clusters
    ssi = ClusterSubspace(clusters[ci.idx])

    num_states_in_p_space = 1
    for sec in FermiCG.possible_spin_focksectors(clusters, init_fspace)
        add_subspace!(ssi, sec[ci.idx], 1:num_states_in_p_space)
    end
    push!(p_spaces, ssi)
end

ci_vector = BSTstate(clusters, p_spaces, cluster_bases, R=nroots) 
ci_vector = FermiCG.add_spin_focksectors(ci_vector)

na = FermiCG.n_elec_a(init_fspace)
nb = FermiCG.n_elec_b(init_fspace)

FermiCG.fill_p_space!(ci_vector, na, nb)
FermiCG.eye!(ci_vector)
e_ci, v = FermiCG.ci_solve(ci_vector, cluster_ops, clustered_ham)

ept = FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=1e-8)

e_var, v_var = block_sparse_tucker(v, cluster_ops, clustered_ham,
                                   max_iter    = 20,
                                   nbody       = 4,
                                   H0          = "Hcmf",
                                   thresh_var  = 1e-1,
                                   thresh_foi  = 1e-6,
                                   thresh_pt   = 1e-3,
                                   ci_conv     = 1e-5,
                                   ci_max_iter = 100,
                                   do_pt       = true,
                                   resolve_ss  = false,
                                   tol_tucker  = 1e-4,
                                   solver      = "davidson")

#ept = FermiCG.compute_pt2_energy(v_var, cluster_ops, clustered_ham, thresh_foi=1e-6, prescreen=false, compress_twice=false)


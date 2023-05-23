using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2


@load "/home/arnab22/FermiCG-data/bimetallics/cr2_morokuma/13__3d_2p_3d/data_cmf_13_cr2.jld2"

M = 200 

init_fspace = FockConfig([(3, 0), (3, 3), (0, 3)])

cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3], init_fspace, max_roots=M, verbose=1);

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=6

# start by defining P/Q spaces
p_spaces = Vector{ClusterSubspace}()

ssi = ClusterSubspace(clusters[1])
add_subspace!(ssi, (3,0), 1:1)
add_subspace!(ssi, (2,1), 1:1)
add_subspace!(ssi, (1,2), 1:1)
add_subspace!(ssi, (0,3), 1:1)
push!(p_spaces, ssi)

ssi = ClusterSubspace(clusters[2])
add_subspace!(ssi, (3,3), 1:1)
push!(p_spaces, ssi)

ssi = ClusterSubspace(clusters[3])
add_subspace!(ssi, (3,0), 1:1)
add_subspace!(ssi, (2,1), 1:1)
add_subspace!(ssi, (1,2), 1:1)
add_subspace!(ssi, (0,3), 1:1)
push!(p_spaces, ssi)


ci_vector = BSTstate(clusters, p_spaces, cluster_bases, R=4) 

na = 6
nb = 6
FermiCG.fill_p_space!(ci_vector, na, nb)
FermiCG.eye!(ci_vector)
ebst, vbst = FermiCG.ci_solve(ci_vector, cluster_ops, clustered_ham)

ebst, vbst = FermiCG.block_sparse_tucker(vbst, cluster_ops, clustered_ham, thresh_var=1e-2, thresh_foi=1e-6, thresh_pt=1e-3)
e2_bst = FermiCG.compute_pt2_energy(vbst, cluster_ops, clustered_ham, thresh_foi=1e-6);
#e2_bst = FermiCG.compute_pt2_energy(vbst, cluster_ops, clustered_ham, thresh_foi=1e-7);

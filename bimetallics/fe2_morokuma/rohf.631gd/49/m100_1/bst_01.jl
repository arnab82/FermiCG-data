using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf

@load "../data_cmf.jld2"

M = 3 


init_fspace =  [(5, 0), (4, 4), (0, 5), (12, 12), (0, 0)]
init_fspace = FermiCG.FockConfig(init_fspace)
cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,1,1], init_fspace, max_roots=M, verbose=1);

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=6
# start by defining P/Q spaces
p_spaces = Vector{ClusterSubspace}()

    
ssi = ClusterSubspace(clusters[1])
add_subspace!(ssi, (5,0), 1:1)
add_subspace!(ssi, (4,1), 1:1)
add_subspace!(ssi, (3,2), 1:1)
add_subspace!(ssi, (2,3), 1:1)
add_subspace!(ssi, (1,4), 1:1)
add_subspace!(ssi, (0,5), 1:1)
push!(p_spaces, ssi)

ssi = ClusterSubspace(clusters[2])
add_subspace!(ssi, init_fspace[2], 1:1)
push!(p_spaces, ssi)

ssi = ClusterSubspace(clusters[3])
add_subspace!(ssi, (5,0), 1:1)
add_subspace!(ssi, (4,1), 1:1)
add_subspace!(ssi, (3,2), 1:1)
add_subspace!(ssi, (2,3), 1:1)
add_subspace!(ssi, (1,4), 1:1)
add_subspace!(ssi, (0,5), 1:1)
push!(p_spaces, ssi)

ssi = ClusterSubspace(clusters[4])
add_subspace!(ssi, init_fspace[4], 1:1)
push!(p_spaces, ssi)

ssi = ClusterSubspace(clusters[5])
add_subspace!(ssi, init_fspace[5], 1:1)
push!(p_spaces, ssi)


ci_vector = BSTstate(clusters, p_spaces, cluster_bases, R=6) 

na = sum([i[1] for i in init_fspace]) 
nb = sum([i[2] for i in init_fspace]) 

FermiCG.fill_p_space!(ci_vector, na, nb)
FermiCG.eye!(ci_vector)
e_ci, v = FermiCG.ci_solve(ci_vector, cluster_ops, clustered_ham);

FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=1e-8)

# v = FermiCG.compress(v, thresh=1e-3)
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham, thresh_var=1e-2, thresh_foi=1e-3);
@save "data_bst.jld2" clusters init_fspace ints cluster_bases v

FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=1e-8)

bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham, thresh_var=1e-2, thresh_foi=1e-4);
@save "data_bst.jld2" clusters init_fspace ints cluster_bases v

FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=1e-8)

bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham, thresh_var=.8e-3, thresh_foi=1e-4);
@save "data_bst.jld2" clusters init_fspace ints cluster_bases v

FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=1e-8)

bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham, thresh_var=.6e-3, thresh_foi=1e-4);
@save "data_bst.jld2" clusters init_fspace ints cluster_bases v

FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=1e-8)

bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham, thresh_var=.4e-3, thresh_foi=1e-4);
@save "data_bst.jld2" clusters init_fspace ints cluster_bases v

FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=1e-8)

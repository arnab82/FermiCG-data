using QCBase
using ClusterMeanField
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2


@load "data_cmf.jld2"

M = 40 

init_fspace = FockConfig([(4,4), (2,5), (4,4), (5,2)])


cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [2,2,2,2], init_fspace, max_roots=M, verbose=1);

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=4
ci_vector = FermiCG.TPSCIstate(clusters, init_fspace, R=nroots)

ci_vector = FermiCG.add_spin_focksectors(ci_vector)

eci, v = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);

e0a, v0a = FermiCG.tpsci_ci(ci_vector, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 1e-3, 
                            thresh_foi   = 1e-5);

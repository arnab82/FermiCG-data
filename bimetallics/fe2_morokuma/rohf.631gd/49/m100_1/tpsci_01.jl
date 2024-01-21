using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf

function print_pt2(e, nroots)
    @printf("*E(PT2): ")
    for i in e
        @printf("%14.12f ", i)
    end
    println()
end

@load "../data_cmf.jld2"

M = 100 


init_fspace =  [(5, 0), (4, 4), (0, 5), (12, 12), (0, 0)]
init_fspace = FermiCG.FockConfig(init_fspace)
cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,1,1], init_fspace, max_roots=M, verbose=1);

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=6


ci_vector = FermiCG.TPSCIstate(clusters, init_fspace, R=nroots)

ci_vector = FermiCG.add_spin_focksectors(ci_vector)

eci, v = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);


e0a, v0a = FermiCG.tpsci_ci(ci_vector, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 1e-2, 
                            thresh_spin  = 1e-3,
                            thresh_foi   = 1e-5,
                            max_mem_ci   = 60);
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)

e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 1e-3, 
                            thresh_spin  = 1e-3,
                            thresh_foi   = 1e-5,
                            max_mem_ci   = 60);
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)

e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 8e-4, 
                            thresh_spin  = 8e-4,
                            thresh_foi   = 1e-5,
                            max_mem_ci   = 60);
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)

e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 6e-4, 
                            thresh_spin  = 6e-4,
                            thresh_foi   = 1e-5,
                            max_mem_ci   = 60);
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)

e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 4e-4, 
                            thresh_spin  = 4e-4,
                            thresh_foi   = 1e-5,
                            max_mem_ci   = 60);
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)

e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 2e-4, 
                            thresh_spin  = 2e-4,
                            thresh_foi   = 1e-5,
                            max_mem_ci   = 60);
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)

e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 1e-4, 
                            thresh_spin  = 1e-4,
                            thresh_foi   = 1e-5,
                            max_mem_ci   = 60);
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)


@save "data_tpsci.jld2" clusters init_fspace ints cluster_bases v0a



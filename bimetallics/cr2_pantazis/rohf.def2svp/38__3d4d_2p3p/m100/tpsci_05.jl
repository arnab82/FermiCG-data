using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using LinearAlgebra
using MKL
using Printf

function run()
	@load "data_cmf.jld2"

	M = 100

	display(BLAS.get_config())

	init_fspace = FockConfig(init_fspace)
	cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,3,3], init_fspace, max_roots=M, verbose=1);

	clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
	cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

	FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

	nroots=4


	ci_vector = FermiCG.TPSCIstate(clusters, init_fspace, R=nroots)

	ci_vector = FermiCG.add_spin_focksectors(ci_vector)

	eci, v0 = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);

	values = [10e-4, 8e-4, 6e-4, 4e-4, 2e-4, 1e-4]
	
	v0a = deepcopy(v0)

	idx = 0

	for v in values

		e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
					    thresh_cipsi      = v,
					    thresh_spin       = v,
					    thresh_foi        = v/10,
					    ci_max_iter       = 150,
					    ci_lindep_thresh  = 1e-12,
					    conv_thresh       = 1e-6,
					    max_iter          = 20,
					    nbody             = 4,
					    max_mem_ci        = 700);
		print(" E(VAR): ")
		print(e0a)
		println()

		@save "out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace v 

		e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)

		print(" E(PT2): ")
		print(e2a)
		println()

		@save @sprintf("data_tpsci_%02i.jld2",idx) ci_vector v0 eci v0a e0a e2a v 

		idx += 1
	end
end

run()



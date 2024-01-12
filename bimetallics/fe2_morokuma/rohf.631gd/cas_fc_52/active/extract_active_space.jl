using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf


@load "../data_cmf.jld2"

#C = npzread("mo_coeffs.npy")
#h0 = npzread("ints_h0.npy")
#h1 = npzread("ints_h1.npy")
#h2 = npzread("ints_h2.npy")

#ints = InCoreInts(h0, h1, h2)

ecore = ints.h0

init_fspace =  [(5, 0), (3, 3), (0, 5), (3, 3), (3, 3), (3, 3), (3, 3), (3, 3), (3, 3)]
clusters    =  [[16, 17, 18, 19, 20], [21, 22, 23, 24, 25, 26], [27, 28, 29, 30, 31], [32, 33, 34, 35, 36, 37], [38, 39, 40, 41, 42, 43], [44, 45, 46, 47, 48, 49], [50, 51, 52, 53, 54, 55], [56, 57, 58, 59, 60, 61], [62, 63, 64, 65, 66, 67]]


active_cluster = []


for i in 1:length(clusters)
    for j in 1:length(clusters[i])
        push!(active_cluster, clusters[i][j])
        clusters[i][j] -= 15
    end
end

frozen_cluster = MOCluster(0, [i for i in 1:15])
active_cluster = MOCluster(0, active_cluster)

display(frozen_cluster)
display(active_cluster)

Nold = size(ints.h1,1)

dcore = zeros(Nold, Nold)
for i in frozen_cluster.orb_list
    dcore[i,i] = 1
end
dcore = RDM1(dcore, dcore)
d1 = subset(d1, active_cluster)
ints = subset(ints, active_cluster, dcore)

ecore += compute_energy(subset(ints, frozen_cluster), subset(d1, frozen_cluster))
ints = InCoreInts(ecore, ints.h1, ints.h2)
@printf(" h0          :    %12.8f\n", ints.h0)
@printf(" Core energy :    %12.8f\n", ecore)
@printf(" Input energy:    %12.8f\n", compute_energy(ints, d1))

clusters = [MOCluster(i, collect(clusters[i])) for i = 1:length(clusters)]

@save "data_cmf.jld2" clusters init_fspace ints d1 e_cmf U 

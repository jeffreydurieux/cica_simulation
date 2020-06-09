# Tue Jun  9 13:47:26 2020
# Author: Jeffrey Durieux, MSc

# checks for CICA simulation 2

# check percentage of removed voxels

library(TwoStepClustering3Way)

setwd("~/Repositories/cica_simulation/Simulation2/CreateMaps/")

S1 <- get(load("cluster1.Rdata"))
S2 <- get(load("cluster2.Rdata"))
S3 <- get(load("cluster3.Rdata"))

bmask <- get(load("bmaskslice.Rdata"))
idx <- which(bmask==1)

S1 <- S1[idx,]
S2 <- S2[idx,]
S3 <- S3[idx,]

### how many voxels deactivated?

# cluster 1 compared to cluster 2

((sum(S1[,1] != 0) - sum(S2[,1] != 0) ) / 4761 * 100 +
(sum(S1[,2] != 0) - sum(S2[,2] != 0) ) / 4761 * 100 +
(sum(S1[,3] != 0) - sum(S2[,3] != 0) ) / 4761 * 100 +
(sum(S1[,4] != 0) - sum(S2[,4] != 0) ) / 4761 * 100) / 4


# cluster 1 compared to cluster 3

((sum(S1[,1] != 0) - sum(S3[,1] != 0) ) / 4761 * 100 +
    (sum(S1[,2] != 0) - sum(S3[,2] != 0) ) / 4761 * 100 +
    (sum(S1[,3] != 0) - sum(S3[,3] != 0) ) / 4761 * 100 +
    (sum(S1[,4] != 0) - sum(S3[,4] != 0) ) / 4761 * 100) / 4


S <- list(S1,S2,S3)
RVmat <- computeRVmat(S, dist = F)
mean(RVmat[upper.tri(RVmat)])
RVmat

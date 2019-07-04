# Thu Jul  4 15:55:21 2019
# Author: Jeffrey Durieux, MSc

# What: R-script to merge all simulation results from computer cluster
# Output: a dataframe to do an ANOVA on

library(gtools)

setwd("~/Desktop/testd/")

files <- dir()
files <- mixedsort(files)
files

V <- c(500,2000)                # voxels
Q <- c(2,5,20)                  # components
R <- c(2,4)                     # clusters
D <- c('square', 'nonsquare')   # dimension mixing matrix
E <- c(0.05, 0.2, 0.4)          # noise level
rep <- 1:10
design <- expand.grid(V=V, Q = Q,
                      R = R, D = D, E=E, replication = rep)

load(files[i])
des <- data.frame(design[i,], files[i],ARI = ResultVector[1],
                  STuck = ResultVector[2], ATuck = ResultVector[3],
                  Time = ResultVector[4])
for(i in 2:720){
  load(files[i])
  des2 <- data.frame(design[i,], files[i],ARI = ResultVector[1],
  STuck = ResultVector[2], ATuck = ResultVector[3],
  Time = ResultVector[4])
  des <- rbind(des,des2)
}


save(des, file = "PATHNAME")
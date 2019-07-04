# Mon Jul  1 15:09:19 2019
# Author: Jeffrey Durieux, MSc
# Simulationscript for C-ICA simulation number 1

# R-version used:

### passed arguments:
# 1: replication
# 2: row (row of factorial simulation design)

args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)

replication <- args[1]
row <- args[2]

# This will generate a unique seed
seed <- as.numeric(paste(replication,0,row,sep = ""))
set.seed(seed)

library(ica)
library(multiway)
library(neuRosim)
library(mclust)

#source(file = "/home/durieuxj/CICACODE/R/CICA_datagen.R")
source(file = "/home/durieuxj/CICACODE/R/CICA_helpers.R")
source(file = "/home/durieuxj/CICACODE/R/CICA_main.R")
source(file = "/home/durieuxj/CICACODE/R/Simulate_CICA.R")




#between design
V <- c(500,2000)                # voxels
Q <- c(2,5,20)                  # components
R <- c(2,4)                     # clusters
D <- c('square', 'nonsquare')   # dimension mixing matrix
E <- c(0.05, 0.2, 0.4)          # noise level

design <- expand.grid(V=V, Q = Q,
                      R = R, D = D, E=E)


rowd <- design[row , ]

### Simulate CICA ###
dim <- ifelse(test = (rowd$D == 'square'), yes = rowd$Q , no = 100 )
SimData <- Simulate_CICA(N = 40, V = rowd$V, D = dim, Q = rowd$Q, R = rowd$R,E = rowd$E)

### Do CICA ###

time <- proc.time()
Result <- CICA(nStarts = 30, DataList = SimData$Xe, nComp = rowd$Q,
               nClus = rowd$R ,show = 'best',startclus = 'random', 
               icaf = 'normal')
Time <- proc.time()-time

##### compute result measures #####
# ARI
# Tucker Qr
# Tucker Air


### compute ARI ###
ARI <- adjustedRandIndex(SimData$P, Result$beststart$P)

### compute Tuckers S ###
SRest <- Result$beststart$Sr
SR <- SimData$Qr
SRtucker <- numeric()
for(i in 1:length(SRest)){
  x <- numeric()
  for(j in 1:length(SR)){
    x[j] <- mean (apply (abs ( congru(x = SR[[j]], y = SRest[[i]]) ), MARGIN = 2, FUN = max) )
  }
  SRtucker[i] <- max(x)
}
SRtucker <- mean(SRtucker)

### compute Tuckers A ####
# this is not an easy problem, here you need to split the Ai matrices and take into account the 
# cluster permutation


# function to split A matrices
Asplit <- function(Air, P, nTS){
  Unique_clus <- sort(unique(P))
  Al <- list()
  for(i in 1:length(Unique_clus)){
   obs <- which(P == Unique_clus[i])
    parts <- split(1:(nTS*length(obs)), ceiling(seq_along(1:(nTS*length(obs))) /nTS ))
    Al[[i]] <- lapply(parts, function(x) Air[[i]][,x])
  }
  return(Al)
}




A <- SimData$Air
AL <- unlist(A, recursive = F)
Aest <- Result$beststart$Air
Aest <- Asplit(Air = Aest, P = Result$beststart$P, nTS = dim )

if(rowd$R ==2 ){
  idx1 <- which(Result$beststart$P == 1)
  idx2 <- which(Result$beststart$P == 2)
  
  names(Aest[[1]]) <- idx1
  names(Aest[[2]]) <- idx2
}else{
  idx1 <- which(Result$beststart$P == 1)
  idx2 <- which(Result$beststart$P == 2)
  idx3 <- which(Result$beststart$P == 3)
  idx4 <- which(Result$beststart$P == 4)
  names(Aest[[1]]) <- idx1
  names(Aest[[2]]) <- idx2
  names(Aest[[3]]) <- idx3
  names(Aest[[4]]) <- idx4
}


AestL <- unlist(Aest, recursive = F)
AestL <- lapply(AestL,t)

order <- as.numeric(names(AestL))

### order the estimated matrices so that a one-to-one matching can be performed

OrderedEstimatedA <- vector(mode = 'list', length = 40)
for(i in 1:length(order)){
  OrderedEstimatedA[ order[i] ] <- AestL[i]
}

x <- numeric()
for(i in 1:length(AL)){
    x[i] <- mean ( apply (abs ( congru(x = AL[[i]], y = OrderedEstimatedA[[i]]) ), MARGIN = 2, FUN = max) )
}
AallTucker <- mean(x)

ResultVector <- c(ARI=ARI,SRTucker = SRtucker, AiTucker = AallTucker, Time = Time[3])


### save data ###
# what to save:
# 1 SimData + Results CICA
# 2 ResultVector

Bulk <- list()
Bulk$SimData <- SimData
Bulk$Results <- Result
Bulk$sessionInfo <- sessionInfo()

# set right path for output folder
setwd()

save(ResultVector, file = paste("rep",replication,"row",row,  "ResultVector.Rdata",sep = "_"))
save(Bulk, file = paste("rep",replication,"row",row,  "BulkData.Rdata",sep = "_"))





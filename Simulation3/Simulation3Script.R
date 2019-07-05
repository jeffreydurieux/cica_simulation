#### Model Selection Beckmann networks
# true model: 3 clusters 4 components
# range of clusters = 1, 2, 3, 4, 5
# range of comp     = 2, 3, 4, 5, 6


args <- commandArgs(TRUE)
args <- as.numeric(args)

clus <- args[1]
comp <- args[2]
start <- args[3]

set.seed(2407)

##### source CICA files
library(ica)
library(neuRosim)
library(mclust)
source("/home/durieuxj/CICACODE/R/CICA_main.R")
source("/home/durieuxj/CICACODE/R/CICA_helpers.R")
source("/home/durieuxj/CICACODE/R/CICA_datagen.R")

setwd("/exports/fsw/durieuxj/Simulations/CICA_SIM_2/Output/")

#### select data set here #####

X <- get(load("SubSampXE_SNR_1_Gaussian.Rdata"))

#### run CICA
#tmp <- proc.time()
cica <- CICA(nStarts = 30, DataList = X, nComp = comp, nClus = clus, show='best')
#time <- proc.time() - tmp


# output
ext1 <- paste("clus", clus, "_comp", comp,"start", start , "_resultSim3.Rdata",sep = "")
save(cica, file = ext1)


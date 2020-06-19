#### Model Selection Beckmann networks
# true model: 3 clusters 4 components
# range of clusters = 2, 3, 4, 5  # do group icas for K = 1
# range of comp     = 2, 3, 4, 5, 6, 7


# Fri Jun 19 14:33:52 2020
# Author: Jeffrey Durieux, MSc
# adjustments made for doing extra modsel analyses on ALICE cluster

# multiple CICAs on 10 simulated datasets of simulation 2
# take two conditions: .3 and .7

args <- commandArgs(TRUE)
#args <- as.numeric(args)

row <- args[1]

# grid with 480 rows
mygrid <- expand.grid(clus = 2:5, comp = 2:7, rep = 1:10, noise = c(30,70))


set.seed(row)

##### source CICA files
library(ica)
library(neuRosim)
library(mclust)
source("/home/durieuxj/CICACODE/R/CICA_main.R")
source("/home/durieuxj/CICACODE/R/CICA_helpers.R")
source("/home/durieuxj/CICACODE/R/CICA_datagen.R")

setwd("/data/durieuxj/")

#### select data set here #####

if(mygrid[row,]$noise == 30){
  datasetname <- paste('rep_', mygrid[row,]$rep,
                       'CICA_simdata_0.30.Rdata', sep = '')
}else{
  datasetname <- paste('rep_', mygrid[row,]$rep,
                       'CICA_simdata_0.70.Rdata', sep = '')
}

X <- get(load(datasetname))


#### run CICA
tmp <- proc.time()
cica <- CICA(nStarts = 30, DataList = X, nComp = mygrid[row,]$comp,
             nClus = mygrid[row,]$clus, show='best')
time <- proc.time() - tmp

# output
setwd("/data/durieuxj/")

ext1 <- paste('rep', mygrid[row,]$rep, 
              "_clus", mygrid[row,]$clus, 
              "_comp", mygrid[row,]$comp, 
              "_noise", mygrid[row,]$noise, 
              "_resultSim3.Rdata",sep = "")
save(cica, file = ext1)


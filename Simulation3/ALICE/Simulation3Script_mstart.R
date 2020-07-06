# Fri Jul  3 12:09:54 2020
# Author: Jeffrey Durieux, MSc

# redo failed files

#### Model Selection Beckmann networks
# true model: 3 clusters 4 components
# range of clusters = 2, 3, 4, 5  # do group icas for K = 1
# range of comp     = 2, 3, 4, 5, 6, 7


# Fri Jun 19 14:33:52 2020
# Author: Jeffrey Durieux, MSc
# adjustments made for doing extra modsel analyses on ALICE cluster

# multiple CICAs on 10 simulated datasets of simulation 2
# take two conditions: .3 and .7

# alice cluster info:
# will try to run it on the short partition.. Max time = 03:00:00
# maxjobsubmit is 100 so need to rewrite the code for --array=[1-100]


args <- commandArgs(TRUE)
#args <- as.numeric(args)

row <- args[1]


setwd("/data/durieuxj/")
load('Sim3GridModSelALICENotDoneYet.Rdata')

mygrid <- mygridnew

set.seed(mygrid$mstart[row])

##### source CICA files
library(ica)
library(neuRosim)
library(mclust)
source("/data/durieuxj/CICACODE/CICA_main.R")
source("/data/durieuxj/CICACODE/CICA_helpers.R")
source("/data/durieuxj/CICACODE/CICA_datagen.R")

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
#tmp <- proc.time()

if(mygrid[row,]$clus == 1){
  cica <- CICA(nStarts = 1, DataList = X, nComp = mygrid[row,]$comp,
               nClus = mygrid[row,]$clus, show='best')	
}else{
  
  cica <- CICA(nStarts = 1, DataList = X, nComp = mygrid[row,]$comp,
               nClus = mygrid[row,]$clus, show='best')
}

#time <- proc.time() - tmp

# output
setwd("/data/durieuxj/CICAsim3Results/mstarts")

ext1 <- paste('rep', mygrid[row,]$rep, 
              "_clus", mygrid[row,]$clus, 
              "_comp", mygrid[row,]$comp, 
              "_noise", mygrid[row,]$noise,
              '_rowid',row, 
              '_mstart',mygrid$mstart,
              "_resultSim3.Rdata",sep = "")
save(cica, file = ext1)

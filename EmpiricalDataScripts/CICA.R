### CICA P1.5 shark
### script used for empirical cica analysis

### note: did a Q = 5 analysis on new shark slurm schedular
#### some adjustments were made

args <- commandArgs(TRUE)
args <- as.numeric(args)

#clus <- args[1]
#comp <- args[2]
#start_quint <- args[3] 

##### slurm grid ######

grid <- expand.grid(clus = 2:5, comp = 5, start_quint = 1:30)
row <- args[1]
mygrid <- grid[row,]


set.seed(mygrid$start_quint)

# sources and libraries
library(ica)
library(neuRosim)
library(mclust)
source("/home/durieuxj/CICACODE/R/CICA_main.R")
source("/home/durieuxj/CICACODE/R/CICA_helpers.R")
source("/home/durieuxj/CICACODE/R/CICA_datagen.R")


# data
setwd("/exports/fsw/durieuxj/data_nifti/DataTVDownSizedTwiceBrainMasked/")
files <- dir()

dataL <- list()
for(i in 1:250){
  load(files[i])
  
  # subtract voxel means 
  cmean <- colMeans(data_tv)
  data <- sweep(data_tv, 2, cmean)
  
  # scale data such that each block has a ssq of 1000
  data <- xscale(data)

  dataL[[i]] <- data
}

# transpose data blocks
dataL <- lapply(dataL, t)

# do cica
tmp <- proc.time()
#cica <- CICA(nStarts = 1, DataList = dataL, nComp = comp, nClus = clus)
cica <- CICAStart(DataList = dataL, exSig = mygrid$comp, nClus = mygrid$clus, 
                  newclus = clusf(length(dataL), clus),icaf = "no")
time <- proc.time() - tmp

# output
setwd( paste("/exports/fsw/durieuxj/P1.5/CICAOUTPUT/clus", mygrid$clus,sep = "") )
ext1 <- paste("CICA_nc", mygrid$comp, "Batch", mygrid$start_quint ,".Rdata" , sep = "")
save(cica, file = ext1)
#ext2 <- paste("TimingCICA_nc",comp, "Batch", start_quint, ".Rdata" , sep = "")
#save(time, file = ext2)




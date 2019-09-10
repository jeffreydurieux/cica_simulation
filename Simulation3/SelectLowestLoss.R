# Tue Sep 10 13:47:40 2019 ------------------------------
# Author: Jeffrey Durieux, MSc

# This script selects the lowest loss cica object from a parallel 



setwd("~/Repositories/cica_simulation/Data/CICA_SIM_3/clus5/")

files <- dir()

#### select per component ####
select <- 'comp6'

id <- grep(select, files)

files <- files[id]

files

Loss <- numeric()

for(i in 1:length(files)){
  load(files[i])
  Loss[i] <- cica$beststart$LossHistory[cica$beststart$Iterations]
}

minloss <- which.min(Loss)

# file to select
files[minloss]

cica <- get(load(files[minloss]))
cica$beststart$P

# change ext name here 
save(cica, file = "../clus5_comp6_resultSim3.Rdata")

# Fri Jul  3 11:38:17 2020
# Author: Jeffrey Durieux, MSc

# check failed files cica mod sel sim

library(gtools)

setwd("~/Repositories/cica_simulation/Data/CICA_SIM_3/ALICE_results/")

files <- dir()

mygrid <- expand.grid(clus = 1:5, comp = 2:7, rep = 1:10, noise = c(30,70))

mygridnames <- vector(mode = 'character', length = 600)
for(row in 1:nrow(mygrid)){
  ext1 <- paste('rep', mygrid[row,]$rep, 
                "_clus", mygrid[row,]$clus, 
                "_comp", mygrid[row,]$comp, 
                "_noise", mygrid[row,]$noise,
                '_rowid',row, 
                "_resultSim3.Rdata",sep = "")
  mygridnames[row] <- ext1
}



files %in% mygridnames
id_not_done <- !(mygridnames %in% files)

mygrid[id_not_done,]

mygrid <- mygrid[id_not_done,]
mygridnew <- do.call(rbind, replicate(50, mygrid, simplify = F))
mygridnew <- cbind(mygridnew, mstart = rep(1:50, each = 50))
save(mygridnew, file = '../../Sim3GridModSelALICENotDoneYet.Rdata')

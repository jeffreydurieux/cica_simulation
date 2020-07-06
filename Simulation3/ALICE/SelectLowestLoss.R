# Sat Jul  4 12:55:39 2020
# Author: Jeffrey Durieux, MSc


# This script selects the lowest loss cica object from a parallel 
setwd("~/Repositories/cica_simulation/Data/CICA_SIM_3/ALICE_results/CICAmstarts/")

files <- dir()

#### select per rep of data frame ####

load('~/Repositories/cica_simulation/Data/Sim3GridModSelALICENotDoneYet.Rdata')
mygrid <- mygridnew[1:50,]

namesgrid <- character()
for(row in 1:50){
  ext1 <- paste('rep', mygrid[row,]$rep, 
                "_clus", mygrid[row,]$clus, 
                "_comp", mygrid[row,]$comp, 
                "_noise", mygrid[row,]$noise,
                '_rowid',row, 
                '_mstart',mygrid[row,]$mstart,
                "_resultSim3.Rdata",sep = "")
  
  namesgrid[row] <- ext1
  
}

namesgrid
namesgrid <- gsub('.{33}$','',namesgrid)
namesgrid


# check if all necessary C-ICA have at least 30 starts
# conclusion: yes
for(i in 1:length(namesgrid)){
  print(length(grep(namesgrid[i], files)))
}
  
for(i in 1:length(namesgrid)){
  id <- grep(namesgrid[i], files)
  ff <- files[id]
  
  LOSS <- numeric()
  for(ds in 1:length(ff)){
    load(ff[ds])
    LOSS[ds] <- cica$beststart$LossHistory[cica$beststart$Iterations]
  }
  
  minloss <- which.min(LOSS)
  ff[minloss]
  
  ext <- gsub('.{26}$','',ff[minloss])
  ext <- paste('../lowestlossmstarts/', ext, sep = '')
  ext <- paste(ext, 'resultSim3.Rdata', sep = '', collapse = '')
  
  # change ext name here 
  save(cica, file = ext)  
}





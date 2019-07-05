args <- commandArgs(TRUE)
args <- as.numeric(args)
rep <- args[1]


##################################################################
# What dataset:                                                  #
name <- paste('rep_',rep,"CICA_simdata_0.10.Rdata",sep = "")     #
#name <- paste('rep_',rep,"CICA_simdata_0.30.Rdata",sep = "")    #
#name <- paste('rep_',rep,"CICA_simdata_0.70.Rdata",sep = "")    #
#                                                                #
#name <- paste('rep_',rep,"AR_10_CICA_simdata.Rdata",sep = "")   #
#name <- paste('rep_',rep,"AR_30_CICA_simdata.Rdata",sep = "")   #
#name <- paste('rep_',rep,"AR_70_CICA_simdata.Rdata",sep = "")   #
##################################################################



set.seed(2407*rep)


source("/home/durieuxj/CICACODE/R/CICA_main.R")
source("/home/durieuxj/CICACODE/R/CICA_helpers.R")
source("/home/durieuxj/CICACODE/R/CICA_datagen.R")

setwd("/exports/fsw/durieuxj/Simulations/CICA_SIM_2/")

X <- get(load(name))

### data is already masked 
#bmask <- get(load("bmaskslice.Rdata"))
#idx <- which(bmask==1)

# apply mask to data

#Xm <- lapply(seq_along(X), function(x) X[[x]][ idx, ])

cica <- CICA(nStarts = 30, DataList = X, nComp = 4, nClus = 3, show = 'best')

setwd("/exports/fsw/durieuxj/Simulations/CICA_SIM_2/Output")

ext <- paste("result_of_",name,sep = "")
save(cica, file = ext)

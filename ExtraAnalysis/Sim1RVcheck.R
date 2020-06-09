# Tue Jun  9 10:32:14 2020
# Author: Jeffrey Durieux, MSc

# check RV of simulated data

library(TwoStepClustering3Way)
library(ica)

source("~/Repositories/cica_simulation/SharkCode/Simulate_CICA.R")


#between design
V <- c(500,2000)                # voxels
Q <- c(2,5,20)                  # components
R <- c(2,4)                     # clusters
D <- c('square', 'nonsquare')   # dimension mixing matrix
E <- c(0.05, 0.2, 0.4)          # noise level

design <- expand.grid(V=V, Q = Q,
                      R = R, D = D, E=E)


#### replication = 1:10
#### row = 1:72

RVlist <- list()
for(i in 1:10){
  RV <- numeric()
  
  time <- proc.time()
  for(j in 1:72){
    replication <- i
    row <- j
    
    # This will generate a unique seed
    seed <- as.numeric(paste(replication,0,row,sep = ""))
    set.seed(seed)
    
    rowd <- design[row , ]
    
    ### Simulate CICA ###
    dim <- ifelse(test = (rowd$D == 'square'), yes = rowd$Q , no = 100 )
    SimData <- Simulate_CICA(N = 40, V = rowd$V, D = dim, Q = rowd$Q, R = rowd$R,E = rowd$E) 
    
    ### check RV of Qr ###
    cat('Replication: ',replication,'Row of design: ',row,'\n')
    RVmat <- computeRVmat(SimData$Qr, dist = F, verbose = T)
    id <- upper.tri(RVmat)
    RV[j] <- mean(RVmat[id])
  }
  proc.time() - time
  RVlist[[i]] <- RV
}

### put everyting in a dataframe for further analyses ###
design <- expand.grid(V=V, Q = Q,
                      R = R, D = D, E=E, Replication = 1:10)

RVall <- unlist(RVlist)
design <- cbind(design, RV = RVall)

save(design, file = '~/Repositories/cica_simulation/Data/Sim1RV.Rdata')

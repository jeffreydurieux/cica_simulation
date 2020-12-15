# Mon Dec 14 12:58:34 2020
# Author: Jeffrey Durieux, MSc

# modified RV computations on cluster for simulation 3 C-ICA

args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)

replication <- args[1]

set.seed(replication)

library(utils)

load('/home/durieuxj/data/SingleICAsShats.Rdata')

modRV <- function(X, Y){

  if(nrow(X) != nrow(Y)){
    stop('Number of rows of input matrices are not equal')
  }

  XXtilde <- ( X %*% t(X) ) - diag (diag( X %*% t(X) ) )
  YYtilde <- ( Y %*% t(Y) ) - diag (diag( Y %*% t(Y) ) )

  res <-  ( t(c(XXtilde)) %*% c(YYtilde) ) /
    ( sqrt( ( t(c(XXtilde)) %*% c(XXtilde)) * ( t(c(YYtilde)) %*% c(YYtilde)) ) )


  return(res)
}

N <- 40
comb <- t(combn(1:N, 2))

splits <- split(1:780, ceiling(seq_along(1:780)/100))

comb <- comb[splits[[replication]],]

RVsS <- matrix(data = NA, nrow = N , ncol = N)
RVS <- numeric()

cat("Computing pairwise modRV statistics: \n")
pb <- txtProgressBar(min = 0, max = nrow(comb), initial = 0)

for(i in 1:nrow(comb)){

  RVS[i] <- modRV( Shats[[ comb[i,1] ]] , Shats[[ comb[i,2] ]])
  res <- c(comb[i , ] , RVS[i] )

  RVsS[res[1]  , res[2] ] <- res[3]

  setTxtProgressBar(pb, i)
}

RVsS[lower.tri(RVsS)] = t(RVsS)[lower.tri(RVsS)]
diag(RVsS) <- 1
cat('\n')

ext <- paste('/home/durieuxj/data/modRV_',replication,'.Rdata', sep = '')
save(RVsS, file = ext)


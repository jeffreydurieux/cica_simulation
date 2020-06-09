### CICA version 1.0
### Jeffrey Durieux
### 22-02-2017
### repo on github

### load libraries and helper functions
#devtools::use_package("ica")
#devtools::use_package("NMFN")
#devtools::use_package("neuRosim")
#devtools::use_package("expm")
#devtools::use_package("MatrixCorrelation")
#devtools::use_package("MASS")
#devtools::use_package("oro.nifti")
#devtools::use_package("fpc")
#devtools::use_package("gtools")
#devtools::use_package("gplots")
#devtools::use_package("lattice")

CICA <- function(nStarts, DataList, nComp, nClus, show="best", startclus="random", icaf = "normal"){

  if(startclus == "random"){
    if(show == "best"){
      LD <- length( DataList )
      rClus <- clusf(LD, nClus)
      beststart <- CICAStart(DataList, nComp, nClus, rClus,icaf)
      Loss <- numeric()
      floss<- beststart$LossHistory[beststart$Iterations]
      for(i in 2:nStarts){
        rClus <- clusf(LD, nClus)
        nextstart <- CICAStart(DataList, nComp, nClus, rClus, icaf)
        Loss[i] <- nextstart$LossHistory[nextstart$Iterations]
        n <- nextstart$LossHistory[nextstart$Iterations]
        b <- beststart$LossHistory[beststart$Iterations]

        if(n <= b){
          beststart <- nextstart
        }
      }#end for
      Loss[1] <- floss
      res <- list(beststart = beststart,Loss=Loss)
    }else{
      LD <- length( DataList )
      res <- list()
      for(i in 1:nStarts){
        #random clustering
        rClus <- clusf(LD, nClus)
        res[[i]] <- CICAStart(DataList, nComp, nClus, rClus,icaf)
      }# end for
    }# end else
  }# end if startclus random
  else{
    rClus <- startclus
    res <- CICAStart(DataList, nComp, nClus, rClus,icaf)
  }# end else rational start
  return(res)
}# end CICA (main function)



#CICAtest <- CICA(nStarts = 10,DataList = test$XE_List, nComp = 2, nClus = 2, show = "best")
# s1 <- icasamp("b", "rnd", 1000)
# s2 <- icasamp("b", "rnd", 1000)
# s3 <- icasamp("b", "rnd", 1000)
# s4 <- icasamp("b", "rnd", 1000)
# s5 <- icasamp("b", "rnd", 1000)
#
#
# x1 <- cbind(s1,s2) %*% matrix(runif(4,-1,1) , 2, 2)
# x2 <- cbind(s1,s2) %*% matrix(runif(6,-1,1) , 2, 3)
# x3 <- cbind(s1,s2) %*% matrix(runif(8,-1,1) , 2, 4)
# x4 <- cbind(s1,s2) %*% matrix(runif(4,-1,1) , 2, 2)
# x5 <- cbind(s1,s2) %*% matrix(runif(6,-1,1) , 2, 3)
# x6 <- cbind(s1,s2) %*% matrix(runif(8,-1,1) , 2, 4)
#
# x7 <- cbind(s3,s4) %*% matrix(runif(4,-1,1) , 2, 2)
# x8 <- cbind(s3,s4) %*% matrix(runif(6,-1,1) , 2, 3)
# x9 <- cbind(s3,s4) %*% matrix(runif(8,-1,1) , 2, 4)
# x10 <- cbind(s3,s4) %*% matrix(runif(4,-1,1) , 2, 2)
# x11 <- cbind(s3,s4) %*% matrix(runif(6,-1,1) , 2, 3)
# x12 <- cbind(s3,s4) %*% matrix(runif(8,-1,1) , 2, 4)
#
# Xdata <- list(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)
#
# set.seed(2407)
# test2 <- CICA(nStarts = 50,DataList = Xdata, nComp = 2,nClus = 2, show="bet")
# opt <- LowestLoss(test2)
# sum(test$Loss==min(test$Loss))

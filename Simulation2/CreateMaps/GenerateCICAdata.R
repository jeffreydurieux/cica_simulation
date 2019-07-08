#### generate cica data
library(CICA) #Gen_TimeCourses
library(plotly)
library(neuRosim)

####### first load spatial maps of the three clusters #######
setwd("/Volumes/LaCie/MyData/CICA/Project1/SIM2/PosterExample/")
S1 <- get(load("cluster1.Rdata"))
S2 <- get(load("cluster2.Rdata"))
S3 <- get(load("cluster3.Rdata"))

bmask <- get(load("bmaskslice.Rdata"))
idx <- which(bmask==1)

set.seed(2407)

####### simulate rs-timecourses ##########

#### loop this script 10 times from here ######

# N per cluster equals 10
for(rep in 1:10){
  
  A1 <- Gen_TimeCourses(N = 30, nTS = 4, nscan = 150)
  
  
  nameext <- paste('rep_', rep, sep = "") 
  
  save(A1, file = paste(nameext,'TimeCoursesA1.Rdata' ,sep = "") )
  
  X1 <- lapply(seq_along(A1), function(x) tcrossprod(S1,A1[[x]]))
  
  ##### remove outer brain ####
  X1 <- lapply(seq_along(X1), function(x) X1[[x]][ idx, ])

  
  addError2 <- function(datablock,error, type = "Gaussian", additiontype = 1)
  {
    if(type == "Gaussian"){
      errorM<-replicate(ncol(datablock),rnorm(nrow(datablock), sd = 1))
    }else if(type == "AR"){
      nscan <- ncol(datablock)
      vdim <- nrow(datablock)^(1/3)
      dim <- rep(vdim, 3)
      errorM <- neuRosim::temporalnoise(dim = dim, nscan = nscan, sigma = 1,rho = 0.2)
      errorM <- matrix(errorM, ncol = nscan)
    }
    
    errorM<-SSequal(errorM,datablock)
    errorlevel<-error/(1-error)
    
    if(additiontype == 1){
      res<-datablock + (errorM * sqrt(errorlevel))
    }else{
      res<-datablock + (errorM * errorlevel)
    }
    return(res)
  }
  
  
  ######### types of Gaussian noise #############
  
  ### .10
  XE <- lapply(X1, FUN = addError2, error = 0.1, additiontype=1)
  
  #X <- c(XE1,XE2,XE3)
  save(XE, file=paste(nameext,"CICA_simdata_0.10.Rdata",sep = ""))
  
  
  ### .30
  XE <- lapply(X1, FUN = addError2, error = 0.3, additiontype=1)
  
  
  #X <- c(XE1,XE2,XE3)
  save(XE, file=paste(nameext,"CICA_simdata_0.30.Rdata",sep = ""))
  
  ### .70
  XE <- lapply(X1, FUN = addError2, error = 0.7, additiontype=1)
  
  
  
  save(XE, file=paste(nameext,"CICA_simdata_0.70.Rdata",sep = ""))
  
  
  ###### add AR(1) noise ########
  
  nvoxels <- nrow(X1[[1]]); nTime <- ncol(X1[[1]])
  
  ColWiseNoise <- function(ncol, nrow, rho, sd){
    x <- replicate(ncol, arima.sim( n = nrow ,
                                    mean = 0, list(order = c(1,0,0),ar = rho),
                                    sd = sd ))
    return(x)
  }
  
  RowWiseNoise <- function(ncol, nrow, rho, sd){
    x <- replicate(nrow, arima.sim( n = ncol ,
                                    mean = 0, list(order = c(1,0,0),ar = rho),
                                    sd = sd ))
    return(x)
  }
  
  #X <- X1
  
  ### rho = .1
  rho <- c(.1); sd <- 1
  Ec <- lapply(seq_along(X1), FUN = function(x) ColWiseNoise(ncol=nTime, nrow = nvoxels, rho=rho, sd = sd))
  Er <- lapply(seq_along(X1), FUN = function(x) ColWiseNoise(ncol=nTime, nrow = nvoxels, rho=rho, sd = sd))
  Xe <- lapply(seq_along(X1), function(x) X1[[x]] + Ec[[x]] + Er[[x]])
  
  save(Xe, file=paste(nameext,"AR_10_CICA_simdata.Rdata",sep = ""))
  
  
  ### rho = .3
  rho <- c(.3); sd <- 1
  Ec <- lapply(seq_along(X1), FUN = function(x) ColWiseNoise(ncol=nTime, nrow = nvoxels, rho=rho, sd = sd))
  Er <- lapply(seq_along(X1), FUN = function(x) ColWiseNoise(ncol=nTime, nrow = nvoxels, rho=rho, sd = sd))
  Xe <- lapply(seq_along(X1), function(x) X1[[x]] + Ec[[x]] + Er[[x]])
  
  save(Xe, file=paste(nameext,"AR_30_CICA_simdata.Rdata",sep = ""))
  
  
  ### rho = .7
  rho <- c(.7); sd <- 1
  Ec <- lapply(seq_along(X1), FUN = function(x) ColWiseNoise(ncol=nTime, nrow = nvoxels, rho=rho, sd = sd))
  Er <- lapply(seq_along(X1), FUN = function(x) ColWiseNoise(ncol=nTime, nrow = nvoxels, rho=rho, sd = sd))
  Xe <- lapply(seq_along(X1), function(x) X1[[x]] + Ec[[x]] + Er[[x]])
  
  save(Xe, file=paste(nameext,"AR_70_CICA_simdata.Rdata",sep = ""))
  
}

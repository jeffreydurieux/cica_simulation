### CICA data generation script

SIM_CICA <- function(NBlocks, NSig, NVoxels, NClus, TCtype=NULL ,TClength, Error
                     ,errortype = "Gaussian",IrrSig=0, CommonSig=0, sd=1 ){

  if(IrrSig >= 1 | IrrSig < 0) {
    stop("IrrSig should be a number between 0 and 1 (Percentage)")
  }
  if(CommonSig >= 1 | CommonSig < 0) {
    stop("IrrSig should be a number between 0 and 1 (Percentage)")
  }

  if(IrrSig!=0){
    NSigg <- ceiling(NSig * (1-IrrSig))
  }else{
    NSigg <- ceiling(NSig * (1-CommonSig))
  }
  NIrr <- floor(NSig * IrrSig)
  NCommon <- floor(NSig * CommonSig)

  # generate relevant signals and put them in a list
  Sgen <- function(NS, NV){
    return(t(replicate(NS,
                       ica::icasamp(dname = "b",query = "rnd", nsamp = NV) ) ) )
  }
  SR <- replicate(NClus, Sgen(NS=NSigg, NV=NVoxels), simplify = F)

  # generate common signals and put them in a list
  if(NCommon!=0){
    SC <- Sgen(NS=NCommon, NVoxels)
  }

  # generate subject specific signals and put them in a list
  SLf <- function(NBlocks,NIrr){
    out <- lapply(1:NBlocks, function(i)
   t(replicate(NIrr,ica::icasamp("b", "rnd", nsamp = NVoxels) ) ) )
    return(out)
  }

  if(NIrr != 0){
    SLirr <- replicate(NClus, SLf(NBlocks, NIrr), simplify = F)
  }


  # generate time courses
  if( is.null(TCtype) ){
    Agen <- function(NS, TCl){
      return(replicate(NS, rnorm(n = TCl,sd = sd)))
    }
    ARgen <- function(NB){
      return(replicate(NB, Agen(NS=NSig, TCl=TClength), simplify = F ))
    }
    Air <- replicate(NClus, ARgen(NB=NBlocks), simplify = F)
  }else{
    Air <- replicate(NClus, Gen_TimeCourses(NBlocks, NSig, TClength),simplify = F)
  }

  # mix data
  Xmix <- function(Air, Sr, Sir, SC){
    # input Air = cluster list of list Ai's
    # input Sr  = cluster list of Sr's
    # input Sir = list of irrelevant Si's
    # input SC  = matrix of common signals
    Lfun <- function(Al, Sr, SSir, SC){
      if(missing(SSir) & missing(SC)){
        out <- lapply(1:length(Al), function(i, a, s){a[[i]]%*%s}, a=Al,s=Sr)
      }else if(missing(SC)){
        out <- lapply(1:length(Al),
                      function(i, a, s, sir){a[[i]]%*%rbind(s,sir[[i]])},
                      a=Al,s=Sr,sir=SSir)
      }else{
        out <- lapply(1:length(Al),
                      function(i, a, s, sc){a[[i]]%*%rbind(s,sc)},
                      a=Al,s=Sr,sc=SC)
      }
      return(out)
    }

    AirLength <- length(Air)
    if(missing(Sir) & missing(SC)){
      out <- lapply(1:AirLength, function(i) Lfun(Al=Air[[i]], Sr=SR[[i]]))
    }else if(missing(SC)){
      out <- lapply(1:AirLength, function(i) Lfun(Al=Air[[i]], Sr=SR[[i]],
                                                  SSir=SLirr[[i]]) )
    }else{
      out <- lapply(1:AirLength, function(i) Lfun(Al=Air[[i]], Sr=SR[[i]],
                                                  SC=SC) )
    }
  return(out)
  }

  if(NIrr == 0 & NCommon == 0){X <- Xmix(Air, SR)}
  if(NIrr !=0 & NCommon ==0){X <- Xmix(Air, SR, SLirr)}
  if(NIrr == 0 & NCommon !=0){X <- Xmix(Air, SR, SLirr, SC)}


  X <- unlist(X, recursive = F)

  # add gaussian error
  Xe <- lapply(1:length(X), function(i) addError(X[[i]] , error = Error, type = errortype))
  # out list
  out <- list()
  out$P <- sort(rep(1:NClus,NBlocks))
  out$SR  <- SR
  if(NIrr !=0){
    out$SLirr <- unlist(SLirr, recursive = F)
  }
  if(NCommon!=0){
    out$SCommon <- SC
  }
  out$AIR <- Air
  out$X   <- lapply(X , function(i) t(i))
  out$XE  <- lapply(Xe, function(i) t(i))
  return(out)
}

#helper functions
addError<-function(datablock,error, type = "Gaussian", additiontype = 2)
{
  if(type == "Gaussian"){
    errorM<-replicate(ncol(datablock),rnorm(nrow(datablock)))
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

SSequal<-function(m1,m2)
{
  #res<-(m1/sqrt(SSsuga(m1)) * sqrt(SSsuga(m2))) #c++
  res<-(m1/sqrt(sum(m1^2)) * sqrt(sum(m2^2))) #R
  return(res)
}

# function added on 22-02-2017
# mixture noise option added 22-11-2017
Gen_TimeCourses <- function(N ,nTS , nscan, SNR = 1, type="none"){
  A <- list()
  for(nobs in 1:N){
    TS <- matrix(data = NA, nrow = nscan, ncol = nTS)
    for(nts in 1:nTS){# add a piece of code such that NaN's do not occur
      options(warn = -1) #temp warning supression
      repeat{
        TS[ ,nts] <- neuRosim::simTSrestingstate(nscan = nscan, TR = 2, SNR = SNR, noise = type, weights = c(.25,.25,.25,.25))
        if(any(is.finite(TS[,nts])) == TRUE){
          break
        }#end break
      }#end repeat
    }# end for nts
    options(warn = 0) #turn warnings on
    A[[nobs]] <- TS
  }# end for nobs
  return(A)
}


# this function simulates a single region
# is used in function: SignalStructure
SimRegion <- function(dim = c(64,64), int = 6){

  blsize <- floor(sqrt(dim[1]))

  szero <- matrix(0, dim[1], dim[2])
  xseed <- sample( (blsize+1) :(dim[1] - (blsize+1)),1)
  yseed <- sample( (blsize+1) :(dim[1] - (blsize+1)),1)

  for(i in 1:blsize){
    for(j in 1:blsize){
      szero[ xseed+i, yseed+j ] <- runif(1, min = int,max = int+1)
    }
  }
  return(szero)
}


# this function simulates a single spatial map with 'active' regions
SignalStructure <- function(dim = c(64,64), regions = 2, vec = TRUE){

  # first randomly generate a double exponential (superGaussian)
  s <- icasamp("b", "rnd", nsamp = dim[1]*dim[2])

  # simulate regions with high intensity
  str <- replicate(regions ,SimRegion(dim,int = max(s) ))
  arr <- matrix(0, nrow = dim(str)[1], ncol = dim(str)[2])
  for(i in 1:(dim(str)[3])){
    arr <- arr + str[ , , i]
  }

  res <- matrix(s,dim[1], dim[2]) + arr
  if(vec == TRUE){
    res <- as.vector(res)
  }
  return(res)
}


SimCICARealistic <- function(NBlocks = 10, NClus = 2,NSig = 2,
                             dim = c(64,64),TClength = 100, Error = 0.05,
                             regions = 2, Xscale = TRUE, Xscalef = 1000){
  # Generate signals
  SR <- lapply(1:NClus, function(x) replicate(NSig,SignalStructure(dim = dim, regions = 2,
                                                                   vec = TRUE)) )

  # Generate time courses per cluster
  AIR <- lapply(1:NClus, function(x) Gen_TimeCourses(N = NBlocks, nTS = NSig, nscan = TClength) )

  # Mix signals in a for loop, each iteration is for one cluster
  X <- list()
  for(i in 1:NClus){
    X[[i]] <- lapply(1:NBlocks, function(x) SR[[i]] %*% t(AIR[[i]][[x]]) )
  }

  # put X in one list instead of #NClus lists
  X <- unlist(X,recursive = FALSE)

  # add error per data block
  Xe <- lapply(X, FUN = addError, error = Error)

  if(Xscale == TRUE){
    fun <- function(Xi, fac){
      f <- sqrt(fac/sum(Xi^2))
      return(f*Xi)
    }
    Xe <- lapply(Xe, FUN = fun, fac = Xscalef)
  }

  # out
  out <- list()
  out$P <- sort(rep(1:NClus,NBlocks))
  out$SR  <- SR
  out$AIR <- AIR
  out$X   <- X
  out$XE  <- Xe
  return(out)
}

#### add SNR noise, either Gaussian or AR1 process
addSNRnoise <- function(ts, SNR, type = "Gaussian"){
  sigma <- sd(ts) / SNR

  if(type == "Gaussian"){
    noise <- rnorm(length(ts), mean = 0, sd = sigma)
  }else if(type == "AR1"){
    # rho = 0.2, model order 1
    noise <- arima.sim(model = list(ar = 0.2, ma = 0), n = length(ts), sd = sigma)
  }


  out <- list()
  out$newts <- ts+noise
  out$noise <- noise
  out$sigma <- sigma
  return(out)
}

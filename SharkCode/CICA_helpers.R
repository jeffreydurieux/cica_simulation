### CICA version 1.0 helper functions
### Jeffrey Durieux

### CICA start -> a single start of the CICA algorithm
CICAStart <- function(DataList, exSig, nClus, newclus,icaf) {

  ## preperations before repeat loop starts

  #save first clustering. Used for search empty cluster
  OriginalClus <- newclus

  ## concatenate all datablocks for first Loss function value
  # fun: concdata
  ConcData <- ConcDataFun(DataList = DataList, ClusVec = rep(1, length(DataList)))
  LossHis <- numeric()
  LossHis[1] <- sum(ConcData[[1]]^2) + 1
  LastIter <- numeric()
  counter <- 1

  #repeat loop for 1 start of CICA
  repeat{

    # concatenate data based on P
    # fun: concdata
    ClusteredDataList <- ConcDataFun(DataList = DataList, ClusVec = newclus)
    # estimate Sr and Ai's
    # fun: ExtractSk
    if(icaf == "custom"){
      ICA <- ExtractSkcust(ClusteredDataList, exSig)
    }else{
      ICA <- ExtractSk(ClusteredDataList, exSig)
    }

    # Update P
    # fun: Reclus
    newclusL <- Reclus(DataList, ICA$Sk)
    # check if there are no emptyclusters, if so put the worst fitting block into that empty cluster (= function call)
    # fun: SearchEmptyCluster
    newclus <- SearchEmptyCluster(OriginalClus, newclusL$newclus, newclusL$SSList)
    #keep track of Loss history
    LossHis <- c(LossHis, newclusL$Loss)

    #eye candy not necessary for Corleone note that you need to add startnumber to function arguments
    #cat("#Starts: ", startnumber, "Iterations: ", counter,"\n")
    counter <- counter + 1

    #break statement ALS step 4
    if (abs(LossHis[counter] - LossHis[counter - 1]) < 0.000001) {
      LastIter <- counter - 1
      break()
    }

  }#repeat

  ResList <-
    list(
      "Sr" = ICA$Sk, "Air" = ICA$Ak, "P" = newclus,
      "LossHistory" = LossHis, "Iterations" = LastIter
    )
  return(ResList)

}#end CICAstart function

### f_ConcData
ConcDataFun <- function(DataList, ClusVec){
  ClusL <- length(unique(ClusVec))
  NewList <- list()
  for(i in 1:ClusL){
    NewList[[i]] <- do.call(cbind, DataList[ClusVec == i])
  }
  return(NewList)
} # end f_ConcData

### f_ExtractSk
ExtractSk <- function(ClusList,exSig){
  icaListCluster <- lapply(ClusList, ica::icafast, nc = exSig) # could use pbapply for progressbar

  #extract s and A results and put them in a list
  ICA_S <- lapply(icaListCluster, function(anom) anom$S)
  ICA_A <- lapply(icaListCluster, function(anom) t(anom$M))
  ListRes <- list("Sk" = ICA_S, "Ak" = ICA_A)
  return(ListRes)
}# ExtractSk

ExtractSkcust <- function(ClusList,exSig){
  icaListCluster <- lapply(ClusList, icafast.jeffrey, nc = exSig) # could use pbapply for progressbar

  #extract s and A results and put them in a list
  ICA_S <- lapply(icaListCluster, function(anom) anom$S)
  ICA_A <- lapply(icaListCluster, function(anom) t(anom$M))
  ListRes <- list("Sk" = ICA_S, "Ak" = ICA_A)
  return(ListRes)
}# ExtractSk


### f_Reclus
# fun: ssvec
Reclus <- function(BlockList, SkList) {
  SSList <- lapply(BlockList, FUN = SSvec, SkL = SkList) # function call
  SSminIndexVec <- sapply(SSList, FUN = which.min)
  SSminVec <- sapply(SSList, FUN = min)
  Loss <- sum(SSminVec)
  ResList <- list("newclus" = SSminIndexVec,"Loss" = Loss,"SSList" = SSList)
  return(ResList)
}#function

## f_ssvec
SSvec <- function(Xi, SkL)
{
  Ah <- lapply(seq_along(SkL), function(i) {
    Ah <- t(t(Xi) %*% SkL[[i]] %*% NMFN::mpinv(t(SkL[[i]]) %*% SkL[[i]]))
  })

  Xhat <- lapply(seq_along(Ah), function(i) {
    Xhat <- SkL[[i]] %*% Ah[[i]]
  })

  SS <- sapply(seq_along(Xhat), function(i) {
    sum((Xi - Xhat[[i]]) ^ 2)
  })
  return(SS)
}#function

### f_SearchEmptycluster
SearchEmptyCluster <- function(oldcluster, newcluster, SSList) {
  OriCluster <- sort(unique(oldcluster))
  test <- sapply(OriCluster, FUN = '%in%', newcluster)
  #test result = no empty clusters so return original newcluster
  if (all(test == TRUE)) {
    newcluster <- newcluster
  }else{
    EmptyClusters <- which(test == FALSE)
    worst <- sort(sapply(SSList, FUN = max), decreasing = TRUE)
    SSArray <- simplify2array(SSList)
    BlIndexL <- lapply(seq_along(EmptyClusters), function(i) FUN = which(SSArray == worst[i], arr.ind = T))
    BlIndex <- sapply(BlIndexL, function(anom) anom[2])
    newcluster <- replace(newcluster, BlIndex, EmptyClusters)
  }# else some emptyclusters
  return(newcluster)
}# SearchEmptyCluster

### f_LowestLoss
LowestLoss<-function(MS_CICA.ResultStarts)
{
  minSS<-numeric()
  for(i in 1:length(MS_CICA.ResultStarts))
  {
    LastIter<-MS_CICA.ResultStarts[[i]]$Iterations
    minSS[i]<-MS_CICA.ResultStarts[[i]]$LossHistory[LastIter]
  }

  resultIndex<-which(minSS==min(minSS))
  result<-MS_CICA.ResultStarts[[resultIndex[1]]] #select just one

  return(result)
}# Last SS function



#------------------------- random Cluster functions (Tom's code)
GenerateRandomClustering <- function(nElement , nClust , Prob)
{
  ####GenerateRandomClustering = for Random Starts

  # nElement: number of elements to be clustered
  # nClust: number of clusters
  # Prob (1 x nClust): proportion of elements in each cluster

  BestClust = NULL
  ErrorEncountered = F

  if (!(length(Prob) == nClust))
  {
    cat('there should be as much probabilities as clusters')
    ErrorEncountered = T
  }

  if ((abs(sum(Prob) - 1) > .000000001) | (any(Prob < 0)))
  {
    cat('probabilities should sum to one (and cannot be negative)')
    ErrorEncountered = T
  }

  if (!(any(nClust == 1:nElement)))
  {
    cat("nClust should be a number between 1 and 'nElement'")
    ErrorEncountered = T
  }

  if (!(ErrorEncountered))
  {
    if (nElement > nClust)
    {
      if (nClust == 1)
      {
        BestClust = rep(1 , times = nElement)
      }
      else
      {
        ProbVV = round(Prob * nElement)
        if (!(sum(ProbVV) == nElement) |
            (any(ProbVV < 1)))
          #not enough elements, or empty clusters
        {
          ProbVV = AdjustProb(ProbVV , nElement)
        }

        tempclus = rep(1:length(ProbVV) , ProbVV)
        BestClust = tempclus[sample(1:nElement,size = nElement,replace =
                                      FALSE)]
      }
    }
    else
    {
      BestClust = 1:nClust
    }
  }

  if (!(length(unique(BestClust)) == nClust))
  {
    BestClust = NULL
  }

  return(BestClust)
}

AdjustProb <- function(v , MaxElem)
{
  # INPUT
  #   v (1 x nElem): vector
  #   MaxElem: number of elements
  #
  # OUTPUT
  #   v (1 x nElem): vector with 'sum(out)=MaxElem' and no 0-cells

  nElem = length(v)

  if (any(v < 1))
    #add 1 to 0-cells
  {
    tempv = rep(0,nElem)
    tempv[v < 1] = rep(1 , length(which(v < 1)))
    v = v + tempv
  }#### replace(v,which(v<1),1)

  while (!(sum(v) == MaxElem))
  {
    diff = sum(v) - MaxElem
    if (diff < 0)
      #add elements
    {
      if (abs(diff) < (nElem - sum(v == 1)))
      {
        for (tel in 1:abs(diff))
        {
          tempcl = ceiling(runif(1) * nElem)
          while (v[tempcl] == 1)
          {
            tempcl = ceiling(runif(1) * nElem)
          }
          v[tempcl] = v[tempcl] + 1
        }
      }
      else
      {
        v[which(!(v == 1))] = v[which(!(v == 1))] + rep(1 , length(which(!(v == 1))))
      }
    }
    else
      #delete elements
    {
      if (abs(diff) < (nElem - sum(v == 1)))
      {
        for (tel in 1:abs(diff))
        {
          tempcl = ceiling(runif(1) * nElem)
          while (v[tempcl] == 1)
          {
            tempcl = ceiling(runif(1) * nElem)
          }
          v[tempcl] = v[tempcl] - 1
        }
      }
      else
      {
        v[which(!(v == 1))] = v[which(!(v == 1))] - rep(1 , length(which(!(v == 1))))
      }
    }
  }

  return(v)
}

clusf <- function(nBlocks, nClus) {
  #simplyfied cluster generation function using an equal probability
  clus <-
    GenerateRandomClustering(nBlocks, nClus, rep(c(1 / nClus), nClus))
}

Asplit <- function(Air, P, nTS){
  Unique_clus <- sort(unique(P))
  Al <- list()
  for(i in 1:length(Unique_clus)){
    obs <- which(P == Unique_clus[i])
    parts <- split(1:(nTS*length(obs)), ceiling(seq_along(1:(nTS*length(obs))) /nTS ))
    Al[[i]] <- lapply(parts, function(x) Air[[i]][,x])
  }
  return(Al)
}


to2d <- function(nif){
  dim_nif <- dim(nif)
  if( length(dim_nif) > 4){
    stop("Nifti file contains more than four dimensions!")
  }
  two_D <- matrix(nif, nrow = dim_nif[4], byrow = TRUE)
  return(two_D)
}

# 2d to 4d Rfunction
to4d <- function(twod, dim){
  four_D <- array(t(twod), dim = dim)
  return(four_D)
}

# Wed Mar 15 19:49:03 2017 ------------------------------
# function to strip empty voxels
StripEmptyVoxels <- function(Data, onlyidx = T) {
  # input Data t x v
  # output <- stripped data and dim dim vec

  EmptyVoxelIndex <- vector("numeric", length = ncol(Data))
  for(i in seq_along(EmptyVoxelIndex) ){
    test <- all(Data[, i ] == 0)
    EmptyVoxelIndex[i] <- test
  }

  out <- list()
  if(onlyidx != TRUE){
    out$Data <- Data[ , !EmptyVoxelIndex ]
  }
  out$Empty <- EmptyVoxelIndex
  out$Dim <- c(nrow(Data) , ncol(Data) )
  return( out )
}

# this function adds empty voxels to data to return to original dimensions
# input is a StripEmptyVoxels object
AddEmptyVoxels <- function(Data, index ,dimvec){
  Mat <- matrix(data = 0, nrow = dimvec[1] ,
                ncol = dimvec[2] )

  Mat[ , which(index==0)] <- Data
  return( Mat )
}

#CommonNonEmpty
StripEmptyVoxelsMultiple <- function(DataList, OnlyCommon = T){

  #change into lapply etc.
  StripEmptyVoxels_objectList <- lapply(DataList, StripEmptyVoxels)

  empties <- sapply(StripEmptyVoxels_objectList, function(anom) anom$Empty)
  check <- apply(empties, 1, sum) # if all 0 then common non-empty voxel
  if(OnlyCommon == T){
    idx <- which(check == 0)
  }else{
    first <- head( which(check == 1) , n = 1 )
    last <- tail( which(check == 1) , n = 1 )
    idxx <- rep(1, length(check) )
    idxxx <- replace(idxx , list =  first:last, values = 0)
    idx <- which(idxxx == 0)
  }

  DataL <- lapply(1:length(DataList), function(x){
    DataList[[x]][,idx]
  })
  out <- list()
  out$Data <- DataL
  out$NonEmptyIdx <- idx
  out$NonEmptyBool <- check

  return(out)
}

# Thu Mar 16 13:34:31 2017 ------------------------------
# subtract means


# Thu Mar  9 13:55:17 2017 ------------------------------

# Intensity normalization
# divide each voxel time course by its mean, comparable to converting
# to percent signal change
IntNormalization <- function(Data){
  # local function
  Normalize <- function(x){
    return( x / mean(x) )
  }

  if ( is.list(Data) ) {
    out <- lapply(Data, function(anom) apply(anom, 2, Normalize) )
  }else{
    out <- apply(Data, 2, Normalize)
  }
  return(out)
}

# Variance normalization
# z-score each voxel time course to have zero mean and variance of one
VarNormalization <- function(Data) {
  if(is.list(Data)){
    out <- lapply(Data, scale)
  }else{
    out <- scale(Data)
  }
  return(out)
}

# Row normalization
# i.e., center rows
RowCent <- function(Data){
  if ( is.list( Data) ) {
    trans <- lapply(Data, t)
    trans <- lapply(trans, scale, center = TRUE, scale = FALSE)
    out <- lapply(trans, t)
  }
  else{
    out <- t( scale(x = t(Data), center = TRUE, scale = FALSE ) )
  }
  return(out)
}


SubtractMean <- function(Data, t=T, v=T){

  if(v == T){
    mv <- colMeans(Data)
    Dat <- sweep(Data, MARGIN = 2, STATS = mv, FUN = '-')
  }
  if(t == T){
    tm <- rowMeans(Data)
    Dat <- sweep(Data,1,STATS = tm,FUN = '-')
  }

  return(Dat)
}


########## smart start section
ComputeRVmat <- function(S){
  l <- length(S)
  mat <- matrix(NA,l,l)
  diag(mat) <- 1
  comb <- combn(1:l, 2)
  for(i in 1:ncol(comb)){
    mat[comb[1,i] , comb[2,i]] <- MatrixCorrelation::RV(S[[ comb[1,i] ]], S[[ comb[2,i] ]])
    #cat(paste(comb[,i],"\n" ))
  }
  return(mat)
}

ComputeRVmatPar <- function(S,Part){
  l <- length(S)
  pidx <- (1:choose(l,2))[Part]
  mat <- matrix(0,l,l)
  #diag(mat) <- 1
  comb <- combn(1:l, 2)
  comb <- comb[,pidx]
  for(i in 1:ncol(comb)){
    mat[comb[1,i] , comb[2,i]] <- MatrixCorrelation::RV(S[[ comb[1,i] ]], S[[ comb[2,i] ]])
    #cat(paste(comb[,i],"\n" ))
  }
  return(mat)
}

HCL <- function(RVmat, maxk){
  d <- as.dist(1-RVmat)
  s.hcl <- hcl <- list()
  for (i in 2:(maxk+1)){
    s.hcl[[i-1]]<-cutree(hclust(d), k=i)
    hcl[[i-1]]<-fpc::cluster.stats(d, s.hcl[[i-1]])
  }

  outcome.h<-matrix(0,2,maxk)
  for (i in 1:maxk){
    outcome.h[,i]=c(hcl[[i]]$avg.silwidth,hcl[[i]]$pearsongamma)
  }

  best <- apply(outcome.h,1,which.max)+1
  val <- c(outcome.h[1,best[1]-1] , outcome.h[2,best[2]-1]  )
  mat <- rbind(best,val)
  colnames(mat) <- c("silhouette" , "Pearson Gamma")
  res <- list()
  res$best <- mat
  res$labels <- s.hcl
  res$stats <- hcl
  res$dist <- d
  return(res)
}

ComputeTuckerWithPermutation <- function( Atrue , A , nComp )
{
  #require(gtools)
  #source("ComputePhiMatrix.R")

  Perms = gtools::permutations(n = nComp, r = nComp, v = 1:nComp)
  nPerms = dim(Perms)[1]

  BestValue = -9999
  BestValuePerComp = matrix( -9999 , 1 , nComp )
  BestPermut = matrix( -9999 , 1 , nComp )

  for( tel in 1:nPerms )
  {
    Phi = ComputePhiMatrix ( Atrue , A[ , Perms[tel,] ] )
    ValuePerComp = abs ( diag ( Phi ) )
    TuckerValue = mean( ValuePerComp )
    if( TuckerValue > BestValue )
    {
      BestValue = TuckerValue
      BestValuePerComp = ValuePerComp
      BestPermut = Perms[tel,]
    }
  }

  Out = list()
  Out$Tucker = BestValue
  Out$TuckerComp = BestValuePerComp
  Out$BestPermut = BestPermut
  return(Out)
}

dataloader <- function(filenames){
  dataL <- list()
  for(i in 1:length(filenames)){
    dataL[[i]] <- get(load(filenames[i]))
  }

  # temp.space <- new.env()
  # bar <- load('saved.file.rda', temp.space)
  # the.object <- get(bar, temp.space)
  # rm(temp.space)
  names(dataL) <- filenames
  return(dataL)
}

icafast.jeffrey <- function (X, nc, center = TRUE, maxit = 100, tol = 1e-06, Rmat = diag(nc),
                             alg = c("par", "def"), fun = c("logcosh", "exp", "kur"),
                             alpha = 1)
{

  # some internal functions
  sdiag <- function (x)
  {
    if (length(x) < 2L)
      matrix(x)
    else diag(x)
  }

  #

  X <- as.matrix(X)
  nobs <- nrow(X)
  nvar <- ncol(X)
  nc <- as.integer(nc[1])
  if (nc < 1) {
    stop("Must set nc>=1 component.")
  }
  maxit <- as.integer(maxit[1])
  if (maxit < 1) {
    stop("Must set maxit>=1 iteration.")
  }
  tol <- tol[1]
  if (tol <= 0) {
    stop("Must set ctol>0.")
  }
  if (nc > min(nobs, nvar)) {
    stop("Too many components. Set nc<=min(dim(X)).")
  }
  alpha <- alpha[1]
  if (alpha < 1 | alpha > 2) {
    stop("Must set 'alpha' between 1 and 2.")
  }
  if (nrow(Rmat) != nc | ncol(Rmat) != nc) {
    stop("Input 'Rmat' must be nc-by-nc rotation matrix.")
  }
  if (center)
    X <- scale(X, scale = FALSE)
  xeig <- eigen(crossprod(X)/nobs, symmetric = TRUE)
  nze <- sum(xeig$val > xeig$val[1] * .Machine$double.eps)
  if (nze < nc) {
    warning("Numerical rank of X is less than requested number of components (nc).\n  Number of components has been redefined as the numerical rank of X.")
    nc <- nze
    Rmat <- diag(nc)
  }
  Dmat <- sdiag(sqrt(xeig$val[1:nc]))
  Mprt <- tcrossprod(Dmat, xeig$vec[, 1:nc])
  diag(Dmat) <- 1/diag(Dmat)
  Pmat <- xeig$vec[, 1:nc] %*% Dmat
  Xw <- X %*% Pmat
  if (nc == 1L) {
    return(list(S = Xw, M = Mprt, W = t(Pmat), Y = Xw, Q = t(Pmat),
                R = matrix(1), vafs = (sum(Mprt^2) * nobs)/sum(X^2),
                iter = NA, alg = alg, fun = fun, alpha = alpha))
  }
  if (fun[1] == "kur") {
    fun1d <- function(x) {
      x^3
    }
    fun2d <- function(x) {
      3 * (x^2)
    }
  }
  else if (fun[1] == "exp") {
    fun1d <- function(x) {
      x * exp(-(x^2)/2)
    }
    fun2d <- function(x) {
      exp(-(x^2)/2) * (1 - x^2)
    }
  }
  else {
    fun1d <- function(x) {
      tanh(alpha * x)
    }
    fun2d <- function(x) {
      alpha * (1 - tanh(alpha * x)^2)
    }
  }
  if (alg[1] == "def") {
    myiters <- rep(NA, nc)
    for (j in 1:nc) {
      if (j < 2) {
        Rmat[, j] <- Rmat[, j]/sqrt(sum(Rmat[, j]^2))
        iter <- 0
        vtol <- 1
        while (vtol > tol && iter < maxit) {
          svec <- Xw %*% Rmat[, j]
          rnew <- colMeans(Xw * matrix(fun1d(svec), nobs,
                                       nc))
          rnew <- rnew - mean(fun2d(svec)) * Rmat[, j]
          rnew <- rnew/sqrt(sum(rnew^2))
          vtol <- 1 - abs(sum(Rmat[, j] * rnew))
          iter <- iter + 1
          Rmat[, j] <- rnew
        }
        myiters[j] <- iter
      }
      else {
        Rmat[, j] <- Rmat[, j]/sqrt(sum(Rmat[, j]^2))
        svec <- matrix(0, nc, 1)
        for (k in 1:(j - 1)) {
          svec <- svec + sum(Rmat[, k] * Rmat[, j]) *
            Rmat[, k]
        }
        Rmat[, j] <- Rmat[, j] - svec
        Rmat[, j] <- Rmat[, j]/sqrt(sum(Rmat[, j]^2))
        iter <- 0
        vtol <- 1
        while (vtol > tol && iter < maxit) {
          svec <- Xw %*% Rmat[, j]
          rnew <- colMeans(Xw * matrix(fun1d(svec), nobs,
                                       nc))
          rnew <- rnew - mean(fun2d(svec)) * Rmat[, j]
          rnew <- rnew/sqrt(sum(rnew^2))
          svec <- matrix(0, nc, 1)
          for (k in 1:(j - 1)) {
            svec <- svec + sum(Rmat[, k] * rnew) * Rmat[,
                                                        k]
          }
          rnew <- rnew - svec
          rnew <- rnew/sqrt(sum(rnew^2))
          vtol <- 1 - abs(sum(Rmat[, j] * rnew))
          iter <- iter + 1
          Rmat[, j] <- rnew
        }
        myiters[j] <- iter
      }
    }
  }
  else {
    rsvd <- svd(Rmat)
    Rmat <- tcrossprod(rsvd$u, rsvd$v)
    iter <- 0
    vtol <- 1
    while (vtol > tol && iter < maxit) {
      smat <- Xw %*% Rmat
      rnew <- crossprod(Xw, fun1d(smat))/nobs
      rnew <- rnew - Rmat %*% sdiag(colMeans(fun2d(smat)))
      rsvd <- svd(rnew)
      rnew <- tcrossprod(rsvd$u, rsvd$v)
      vtol <- 1 - min(abs(colSums(Rmat * rnew)))
      iter <- iter + 1
      Rmat <- rnew
    }
    myiters <- iter
  }
  M <- crossprod(Rmat, Mprt)
  vafs <- rowSums(M^2)
  ix <- sort(vafs, decreasing = TRUE, index.return = TRUE)$ix
  M <- M[ix, ]
  Rmat <- Rmat[, ix]
  vafs <- (vafs[ix] * nobs)/sum(X^2)

  # put addition here:
  #return(list(S = Xw %*% Rmat, M = t(M), W = t(Pmat %*% Rmat),
  #           Y = Xw, Q = t(Pmat), R = Rmat, vafs = vafs, iter = myiters,
  #          alg = alg[1], fun = fun, alpha = alpha))
  res <- list(S = Xw %*% Rmat, M = t(M), W = t(Pmat %*% Rmat),
              Y = Xw, Q = t(Pmat), R = Rmat, vafs = vafs, iter = myiters,
              alg = alg[1], fun = fun, alpha = alpha)

  s <- res$S; m <- res$M
  var <- apply(m,2,var);var <- var
  sn <- s %*% diag(var); mn <- m %*% diag(1/var)
  ress <- list(S = sn, M = mn)
  return(ress)

}

xscale <- function(Xi){
  f <- sqrt(1000/sum(Xi^2))
  return(f*Xi)
}

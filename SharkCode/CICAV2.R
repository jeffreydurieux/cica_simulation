#### functions #####
# Tucker diag
TuckerDiag <- function(Xi, Xcent, type = "diag"){
  con <- congru( Xi, Xcent)
  if(type == 'diag'){
    TC <- mean( abs( diag( con ) ) )  
  }
  else{
    TC <- mean( apply( abs(con), 2, max) )
  }
  
  return(TC)
}

TuckerDiagList<- function(XL, Xcent, type = 'diag'){
  tucks <- sapply( seq_along(XL), function(x) TuckerDiag(XL[[x]], Xcent, type = type) ) 
  return( tucks )
}

# reconstruction and centroid functions
Centroid <- function(Sr,ALr){
  Xl <- lapply(seq_along(ALr) , function(x) Sr %*% t(ALr[[x]])) 
  Cent <- Reduce('+', Xl) / length(ALr)
  return(Cent)
}

Xrecon <- function(Sr, ALr){
  Xl <- lapply(seq_along(ALr) , function(x) Sr %*% t(ALr[[x]])) 
}

XRecon_Centroid <- function(Sr,ALr){
  Xl <- lapply(seq_along(ALr) , function(x) Sr %*% t(ALr[[x]])) 
  Cent <- Reduce('+', Xl) / length(ALr)
  Res <- list()
  Res$Xl <- Xl
  Res$Cent <- Cent
  return(Res)
}

# concatenate and split functions
ConcData <- function(DataList, ClusVec){
  ClusL <- length(unique(ClusVec))
  NewList <- list()
  for(i in 1:ClusL){
    NewList[[i]] <- do.call(cbind, DataList[ClusVec == i])
  }
  return(NewList)
}

Asplitter <- function(AL, nTS, nV){
  len <- length(AL)
  Al <- list()
  for(i in 1:len){
    obs <- which(P == Unique_clus[i])
    parts <- split(1:(nTS*length(obs)), ceiling(seq_along(1:(nTS*length(obs))) /nTS ))
    Al[[i]] <- lapply(parts, function(x) Air[[i]][,x])
  }
  return(Al)
}

# split group ica A result 
Asplitter <- function(A, nTS, nV){
  obs <- nrow(A)/nTS
  parts <- split(1:(nTS*obs), ceiling(seq_along(1:(nTS*obs)) /nTS ))
  
  As <- lapply(seq_along(parts), function(x) A[ parts[[x]] , ]  )
  return(As)
}

### GROUPICA simplified output ####
GICA <- function(XL, nc = 2){
  GICAS <- lapply(XL, FUN = icafast, nc = nc)
  Res <- list()
  Res$SR <- lapply(seq_along(GICAS), function(x) GICAS[[x]]$S)
  Res$AL <- lapply(seq_along(GICAS), function(x) GICAS[[x]]$M)
  return(Res)
}

#### Random P ####
RandomP <- function(nBlocks, nClus) {
  #simplyfied cluster generation function using an equal probability
  clus <-GenerateRandomClustering(nBlocks, nClus, rep(c(1 / nClus), nClus))
}

GenerateRandomClustering <- function(nElement , nClust , Prob){
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


#### SSQ criterion #####
SSQ_XE_Xhat <- function(XE, Xhat){
  res <- sapply(seq_along(XE), function(x) sum( (XE[[x]] - Xhat[[x]])^2) )
  return(sum(res))
}

SSQ_Xhat_Xcent <- function(Xhat, CentroidL){
  
  X <- unlist(Xhat, recursive = F)
  TotalSS <- sum(sapply(seq_along(X), function(x) sum(X[[x]]^2) ))
  
  len <- length(CentroidL)
  res <- numeric()
  for(i in 1:len){
    lenn <- length(Xhat[[i]])
    res[i] <- sum(sapply(1:lenn, function(x) sum( (Xhat[[i]][[x]] - CentroidL[[i]])^2) ))
  }
  # comment out for more info
  wbss <- res/TotalSS
  res <- sum(wbss) 
  return(res)
}

Xpredicted <- function(Xi, SR){
  Ahat <- t(Xi) %*% SR %*% NMFN::mpinv( t(SR) %*% SR)
  Xhat <- SR %*% t(Ahat)
  return(Xhat)
}

XhatsL <- function(Xl, SR){
  res <- lapply(seq_along(Xl), function(x) 
    Xpredicted(Xl[[x]], SR = SR))
  return(res)
}

SearchEmpty <- function(oldcluster, newcluster, TUCKCOLL) {
  OriCluster <- sort(unique(oldcluster))
  test <- sapply(OriCluster, FUN = '%in%', newcluster)
  #test result = no empty clusters so return original newcluster
  if (all(test == TRUE)) {
    newcluster <- newcluster
  }else{
    while(any(test == FALSE))
      EmptyCluster <- which(test == FALSE)
    mins <- apply(TUCKCOLL[,-EmptyCluster], MARGIN = 2, min)
    idcol <- which.min(mins)
    rowid <- which.min(TUCKCOLL[,-EmptyCluster][,2])
    newcluster[rowid] <- EmptyCluster
    test <- sapply(OriCluster, FUN = '%in%', newcluster)
  }
  return(newcluster)
}


########## main function ##########

CICAV2_1 <- function(XE, nc = 2, nclus = 2, starts = 1, maxiter = 50,
                     criterion = 2, losstest=FALSE){
  # args 
  N <- length(XE)
  nTS <- ncol(XE[[1]])
  nV <- nrow(XE[[1]])
  
  STARTS <- list()
  
  for(st in 1:starts){
    
    LOSS <- numeric()
    
    # first iter using random P
    Pr <- RandomP(nBlocks = N, nClus = nclus)
    Pmat <- do.call(cbind, list(Pr))
    iteration <- 1
    #cat(iteration, "\n")
    repeat{
      #step 2: condata based on P and do GICA
      X <- ConcData(XE, Pmat[ ,iteration])
      GICAS <- GICA(X,nc = nc)
      
      # step 3 centroid calc
      #ALs <- list()
      Centroids <- list()
      Xrecon <- list()
      for(i in 1:nclus){
        ALs <- Asplitter(GICAS$AL[[i]], nTS = nTS, nV = nV)
        Xr <- XRecon_Centroid(GICAS$SR[[i]], ALr = ALs)
        Centroids[[i]] <- Xr$Cent
        Xrecon[[i]] <- Xr$Xl
      }
      
      #LOSS[iteration] <- SSQ_Xhat_Xcent(Xrecon, Centroids)
      
      # step 4 reallocate two types
      # criterion 1) tucker between Xhat and cent
      # criterion 2) tucker between XE and cent
      if(criterion == 1){
        
        TUCKCOL <- matrix(data = NA, nrow = N, ncol = nclus)
        for(i in 1:nclus){
          # compute Xhats based on cluster Si
          # compute tuck centi
          Xhats <- XhatsL(Xl = XE, SR = GICAS$SR[[i]])
          TUCKCOL[ , i] <- TuckerDiagList(Xhats, Centroids[[i]])
        }
      }else{
        TUCKCOL <- matrix(data = NA, nrow = N, ncol = nclus)
        for(i in 1:nclus){
          TUCKCOL[ , i] <- TuckerDiagList(XE, Centroids[[i]])
        }
      }
      
      Pn <- apply(TUCKCOL,1,which.max)
      ### prevent empty cluster
      ### TODO
      #Pn <- SearchEmpty(Pr, Pn, TUCKCOL)
      Pmat <- cbind(Pmat,Pn)
      #Pmat
      #sum(apply(TUCKCOL,1,max))
      Pmat
      iteration <- iteration + 1
      cat('start' , st,'iteration',iteration, "\n")
      ch <- all(Pmat[, ncol(Pmat)] == Pmat[, ncol(Pmat) -1 ])
      ch
      if(ch == TRUE | iteration == maxiter){
        SUMMAX <- sum(apply(TUCKCOL,1,max))
        
        # compute loss based on final P
        if(losstest == TRUE){
          X <- ConcData(XE, Pmat[ ,iteration])
          GICAS <- GICA(X,nc = nc)
          Centroids <- list()
          Xrecon <- list()
          for(i in 1:nclus){
            ALs <- Asplitter(GICAS$AL[[i]], nTS = nTS, nV = nV)
            Xr <- XRecon_Centroid(GICAS$SR[[i]], ALr = ALs)
            Centroids[[i]] <- Xr$Cent
            Xrecon[[i]] <- Xr$Xl
            LOSS[iteration] <- SSQ_Xhat_Xcent(Xrecon, Centroids)
          }  
        }
        
        break
      }
    } # end repeat
    Res <- list()
    Res$P <- Pmat[, ncol(Pmat)]
    Res$SUMMAX <- SUMMAX
    Res$LOSS <- LOSS 
    
    STARTS[[st]] <- Res
  } # end starts
  return(STARTS)
}


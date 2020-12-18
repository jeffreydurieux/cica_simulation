# con 1 results in a low correlation

Sim_A <- function(nscan, n_signals){
  TS <- matrix(data = NA, nrow = nscan, ncol = n_signals)
  for(nts in 1:n_signals){# add a piece of code such that NaN's do not occur
    options(warn = -1) #temp warning supression
    repeat{
      TS[ ,nts] <- neuRosim::simTSrestingstate(nscan = nscan, TR = 2, noise = "none")
      if(any(is.finite(TS[,nts])) == TRUE){
        break
      }#end break
    }#end repeat
  }# end for nts
  options(warn = 0) #turn warnings on
  return(TS)
}

SimData <- function(Q, nscan, Nsamp, n_clusters, con, error, withperm = F){
  
  Sbase <- replicate(n = Q, expr = runif(n = 1000, min = -1, max = 1))
  
  Alist <- list()
  Slist <- list()
  Xlist <- list()
  for(i in 1:n_clusters){
      Alist[[i]] <- replicate( (Nsamp / n_clusters) , 
                               Sim_A(nscan = nscan,n_signals = Q), simplify = F) 
      
      Slist[[i]] <- replicate(n = Q, 
                              runif(n = 1000, min = -con, max = con))
      
    
    
    S <- lapply(seq_along(Slist), function(x) Sbase + Slist[[x]])
    Xlist[[i]] <- lapply(seq_along(Alist[[i]]) , function(x) S[[i]] %*% t(Alist[[i]][[x]]))
    
  }
  
  X <- unlist(Xlist, recursive = F)
  XE <- lapply(X, FUN = addError, error = error)
  
  data <- list()
  data$Xe <- XE
  data$X <- X
  data$S <- S
  data$A <- Alist
  pid <- sapply(Alist, length)
  if(n_clusters == 2){
    data$P <- c(rep(1,pid[1]), rep(2,pid[2]))  
  }else{
    data$P <- c(rep(1,pid[1]), rep(2,pid[2]), rep(3,pid[3]), rep(4,pid[4]))
  }
  
  return(data)
}


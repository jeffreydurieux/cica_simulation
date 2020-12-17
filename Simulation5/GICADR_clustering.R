# Thu Dec 17 13:49:14 2020
# Author: Jeffrey Durieux, MSc

# Script to do GICA + Dual regression + clustering

DRCLUS <- function(simdata, Q, R, hcl = T, pam = T){
  
  Xe <- lapply(simdata$Xe, FUN = CICA:::xscale)
  X <- do.call(Xe, what = cbind)
  
  GICA <- icafast(X, nc = Q)
  
  XX <- Xe
  Ahats <- lapply(seq_along(XX),
                  function(lam) mpinv(GICA$S) %*% XX[[lam]])
  
  ### dual regression step 2:
  
  Shats <- lapply(seq_along(Ahats),
                  function(lam) XX[[lam]] %*% mpinv(Ahats[[lam]]))
  
  #### clustering ####
  N <- length(Xe)
  comb <- t(utils::combn(1:N, 2))
  
  RVsS <- matrix(data = NA, nrow = N , ncol = N)
  RVS <- numeric()
  
  cat("Computing pairwise Tucker statistics: \n")
  pb <- txtProgressBar(min = 0, max = nrow(comb), initial = 0)
  
  for(i in 1:nrow(comb)){
    
    RVS[i] <- mean(diag(Tucker( Shats[[ comb[i,1] ]] , Shats[[ comb[i,2] ]])))
    res <- c(comb[i , ] , RVS[i] )
    
    RVsS[res[1]  , res[2] ] <- res[3]
    
    
    
    setTxtProgressBar(pb, i)
  }
  
  RVsS[lower.tri(RVsS)] = t(RVsS)[lower.tri(RVsS)]
  diag(RVsS) <- 1
  cat('\n')
  RVmat <- as.dist(1-RVsS)
  # hcl 
  if(hcl == T){
    hcl <- hclust(d = RVmat, method = 'ward.D2')
    khcl <- cutree(hcl, k = R)  
  }
  
  # hcl 
  if(pam == T){
    pam <- cluster::pam(RVmat, k = R)
    kpam <- pam$clustering
  }
  
  output <- list()
  output$RV <- RVsS
  output$clustering <- data.frame(khcl = khcl, kpam = kpam)
  return(output)
}
# Fri Dec 11 09:31:48 2020
# Author: Jeffrey Durieux, MSc

# What:
# Extra simulations C-ICA
# Group-ICA + dual regression and clustering:
# - hclust ward
# - PAM
# - kmeans + tandem analysis
# computation of modRV takes a long time
# put design (1:72, parallel )

library(ica)
library(NMFN) #moore penrose
library(mclust) # ARI
library(cluster) # pam


##### functions to use #######

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


Simulate_CICA <- function(N = 40, V = 1000, D = 100, Q = 2, R = 2, E = 0.1){
  # Note input args; D = time Dimension so that it is not confused with T (True)

  # generate cluster specific components (FC patterns)
  # use icasamp from ica package
  # double exponential = b
  # Dimensions: R lists of V X Q matrices
  Qr <- lapply(seq_along(1:R), function(anom) {
    replicate(n = Q, expr = icasamp(dname = 'b', query = 'rnd', nsamp = V))
  })

  # generate random time courses
  # dimension: R lists of T X Q matrices
  NperClus <- N / R
  Air <- lapply(seq_along(1:R), function(anom1) {

    lapply(seq_along(1:NperClus), function(anom2) {
      replicate(n = Q, runif(n = D, min = -2, max = 2))
    })

  })

  # Mix data -> Qr %*% t(Air)
  X <- lapply(seq_along(Qr), function(anom1){
    lapply( seq_along( Air[[anom1]]), function(anom2){
      Qr[[anom1]] %*% t(Air[[anom1]][[anom2]] )
    })
  })

  # concatenate lists into one large list
  X <- unlist(X, recursive = F)

  # add error
  Xe <- lapply(X, FUN = addError, error = E)

  Out <- list()
  Out$P <- rep(1:R, each = NperClus)
  Out$Qr <- Qr
  Out$Air <- Air
  Out$X <- X
  Out$Xe <- Xe
  return(Out)
}


### helper functions for generating gaussian
addError<-function(datablock,error)
{

  errorM<-replicate(ncol(datablock),rnorm(nrow(datablock)))
  errorM<-SSequal(errorM,datablock)
  errorlevel<-error/(1-error)

  res<-datablock + (errorM * sqrt(errorlevel))
  return(res)
}

SSequal<-function(m1,m2)
{
  res<-(m1/sqrt(sum(m1^2)) * sqrt(sum(m2^2))) #R
  return(res)
}

V <- c(500,2000)                # voxels
Q <- c(2,5,20)                  # components
R <- c(2,4)                     # clusters
D <- c('square', 'nonsquare')   # dimension mixing matrix
E <- c(0.05, 0.2, 0.4)          # noise level

design <- expand.grid(V=V, Q = Q,
                      R = R, D = D, E=E)



######### analysis ########
args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)

replication <- args[1]
#row <- args[2]

replication <- 1
#row <- 66

tmp <- proc.time()

for(row in 1:12){
  # This will generate a unique seed

  seed <- as.numeric(paste(replication,0,row,sep = ""))
  set.seed(seed)



  rowd <- design[row , ]

  ### Simulate CICA ###
  dim <- ifelse(test = (rowd$D == 'square'), yes = rowd$Q , no = 100 )
  SimData <- Simulate_CICA(N = 40, V = rowd$V, D = dim, Q = rowd$Q, R = rowd$R,E = rowd$E)



  ##### Group ICA #####

  X <- do.call(cbind, SimData$Xe)
  ica <- icafast(X = X, nc = rowd$Q)

  ### dual regression step 1:
  XX <- SimData$Xe
  Ahats <- lapply(seq_along(XX),
                  function(lam) mpinv(ica$S) %*% XX[[lam]])

  ### dual regression step 2:

  Shats <- lapply(seq_along(Ahats),
                  function(lam) XX[[lam]] %*% mpinv(Ahats[[lam]]))



  ##### Group ICA 2 #####

  if(rowd$R == 2){
    X <- do.call(cbind, SimData$Xe[SimData$P==1])
    ica1 <- icafast(X = X, nc = rowd$Q)

    XX <- SimData$Xe[SimData$P==1]
    Ahats <- lapply(seq_along(XX),
                    function(lam) mpinv(ica1$S) %*% XX[[lam]])

    Shats1 <- lapply(seq_along(Ahats),
                     function(lam) XX[[lam]] %*% mpinv(Ahats[[lam]]))


    X <- do.call(cbind, SimData$Xe[SimData$P==2])
    ica2 <- icafast(X = X, nc = rowd$Q)

    XX <- SimData$Xe[SimData$P==2]
    Ahats <- lapply(seq_along(XX),
                    function(lam) mpinv(ica2$S) %*% XX[[lam]])

    Shats2 <- lapply(seq_along(Ahats),
                     function(lam) XX[[lam]] %*% mpinv(Ahats[[lam]]))
    ShatsG2 <- c(Shats1,Shats2)
  }else{
    X <- do.call(cbind, SimData$Xe[SimData$P==1])
    ica1 <- icafast(X = X, nc = rowd$Q)

    XX <- SimData$Xe[SimData$P==1]
    Ahats <- lapply(seq_along(XX),
                    function(lam) mpinv(ica1$S) %*% XX[[lam]])

    Shats1 <- lapply(seq_along(Ahats),
                     function(lam) XX[[lam]] %*% mpinv(Ahats[[lam]]))

    X <- do.call(cbind, SimData$Xe[SimData$P==2])
    ica2 <- icafast(X = X, nc = rowd$Q)

    XX <- SimData$Xe[SimData$P==2]
    Ahats <- lapply(seq_along(XX),
                    function(lam) mpinv(ica2$S) %*% XX[[lam]])

    Shats2 <- lapply(seq_along(Ahats),
                     function(lam) XX[[lam]] %*% mpinv(Ahats[[lam]]))

    X <- do.call(cbind, SimData$Xe[SimData$P==3])
    ica3 <- icafast(X = X, nc = rowd$Q)

    XX <- SimData$Xe[SimData$P==3]
    Ahats <- lapply(seq_along(XX),
                    function(lam) mpinv(ica3$S) %*% XX[[lam]])

    Shats3 <- lapply(seq_along(Ahats),
                     function(lam) XX[[lam]] %*% mpinv(Ahats[[lam]]))

    X <- do.call(cbind, SimData$Xe[SimData$P==4])
    ica4 <- icafast(X = X, nc = rowd$Q)

    XX <- SimData$Xe[SimData$P==4]
    Ahats <- lapply(seq_along(XX),
                    function(lam) mpinv(ica4$S) %*% XX[[lam]])

    Shats4 <- lapply(seq_along(Ahats),
                     function(lam) XX[[lam]] %*% mpinv(Ahats[[lam]]))

    ShatsG2 <- c(Shats1, Shats2, Shats3, Shats4)
  }


  ##### Compute similarity matrices ######
  N <- 40
  comb <- t(utils::combn(1:N, 2))

  RVsS <- matrix(data = NA, nrow = N , ncol = N)
  RVS <- numeric()

  RVsSG2 <- matrix(data = NA, nrow = N , ncol = N)
  RVSG2 <- numeric()


  cat("Computing pairwise modified-RV statistics: \n")
  pb <- txtProgressBar(min = 0, max = nrow(comb), initial = 0)

  for(i in 1:nrow(comb)){

    RVS[i] <- modRV( Shats[[ comb[i,1] ]] , Shats[[ comb[i,2] ]])
    res <- c(comb[i , ] , RVS[i] )

    RVsS[res[1]  , res[2] ] <- res[3]

    RVSG2[i] <- modRV( ShatsG2[[ comb[i,1] ]] , ShatsG2[[ comb[i,2] ]])
    res <- c(comb[i , ] , RVSG2[i] )

    RVsSG2[res[1]  , res[2] ] <- res[3]


    setTxtProgressBar(pb, i)
  }

  RVsS[lower.tri(RVsS)] = t(RVsS)[lower.tri(RVsS)]
  diag(RVsS) <- 1

  RVsSG2[lower.tri(RVsSG2)] = t(RVsSG2)[lower.tri(RVsSG2)]
  diag(RVsSG2) <- 1

  cat('\n')

  ######## Custering procedures #######
  #### hclust ward ####
  RVmat <- as.dist(1-RVsS)
  RVmatG2 <- as.dist(1-RVsSG2)

  hcl <- hclust(d = RVmat, method = 'ward.D2')
  hclG2 <- hclust(d = RVmatG2, method = 'ward.D2')

  #plot(hcl)
  #plot(hclG2)
  k <- cutree(hcl, k = rowd$R)
  kG2 <- cutree(hclG2, k = rowd$R)


  hclARI <- adjustedRandIndex(SimData$P, k)
  hclG2ARI <- adjustedRandIndex(SimData$P, kG2)


  #### PAM ####
  pam <- pam(RVmat, k = rowd$R)
  pamG2 <- pam(RVmatG2, k = rowd$R)
  pamARI <- adjustedRandIndex(SimData$P, pam$clustering)
  pamARIG2 <- adjustedRandIndex(SimData$P, pamG2$clustering)


  #### kmeans #####
  XkmL <- lapply(seq_along(Shats),
                 function(lam) as.vector(Shats[[lam]]))

  Xkm <- do.call(rbind, XkmL)

  XkmLG2 <- lapply(seq_along(ShatsG2),
                   function(lam) as.vector(ShatsG2[[lam]]))

  XkmG2 <- do.call(rbind, XkmLG2)


  km <- kmeans(Xkm, centers = rowd$R)
  kmARI <- adjustedRandIndex(SimData$P, km$cluster)

  kmG2 <- kmeans(XkmG2, centers = rowd$R)
  kmG2ARI <- adjustedRandIndex(SimData$P, kmG2$cluster)

  ##### tandem analysis #####
  pca <- prcomp(x = Xkm, retx = T, rank. = rowd$Q, scale. = T)
  pcakm <- kmeans(pca$x, centers = rowd$R)
  pcakmARI <- adjustedRandIndex(SimData$P, pcakm$cluster)

  pcaG2 <- prcomp(x = XkmG2, retx = T, rank. = rowd$Q, scale. = T)
  pcakmG2 <- kmeans(pcaG2$x, centers = rowd$R)
  pcakmG2ARI <- adjustedRandIndex(SimData$P, pcakmG2$cluster)

  ###### Output #######
  clusterings <- data.frame(hcl = k, pam = pam$clustering,
                            km = km$cluster, pcakm = pcakm$cluster,
                            hclG2 = kG2, pamG2 = pamG2$clustering,
                            kmG2 = kmG2$cluster, pcakmG2 = pcakmG2$cluster)

  ARIs <- c(hclARI, pamARI, kmARI, pcakmARI,
            hclG2ARI, pamARIG2, kmG2ARI, pcakmG2ARI)
  names(ARIs) <- c('hcl','pam','km','pcakm',
                   'hclG2','pamG2','kmG2','pcakmG2')

  output <- list()
  output$clusterings <- clusterings
  output$ARIs <- ARIs



  setwd('~/Downloads/')

  save(output, file = paste("DualRegRV_rep",replication,"row",row,  "BulkData.Rdata",sep = "_"))
  cat('row ', row, '\n')

}
time <- proc.time() - tmp
save(time, file = 'timingfirst12.Rdata')


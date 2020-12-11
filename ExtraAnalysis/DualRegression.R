# Thu Dec 10 15:24:12 2020
# Author: Jeffrey Durieux, MSc


# What: script to do dual regression and clustering

library(ica)
library(CICA)
library(NMFN)
library(mclust)

source('~/Repositories/cica_simulation/Simulation1/Simulate_CICA.R')

### passed arguments:
# 1: replication
# 2: row (row of factorial simulation design)

args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)

replication <- args[1]
row <- args[2]

replication <- 3
row <- 55

# This will generate a unique seed
seed <- as.numeric(paste(replication,0,row,sep = ""))
set.seed(seed)

V <- c(500,2000)                # voxels
Q <- c(2,5,20)                  # components
R <- c(2,4)                     # clusters
D <- c('square', 'nonsquare')   # dimension mixing matrix
E <- c(0.05, 0.2, 0.4)          # noise level

design <- expand.grid(V=V, Q = Q,
                      R = R, D = D, E=E)


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

multiway::congru(Shats[[1]], Shats[[2]])

#RVmat <- CICA:::computeRVmat(DataList = Shats)
#CICA:::computeRVmat

##### Group ICA 2 #####
SimData$P
X <- do.call(cbind, SimData$Xe[SimData$P==1])
ica1 <- icafast(X = X, nc = rowd$Q)
X <- do.call(cbind, SimData$Xe[SimData$P==2])
ica2 <- icafast(X = X, nc = rowd$Q)


### dual regression step 1:
XX <- SimData$Xe[SimData$P==1]
Ahats <- lapply(seq_along(XX),
                function(lam) mpinv(ica1$S) %*% XX[[lam]])

### dual regression step 2:

Shats1 <- lapply(seq_along(Ahats),
                function(lam) XX[[lam]] %*% mpinv(Ahats[[lam]]))

### dual regression step 1:
XX <- SimData$Xe[SimData$P==2]
Ahats <- lapply(seq_along(XX),
                function(lam) mpinv(ica2$S) %*% XX[[lam]])

### dual regression step 2:

Shats2 <- lapply(seq_along(Ahats),
                 function(lam) XX[[lam]] %*% mpinv(Ahats[[lam]]))

Shats <- c(Shats1,Shats2)


###### hcl on dist #####
N <- length(Shats)
comb <- t(utils::combn(1:N, 2))
RVsS <- matrix(data = NA, nrow = N , ncol = N)
RVS <- numeric()

cat("Computing pairwise modified-RV statistics: \n")
pb <- txtProgressBar(min = 0, max = nrow(comb), initial = 0)

for(i in 1:nrow(comb)){
  RVS[i] <- mean(diag(Tucker( Shats[[ comb[i,1] ]] , Shats[[ comb[i,2] ]])))
  #RVS[i] <- mean(diag(cor( Shats[[ comb[i,1] ]] , Shats[[ comb[i,2] ]])))
  res <- c(comb[i , ] , RVS[i] )

  RVsS[res[1]  , res[2] ] <- res[3]
  setTxtProgressBar(pb, i)

}

RVsS[lower.tri(RVsS)] = t(RVsS)[lower.tri(RVsS)]
diag(RVsS) <- 1
cat('\n')

RVmat <- as.dist(1-RVsS)
hcl <- hclust(d = RVmat, method = 'ward.D2')
plot(hcl)
k <- cutree(hcl, k = rowd$R)
hclARI <- adjustedRandIndex(SimData$P, k)

#### pam #####
library(cluster)
pam <- pam(RVmat, k = rowd$R)
pam$clustering
pamARI <- adjustedRandIndex(SimData$P, pam$clustering)

#### kmeans ####
XkmL <- lapply(seq_along(Shats),
                function(lam) as.vector(Shats[[lam]]))

Xkm <- do.call(rbind, XkmL)

km <- kmeans(Xkm, centers = rowd$R)
kmARI <- adjustedRandIndex(SimData$P, km$cluster)


##### tandem analysis #####
pca <- prcomp(x = Xkm, retx = T, rank. = rowd$Q, scale. = T)
pcakm <- kmeans(pca$x, centers = rowd$R)
pcakmARI <- adjustedRandIndex(SimData$P, pcakm$cluster)

#### reduced k-means #####
#library(clustrd)
### very slow
#rdkm <- cluspca(data = Xkm, nclus = rowd$R, ndim = rowd$Q, method = 'RKM')
print(hclARI)
print(pamARI)
print(kmARI)
print(pcakmARI)
des$ARI[1+72+72]
rowd
id <- which(des$V==2000&des$Q==20&des$R==2&des$D=='nonsquare'&des$E==.40)
des$ARI[id]



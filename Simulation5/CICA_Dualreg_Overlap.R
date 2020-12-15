# Tue Dec 15 14:27:35 2020
# Author: Jeffrey Durieux, MSc


## C-ICA extra simulation, compare C-ICA with dual regression for difficult conditions

## datasets are generated with cluster overlap. Such as done in paper 1

library(CICA)
library(ica)
library(mclust)
library(cluster)
library(NMFN)

source('~/Repositories/PhDthesis/P1/SpecialIssue/SPECIALISSUE/Files/Scripts/SimulateData.R')
source('~/Repositories/cica_simulation/SharkCode/CICA_datagen.R')

# sim design paper 1
n_clusters <- c(2,4)
c_size <- c("equal", "unequal")
con <- c(.08, .12, .15, .23, .395) # smaller means more overlap (.99 versus .75)
error <- c(0.1, 0.3, 0.6, 0.8)

design <- expand.grid(n_clusters = n_clusters, c_size = c_size,
                      con = con, error = error)

##### sim testdata #####
set.seed(2407)
set.seed(240702)
testdata <- SimData(n_clusters = 2, c_size = 'equal', con = .08, error = 0.6)



set.seed(2407)

testdata$P
rat <- list()
rat$rationalstarts <- data.frame(testdata$P, k)
class(rat) <- 'rationalstarts'

set.seed(2407)
CICA <- CICA(DataList = testdata$Xe, nStarts = 1, nComp = 20,nClus = 2, scale = T, center = T, rational = rat)

save(CICA, file = '~/Downloads/CICA_clusteroverlap_seed2407_2R_equal_con15_E06.Rdata')
CICA$P

adjustedRandIndex(testdata$P, CICA$P)
Tucker(testdata$S[[1]], CICA$Sr[[1]]) %>% abs %>% apply(MARGIN = 2, FUN = max) %>% mean
Tucker(testdata$S[[2]], CICA$Sr[[2]]) %>% abs %>% apply(MARGIN = 2, FUN = max) %>% mean


### GICA ####
set.seed(2407)
X <- do.call(testdata$Xe, what = cbind)

GICA <- icafast(X, nc = 20)

XX <- testdata$Xe
Ahats <- lapply(seq_along(XX),
                function(lam) mpinv(GICA$S) %*% XX[[lam]])

### dual regression step 2:

Shats <- lapply(seq_along(Ahats),
                function(lam) XX[[lam]] %*% mpinv(Ahats[[lam]]))

Tucker(testdata$S[[2]], Shats[[1]]) %>% abs %>% apply(MARGIN = 2, FUN = max) %>% mean

Tucker(testdata$S[[1]], Shats[[31]]) %>% abs %>% apply(MARGIN = 2, FUN = max) %>% mean


####### pairwise tuckers #####
N <- 60
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

library(plotly)

######## Custering procedures #######
#### hclust ward ####

RVmat <- as.dist(1-RVsS)

mds <- cmdscale(RVmat)
plot_ly(data.frame(mds), x= ~X1, y= ~X2, color = as.factor(testdata$P))

hcl <- hclust(d = RVmat, method = 'ward.D2')
plot(hcl)
k <- cutree(hcl, k = 2)
adjustedRandIndex(testdata$P, k)

plot_ly(data.frame(mds), x= ~X1, y= ~X2, color = as.factor(k))


pam <- pam(RVmat, k = 2)
adjustedRandIndex(testdata$P, pam$clustering)

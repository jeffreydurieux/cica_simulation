# Tue Dec 15 14:27:35 2020
# Author: Jeffrey Durieux, MSc


## C-ICA extra simulation, compare C-ICA with dual regression for difficult conditions

## datasets are generated with cluster overlap. Such as done in paper 1

library(CICA)
library(ica)
library(mclust)
library(cluster)
library(NMFN)
library(plotly)

source('~/Repositories/cica_simulation/Simulation5/SimData.R')
source('~/Repositories/cica_simulation/SharkCode/CICA_datagen.R')
source('~/Repositories/cica_simulation/Simulation5/GICADR_clustering.R')
source('~/Repositories/cica_simulation/Simulation5/ShuffleP4.R')
# sim design paper 1

n_clusters <- c(2,4)
con <- c(.08, .12, .15, .23, .395) # smaller means more overlap (.99 versus .75)
error <- c(0.6, 0.8)

design <- expand.grid(n_clusters = n_clusters,
                      con = con, error = error)

##### sim testdata #####
set.seed(2407)
testdata <- SimData(Nsamp = 40, n_clusters = 4, con = .395, error = 0.7)

shuffle <- shuffleR4(testdata$Xe, R = 4)

testdata$Xe <- shuffle$newX
testdata$P <- shuffle$newP
#### permute data ####

### GICA dual reg clustering ####
DR <- DRCLUS(testdata, Q = 5, R = 4)

mds <- cmdscale(1-DR$RV,k = 3)
cc <- as.factor(testdata$P)
#cc <- as.factor(CICA$rationalstart$P)
plot_ly(data = data.frame(mds), x = ~X1, y = ~X2, color = cc)

DR$clustering$khcl
DR$clustering$kpam

DR$clustering$khcl %>% adjustedRandIndex(testdata$P)
DR$clustering$kpam %>% adjustedRandIndex(testdata$P)

set.seed(2407)
testdata$P
rat <- list()
rat$rationalstarts <- data.frame(testdata$P,CICA:::perturbation(testdata$P), DR$clustering$khcl,DR$clustering$kpam)

rat$rationalstarts <- data.frame(testdata$P,
                                 CICA:::perturbation(testdata$P,percentage = .05),
                                 CICA:::perturbation(testdata$P,percentage = .1),
                                 CICA:::perturbation(testdata$P,percentage = .15),
                                 CICA:::perturbation(testdata$P,percentage = .20),
                                 DR$clustering$khcl,DR$clustering$kpam)


#rat$rationalstarts <- data.frame(testdata$P,
#                                 CICA:::perturbation(testdata$P))

class(rat) <- 'rationalstarts'

set.seed(2407)
CICA <- CICA(DataList = testdata$Xe, nStarts = 1, nComp = 5,nClus =4, scale = T, center = T, rational = rat)


#save(CICA, file = '~/Downloads/CICA_clusteroverlap_seed2407_2R_equal_con15_E06.Rdata')
CICA$rationalstart$P
CICA$rationalstart$P %>% table(testdata$P)
CICA$rationalstart$P %>% adjustedRandIndex(testdata$P)

CICA$rationalstart$P %>% table(DR$clustering$khcl)
CICA$rationalstart$P %>% adjustedRandIndex(DR$clustering$khcl)

CICA$rationalstart$P %>% table(DR$clustering$kpam)
CICA$rationalstart$P %>% adjustedRandIndex(DR$clustering$kpam)

CICA$P
CICA$P %>% table(testdata$P)
CICA$P %>% adjustedRandIndex(testdata$P)

CICA$P %>% table(DR$clustering$khcl)
CICA$P %>% adjustedRandIndex(DR$clustering$khcl)

CICA$P %>% table(DR$clustering$kpam)
CICA$P %>% adjustedRandIndex(DR$clustering$kpam)




CICA$rationalstart$P
Tucker(testdata$S[[1]], CICA$rationalstart$Sr[[1]]) %>% abs %>% apply(MARGIN = 2, FUN = max) %>% mean
Tucker(testdata$S[[2]], CICA$rationalstart$Sr[[2]]) %>% abs %>% apply(MARGIN = 2, FUN = max) %>% mean
Tucker(testdata$S[[3]], CICA$rationalstart$Sr[[3]]) %>% abs %>% apply(MARGIN = 2, FUN = max) %>% mean
Tucker(testdata$S[[4]], CICA$rationalstart$Sr[[4]]) %>% abs %>% apply(MARGIN = 2, FUN = max) %>% mean




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

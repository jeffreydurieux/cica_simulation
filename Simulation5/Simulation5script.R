# Fri Dec 18 16:25:32 2020
# Author: Jeffrey Durieux, MSc

# Simulation study C-ICA vs Dual Regression vs Two Step.
# Using overlapping clusters

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


#n_clusters <- c(2,4)
#con <- c(.08, .12, .15, .23, .395) # smaller means more overlap (.99 versus .75)
con <- c(.08, .12, .15, .23, .395) # smaller means more overlap (.99 versus .75)
error <- c(0.6, 0.8)
rep <- 1:10

design <- expand.grid(con = con, error = error)

### expecting interaction effect between overlap and error
### error = .50 vs .75
### con = .08 vs. .15 


set.seed(2)
testdata <- SimData(Q = 5,nscan = 50, Nsamp = 40, n_clusters = 4, con = .15, error = 0.5)

shuffle <- shuffleR4(testdata$Xe, R = 4)

testdata$Xe <- shuffle$newX
testdata$P <- shuffle$newP

DR <- DRCLUS(testdata, Q = 5, R = 4)
DR$clustering

mds <- cmdscale(1-DR$RV, k = 3)
cc <- as.factor(testdata$P)
#cc <- as.factor(CICA$rationalstart$P)
plot_ly(data = data.frame(mds), x = ~X1, y = ~X2, color = cc)
plot_ly(data = data.frame(mds), x = ~X1, y = ~X2, color = as.factor(DR$clustering$khcl))
plot_ly(data = data.frame(mds), x = ~X1, y = ~X2, color = as.factor(DR$clustering$kpam))



rat <- list()
rat$rationalstarts <- data.frame(testdata$P,CICA:::perturbation(testdata$P), DR$clustering$khcl,DR$clustering$kpam)

rat$rationalstarts <- data.frame(testdata$P,
                                 #DR$clustering$khcl,DR$clustering$kpam,
                                 #CICA:::perturbation(testdata$P,percentage = .05),
                                 #CICA:::perturbation(testdata$P,percentage = .1),
                                 CICA:::perturbation(testdata$P,percentage = .1),
                                 CICA:::perturbation(testdata$P,percentage = .33))
                                 


#rat$rationalstarts <- data.frame(testdata$P,
#                                 CICA:::perturbation(testdata$P))

class(rat) <- 'rationalstarts'

CICA <- CICA(DataList = testdata$Xe, nStarts = 1, nComp = 5,nClus =4, scale = T, center = T, rational = rat)

# check if rational start is perfect, then other comparisons are ok
CICA$rationalstart$P
CICA$rationalstart$P %>% table(testdata$P)
CICA$rationalstart$P %>% adjustedRandIndex(testdata$P)

CICA$rationalstart$P %>% table(DR$clustering$khcl)
CICA$rationalstart$P %>% adjustedRandIndex(DR$clustering$khcl)

CICA$rationalstart$P %>% table(DR$clustering$kpam)
CICA$rationalstart$P %>% adjustedRandIndex(DR$clustering$kpam)


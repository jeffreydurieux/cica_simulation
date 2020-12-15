# Tue Dec 15 10:00:40 2020
# Author: Jeffrey Durieux, MSc

# What: script for adding modRV matrices from shark
# modRV computed between Group-ICA on subselection dataset simulation 3/4
# C-ICA paper

setwd('~/Repositories/cica_simulation/Data/CICA_SIM_4/modRV/')

library(gtools)

files <- dir()
files <- mixedsort(files)

RV1 <- get(load(files[1]))
RV2 <- get(load(files[2]))
RV3 <- get(load(files[3]))
RV4 <- get(load(files[4]))
RV5 <- get(load(files[5]))
RV6 <- get(load(files[6]))
RV7 <- get(load(files[7]))
RV8 <- get(load(files[8]))

RV1[is.na(RV1)] <- 0
RV2[is.na(RV2)] <- 0
RV3[is.na(RV3)] <- 0
RV4[is.na(RV4)] <- 0
RV5[is.na(RV5)] <- 0
RV6[is.na(RV6)] <- 0
RV7[is.na(RV7)] <- 0
RV8[is.na(RV8)] <- 0

RV <- RV1 + RV2 + RV3 + RV4 + RV5 + RV6 + RV7 + RV8
diag(RV) <- 1

heatmap(RV)

labels <- c(rep('EC',20), rep('AD',20))

RVmat <- as.dist(1-RV)
RVmat <- as.dist(RV)

# change to centroid for inversion example
# change to single to see trailing
# visually inspect complete--> interesting separation

col <- ifelse(labels=='EC', 'blue', 'red')
mds <- cmdscale(RVmat, k = 3)
plotly::plot_ly(data = data.frame(mds), 
                x = ~X1, y = ~X2, z = ~X3, color = labels)

plotly::plot_ly(data = data.frame(mds), 
                x = ~X1, y = ~X2, color = labels)


hcl <- hclust(d = RVmat, method = 'ward.D2')
plot(hcl)
rect.hclust(hcl, k = 2)
k <- cutree(hcl, k = 2)

plotly::plot_ly(data = data.frame(mds), 
                x = ~X1, y = ~X2, color = as.factor(k), text = 1:40)

plotly::plot_ly(data = data.frame(mds), 
                x = ~X1, y = ~X2, z = ~X3,color = as.factor(k), text = 1:40)

mclust::adjustedRandIndex(labels,k)
addmargins(table(labels,k))

pam <- cluster::pam(RVmat, k = 2)
mclust::adjustedRandIndex(labels,pam$clustering)
addmargins(table(labels,pam$clustering))


#### compare to CICA #####
load('../CICA_nc20Clus2.Rdata')
cica$P
addmargins(table(labels,cica$P))


cbind(labels,cica$P,k)

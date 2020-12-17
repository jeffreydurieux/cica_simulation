# Tue Dec 15 10:00:40 2020
# Author: Jeffrey Durieux, MSc

# What: script for adding modRV matrices from shark
# modRV computed between Group-ICA on subselection dataset simulation 3/4
# C-ICA paper

### Contents ####
# First part is modRV single subjects ICAs #
# Second part is modRV based on GICA       #
# Comparison with CICA is made #

library(NMFN) #mpinv

#### Part 1: single subject ICA #####

setwd('~/Repositories/cica_simulation/Data/CICA_SIM_4/modRV/')

library(gtools)

filess <- dir()
filess <- mixedsort(filess)

files <- filess[1:8]
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
#RVmat <- as.dist(RV)

# change to centroid for inversion example
# change to single to see trailing
# visually inspect complete--> interesting separation

col <- ifelse(labels=='EC', 'blue', 'red')
mds <- cmdscale(RVmat, k = 3)

#### HCL ####
hcl <- hclust(d = RVmat, method = 'ward.D2')
plot(hcl)
rect.hclust(hcl, k = 2)
k <- cutree(hcl, k = 2)
khcl1 <- k
mclust::adjustedRandIndex(labels,k)
addmargins(table(labels,k))

hclp1.1 <- plotly::plot_ly(data = data.frame(mds), 
                x = ~X1, y = ~X2, color = labels,
                colors = 'Pastel1', text = 1:40)

hclp2.1 <- plotly::plot_ly(data = data.frame(mds), 
                         x = ~X1, y = ~X2, 
                         color = as.factor(k),
                         colors = 'Dark2', text = 1:40)

phcl_singleICA_RV <- plotly::subplot(hclp1.1,hclp2.1)

phcl_singleICA_RV

#### PAM ####
pam <- cluster::pam(RVmat, k = 2)
mclust::adjustedRandIndex(labels,pam$clustering)
kpam1 <- pam$clustering
addmargins(table(labels,pam$clustering))

pamp1.1 <- plotly::plot_ly(data = data.frame(mds), 
                         x = ~X1, y = ~X2, color = labels,
                         colors = 'Pastel1', text = 1:40)

pamp2.1 <- plotly::plot_ly(data = data.frame(mds), 
                         x = ~X1, y = ~X2, 
                         color = as.factor(pam$clustering),
                         colors = 'Dark2', text = 1:40)

ppam_singleICA_RV <- plotly::subplot(pamp1.1,pamp2.1)

ppam_singleICA_RV


#### Part 2: Dual Reg #####

filess <- dir()
filess <- mixedsort(filess)
filess
files <- filess[9:16]
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
#RVmat <- as.dist(RV)

# change to centroid for inversion example
# change to single to see trailing
# visually inspect complete--> interesting separation

col <- ifelse(labels=='EC', 'blue', 'red')
mds <- cmdscale(RVmattuck, k = 3)
plotly::plot_ly(data = data.frame(mds), 
                x = ~X1, y = ~X2, z = ~X3, color = labels)

plotly::plot_ly(data = data.frame(mds), 
                x = ~X1, y = ~X2, color = labels)

#### HCL ####
hcl <- hclust(d = RVmat, method = 'ward.D2')
plot(hcl)
rect.hclust(hcl, k = 2)
k <- cutree(hcl, k = 2)
khcl2 <- k
mclust::adjustedRandIndex(labels,k)
addmargins(table(labels,k))


hclp1.2 <- plotly::plot_ly(data = data.frame(mds), 
                         x = ~X1, y = ~X2, color = labels,
                         colors = 'Pastel1', text = 1:40)

hclp2.2 <- plotly::plot_ly(data = data.frame(mds), 
                         x = ~X1, y = ~X2, 
                         color = as.factor(k),
                         colors = 'Dark2', text = 1:40)

phcl_GICA_RV <- plotly::subplot(hclp1.2,hclp2.2)

phcl_GICA_RV

#### PAM ####
pam <- cluster::pam(RVmat, k = 2)
mclust::adjustedRandIndex(labels,pam$clustering)
kpam2 <- pam$clustering
addmargins(table(labels,pam$clustering))

pamp1.2 <- plotly::plot_ly(data = data.frame(mds), 
                         x = ~X1, y = ~X2, color = labels,
                         colors = 'Pastel1', text = 1:40)

pamp2.2 <- plotly::plot_ly(data = data.frame(mds), 
                         x = ~X1, y = ~X2, 
                         color = as.factor(pam$clustering),
                         colors = 'Dark2', text = 1:40)

ppam_GICA_RV <- plotly::subplot(pamp1.2,pamp2.2)

ppam_GICA_RV

#### compare to CICA #####
load('~/Repositories/cica_simulation/Data/CICA_SIM_4/CICA_nc20Clus2.Rdata')
cica$P
addmargins(table(labels,cica$P))


#### CICA DUAL REG ####
setwd('~/Repositories/cica_simulation/Data/CICA_SIM_4/dataTV/')

files <- dir()
files <- files[-41]
files

dataL <- list()
for(i in 1:40){
  x <- get(load(files[i]))
  dataL[[i]] <- x
}

XX <- lapply(dataL, t)

XX1 <- XX[cica$P==1]
Ahats1 <- lapply(seq_along(XX1),
                function(lam) mpinv(cica$Sr[[1]]) %*% XX1[[lam]])

XX2 <- XX[cica$P==2]
Ahats2 <- lapply(seq_along(XX2),
                 function(lam) mpinv(cica$Sr[[2]]) %*% XX2[[lam]])

### dual regression step 2:

Shats1 <- lapply(seq_along(Ahats1),
                function(lam) XX1[[lam]] %*% mpinv(Ahats1[[lam]]))

Shats2 <- lapply(seq_along(Ahats2),
                 function(lam) XX2[[lam]] %*% mpinv(Ahats2[[lam]]))

#### CICA P is already in order in this specific case ###
# Therefore Shats1 and Shats2 can be concatenated

ShatsCICA <- c(Shats1, Shats2)

N <- 40
comb <- t(utils::combn(1:N, 2))


RVsS <- matrix(data = NA, nrow = N , ncol = N)
RVS <- numeric()

cat("Computing pairwise Tucker statistics: \n")
pb <- txtProgressBar(min = 0, max = nrow(comb), initial = 0)

for(i in 1:nrow(comb)){
  
  if( (comb[i,1] < 21 & comb[i,2] < 21) == TRUE ){
    RVS[i] <- mean(diag(Tucker( ShatsCICA[[ comb[i,1] ]] , ShatsCICA[[ comb[i,2] ]])))
  }
  
  else if( (comb[i,1] >21 & comb[i,2] > 21) == TRUE ){
    RVS[i] <- mean(diag(Tucker( ShatsCICA[[ comb[i,1] ]] , ShatsCICA[[ comb[i,2] ]])))
  }else{
    RVS[i] <- mean(apply(abs(Tucker( ShatsCICA[[ comb[i,1] ]] , ShatsCICA[[ comb[i,2] ]])),MARGIN = 2, max))  
  }

  
  res <- c(comb[i , ] , RVS[i] )
  
  RVsS[res[1]  , res[2] ] <- res[3]

  setTxtProgressBar(pb, i)
}


RVsS[lower.tri(RVsS)] = t(RVsS)[lower.tri(RVsS)]
diag(RVsS) <- 1

RVmat <- as.dist(1-RVsS)

mds <- cmdscale(RVmat, k = 3)
plotly::plot_ly(data = data.frame(mds), 
                x = ~X1, y = ~X2, z = ~X3, color = labels)

cicap1 <- plotly::plot_ly(data = data.frame(mds), 
                x = ~X1, y = ~X2, color = labels, 
                colors = 'Pastel1',  text = ~1:40)

cicap2 <- plotly::plot_ly(data = data.frame(mds), 
                x = ~X1, y = ~X2, color = as.factor(cica$P),
                colors = 'Dark2' , text = ~1:40)

CICATuckpairwise <- plotly::subplot(cicap1,cicap2)

CICATuckpairwise 

#### output to make markdown report
reportoutput <- list()

reportoutput$figures$SingleICAs$hcl <- phcl_singleICA_RV

reportoutput$figures$SingleICAs$pam <- ppam_singleICA_RV

reportoutput$figures$GICA$hcl <- phcl_GICA_RV
reportoutput$figures$GICA$pam <- ppam_GICA_RV

reportoutput$figures$CICA <- CICATuckpairwise
reportoutput$figures$SingleICAs$hcl

clustering <- data.frame(labels = labels,
                         cica = cica$P,
                         khcl1 = khcl1,
                         kpam1 = kpam1,
                         khcl2 = khcl2,
                         kpam2 = kpam2)

reportoutput$clustering <- clustering

##### nice dendrograms #####
# http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning#plot.dendrogram-function

library(ape) #phlyo trees

colors = c("red", "blue")
hcl <- reportoutput$hcl$hcl2
clus2 = cutree(hcl, 2)
hcl2 <- plot(as.phylo(hcl), type = "fan", tip.color = colors[clus2],
     label.offset = .0001, cex = 1)

# manually add hcl objects 
reportoutput$hcl$hcl2 <- hcl
plot(reportoutput$hcl$hcl2)
rect.hclust(reportoutput$hcl$hcl2, k = 2)


### save report output ###
save(reportoutput, file = '../../Reports/reportoutput_subselectionGICAvsCICAvsTwoStep.Rdata')

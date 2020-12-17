# Tue Dec 15 10:00:40 2020
# Author: Jeffrey Durieux, MSc

# What: script for adding Tuck matrices from shark
# Tuck computed between Group-ICA on subselection dataset simulation 3/4
# C-ICA paper
# This does not work properly, taking the column max abs tucker value is arbitrary, see it by changing the columns of comb and selecting the lower.tri vs the upper.tri


### Contents ####
# First part is Tuck single subjects ICAs #
# Second part is Tuck based on GICA       #
# Comparison with CICA is made #

library(NMFN) #mpinv
library(CICA)
library(gtools)
#### Part 1: single subject ICA #####

setwd('~/Repositories/cica_simulation/Data/CICA_SIM_4/')
load(file = 'SingleICAsShats.Rdata')

Shats <- c(Shats[21:40], Shats[1:20])

N <- 40
comb <- t(utils::combn(1:N, 2)) # take this with lower tri line

comb <- cbind(comb[,2],comb[,1]) # take this with upper tri

RVsS <- matrix(data = NA, nrow = N , ncol = N)
RVS <- numeric()

cat("Computing pairwise Tucker statistics: \n")
pb <- txtProgressBar(min = 0, max = nrow(comb), initial = 0)


Tucker( Shats[[ comb[i,1] ]] , Shats[[ comb[i,2] ]]) %>% abs %>% apply(MARGIN = 2,FUN = max) %>% mean

signmax <- function(x){
  a <- range(x)
  index <- which.max(abs(a))
  return(a[index])
}


for(i in 1:nrow(comb)){
  
  #RVS[i] <- mean(diag(Tucker( Shats[[ comb[i,1] ]] , Shats[[ comb[i,2] ]])))
  RVS[i] <- mean(apply(abs(Tucker( Shats[[ comb[i,1] ]] , Shats[[ comb[i,2] ]])),MARGIN = 2, max))
  
  res <- c(comb[i , ] , RVS[i] )
  
  RVsS[res[1]  , res[2] ] <- res[3]
  
  setTxtProgressBar(pb, i)
}

#RVsS[upper.tri(RVsS)] = t(RVsS)[upper.tri(RVsS)]
RVsS[lower.tri(RVsS)] = t(RVsS)[lower.tri(RVsS)]
diag(RVsS) <- 1

RVmat <- as.dist(1-RVsS)

labels <- c(rep('EC',20), rep('AD',20))

#RVmat <- as.dist(RV)

# change to centroid for inversion example
# change to single to see trailing
# visually inspect complete--> interesting separation

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

rm(Shats)
load(file = 'ShatsGICADUAL.Rdata')

N <- 40
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

RVmat <- as.dist(1-RVsS)



# change to centroid for inversion example
# change to single to see trailing
# visually inspect complete--> interesting separation

mds <- cmdscale(RVmat, k = 3)

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
  
  RVS[i] <- mean(diag(Tucker( ShatsCICA[[ comb[i,1] ]] , ShatsCICA[[ comb[i,2] ]])))
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
reportoutputtuck <- list()

reportoutputtuck$figures$SingleICAs$hcl <- phcl_singleICA_RV

reportoutputtuck$figures$SingleICAs$pam <- ppam_singleICA_RV

reportoutputtuck$figures$GICA$hcl <- phcl_GICA_RV
reportoutputtuck$figures$GICA$pam <- ppam_GICA_RV

reportoutputtuck$figures$CICA <- CICATuckpairwise
reportoutputtuck$figures$SingleICAs$hcl

clusteringtuck <- data.frame(labels = labels,
                         cica = cica$P,
                         khcl1 = khcl1,
                         kpam1 = kpam1,
                         khcl2 = khcl2,
                         kpam2 = kpam2)

reportoutputtuck$clustering <- clusteringtuck

##### nice dendrograms #####
# http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning#plot.dendrogram-function

library(ape) #phlyo trees

colors = c("red", "blue")
hcl <- hcl
clus2 = cutree(hcl, 2)
hcl2 <- plot(as.phylo(hcl), type = "fan", tip.color = colors[clus2],
             label.offset = .0001, cex = 1)

# manually add hcl objects 
reportoutputtuck$hcl$hcl2 <- hcl
plot(reportoutputtuck$hcl$hcl2)
rect.hclust(reportoutputtuck$hcl$hcl2, k = 2)


### save report output ###
save(reportoutputtuck, file = '~/Repositories/cica_simulation/Data/Reports/reportoutputtuck_subselectionGICAvsCICAvsTwoStep.Rdata')

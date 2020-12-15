# Mon Dec 14 11:12:16 2020
# Author: Jeffrey Durieux, MSc

# G_ICA dual regression on empirical data

setwd('~/Repositories/cica_simulation/Data/CICA_SIM_4/dataTV/')

files <- dir()
files <- files[-41]
files

labels <- c(rep('EC',20), rep('AD',20))


dataL <- list()
for(i in 1:40){
  x <- get(load(files[i]))
  dataL[[i]] <- x
}

dataL <- lapply(dataL, t)

#### GICA dual regression selecting for 20 components #####
library(ica)
library(CICA)
library(NMFN)
library(mclust)
library(cluster)

# Concatenate data
X <- do.call(cbind, dataL)
dim(X)

set.seed(2407)
GICA <- icafast(X = X, nc = 20)

load('~/Repositories/cica_simulation/Data/CICA_SIM_4/GICA20.Rdata')

### dual regression step 1:
XX <- dataL
Ahats <- lapply(seq_along(XX),
                function(lam) mpinv(GICA$S) %*% XX[[lam]])

### dual regression step 2:

Shats <- lapply(seq_along(Ahats),
                function(lam) XX[[lam]] %*% mpinv(Ahats[[lam]]))

save(Shats, file = '../ShatsGICADUAL.Rdata')

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
cat('\n')

library(plotly)

######## Custering procedures #######
#### hclust ward ####

RVmat <- as.dist(1-RVsS)

#### MDS####
mds <- cmdscale(RVmat)
plot_ly(data.frame(mds), x= ~X1, y= ~X2, color = labels)

hcl <- hclust(d = RVmat, method = 'ward.D2')
plot(hcl)

k <- cutree(hcl, k = 2)
p1 <- plot_ly(data.frame(mds), x= ~X1, y= ~X2, color = as.factor(k), colors = c('blue','red'))

labels
p2 <- plot_ly(data.frame(mds), x= ~X1, y= ~X2, color = labels)
subplot(p1,p2)

hclARI <- adjustedRandIndex(labels, k)

#### PAM ####
pam <- pam(RVmat, k = 2)
pamARI <- adjustedRandIndex(labels, pam$clustering)

addmargins(table(labels,k))
addmargins(table(labels,pam$clustering))

load('~/Repositories/cica_simulation/Data/CICA_SIM_4/CICA_nc20Clus2.Rdata')

cbind(cica$P, k,pam$clustering)


# modRV procedure paper 1
# fails due to vector memory exhaustion
# done on shark, scripts located in simulation 4 folder

#rat <- FindRationalStarts(DataList = dataL, nComp = 20, nClus = 2,pseudo = F)



Shats <- lapply(dataL, FUN = icafast, nc = 20)
Shats <- lapply(seq_along(Shats), function(lam) Shats[[lam]]$S)

#save Shats for cluster
save(Shats, file = '../SingleICAsShats.Rdata')


##### MDS picture of unreduced data #####

RVsX <- matrix(data = NA, nrow = N , ncol = N)
RVS <- numeric()

for(i in 1:nrow(comb)){

  RVS[i] <- mean(diag(Tucker( XX[[ comb[i,1] ]] , XX[[ comb[i,2] ]])))
  res <- c(comb[i , ] , RVS[i] )

  RVsX[res[1]  , res[2] ] <- res[3]

  setTxtProgressBar(pb, i)
}

RVsX[lower.tri(RVsX)] = t(RVsX)[lower.tri(RVsX)]
diag(RVsX) <- 1
cat('\n')

save(RVsX, file = '../RVX_datasubselection.Rdata')

RVXmat <- as.dist(1-RVsX)

#### MDS####
mds <- cmdscale(RVXmat)
p3 <- plot_ly(data.frame(mds), x= ~X1, y= ~X2, color = labels, text = 1:40)


#### final MDS plots #####
ax <- list(
  title = "",
  zeroline = TRUE,
  showline = TRUE,
  showticklabels = F,
  showgrid = TRUE
)

p3 %>% layout(xaxis = ax, yaxis = ax)
subplot(p2 %>% layout(xaxis = ax, yaxis = ax),
        p1 %>% layout(xaxis = ax, yaxis = ax))

p4 <- plot_ly(data.frame(mds), x= ~X1, y= ~X2, color = as.factor(cica$P), colors = c('red','blue'), text = 1:40)

subplot(p2 %>% layout(xaxis = ax, yaxis = ax),
        p4 %>% layout(xaxis = ax, yaxis = ax))



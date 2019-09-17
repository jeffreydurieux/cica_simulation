# Wed Jul 10 10:57:09 2019
# Author: Jeffrey Durieux, MSc

# What: R-script for the sequential scree model selection
# this script is not generally applicable, it is interactive and dependent on the data that you need to check. If you want to use it you need to change some input parameters like: totalSSQ, ncomp and nclus

######## X #######
setwd("/Volumes/LaCie/MyData/CICA/Project1/AllDataResults/dataTV/")
files <- dir()
X <- list()
for(i in 1:250){
  X[[i]] <- get(load(files[i]))
}
X <- lapply(X,t)
X <- lapply(X =X, FUN = CICA::xscale)
X <- CICA::ConcDataFun(DataList = X, ClusVec = rep(1,250))


###### compute loss of GICA ######
gicaloss <- numeric()
nc <- c(20,25,30,35)
for(i in 1:4){
  load(paste("/Volumes/LaCie/MyData/CICA/Project1/AllDataResults/ResultsFromShark/GICA_nc",nc[i],".Rdata", sep = ""))
  Xhat <- gica$S %*% t(gica$M)
  gicaloss[i] <- sum((X[[1]]-Xhat)^2)
}
options(scipen=999)
gicaloss


library(plotly)

## sub selection 40 * 1000
## change this value depending on the data
totalSSQ <- 250000




############# Functions ##############

# convert to VAF
VAF <- function(Loss, totalSSQ){
  vaf <- (totalSSQ - Loss) / totalSSQ
  return(vaf*100)
}


############# Model Selection Prep ##############
ncomp <- c(20,25,30,35)
nclus <- c(2,3,4,5)

grid <- expand.grid(n.comp = ncomp, n.clus = nclus)

setwd("/Volumes/LaCie/MyData/CICA/Project1/AllDataResults/ResultsFromShark/")

files <- dir(pattern = "CICA")
files <- files[1:16]

loss <- numeric()
for(i in 1:16){
  vals <- grid[i,]
  nam <- paste("full_CICA_nc", vals$n.comp, "Clus", vals$n.clus, sep = "") 
  nam <- paste(nam,".Rdata", sep="")
  cica <- get( load(nam))
  loss[i] <- cica$LossHistory[cica$Iterations]
}
loss <- unlist(loss)
loss <- c(gicaloss,loss)


VAFS <- sapply(loss, VAF, totalSSQ=totalSSQ)
VAFS

VAFdata <- data.frame(grid,loss=loss,VAF=VAFS)



#step 1 of sequentials: calculate scree ratio for each value of K (cluster) given different diff values of Q (components)
#Then average all the SR_K over all components. and select the K that yields the highest average

########## step 1 function #######
# compute scree ratio for a certain number of clusters r fixing for q components

# ( Lr-1,q - Lr,q)  /  (Lr,q - Lr+1, q )

# Lq is the loss vector of length K, screes conditional on Q
SR_rq <- function(Lq){
  
  Screes <- numeric()
  
  for(i in 1:length(Lq)){
    if(i == 1){
      Screes[i] <- NA
    }
    else if(i == length(Lq) ){
      Screes[i] <- NA
    }else{
      Screes[i] <- (Lq[i-1] - Lq[i]) / (Lq[i] - Lq[i+1])
    }
  }
  return(Screes)
}


#### compute step 1: calculate first step formula for all compentens (i.e., screeratio_r | q)
#### average over this matrix, the position of the highest average is the optimal number of clusters

Screes <- rbind(SR_rq(VAFdata[VAFdata$n.comp==20,]$loss) ,
                SR_rq(VAFdata[VAFdata$n.comp==25,]$loss) ,
                SR_rq(VAFdata[VAFdata$n.comp==30,]$loss) ,
                SR_rq(VAFdata[VAFdata$n.comp==35,]$loss) )
                
                

Screes
colMeans(Screes, na.rm = TRUE)


########### step 2 function #############
# after computing the optimal number of clusters (see above)
# compute the optimal number of components

SR_qR <- function(LqR){
  Screes <- numeric()
  
  for(i in 1:length(LqR)){
    if(i == 1){
      Screes[i] <- NA
    }
    else if(i == length(LqR) ){
      Screes[i] <- NA
    }else{
      Screes[i] <- (LqR[i-1] - LqR[i]) / (LqR[i] - LqR[i+1])
    }
  }
  return(Screes)
}

# optimal number of clusters = 3
# the highest number indicates the position of the optimal Q. So if the fourth number is the highest then the fourth integer of the numbered vector of Q (e.g., c(2,3,4,5,6,7)) is the optimal number of Q.

VAFdata[VAFdata$n.clus==3,]
SR_qR(VAFdata[VAFdata$n.clus==3,]$loss)


p <- plot_ly(data = VAFdata, x = unique(VAFdata$n.comp), y = ~loss[VAFdata$n.clus==2], name = "Clus 2", type = 'scatter', mode = 'lines+markers') %>% 
  add_trace(y = ~loss[VAFdata$n.clus==3], name = "Clus 3") %>% 
  add_trace(y = ~loss[VAFdata$n.clus==4], name = "Clus 4") %>% 
  add_trace(y = ~loss[VAFdata$n.clus==5], name = "Clus 5") %>% 
  layout(title = "Sequential scree plot",
         yaxis = list(title = "Variance accounted for",
                      range=c(160000,180000)),
         xaxis = list (title = "Components")
  )
p


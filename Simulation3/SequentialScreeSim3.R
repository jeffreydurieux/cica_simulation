# Wed Sep 10 10:57:09 2019
# Author: Jeffrey Durieux, MSc

# What: R-script for the sequential scree model selection for simulation 3
# this script is not generally applicable, it is interactive and dependent on the data that you need to check. If you want to use it you need to change some input parameters like: totalSSQ, ncomp and nclus

#### VAF model selection
####


library(plotly)

## sub selection 40 * 1000
## change this value depending on the data

data <- get(load("~/Repositories/cica_simulation/Data/CICA_SIM_3/rep_1CICA_simdata_0.30.Rdata"))
totalSSQ <- sum( sapply(data, function(i) sum(i^2) ) )


############# Model Selection Prep ##############
ncomp <- c(2,3,4,5,6,7)
nclus <- c(1,2,3,4,5)

grid <- expand.grid(n.comp = ncomp, n.clus = nclus)

setwd("~/Repositories/cica_simulation/Data/CICA_SIM_3/")

files <- dir(pattern = "resultSim3")

# convert to VAF
VAF <- function(Loss, totalSSQ){
  vaf <- (totalSSQ - Loss) / totalSSQ
  return(vaf*100)
}


loss <- numeric()
for(i in 1:length(files)){
  vals <- grid[i,]
  nam <- paste("clus", vals$n.clus, "_comp", vals$n.comp, "_resultSim3.Rdata", sep = "") 
  cica <- get( load(nam))
  loss[i] <- cica$beststart$LossHistory[cica$beststart$Iterations]
}
loss

VAFS <- sapply(loss, VAF, totalSSQ=totalSSQ)
VAFS

VAFdata <- data.frame(grid,loss=loss,VAF=VAFS)

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

Screes <- rbind(SR_rq(VAFdata[VAFdata$n.comp==2,]$loss) ,
SR_rq(VAFdata[VAFdata$n.comp==3,]$loss) ,
SR_rq(VAFdata[VAFdata$n.comp==4,]$loss) ,
SR_rq(VAFdata[VAFdata$n.comp==5,]$loss) ,
SR_rq(VAFdata[VAFdata$n.comp==6,]$loss) ,
SR_rq(VAFdata[VAFdata$n.comp==7,]$loss) )

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


#### plotting options ######
# change VAF with loss if necessary
p <- plot_ly(data = VAFdata, x = unique(VAFdata$n.comp), y = ~VAF[VAFdata$n.clus==2], name = "Clus 2", type = 'scatter', mode = 'lines+markers') %>% 
  add_trace(y = ~VAF[VAFdata$n.clus==3], name = "Clus 3") %>% 
  add_trace(y = ~VAF[VAFdata$n.clus==4], name = "Clus 4") %>% 
  add_trace(y = ~VAF[VAFdata$n.clus==5], name = "Clus 5") %>% 
  layout(title = "Sequential scree plot",
         yaxis = list(title = "Variance accounted for",
                      range=c(40,70)),
         xaxis = list (title = "Components")
  )
p


# Wed Jul 10 10:57:09 2019
# Author: Jeffrey Durieux, MSc

# What: R-script for the sequential scree model selection
# this script is not generally applicable, it is interactive and dependent on the data that you need to check. If you want to use it you need to change some input parameters like: totalSSQ, ncomp and nclus

#### VAF model selection
####


library(plotly)

## sub selection 40 * 1000
## change this value depending on the data
totalSSQ <- 40000


############# Functions ##############

# convert to VAF
VAF <- function(Loss, totalSSQ){
  vaf <- (totalSSQ - Loss) / totalSSQ
  return(vaf*100)
}

#step 1 of sequentials: calculate scree ratio for each value of K (cluster) given different diff values of Q (components)
#Then average all the SR_K over all components. and select the K that yields the highest average

SR_K <- function(vafKvec){
  res <- numeric()
  for(i in 2:(length( vafKvec)-1 ) ){
    res[i] <- (vafKvec[i] - vafKvec[i-1]) / (vafKvec[i+1] - vafKvec[i])  
  }
  res <- c(res, NA)
  return(res)
}

#put SR_K in matrix (each column is a Q) -> compute rowmeans and select highest row = Kbest


#step 2 of sequentials: given the optimal K, for each number of Q compute scree ratios
# the best number of Q is the number of comp for which the scree ratio is maximal


SR_Q <- function(vafQvec){
  res <- numeric()
  for(i in 2:(length(vafQvec)-1 ) ){
    res[i] <- (vafQvec[i] - vafQvec[i-1] ) / (vafQvec[i+1] - vafQvec[i] )
  }
  res <- c(res, NA)
  return(res)
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
loss

VAFS <- sapply(loss, VAF, totalSSQ=totalSSQ)
VAFS

VAFdata <- data.frame(grid,loss=loss,VAF=VAFS)




############# Model Selection ##############

Kscree <- matrix(data = NA, nrow = 4,ncol = 4)
row.names(Kscree) <- 2:5

for(i in 2:5){
  vafKvec <- VAFdata[VAFdata$n.clus==i,]  
  res <- SR_K(vafKvec = vafKvec$VAF)
  Kscree[i-1,] <- res
}

#rowsums 
rm <- rowMeans(Kscree, na.rm = TRUE)
max(rm, na.rm = TRUE)
which.max(rm) ### max rm is cluster to select

Qvec <- VAFdata[ VAFdata$n.clus==5,]
Qscree <- SR_Q(Qvec$VAF)
which.max(Qscree)
max(Qscree, na.rm = TRUE)
Qvec


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


# Sat Jul  4 12:56:33 2020
# Author: Jeffrey Durieux, MSc

# What: R-script for the sequential scree model selection for simulation 3
# this script is not generally applicable, it is interactive and dependent on the data that you need to check. If you want to use it you need to change some input parameters like: totalSSQ, ncomp and nclus

##### first part is to collect vafdata, object is saved on repository so load object to do model selection halfway the script


#### VAF model selection
####

# convert to VAF
VAF <- function(Loss, totalSSQ){
  vaf <- (totalSSQ - Loss) / totalSSQ
  return(vaf*100)
}


library(plotly)


## change this value depending on the data

setwd('~/Repositories/cica_simulation/Data/CICA_SIM_3/ALICEsimdata/')
simfiles <- dir()
simfiles <- gtools::mixedsort(simfiles)
simfiles


# 600 files with cica results
# 5 clusters X 6 components X 2 noise structures X 10 replications
setwd('~/Repositories/cica_simulation/Data/CICA_SIM_3/ALICE_results/')

files <- dir()
files <- files[grep('.Rdata', files)]


# this code will first load simulated dataset (20 files), per rep
# then load the associated cica results of the model selection

# subf contains the files names for each sim iteration. Manually checked it and everyting is in the right order: first check sim = 1 and clusf for order

ListVAFDATA <- list()
simm <- rep(1:20, each = 2)
for(sim in 1:length(simfiles)){
  setwd('~/Repositories/cica_simulation/Data/CICA_SIM_3/ALICEsimdata/')
  
  simfiles[sim]
  data <- get(load(simfiles[sim]))
  totalSSQ <- sum( sapply(data, function(i) sum(i^2) ) )
  
  setwd('~/Repositories/cica_simulation/Data/CICA_SIM_3/ALICE_results/')
  
  
  patt <- paste('rep',simm[sim],'_', sep = '')
  subf <- files[grep(pattern = patt, files)]
  
  if( (sim %% 2) == 1){
    # noise 30
    subf <- subf[grep(pattern = 'noise30', subf)]
  }else{
    # noise 70
    subf <- subf[grep(pattern = 'noise70', subf)]
  }
  subf
  
  
  ### extract loss value from dataset
  ncomp <- c(2,3,4,5,6,7)
  nclus <- c(1,2,3,4,5)
  
  grid <- expand.grid(n.comp = ncomp, n.clus = nclus)
  
  loss <- numeric()
  for(i in 1:length(subf)){
   load(subf[i])
   loss[i] <- cica$beststart$LossHistory[cica$beststart$Iterations]
   rm(cica)
  }
  
  VAFS <- sapply(loss, VAF, totalSSQ=totalSSQ)
  VAFdata <- data.frame(grid,loss=loss,VAF=VAFS)
  
  rm(totalSSQ, data)
  
  ListVAFDATA[[sim]] <- VAFdata
}

names <- gsub('.{6}$','',simfiles)
names(ListVAFDATA) <- names

save(ListVAFDATA, file = '~/Repositories/cica_simulation/Data/ListVafdata_Sim3.Rdata')


##### load list of vaf data #####
load('~/Repositories/cica_simulation/Data/ListVafdata_Sim3.Rdata')



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




##### select vaf data from list #####

OptimalK <- numeric()
OptimalQ <- numeric()
for(modsel in 1:20){
  VAFdata <- ListVAFDATA[[modsel]]  
  
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
  which.max(colMeans(Screes, na.rm = TRUE))
  
  Kopt <- which.max(colMeans(Screes, na.rm = TRUE))
  OptimalK[modsel] <- Kopt
  
  #### compute step 2: given K calculate opt Q. Is max value
  VAFdata[VAFdata$n.clus==Kopt,]
  SR_qR(VAFdata[VAFdata$n.clus==Kopt,]$loss)
  # + 1 since Q starts at 2 and not at 1
  which.max( SR_qR(VAFdata[VAFdata$n.clus==Kopt,]$loss) ) + 1
  
  Qopt <- which.max( SR_qR(VAFdata[VAFdata$n.clus==Kopt,]$loss) ) + 1
  OptimalQ[modsel] <- Qopt
  rm(VAFdata)
}

table(OptimalK)
table(OptimalQ)

names(ListVAFDATA)
noise <- as.factor(rep(c(.3,.7), 10))
noise
result_modsel <- data.frame(noise, K = OptimalK, Q = OptimalQ )



result_modsel

precK <- abs( result_modsel$K - 3)
precQ <- abs( result_modsel$Q - 4)

result_modsel$precK <- precK
result_modsel$precQ <- precQ

doBy::summaryBy(precK~noise, data=result_modsel)
doBy::summaryBy(precQ~noise, data=result_modsel)

table(result_modsel$K)
u <- sort(union(result_modsel$K, 1:5))
tabK <- table( factor(result_modsel$K,u))
tabK



table(result_modsel$Q)
u <- sort(union(result_modsel$Q, 2:7))
tabQ <- table( factor(result_modsel$Q,u))
tabQ



cbind(result_modsel, HIT = (result_modsel$K == 3 & result_modsel$Q == 4))

sum ( (result_modsel$K == 3 & result_modsel$Q == 4) )
7/20

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


#### table original loss values  CICA paper ####
# Q = 2 ..7 in rows, R =1 to 5 in cols

VAFdata









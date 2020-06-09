# Tue Jul  9 14:16:57 2019
# Author: Jeffrey Durieux, MSc

# What:R-script to check results of simulation 2 of C-ICA paper

library(mclust)    # Adjusted Rand Index
library(multiway)  # Tucker congruence
library(gtools)    # mixed sort
library(doBy)      # summaryBy
library(ez)        # ez anova
library(CICA)

setwd("~/Repositories/cica_simulation/Simulation2/CreateMaps/")

bmask <- get(load("bmaskslice.Rdata"))
idx <- which(bmask==1)

S1 <- get(load("cluster1.Rdata"))
S2 <- get(load("cluster2.Rdata"))
S3 <- get(load("cluster3.Rdata"))

S1 <- S1[ idx, ]
S2 <- S2[ idx, ]
S3 <- S3[ idx, ]

TRUEP <- c(rep(1,10),rep(2,10),rep(3,10))

#### simulation is run on Shark, download output to local desktop ####
setwd("~/Desktop/CICA_Sim2/Output/")
results <- mixedsort(dir())

ARI <- numeric()
Tuck <- numeric()
for(i in 1:60){
  #### c-ica result files
  resdatafile <- get(load(results[i]))
  
  estP <- resdatafile$beststart$P
  ARI[i] <- adjustedRandIndex(TRUEP, estP)
  
  idp <- unique(estP) # only works when ARI is 1 (which is the case for this Simulation)
  
  # ugly piece of code according to some R wizards I guess. Should I use the forward piping operator (i.e., %>% ) from magrittr/tidyverse? 
  S1m <- mean (apply( abs( congru(S1, resdatafile$beststart$Sr[[ idp[1] ]])), MARGIN = 2, max ) )
  S2m <- mean (apply( abs( congru(S2, resdatafile$beststart$Sr[[ idp[2] ]])), MARGIN = 2, max ) )
  S3m <- mean (apply( abs( congru(S3, resdatafile$beststart$Sr[[ idp[3] ]])), MARGIN = 2, max ) )
  
  Tuck[i] <- Smeans <- mean(S1m,S2m,S3m)
  
}

####### extract Tucker values for time courses Air #########

setwd("~/Desktop/CICA_Sim2/")
TuckA <- numeric()

repidx <- rep(1:10, each = 6)

setwd("Output/")
for(i in 1:60){
  
  if( i > 1){
    setwd("../Output/")  
  }
  #### c-ica result files
  resdatafile <- get(load(results[i]))
  
  estP <- resdatafile$beststart$P
  idp <- unique(estP)
  
  
  ## estimated time courses
  EstA <- Asplit(resdatafile$beststart$Air, resdatafile$beststart$P, nTS = 150)
  EstA <- EstA[idp]
  
  EstA <- unlist(EstA, recursive = F)
  EstA <- lapply(EstA, t)
  
  #### get timecourses 
  setwd("../timecourses/")
  files <- dir()
  files <- mixedsort(files)
  timecoursesdata <- get(load(files[repidx[i]]) )
  
  AperSubject <- numeric()
  for(j in 1:30){
    AperSubject[j] <- mean(apply( abs( congru(timecoursesdata[[j]], EstA[[j]])), MARGIN = 2, max ) )  
  }  
  TuckA[i] <- mean(AperSubject)
}



Noise <- c("AR10","AR30","AR70","G10","G30","G70")

Noise <- rep(Noise,10)
Noise <- as.factor(Noise)

Results <- data.frame(Case = 1:60, Noise = Noise, ARI = ARI, Tucker = Tuck, TuckerA = TuckA)

summaryBy(formula = ARI~Noise, data = Results, FUN = list(mean, sd))
summaryBy(formula = Tucker~Noise, data = Results, FUN = list(mean, sd))
summaryBy(formula = TuckerA~Noise, data = Results, FUN = list(mean, sd))

### overall results

#### change Tucker with Tucker A to see either spatial maps or time courses
anova <- ezANOVA(
  data = Results
  , dv = TuckerA
  , wid = Case
  , within = NULL
  , within_full = NULL
  , within_covariates = NULL
  , between = Noise
  , between_covariates = NULL
  , observed = NULL
  , diff = NULL
  , reverse_diff = FALSE
  , type = 3
  , white.adjust = FALSE
  , detailed = FALSE
  , return_aov = TRUE
)
anova

ezPlot(
  data = Results
  , dv = TuckerA
  , wid = Case
  , within = NULL
  , within_full = NULL
  , within_covariates = NULL
  , between = Noise
  , between_full = NULL
  , between_covariates = NULL
  , x = .(Noise)
  , do_lines = TRUE
  , do_bars = TRUE
  , bar_width = NULL
  , bar_size = NULL
  , split = NULL
  , row = NULL
  , col = NULL
  , to_numeric = NULL
  , x_lab = NULL
  , y_lab = NULL
  , split_lab = NULL
  , levels = NULL
  , diff = NULL
  , reverse_diff = FALSE
  , type = 3
  , dv_levs = NULL
  , dv_labs = NULL
  , y_free = FALSE
  , print_code = FALSE
)

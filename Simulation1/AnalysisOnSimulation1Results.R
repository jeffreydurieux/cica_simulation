# Fri Jul  5 10:29:45 2019
# Author: Jeffrey Durieux, MSc

# What: Data analysis script for anova data

# marginal means for three outcomes: ARI, S tucker and A tucker


library(doBy)
library(repmis)
library(ez)

##### Load anova input dataframe #######
# this does not work yet since repository is private
source_data("https://github.com/jeffreydurieux/cica_simulation/blob/master/Simulation1ANOVAdata.Rdata?raw=true")

load("~/Repositories/cica_simulation/Data/Simulation1ANOVAdata.Rdata")


####### Check means #########

dat <- data.frame(V=as.factor(des$V),Q = as.factor(des$Q), R = as.factor(des$R), 
                  D = as.factor(des$D), E = as.factor(des$E), 
                  ARI = des$ARI, Stuck = des$STuck, Atuck = des$ATuck, Time = des$Time)


###### ARI ########

mean(dat$ARI);sd(dat$ARI)

summaryBy(ARI~V ,data = dat, FUN = c(mean,sd))
summaryBy(ARI~Q ,data = dat, FUN = c(mean,sd))
summaryBy(ARI~R ,data = dat, FUN = c(mean,sd))
summaryBy(ARI~D ,data = dat, FUN = c(mean,sd))
summaryBy(ARI~E ,data = dat, FUN = c(mean,sd))


###### S Tuck ########

mean(dat$Stuck);sd(dat$Stuck)

summaryBy(Stuck~V ,data = dat, FUN = c(mean,sd))
summaryBy(Stuck~Q ,data = dat, FUN = c(mean,sd))
summaryBy(Stuck~R ,data = dat, FUN = c(mean,sd))
summaryBy(Stuck~D ,data = dat, FUN = c(mean,sd))
summaryBy(Stuck~E ,data = dat, FUN = c(mean,sd))


###### A Tuck ########

mean(dat$Atuck);sd(dat$Atuck)

summaryBy(Atuck~V ,data = dat, FUN = c(mean,sd))
summaryBy(Atuck~Q ,data = dat, FUN = c(mean,sd))
summaryBy(Atuck~R ,data = dat, FUN = c(mean,sd))
summaryBy(Atuck~D ,data = dat, FUN = c(mean,sd))
summaryBy(Atuck~E ,data = dat, FUN = c(mean,sd))



####### ANOVA ##########

# DO THREE ANOVA'S (ARI, STUCK AND ATUCK)
# 5 BETWEEN FACTORS



dat <- data.frame(id = 1:720, dat)

#### ARI
#note: nothing is significant
ari_anova_res <- ezANOVA( data = dat , ARI , wid = id, within = NULL , within_full = NULL , within_covariates = NULL , between = .('V','Q','R','D','E') , between_covariates = NULL , observed = NULL , diff = NULL , reverse_diff = FALSE , type = 3 , white.adjust = FALSE , detailed = FALSE , return_aov = TRUE)

print(ari_anova_res)

#### Stuck
# given etasq > .15 at most two-way interactions:
# V:Q = .58
# R:D = .17
# R:E = .17
# D:E = .47
Stuck_anova_res <- ezANOVA( data = dat , Stuck , wid = id, within = NULL , within_full = NULL , within_covariates = NULL , between = .('V','Q','R','D','E') , between_covariates = NULL , observed = NULL , diff = NULL , reverse_diff = FALSE , type = 3 , white.adjust = FALSE , detailed = FALSE , return_aov = TRUE)

print(Stuck_anova_res)
Stuck_anova_res$aov


###### check how effects are going 

summaryBy(Stuck~V ,data = dat, FUN = c(mean))
summaryBy(Stuck~Q ,data = dat, FUN = c(mean))
summaryBy(Stuck~R ,data = dat, FUN = c(mean))
summaryBy(Stuck~D ,data = dat, FUN = c(mean))
summaryBy(Stuck~E ,data = dat, FUN = c(mean)) 

summaryBy(Stuck~V:Q ,data = dat, FUN = c(mean))
summaryBy(Stuck~D:E ,data = dat, FUN = c(mean))
summaryBy(Stuck~R:E ,data = dat, FUN = c(mean))
summaryBy(Stuck~D:R ,data = dat, FUN = c(mean))


#### Atuck
# given etasq > .15 at most two-way interactions:
# V:Q = .80
# Q:E = .21

Atuck_anova_res <- ezANOVA( data = dat , Atuck , wid = id, within = NULL , within_full = NULL , within_covariates = NULL , between = .('V','Q','R','D','E') , between_covariates = NULL , observed = NULL , diff = NULL , reverse_diff = FALSE , type = 3 , white.adjust = FALSE , detailed = FALSE , return_aov = TRUE)

print(Atuck_anova_res)
Atuck_anova_res$aov

###### check how effects are going 

summaryBy(Atuck~V ,data = dat, FUN = c(mean))
summaryBy(Atuck~Q ,data = dat, FUN = c(mean))
summaryBy(Atuck~E ,data = dat, FUN = c(mean))


summaryBy(Atuck~V:Q ,data = dat, FUN = c(mean))
summaryBy(Atuck~Q:E ,data = dat, FUN = c(mean))




######## Time Analysis #######
# note: this has to be tested and is not implemented yet
mean(dat$Time);sd(dat$Time)

summaryBy(Time~V ,data = dat, FUN = c(mean,sd))
summaryBy(Time~Q ,data = dat, FUN = c(mean,sd))
summaryBy(Time~R ,data = dat, FUN = c(mean,sd))
summaryBy(Time~D ,data = dat, FUN = c(mean,sd))
summaryBy(Time~E ,data = dat, FUN = c(mean,sd))


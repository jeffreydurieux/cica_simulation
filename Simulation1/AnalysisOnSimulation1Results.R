# Fri Jul  5 10:29:45 2019
# Author: Jeffrey Durieux, MSc

# What: Data analysis script for anova data

# marginal means for three outcomes: ARI, S tucker and A tucker

# Tue Jun 30 15:11:18 2020
# added plotly figures for anova interactions

library(doBy)
library(repmis)
library(ez)
library(plotly)

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

summaryBy(Stuck~V:Q ,data = dat, FUN = c(mean,sd))
summaryBy(Stuck~D:E ,data = dat, FUN = c(mean,sd))
summaryBy(Stuck~R:E ,data = dat, FUN = c(mean,sd))
summaryBy(Stuck~D:R ,data = dat, FUN = c(mean,sd))


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


#### plots #####

## S
# given etasq > .15 at most two-way interactions:
# V:Q = .58
# D:E = .47
# R:E = .17
# R:D = .17


### fig1a ####
D_VQ <- summaryBy(Stuck~V:Q ,data = dat, FUN = c(mean,plotrix::std.error))
D_VQ
D_VQ <- rename(D_VQ, c("Stuck.FUN1" = "mean"))
D_VQ <- rename(D_VQ, c("Stuck.FUN2" = "se"))
D_VQ$se2 <- D_VQ$se*2
D_VQ
data <- D_VQ

fig1a <- plot_ly(data = data[which(data$V == 500),], x = ~Q, y = ~mean, type = 'scatter', mode = 'lines+markers',
               name = '500 voxels',
               error_y = ~list(array = se2,
                               color = '#000000'))
fig1a <- fig1a %>% add_trace(data = data[which(data$V == 2000),], name = '2000 voxels')
fig1a <- fig1a %>% layout(yaxis = list(title = 'Tucker congruence',
                                       range = c(0.90,1.01) ),
                          xaxis = list(title = 'FC patterns Q')
                          )
fig1a


### fig1b ####
D_DE <- summaryBy(Stuck~D:E ,data = dat, FUN = c(mean,plotrix::std.error))
D_DE
D_DE <- rename(D_DE, c("Stuck.FUN1" = "mean"))
D_DE <- rename(D_DE, c("Stuck.FUN2" = "se"))
D_DE$E <- plyr::revalue(D_DE$E, c('0.05'='5','0.2'='20','0.4'='40'))
D_DE$se2 <- D_DE$se*2
D_DE
data <- D_DE

fig1b <- plot_ly(data = data[which(data$D == 'square'),], x = ~E, y = ~mean, type = 'scatter', mode = 'lines+markers',
                 name = 'Square',
                 error_y = ~list(array = se2,
                                 color = '#000000'))
fig1b <- fig1b %>% add_trace(data = data[which(data$D == 'nonsquare'),], name = 'Non-square')
fig1b <- fig1b %>% layout(yaxis = list(title = 'Tucker congruence',
                                       range = c(0.90,1.01) ),
                          xaxis = list(title = 'Noise %'))
fig1b

### fig1c ####
D_RE <- summaryBy(Stuck~R:E ,data = dat, FUN = c(mean,plotrix::std.error))
D_RE
D_RE <- rename(D_RE, c("Stuck.FUN1" = "mean"))
D_RE <- rename(D_RE, c("Stuck.FUN2" = "se"))
D_RE$E <- plyr::revalue(D_RE$E, c('0.05'='5','0.2'='20','0.4'='40'))
D_RE$se2 <- D_RE$se*2
D_RE
data <- D_RE

fig1c <- plot_ly(data = data[which(data$R == 2),], x = ~E, y = ~mean, type = 'scatter', mode = 'lines+markers',
                 name = '2 Clusters',
                 line = list(color = 'rgb(140, 86, 75)'),
                 marker = list(color = 'rgb(140, 86, 75)'),
                 error_y = ~list(array = se2,
                                 color = '#000000'))
fig1c <- fig1c %>% add_trace(data = data[which(data$R == 4),], name = '4 Clusters',
                             line = list(color = 'rgb(148, 103, 189)'),
                             marker = list(color = 'rgb(148, 103, 189)'))
fig1c <- fig1c %>% layout(yaxis = list(title = 'Tucker congruence',
                                       range = c(0.90,1.01) ),
                          xaxis = list(title = 'Noise %'))
fig1c


### fig1d ####

D_DR <- summaryBy(Stuck~D:R ,data = dat, FUN = c(mean,plotrix::std.error))
D_DR
D_DR <- rename(D_DR, c("Stuck.FUN1" = "mean"))
D_DR <- rename(D_DR, c("Stuck.FUN2" = "se"))
#D_DR$E <- plyr::revalue(D_RE$E, c('0.05'='5','0.2'='20','0.4'='40'))
D_DR$se2 <- D_DR$se*2
D_DR
data <- D_DR

fig1d <- plot_ly(data = data[which(data$R == 2),], x = ~D, y = ~mean, type = 'scatter', mode = 'lines+markers',
                 name = '2 clusters', 
                 line = list(color = 'rgb(140, 86, 75)'),
                 marker = list(color = 'rgb(140, 86, 75)'),
                 error_y = ~list(array = se2,
                                 color = '#000000'))
fig1d <- fig1d %>% add_trace(data = data[which(data$R == 4),], 
                             name = '4 clusters',
                             line = list(color = 'rgb(148, 103, 189)'),
                             marker = list(color = 'rgb(148, 103, 189)'))
fig1d <- fig1d %>% layout(yaxis = list(title = 'Tucker congruence',
                                       range = c(0.90,1.01) ),
                          xaxis = list(title = 'Dimension A'))
fig1d

# V:Q = .58
# D:E = .47
# R:E = .17
# R:D = .17

subplot( style(fig1a, showlegend = T),
         style(fig1b, showlegend = T),
         style(fig1c, showlegend = F),
         style(fig1d, showlegend = T),
         nrows = 2, shareY = T, margin = 0.10,titleX = T) 
  
subplot( style(fig1a, showlegend = T),
         style(fig1b, showlegend = T),
         nrows = 1, shareY = T, margin = 0.10,titleX = T)
  
subplot( style(fig1c, showlegend = F),
         style(fig1d, showlegend = T),
         nrows = 1, shareY = T, margin = 0.10,titleX = T)


  # layout(annotations = list(
  #   list(x = 0.13 , y = 0.45, text = "<b>Clusters = 3<b>", showarrow = F, xref='paper', yref='paper', font = list(size = 15)),
  #   list(x = 0.88 , y = 0.45, text = "<b>Clusters = 5<b>", showarrow = F, xref='paper', yref='paper', font = list(size = 15))
  # )
  # )



## A
# V:Q = .80
#### fig2a #####
Dd_DR <- summaryBy(Atuck~V:Q ,data = dat, FUN = c(mean,plotrix::std.error))
Dd_DR
Dd_DR <- rename(Dd_DR, c("Atuck.FUN1" = "mean"))
Dd_DR <- rename(Dd_DR, c("Atuck.FUN2" = "se"))
#D_DR$E <- plyr::revalue(D_RE$E, c('0.05'='5','0.2'='20','0.4'='40'))
Dd_DR$se2 <- Dd_DR$se*2
Dd_DR
data <- Dd_DR

fig2a <- plot_ly(data = data[which(data$Q == 2),], x = ~V, y = ~mean, type = 'scatter', mode = 'lines+markers',
                 name = 'Q = 2',
                 error_y = ~list(array = se2,
                                 color = '#000000'))
fig2a <- fig2a %>% add_trace(data = data[which(data$Q == 5),], name = 'Q = 5')
fig2a <- fig2a %>% add_trace(data = data[which(data$Q == 20),], name = 'Q = 20')
fig2a <- fig2a %>% layout(yaxis = list(title = 'Tucker congruence',
                                       range = c(0.90,1.01) ),
                          xaxis = list(title = 'Voxels'))
fig2a


# Q:E = .21
#### fig2b #####
Dd_QE <- summaryBy(Atuck~Q:E ,data = dat, FUN = c(mean,plotrix::std.error))
Dd_QE
Dd_QE <- rename(Dd_QE, c("Atuck.FUN1" = "mean"))
Dd_QE <- rename(Dd_QE, c("Atuck.FUN2" = "se"))
Dd_QE$E <- plyr::revalue(Dd_QE$E, c('0.05'='5','0.2'='20','0.4'='40'))
Dd_QE$se2 <- Dd_QE$se*2
Dd_QE
data <- Dd_QE

fig2b <- plot_ly(data = data[which(data$Q == 2),], x = ~E, y = ~mean, type = 'scatter', mode = 'lines+markers',
                 name = 'Q = 2',
                 line = list(color = 'rgb(31, 119, 180)'),
                 marker = list(color = 'rgb(31, 119, 180)'),
                 error_y = ~list(array = se2,
                                 color = '#000000'))
fig2b <- fig2b %>% add_trace(data = data[which(data$Q == 5),], 
                             line = list(color = 'rgb(255, 127, 14)'),
                             marker = list(color = 'rgb(255, 127, 14)'),
                             name = 'Q = 5')
fig2b <- fig2b %>% add_trace(data = data[which(data$Q == 20),], 
                             line = list(color = 'rgb(44, 160, 44)'),
                             marker = list(color = 'rgb(44, 160, 44)'),
                             name = 'Q = 20')
fig2b <- fig2b %>% layout(yaxis = list(title = 'Tucker congruence',
                                       range = c(0.90,1.01) ),
                          xaxis = list(title = 'Noise %'))
fig2b

subplot( style(fig2a, showlegend = F),
         style(fig2b, showlegend = T),
         nrows = 1, shareY = T, margin = 0.10,titleX = T)

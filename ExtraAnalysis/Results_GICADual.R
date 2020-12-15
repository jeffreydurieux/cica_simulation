setwd('~/Repositories/cica_simulation/Data/CICA_SIM_4/CICA_GICAdualResults/')
library(gtools)

files <- dir()

files <- mixedsort(files)


load(files[1])
output$ARIs
ARI <- output$ARIs

for(i in 2:length(files)){
  load(files[i])
  ARI <- rbind(ARI,output$ARIs)
}
ARI
apply(ARI, 2, mean)
apply(ARI, 2, sd)


V <- c(500,2000)                # voxels
Q <- c(2,5,20)                  # components
R <- c(2,4)                     # clusters
D <- c(1, 2)                    # 1 = square 2 = nonsquare dimension mixing matrix
E <- c(0.05, 0.2, 0.4)          # noise level
rep <- 1:10
design <- expand.grid(V=V, Q = Q,
                      R = R, D = D, E=E,Rep = rep)

res <- cbind(design, ARI)

library(doBy)
summary_by(data = res, formula = hcl+pam+km+pcakm~.,FUN = c(mean,sd))
summary_by(data = res, formula = hcl+pam+km+pcakm~V)
summary_by(data = res, formula = hcl+pam+km+pcakm~Q)
summary_by(data = res, formula = hcl+pam+km+pcakm~D)
summary_by(data = res, formula = hcl+pam+km+pcakm~R)
summary_by(data = res, formula = hcl+pam+km+pcakm~E)



# G2
summary_by(data = res, formula = hclG2+pamG2+kmG2+pcakmG2~.)
summary_by(data = res, formula = hclG2+pamG2+kmG2+pcakmG2~V)
summary_by(data = res, formula = hclG2+pamG2+kmG2+pcakmG2~Q)
summary_by(data = res, formula = hclG2+pamG2+kmG2+pcakmG2~R)
summary_by(data = res, formula = hclG2+pamG2+kmG2+pcakmG2~E)


#### add CICA results #####
load('~/Repositories/cica_simulation/Data/Simulation1ANOVAdata.Rdata')

library(plotly)

res <- cbind(res, CICA=des$ARI)

summary_by(data = res, formula = CICA+hcl+pam+km+pcakm~.,FUN = c(mean,sd))

library(plotly)


fig <- plot_ly(data = res, y = res$CICA,boxmean= 'sd', type = "box", name="CICA")%>%
  add_trace( y = res$hcl,  name="HCL") %>%
  add_trace( y = res$pam,  name="PAM") %>%
  add_trace( y = res$km,  name="KM") %>%
  add_trace( y = res$pcakm,  name="PCAKM") %>%
  layout(title = "GICA dual regression and clustering")

fig

mcica <- mean(res$CICA)
sdcica <- sd(res$CICA)

mhcl <- mean(res$hcl)
sdhcl <- sd(res$hcl)

mpam <- mean(res$pam)
sdpam <- sd(res$pam)

mkm <- mean(res$km)
sdkm <- sd(res$km)

mpcakm <- mean(res$pcakm)
sdpcakm <- sd(res$pcakm)

means <- c(mcica,mpam,mhcl, mpcakm,mkm)
sd <- c(sdcica,sdpam,sdhcl,sdpcakm,sdkm)
type <- c('1CICA', '2PAM', '3HCL', '4PCAKM', '5KM')
figdata <- data.frame(means,sd,type)

fig <- plot_ly(data = figdata, split=~type,x= ~type, y = mcica, type = "box", name="CICA")%>%
  add_trace( y = mpam,  name="PAM") %>%
  add_trace( y = mhcl,  name="HCL") %>%
  add_trace( y = mpcakm,  name="PCAKM") %>%
  add_trace( y = mkm,  name="KM") %>%
  layout(title = "GICA dual regression and clustering")

fig

fig1 <- plot_ly(data = figdata,split=~type, x= 1:5, y = ~means, type = "scatter", name="ARI", showlegende = F) %>%
  add_trace(error_y =list(array = ~sd), showlegend = F) %>%
  layout(xaxis = list(title = '',
         ticktext = c('CICA', 'PAM', 'HCL', 'PCAKM', 'KM'),
         tickvals = c(1, 2, 3, 4, 5),
         tickmode = "array"),
         yaxis = list(title='ARI'),
         title ='G-ICA on all data')

fig1

#### group ICA on sep clusters
mcica <- mean(res$CICA)
sdcica <- sd(res$CICA)

mhcl <- mean(res$hclG2)
sdhcl <- sd(res$hclG2)

mpam <- mean(res$pamG2)
sdpam <- sd(res$pamG2)

mkm <- mean(res$kmG2)
sdkm <- sd(res$kmG2)

mpcakm <- mean(res$pcakmG2)
sdpcakm <- sd(res$pcakmG2)

means <- c(mcica,mpam,mhcl, mpcakm,mkm)
sd <- c(sdcica,sdpam,sdhcl,sdpcakm,sdkm)
type <- c('1CICA', '2PAM', '3HCL', '4PCAKM', '5KM')
figdata2 <- data.frame(means,sd,type)

fig2 <- plot_ly(data = figdata2,split=~type, x= 1:5, y = ~means, type = "scatter", name="ARI", showlegende = F) %>%
  add_trace(error_y =list(array = ~sd), showlegend = F) %>%
  layout(xaxis = list(title = '',
                      ticktext = c('CICA', 'PAM', 'HCL', 'PCAKM', 'KM'),
                      tickvals = c(1, 2, 3, 4, 5),
                      tickmode = "array"),
         yaxis = list(title='ARI'),
         title ='')

fig2

p <- subplot(fig1,fig2, shareX = T, shareY = T)
p %>% layout(annotations = list(
  list(x = 0.2 , y = 1.05, text = "G-ICA on all data", showarrow = F, xref='paper', yref='paper'),
  list(x = .9 , y = 1.05, text = "G-ICA per cluster", showarrow = F, xref='paper', yref='paper'))
)


#### same plots but only square vs non square #####
### toggle id to 1 or 2
id <- which(res$D==2)
ress <- res[id,]

mcica <- mean(ress$CICA)
sdcica <- sd(ress$CICA)

mhcl <- mean(ress$hcl)
sdhcl <- sd(ress$hcl)

mpam <- mean(ress$pam)
sdpam <- sd(ress$pam)

mkm <- mean(ress$km)
sdkm <- sd(ress$km)

mpcakm <- mean(ress$pcakm)
sdpcakm <- sd(ress$pcakm)

means <- c(mcica,mpam,mhcl, mpcakm,mkm)
sd <- c(sdcica,sdpam,sdhcl,sdpcakm,sdkm)
type <- c('1CICA', '2PAM', '3HCL', '4PCAKM', '5KM')
figdata3 <- data.frame(means,sd,type)

mcica <- mean(ress$CICA)
sdcica <- sd(ress$CICA)

mhcl <- mean(ress$hclG2)
sdhcl <- sd(ress$hclG2)

mpam <- mean(ress$pamG2)
sdpam <- sd(ress$pamG2)

mkm <- mean(ress$kmG2)
sdkm <- sd(ress$kmG2)

mpcakm <- mean(ress$pcakmG2)
sdpcakm <- sd(ress$pcakmG2)

means <- c(mcica,mpam,mhcl, mpcakm,mkm)
sd <- c(sdcica,sdpam,sdhcl,sdpcakm,sdkm)
type <- c('1CICA', '2PAM', '3HCL', '4PCAKM', '5KM')
figdata4 <- data.frame(means,sd,type)

fig5 <- plot_ly(data = figdata3,split=~type, x= 1:5, y = ~means, type = "scatter", name="ARI", showlegende = F) %>%
  add_trace(error_y =list(array = ~sd), showlegend = F) %>%
  layout(xaxis = list(title = '',
                      ticktext = c('CICA', 'PAM', 'HCL', 'PCAKM', 'KM'),
                      tickvals = c(1, 2, 3, 4, 5),
                      tickmode = "array"),
         yaxis = list(title='ARI'),
         title ='')

fig5


fig6 <- plot_ly(data = figdata4,split=~type, x= 1:5, y = ~means, type = "scatter", name="ARI", showlegende = F) %>%
  add_trace(error_y =list(array = ~sd), showlegend = F) %>%
  layout(xaxis = list(title = '',
                      ticktext = c('CICA', 'PAM', 'HCL', 'PCAKM', 'KM'),
                      tickvals = c(1, 2, 3, 4, 5),
                      tickmode = "array"),
         yaxis = list(title='ARI'),
         title ='')

fig6

p <- subplot(fig3,fig4,fig5,fig6,nrows = 2, shareX = T, shareY = T)
p %>% layout(annotations = list(
  list(x = 0.2 , y = 1.05, text = "G-ICA on all data", showarrow = F, xref='paper', yref='paper'),
  list(x = .9 , y = 1.05, text = "G-ICA per cluster", showarrow = F, xref='paper', yref='paper'))
)


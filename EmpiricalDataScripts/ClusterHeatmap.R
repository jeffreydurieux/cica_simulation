# Tue Sep 17 11:37:04 2019
# Author: Jeffrey Durieux, MSc

# Code for a cluster heatmap of the results


library(plotly)


load("~/Repositories/cica_simulation/Data/PandLShiny.Rdata")



clusters <- data[[1]][,3:18]
clustersmat <- apply(clusters, MARGIN = 2, as.numeric)

clustersmat <- cbind(label = c(rep(-1, 173), rep(0,77)), clustersmat[,1:4])

ax <- list( ticktext = list('Label', 'R1', 'R2', 'R3', 'R4', 'R5'),
            tickvals = list(0, 1, 2, 3, 4),            
            tickmode = "array"
)

leg <- list()
plot_ly(z= clustersmat, type = 'heatmap',
        colors = colorRamp(c('grey', 'black','red', 'blue', 'darkgreen', 'orange', 'purple')), 
        xgap=50,
        ygap=0) %>% layout(xaxis = ax)


##### confusion matrices #########


addmargins(table(data[[1]]$labels, data[[1]]$Comp25_Clus2))
addmargins(table(data[[1]]$labels, data[[1]]$Comp25_Clus3))
addmargins(table(data[[1]]$labels, data[[1]]$Comp25_Clus4))
addmargins(table(data[[1]]$labels, data[[1]]$Comp25_Clus5))


######### measures #########

library(mclust)
library(caret)

adjustedRandIndex(data[[1]]$labels, data[[1]]$Comp25_Clus2)

lab <- data[[1]]$labels
lab <- ifelse((lab == 'AD'), 2, 1)
res <- data[[1]]$Comp20_Clus2
table(lab,res)

confusionMatrix(table(lab,res))

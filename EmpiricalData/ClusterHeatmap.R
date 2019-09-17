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

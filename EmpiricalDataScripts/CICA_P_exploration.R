library(mclust)
library(dplyr)

load('~/Repositories/cica_simulation/Data/CICAdashboard_P.Rdata')

adjustedRandIndex(data$DFlab$Comp20_Clus2, data$DFlab$Comp25_Clus2) %>% round(digits = 2)
adjustedRandIndex(data$DFlab$Comp30_Clus2, data$DFlab$Comp25_Clus2) %>% round(digits = 2)
adjustedRandIndex(data$DFlab$Comp35_Clus2, data$DFlab$Comp25_Clus2) %>% round(digits = 2)

# comp 20 vs comp 25
adjustedRandIndex(data$DFlab$Comp20_Clus2, data$DFlab$Comp25_Clus2) %>% round(digits = 2)
adjustedRandIndex(data$DFlab$Comp20_Clus3, data$DFlab$Comp25_Clus3) %>% round(digits = 2)
adjustedRandIndex(data$DFlab$Comp20_Clus4, data$DFlab$Comp25_Clus4) %>% round(digits = 2)
adjustedRandIndex(data$DFlab$Comp20_Clus5, data$DFlab$Comp25_Clus5) %>% round(digits = 2)

# comp 25 vs comp 30
adjustedRandIndex(data$DFlab$Comp25_Clus2, data$DFlab$Comp30_Clus2) %>% round(digits = 2)
adjustedRandIndex(data$DFlab$Comp25_Clus3, data$DFlab$Comp30_Clus3) %>% round(digits = 2)
adjustedRandIndex(data$DFlab$Comp25_Clus4, data$DFlab$Comp30_Clus4) %>% round(digits = 2)
adjustedRandIndex(data$DFlab$Comp25_Clus5, data$DFlab$Comp30_Clus5) %>% round(digits = 2)

# comp 30 vs comp 35
adjustedRandIndex(data$DFlab$Comp30_Clus2, data$DFlab$Comp35_Clus2) %>% round(digits = 2)
adjustedRandIndex(data$DFlab$Comp30_Clus3, data$DFlab$Comp35_Clus3) %>% round(digits = 2)
adjustedRandIndex(data$DFlab$Comp30_Clus4, data$DFlab$Comp35_Clus4) %>% round(digits = 2)
adjustedRandIndex(data$DFlab$Comp30_Clus5, data$DFlab$Comp35_Clus5) %>% round(digits = 2)


### means components comparison

# comp 20 vs comp 25
c( adjustedRandIndex(data$DFlab$Comp20_Clus2, data$DFlab$Comp25_Clus2) %>% round(digits = 2),
adjustedRandIndex(data$DFlab$Comp20_Clus3, data$DFlab$Comp25_Clus3) %>% round(digits = 2),
adjustedRandIndex(data$DFlab$Comp20_Clus4, data$DFlab$Comp25_Clus4) %>% round(digits = 2),
adjustedRandIndex(data$DFlab$Comp20_Clus5, data$DFlab$Comp25_Clus5) %>% round(digits = 2) ) %>% mean()

# comp 25 vs comp 30
c( adjustedRandIndex(data$DFlab$Comp25_Clus2, data$DFlab$Comp30_Clus2) %>% round(digits = 2),
adjustedRandIndex(data$DFlab$Comp25_Clus3, data$DFlab$Comp30_Clus3) %>% round(digits = 2),
adjustedRandIndex(data$DFlab$Comp25_Clus4, data$DFlab$Comp30_Clus4) %>% round(digits = 2),
adjustedRandIndex(data$DFlab$Comp25_Clus5, data$DFlab$Comp30_Clus5) %>% round(digits = 2)) %>% mean()

# comp 30 vs comp 35
c( adjustedRandIndex(data$DFlab$Comp30_Clus2, data$DFlab$Comp35_Clus2) %>% round(digits = 2),
adjustedRandIndex(data$DFlab$Comp30_Clus3, data$DFlab$Comp35_Clus3) %>% round(digits = 2),
adjustedRandIndex(data$DFlab$Comp30_Clus4, data$DFlab$Comp35_Clus4) %>% round(digits = 2),
adjustedRandIndex(data$DFlab$Comp30_Clus5, data$DFlab$Comp35_Clus5) %>% round(digits = 2)) %>% mean()


### means cluster comparison
# clus 2
c( adjustedRandIndex(data$DFlab$Comp20_Clus2, data$DFlab$Comp25_Clus2) %>% round(digits = 2),
   adjustedRandIndex(data$DFlab$Comp25_Clus2, data$DFlab$Comp30_Clus2) %>% round(digits = 2),
   adjustedRandIndex(data$DFlab$Comp30_Clus2, data$DFlab$Comp35_Clus2) %>% round(digits = 2)) %>% mean()

# clus 3
c( adjustedRandIndex(data$DFlab$Comp20_Clus3, data$DFlab$Comp25_Clus3) %>% round(digits = 2),
   adjustedRandIndex(data$DFlab$Comp25_Clus3, data$DFlab$Comp30_Clus3) %>% round(digits = 2),
   adjustedRandIndex(data$DFlab$Comp30_Clus3, data$DFlab$Comp35_Clus3) %>% round(digits = 2)) %>% mean()

# clus 4
c( adjustedRandIndex(data$DFlab$Comp20_Clus4, data$DFlab$Comp25_Clus4) %>% round(digits = 2),
   adjustedRandIndex(data$DFlab$Comp25_Clus4, data$DFlab$Comp30_Clus4) %>% round(digits = 2),
   adjustedRandIndex(data$DFlab$Comp30_Clus4, data$DFlab$Comp35_Clus4) %>% round(digits = 2)) %>% mean()

# clus 5
c( adjustedRandIndex(data$DFlab$Comp20_Clus5, data$DFlab$Comp25_Clus5) %>% round(digits = 2),
   adjustedRandIndex(data$DFlab$Comp25_Clus5, data$DFlab$Comp30_Clus5) %>% round(digits = 2),
   adjustedRandIndex(data$DFlab$Comp30_Clus5, data$DFlab$Comp35_Clus5) %>% round(digits = 2)) %>% mean()



adjustedRandIndex(data$DFlab$Comp20_Clus5, data$DFlab$Comp25_Clus5) %>% round(digits = 2)
adjustedRandIndex(data$DFlab$Comp25_Clus5, data$DFlab$Comp30_Clus5) %>% round(digits = 2)
adjustedRandIndex(data$DFlab$Comp30_Clus5, data$DFlab$Comp35_Clus5) %>% round(digits = 2)

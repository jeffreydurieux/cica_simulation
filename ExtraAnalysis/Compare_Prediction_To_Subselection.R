# Mon Sep  7 13:36:09 2020
# Author: Jeffrey Durieux, MSc 


# Compare prediction scores subject labels to subselection of simulation 3

### notes ###

# Bijgevoegd een R data file met 250 rijen (subjects) en vijf kolommen (vijf maal kruisvalidatieprocedure herhaald). Je kunt bijvoorbeeld de kolom nemen met de gemiddelde accuratesse/AUC/ARI waarde. De volgorde van de subjecten is alfabetisch ASPS004, ASPS006, … PROD157.
# 
# Dit zijn de uitkomsten als resultaat van een logistische regressie met een group lasso penalty. Het is een model met alle type resting state fMRI maten die ik heb berekend (31 type maten, ~3.000.000 predictoren).
# 
# Laat maar weten als er iets niet duidelijk is.
# 
# Groeten, Frank


load('~/Repositories/cica_simulation/Data/outcome_FrankdeVos_PredictionScores.RData')


ave_pred <- rowMeans(outcome)

# sub selection names:

load('~/Repositories/cica_simulation/Data/CICA_SIM_4/dataTV/subselectionnames.Rdata')

# all data names:
data <- read.csv('~/Desktop/Graz_clinical_data_N250.csv')
data <- data.frame(name = data$Subject, pred = ave_pred, 
                   class = c(rep(0,173), rep(1, 77)) )
data

# select 20 lowest pred healthy and 20 highest pred AD

dfhc <- data[data$class==0,]
dfhc_ord <- dfhc[order(dfhc$pred),]
dfhc_ord <- dfhc_ord[1:20,]


dfad <- data[data$class==1,]
dfad_ord <- dfad[order(dfad$pred, decreasing = T),]
dfad_ord
dfad_ord <- dfad_ord[1:20,]

dff_ord <- data[order(data$pred),]
dff_ord

low <- dff_ord[1:20,]
high <- dff_ord[231:250,]
low
high

#### compare ####

files <- substr(files, 1, 7)

### select line 58 or 59 
pred <- c(dfhc_ord$name,dfad_ord$name)
pred <- c(low$name, high$name)

cbind(files, pred)

check <- files %in% pred
pred%in%files
match(files, pred)
match(pred, files)

sum(files[21:40] %in% pred[21:40])
pred
dfad_ord
files


#check score
id <- data$name %in% files

checkdata <- data[id,]
boxplot(checkdata$pred~checkdata$class)
abline(h = 0.5)

boxplot(data$pred~data$class)
abline(h = 0.5)

library(plotly)

#prod010 wrongly allocated

mis <- which(checkdata$name == 'prod010')

data$class <- as.factor(data$class)
p1 <- ggplot(data, aes(x=class, y=pred)) + 
  geom_hline(yintercept = 0.5, linetype = 'dashed')+
  geom_boxplot()
p1

p1 <- ggplotly(p1)
p1

checkdata$class <- as.factor(checkdata$class)
p2 <- ggplot(checkdata, aes(x=class, y=pred)) + 
  geom_hline(yintercept = 0.5, linetype = 'dashed')+
  geom_boxplot() +
  geom_point(aes( x = '1', y = pred[mis]), shape = 2, color = 'red', size = 3)
p2

#without point/triangle
p2 <- ggplot(checkdata, aes(x=class, y=pred)) + 
  geom_hline(yintercept = 0.5, linetype = 'dashed')+
  geom_boxplot() 
p2

p2 <- ggplotly(p2)
p2

subplot(p1,p2)

predlabel <- ifelse(data$pred>0.5, 1, 0)
tab <- table(predlabel,data$class)
caret::confusionMatrix(tab)

ROCit::rocit(predlabel, data$class)

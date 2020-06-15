### script for checking G-ICA loss per subject

# Fri Jul  7 12:03:44 2017 ------------------------------

# Tue Jul 11 16:03:03 2017 ------------------------------



setwd("~/surfdrive/Documents/Projects/CICA/Project1/P1.4_Application/P1.4.1_Shark/June27/XPCADATA/")

files <- dir()

#ad vs prod
id_hc <- grep("asps", files)
id_ad <- grep("prod", files)

hcL <- list(); adL <- list()

for(i in 1:length(id_hc)){
  fhcl <- files[id_hc]
  load(fhcl[i])
  hcL[[i]] <- t(dataTV)
}
names(hcL) <- fhcl

for(i in 1:length(id_ad)){
  fadl <- files[id_ad]
  load(fadl[i])
  adL[[i]] <- t(dataTV)
}
names(adL) <- fadl

dataL <- c(hcL, adL)

setwd("~/surfdrive/Documents/Projects/CICA/Project1/P1.4_Application/P1.4.1_Shark/June27/GROUPICA/")

load("asps.Rdata")
gica_hc <- ica

load("prod.Rdata")
gica_ad <- ica

rm(dataTV,files,i,ica, id_ad, id_hc)


############ compute SSQ ########
# compute SSQ ||X - SA||^2

#split matrices
splitidx <- rep(1:77,each=50)
adM <- list()
for(i in 1:77){
  adM[[i]] <- gica_ad$M[ which( splitidx==i ) , ] 
}

splitidx <- rep(1:173,each=50)
hcM <- list()
for(i in 1:173){
  hcM[[i]] <- gica_hc$M[ which( splitidx==i ) , ] 
}

#list is first ad (77) and then the hc (173)
RecontrueAD <- lapply(1:77, function(x) gica_ad$S%*% t(adM[[x]]) )
RecontrueHC <- lapply(1:173, function(x) gica_hc$S%*% t(hcM[[x]]) )
Recontrue <- c(RecontrueHC,RecontrueAD)

ReconotherAD <- lapply(1:77, function(x) gica_hc$S%*% t(adM[[x]]) )
ReconotherHC <- lapply(1:173, function(x) gica_ad$S%*% t(hcM[[x]]) )
Reconother <- c(ReconotherHC,ReconotherAD)


rm(ReconotherAD, ReconotherHC,RecontrueAD,RecontrueHC,splitidx)

SSQtrue <- sapply(1:250, function(x) sum( (dataL[[x]] - Recontrue[[x]])^2 ) )
SSQother <- sapply(1:250, function(x) sum( (dataL[[x]] - Reconother[[x]])^2 ) )

SSQ <- cbind(SSQtrue,SSQtrue)
apply(SSQ,MARGIN = 1, which.min)

############ compute percdiff ########

percd=sapply(1:250, function(x) ((SSQother[[x]] - SSQtrue)[[x]]) / SSQtrue[[x]]  )
#percd*100

lab <- c(rep(1,173),rep(2,77))
diff <- SSQtrue-SSQother
absdiff <- abs(SSQtrue-SSQother)

per2 <- absdiff/ ((SSQtrue+SSQother)/2)

files <- c(fhcl,fadl)
file <- substr(files,1,7)
filenum <- substr(file, 5,7)
data <- cbind(lab, SSQtrue, SSQother, diff,absdiff, percd=percd*100)
data <- data.frame(file,filenum,data)
View(data)


################# checks #############

ADdata <- data[ data$lab==2 , ]
HCdata <- data[ data$lab==1 , ]

ADdata[order(ADdata$diff, decreasing = T), ] #ADdata$V6
HCdata[order(HCdata$diff, decreasing = T), ]


#### check misclus cl2cp5 

setwd("~/surfdrive/Documents/Projects/CICA/Project1/P1.4_Application/P1.4.1_Shark/June27/CICA/nc5/")
load("start_21_ncomp_5_nclus_2.RData")
cica <- cica_res[[1]]
cica$P
Phc <- cica$P[1:173]
Pad <- cica$P[174:250]



# Ranking -----------------------------------------------------------------
checkrankingHC <- cbind(1:173,HCdata[order(HCdata$diff, decreasing = T), ])
checkrankingAD <- cbind(1:77,ADdata[order(ADdata$diff, decreasing = T), ])


# Comparison 5cp2cl solution ----------------------------------------------
# 14 hc mis
sum((Phc==1))
which(Phc==1)
mishc <- file[which(Phc==1)]
match(mishc, checkrankingHC$file)

#50 AD mis
sum(Pad==2)
which(Pad==2)
misad <- file[which(Pad==2)+173]
match(misad, checkrankingAD$file)

#27 ad o
sum((Pad==1))
which(Pad==1)
okad <- file[which(Pad==1)]
match(misad, checkrankingAD$file)


# Selection of files ------------------------------------------------------
adsel <- checkrankingAD[1:20,]$file
hcsel <- checkrankingHC[1:20,]$file

adsel <- paste(adsel, collapse = ",")
hcsel <- paste(hcsel, collapse = ",")

hcselidx <- match(hcsel, substr(files,1,7))
adselidx <- match(adsel, substr(files,1,7))
hcselidx

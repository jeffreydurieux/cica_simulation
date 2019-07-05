library(oro.nifti)
library(fslr)
library(CICA)
library(arf3DS4)
library(ica)
library(MatrixCorrelation)

sall <- readNIfTI("/Users/jeffreydurieux/surfdrive/Data/Parkinson/8_netwerken/all_beckmann_8_and_csf_and_wm.nii.gz")

s1 <- sall[,,45,5]
s2 <- sall[,,45,1]
s3 <- sall[,,45,7]
s4 <- sall[,,45,8]

s1[s1<0] <- 0
s2[s2<0] <- 0
s3[s3<0] <- 0
s4[s4<0] <- 0

image(s1)
image(s2)
image(s3)
image(s4)

S <- cbind(as.vector(s1),as.vector(s2)
           ,as.vector(s3),as.vector(s4))


##### make different cluster maps #######

s11 <- s1
s12 <- s2
s13 <- s3
s14 <- s4

####### cluster 2 #########
s21 <- s1
id <- which(s1!=0, arr.ind = TRUE)
min(id[,1])
max(id[,1])
min(id[,2])
max(id[,2])

s21[ 13:20 ,  ] <- 0 
s21[ 69:76 ,  ] <- 0 
s21[       , 95:100  ] <- 0 

image(s21)

s22 <- s2
image(s22)
id <- which(s2!=0, arr.ind = TRUE)
min(id[,1])
max(id[,1])
min(id[,2])
max(id[,2])

s22[ 31:40 ,  ] <- 0 
s22[ 60:69 ,  ] <- 0 
s22[       , 94:99  ] <- 0 

image(s22)


s23 <- s3
image(s23)
id <- which(s3!=0, arr.ind = TRUE)
min(id[,1])
max(id[,1])
min(id[,2])
max(id[,2])

s23[ 10:15 ,  ] <- 0 
s23[ 39:44 ,  ] <- 0 
s23[       , 38:45  ] <- 0 

image(s23)


s24 <- s4
image(s24)
id <- which(s4!=0, arr.ind = TRUE)
min(id[,1])
max(id[,1])
min(id[,2])
max(id[,2])

s24[ 19:30 ,  ] <- 0 
s24[ 75:80 ,  ] <- 0 
s24[       , 90:98  ] <- 0 

image(s24)


############# cluster 3 #########

s31 <- s1
id <- which(s1!=0, arr.ind = TRUE)
min(id[,1])
max(id[,1])
min(id[,2])
max(id[,2])

s31[ 13:25 ,  ] <- 0 
s31[ 64:76 ,  ] <- 0 
s31[       , 90:100  ] <- 0 

image(s31)

s32 <- s2
image(s22)
id <- which(s2!=0, arr.ind = TRUE)
min(id[,1])
max(id[,1])
min(id[,2])
max(id[,2])

s32[ 32:48 ,  ] <- 0 
s32[ 58:68 ,  ] <- 0 
s32[       , 90:99  ] <- 0 

image(s32)


s33 <- s3
image(s23)
id <- which(s3!=0, arr.ind = TRUE)
min(id[,1])
max(id[,1])
min(id[,2])
max(id[,2])

s33[ 11:20 ,  ] <- 0 
s33[ 30:44 ,  ] <- 0 
s33[       , 30:45  ] <- 0 

image(s33)


s34 <- s4
image(s24)
id <- which(s4!=0, arr.ind = TRUE)
min(id[,1])
max(id[,1])
min(id[,2])
max(id[,2])

s34[ 19:35 ,  ] <- 0 
s34[ 70:80 ,  ] <- 0 
s34[       , 85:98  ] <- 0 

image(s34)

################ save clusters/ nifti's #############
S1 <- cbind(as.vector(s11),as.vector(s12)
           ,as.vector(s13),as.vector(s14))

S2 <- cbind(as.vector(s21),as.vector(s22)
            ,as.vector(s23),as.vector(s24))

S3 <- cbind(as.vector(s31),as.vector(s32)
            ,as.vector(s33),as.vector(s34))

setwd("/Volumes/LaCie/MyData/CICA/Project1/SIM2/PosterExample/")
save(S1,file = "cluster1.Rdata")
save(S2,file = "cluster2.Rdata")
save(S3,file = "cluster3.Rdata")

congru(S1,S2)
congru(S1,S3)
congru(S2,S3)
RV2(S1,S2)
RV2(S1,S3)
RV2(S2,S3)

bmask <- readNIfTI("/Volumes/LaCie/MyData/Masks/MNI/MNI152_T1_2mm_brain_mask.nii.gz")

bmask<- bmask[,,45]

save(bmask, file = "bmaskslice.Rdata")

######## view in fsleyes ########
nif1 <- array(data = 0,dim = c(91,109,91,4))
nif1[,,45,] <-S1
file <- readData("/Users/jeffreydurieux/Documents/MNI152_T1_2mm_brain.nii.gz")

file@fullpath = "/Volumes/LaCie/MyData/CICA/Project1/SIM2/PosterExample/NIFTI/"
file@filename = "ORIS1"
file@dims= c(4,91,109,91,4,1,1,1)
writeData(headinf = file, datavec = as.vector(nif1))

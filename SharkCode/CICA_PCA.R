# Wed Mar 08 11:04:34 2017 ------------------------------
# CICA PCA methods

# normalization 1) none
#               2) intensity nor DONE
#               3) variance normalization (z-scoring) DONE

# Four types:   1) Subject specific DONE -> OK
#               2) Spatial concatenation DONE -> ~
#               3) group data mean DONE -> ~
#               4) group level large PCA DONE
#                   Smith MIGP: DONE -> check it
#                   Smith SMIG: DONE -> OK


#library(CICA)
#library(expm) # matrix square root, complex numbers


# Sat Mar 11 12:27:26 2017 ------------------------------
MatSplit <- function(Data, mats){
  # input: Data = a matrix. mats number of matrices to return
  nr <- nrow(Data)
  nc <- ncol(Data)
  Data <- cbind(Data, rep(1:mats, each = nr/mats ) )
  Data <- lapply( split( Data[, 1:nc], Data[ ,nc+1] ), matrix, ncol=nc)
  return(Data)
}


# Wed Mar 08 11:10:03 2017 ------------------------------
PCA_single_subject <- function(Yi, pc) {
  # Input: Y = a Time X Voxel matrix, pc = number of components
  # Output: A data reduced Time_pc x Voxel data matrix
  # Note that this function also demeans (row-centering) the data

  if( pc > nrow(Yi) ){
    stop("Set pc <= nrow(Yi)")
  }

  Yic <- RowCent(Yi)

  #Covariance mat TxT
  Cov <- Yic%*%t(Yic)

  SVD <- svd(Cov)
  Yred <- t(SVD$v[,1:pc]) %*% Yic

  # calculate cum var
  vaf <- cumsum( SVD$d^2) / sum( SVD$d^2 ) *100

  out <- list()
  out$Yred <- Yred
  out$SVD  <- SVD
  out$VAF <- vaf
  return(out)
}



# Thu Mar  9 13:48:01 2017 ------------------------------
PCA_spatially_concatenated <- function(Data, pc, List = TRUE){
  #input: Data = list of datablocks (TxV), pc= number of reduced t
  #output:
  # row center data
  Data <- RowCent(Data)

  # concatenated voxelwise
  ConcVoxels <- do.call(cbind, Data)

  #SVD on conc. data
  Cov <- ConcVoxels %*% t( ConcVoxels )
  SVD <- svd(Cov)

  # common reducing matrix
  F0bar <- t(SVD$v[,1:pc])

  # temporally concatenate data into Ystar using commo reducing matrix
  YLs<- lapply(Data, function(anom) F0bar%*%anom)
  Yred <- do.call(rbind, YLs)

  out <- list()
  if(List == TRUE){
   out$Yred <- YLs
  }else{
    out$Yred <- Yred
  }
  out$SVD <- SVD
  out$F0bar <- F0bar

  return(out)
}

# Thu Mar  9 16:28:00 2017 ------------------------------
# group data mean

PCA_groupmean <- function(Data, pc) {
  #input : Data = list of data, pc = dim reduced data

  Data <- RowCent(Data)
  MeanData <- Reduce(`+`, Data) / length(Data)

  Cov <- MeanData %*% t(MeanData)
  SVD <- svd(Cov)

  # group mean reducing matrix
  F0bar <- t( SVD$v[, 1:pc] )


  Single_sub_PCA <- lapply(Data, PCA_single_subject, pc = pc)
  FiL <- lapply(seq_along(Single_sub_PCA),
                function(i) t( Single_sub_PCA[[i]]$SVD$v[, 1:pc]) )



  YRed <- lapply( seq_along(FiL),
                  function(i) FiL[[i]] %*% t(SVD$v) %*% Data[[i]] )


  out <- list()
  out$Yred <- YRed
  out$SVDgMean <- SVD
  out$FiL <- FiL
  return(out)
}


# Fri Mar 10 12:08:49 2017 ------------------------------
# Group level pca
PCA_grouplevel <- function(RedData, T2){

  # Concatenate temporally
  if( is.list(RedData) ){
    RedData <- do.call(rbind, RedData)
  }

  #call single subject data
  out <- PCA_single_subject(RedData, T2)

  #subject predicted compressed data
 # Ycompress <- t(out$SVD$v[,1:T2]) %*% out$SVD$v[,1:T2] %*% RedData
  #split into T2 parts
  #parts <- split(data.frame(Ycompress), ceiling(seq_along(1:400)/100))
  #out <- list()
  #out <- compressed()

  return(out)
}

# Sat Mar 11 12:27:45 2017 ------------------------------
# Smith et al large group PCA
PCA_LargeGroup_MIGP <- function(Data, pc) {
  t <- nrow(Data[[1]])
  r <- sample(1:length(Data), size = length(Data))
  Y <- Data[[ r[1] ]]

  for(i in 2:length(Data)){
    Y <- rbind(Y, Data[[ r[i] ]])

    SVD <- svd(Y %*% t(Y),nv = t * 2-1)
    Y <- t(SVD$v) %*% Y
  }

  return( Y[1:pc, ] )
}

# Sat Mar 11 14:47:27 2017 ------------------------------


# SMIG
# get mean of Y

#library used: expm
PCA_LargeGroup_SMIG <- function(Data, pc, iter = 50){
  #meanY can only be computed if all Datablocks have the same dimensions!
  meanY <- Reduce("+", Data) / length(Data)

  if (pc < nrow(Data[[1]]) ){
    SVD <- svd(meanY)
    W <- t( SVD$v[ ,1:pc ] )
  }else{
    W = meanY
  }

  for(i in 1:iter){
    TT <- matrix(data = 0, pc, ncol( Data[[1]] ) )
    RR <- matrix(data = 0, pc, pc)

    Mi <- list()
    for (j in 1:length(Data)) {
      Mi[[j]] <- Data[[j]] %*% t(W)
      TT <- TT + t(Mi[[j]]) %*% Data[[j]]
      RR <- RR + t(Mi[[j]]) %*% Mi[[j]]
    }
    W <- Re( solve( sqrtm (RR + 1e-3) %*% diag(1,pc) ) ) %*% TT
  }

  W <- svd(W)
  W <- t(W$v[,1:pc])

  out <- list()
  out$W <- W
  out$M <- Mi

  return(out)
}


ComputePhiMatrix <- function ( Ytarget , Yhat ) # computes Tucker's phi among columns of input matrix
{
  Out = solve( diag( diag( t(Ytarget) %*% Ytarget ) ) ^ .5 ) %*% t(Ytarget) %*% Yhat %*% solve( diag( diag( t(Yhat) %*% Yhat ) ) ^ .5 )
  return( Out )
}

procr <- function ( X , Ytarget )
{
  # find orthonormal rotation of X that fits Ytarget best in least squares sense

  # INPUT
  #   X (n x p1): matrix that will be rotated (orthogonally)
  #   Ytarget (n x p2): target matrix

  # if p1 >= p2: Tmat is orthonormal
  # if p1 < p2: t(Tmat) is orthonormal (Tmat is not orthonormal)



  if ( dim(X)[1] == dim(Ytarget)[1]  )
  {
    #source("ComputePhiMatrix.R")
    temp = svd( t(X) %*% Ytarget )
    #LowestDim = min( dim( t(X) %*% Ytarget ) )
    #p = temp$u[ , 1:LowestDim ]
    #q = temp$v[ , 1:LowestDim ]
    #Tmat = p %*% t(q) # rotation matrix (p1 x p2)
    Tmat = temp$u %*% t(temp$v)
    Yhat = X %*% Tmat # (n x p2)
    RotationMatrix = Tmat
    Fit = sum( ( Yhat - Ytarget ) ^ 2 )
    PhiMatrix = ComputePhiMatrix( Ytarget , Yhat )

    RotationMatrixNonorthogonal = ginv( t(X) %*% X ) %*% t(X) %*% Ytarget
    ###RotationMatrixNonorthogonal = solve( t(X) %*% X ) %*% t(X) %*% Ytarget
    YhatNonorthogonal = X %*% RotationMatrixNonorthogonal
    FitNonorthogonal = sum( ( YhatNonorthogonal - Ytarget ) ^ 2 )
    PhiMatrixNonorthogonal = ComputePhiMatrix( Ytarget , YhatNonorthogonal )

    Out = list()
    Out$Yhat = Yhat #orthogonally rotated X matrix
    Out$RotationMatrix = RotationMatrix
    Out$Fit = Fit
    Out$PhiMatrix = PhiMatrix
    Out$TuckerCongruence = mean( abs( diag( PhiMatrix ) ) )

    Out$NonOrthogonal$Yhat = YhatNonorthogonal
    Out$NonOrthogonal$RotationMatrix = RotationMatrixNonorthogonal
    Out$NonOrthogonal$Fit = FitNonorthogonal
    Out$NonOrthogonal$PhiMatrix = PhiMatrixNonorthogonal
    Out$NonOrthogonal$TuckerCongruence = mean( abs( diag( PhiMatrixNonorthogonal ) ) )
  }
  else
  {
    cat(" ",fill=TRUE)
    cat("X and Ytarget should have the same number of rows", fill=TRUE )
    cat(" ",fill=TRUE)
    Out = list()
  }
  return(Out)
}

# ###############
# # Sat Mar 11 12:28:12 2017 ------------------TESTS
# library(CICA)
#
# #generate 5 datablocks with 10 spatial maps, time course length =100 -> 100x1000 time x voxel data
# # percentage of gaussian error = 5%
#
# D <- SIM_CICA(5,10, 1000,1,NULL,100,0.05, IrrSig = 0)
# Data <- lapply(D$XE, t) # SIM_CICA gives voxel x time -> transpose datablocks first to time x voxel
# Data <- RowCent(Data) # row center the data. center T
#
#
#
#
# ### Predicted data after dimension reduction single subject
# Test=PCA_single_subject(Data[[1]], 10) #reduce to 10 pc
# Pred <- Test$SVD$v%*%Test$SVD$v[,1:10] %*% Test$Yred
# RV2(t(Pred), t(Data[[1]]))
#
# ic <- icafast( t(Test$Yred), 10)
# apply( abs( congru( t(D$SR[[1]]) , ic$S) ), 2, max )
#
#
# ### Predicted data after dimension reduction spatially concatenated
# Test = PCA_spatially_concatenated(Data, 10 , T)
# Pred1 <- t(Test$F0bar)%*%Test$F0bar%*%RowCent(Data[[1]])
# Pred2 <- t(Test$F0bar)%*%Test$F0bar%*%RowCent(Data[[2]])
# RV2( t(Pred1), t(Data[[1]]))
# RV2( t(Pred1), t(RowCent( Data[[1]] )))
# RV2( t(Pred2), t(Data[[2]]))
# RV2( t(Pred2), t(RowCent( Data[[2]] )))
#
# ic <- icafast( t(Test$Yred[[1]]) , 10)
# apply( abs( congru( t(D$SR[[1]]) , ic$S) ), 2, max )
#
# #for all subs
# conM <- vector("numeric", 5)
# for(i in 1:5){
#   ic <- icafast( t(Test$Yred[[i]]) , 10)
#   conM[i] <- mean(apply( abs( congru( t(D$SR[[1]]) , ic$S) ), 2, max ) )
# }
# conM
#
#
# ### Predicted data after dimension reduction group mean
# # not sure about predicted data after dimension reduction part
# Data <- RowCent(Data)
# Test <- PCA_groupmean(Data,10)
# Pred1 <-   Test$SVDgMean$v[, 1:10] %*% Test$FiL[[1]] %*% t(Test$FiL[[1]]) %*% t(Test$SVDgMean$v[,1:10]) %*% Data[[1]]
# Pred2 <-   Test$SVDgMean$v[, 1:10] %*% Test$FiL[[2]] %*% t(Test$FiL[[2]]) %*% t(Test$SVDgMean$v[,1:10]) %*% Data[[2]]
# RV2(t(Pred1) , t(Data[[1]]))
# RV2(t(Pred2) , t(Data[[2]]))
# RV2(t(Pred2) , t(Data[[3]]))
# RV2(t(Pred2) , t(Data[[4]]))
# RV2(t(Pred2) , t(Data[[5]]))
#
# ic <- icafast( t(Test$Yred[[1]]), 10 )
# apply( abs( congru( t(D$SR[[1]]) , ic$S) ), 2, max )
#
# conM <- vector("numeric", 5)
# for(i in 1:5){
#   ic <- icafast( t(Test$Yred[[i]]) , 10)
#   conM[i] <- mean(apply( abs( congru( t(D$SR[[1]]) , ic$S) ), 2, max ) )
# }
# conM
#
#
#
# ### MIGP
#
#
# ###SMIG
#
# #what about the iteration within the code?!?
# # note: if only one iteration -> method equal to hyvarinen & smith 2012
# Test <- PCA_LargeGroup_SMIG(Data,10,iter = 1)
# ic1 <- icafast(t(Reduce("+", Data) / length(Data)), 10)
# apply(abs(congru(t(D$SR[[1]]),ic1$S)),2,max)
#
# ic2 <- icafast(t(Test),10)
# apply(abs(congru(t(D$SR[[1]]),ic2$S)),2,max)
#

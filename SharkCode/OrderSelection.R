# singular values ^ 2 = eigenvalues


# Log likelihood formula for PCA model

LLfun <- function(singularvalues, samplesize, pc){
  # input: singular values from svd output, pc = desired num components
  # output maximum log likelihood value
  lambda <- singularvalues^2
  A <- samplesize/2
  B <- tail(cumprod(lambda[(pc+1):length(lambda)] ^
                      (1 / (length(lambda) - pc ) ) ), n = 1)

  C <- (1 / (length(lambda-pc)) ) *
    tail(cumsum( lambda[(pc+1):length(lambda)] ) , n = 1)

  LL <- A * log( (B/C) ^ (length(lambda) - pc) )
  return(LL)
}

FreeParameters <- function(Dim, pc){
  out <- 1 + (Dim*pc) - (0.5 * pc * (pc - 1))
  return(out)
}

# AIC
AIC <- function(Dim, pc, LL) {
  return( (-2*LL) + (2 * FreeParameters(Dim, pc) ) )
}

BIC <- function(Dim, samplesize, pc, LL) {
  return( (-2*LL) + (log(samplesize) * FreeParameters(Dim, pc) ) )
}

# KIC (kullback leibler info criteria)
# MDL

# library(CICA)
# D <- SIM_CICA(1,10,1000,1,NULL,100,0.50)
# X <- D$XE
# svd <- svd(X[[1]])
# val <- svd$d
#
# lltest <- numeric()
# aic <- numeric()
# bic <- numeric()
# for(i in 2:50){
#   lltest[i] <- LLfun(singularvalues = val, samplesize = 2000,i)
#   aic[i] <- AIC(100,10,lltest[i])
#   bic[i] <- BIC(Dim = 100,samplesize = 2000,pc = 10,LL = lltest[i])
# }
# plot(lltest,type="l")
# plot(aic, type="l")
# plot(bic, type="l")
#

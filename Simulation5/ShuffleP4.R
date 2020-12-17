shuffleR4 <- function(X, R){
  
  X1 <- list()
  X2 <- list()
  X3 <- list()
  X4 <- list()
  a <- seq(from = 1, to = 40, by = 4)
  a <- a[1:10]
  b <- seq(from = 2, to = 40, by = 4)
  b <- b[1:10]
  c <- seq(from = 3, to = 40, by = 4)
  c <- c[1:10]
  d <- seq(from = 4, to = 40, by = 4)
  d <- d[1:10]
  
  iter <- 1
  for(i in 1:length(a)){
    X1[[i]] <- X[[a[iter]]]
    X2[[i]] <- X[[b[iter]]]
    X3[[i]] <- X[[c[iter]]]
    X4[[i]] <- X[[d[iter]]]
    iter <- iter + 1
  }

  newX <- c(X1,X2,X3,X4)
  
  P <- c(a,b,c,d)
  newP <- numeric()
  for(i in 1:40){
    if(P[i] <= 10){
      newP[i] <- 1
    }else if(P[i]>10 & P[i]<=20){
      newP[i] <- 2
    }else if(P[i]>20 & P[i]<=30){
      newP[i] <- 3
    }else{
      newP[i] <- 4
    }
  }

  out <- list()  
  out$newX <- newX
  out$newP <- newP 
  return(out)
}

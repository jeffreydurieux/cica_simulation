# Thu Jul  4 10:02:34 2019
# Author: Jeffrey Durieux, MSc

# What: functions for generating C-ICA data

### safeguard if ica package is not installed
list.of.packages <- c("ica")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) stop(list.of.packages,' package not installed \n')# "uncomment this for auto install" install.packages(new.packages)

Simulate_CICA <- function(N = 40, V = 1000, D = 100, Q = 2, R = 2, E = 0.1){
  # Note input args; D = time Dimension so that it is not confused with T (True)
  
  # generate cluster specific components (FC patterns)
  # use icasamp from ica package
  # double exponential = b
  # Dimensions: R lists of V X Q matrices
  Qr <- lapply(seq_along(1:R), function(anom) {
    replicate(n = Q, expr = icasamp(dname = 'b', query = 'rnd', nsamp = V))
  })
  
  # generate random time courses
  # dimension: R lists of T X Q matrices
  NperClus <- N / R
  Air <- lapply(seq_along(1:R), function(anom1) {
    
    lapply(seq_along(1:NperClus), function(anom2) {
      replicate(n = Q, runif(n = D, min = -2, max = 2))
    }) 
    
  }) 
  
  # Mix data -> Qr %*% t(Air)
  X <- lapply(seq_along(Qr), function(anom1){
    lapply( seq_along( Air[[anom1]]), function(anom2){
      Qr[[anom1]] %*% t(Air[[anom1]][[anom2]] )
    })  
  })
  
  # concatenate lists into one large list
  X <- unlist(X, recursive = F)
  
  # add error
  Xe <- lapply(X, FUN = addError, error = E)
  
  Out <- list() 
  Out$P <- rep(1:R, each = NperClus)
  Out$Qr <- Qr
  Out$Air <- Air
  Out$X <- X
  Out$Xe <- Xe
  return(Out)
}


### helper functions for generating gaussian
addError<-function(datablock,error)
{
  
  errorM<-replicate(ncol(datablock),rnorm(nrow(datablock)))
  errorM<-SSequal(errorM,datablock)
  errorlevel<-error/(1-error)
  
  res<-datablock + (errorM * sqrt(errorlevel))
  return(res)
}

SSequal<-function(m1,m2)
{
  res<-(m1/sqrt(sum(m1^2)) * sqrt(sum(m2^2))) #R
  return(res)
}


###### Tests  ########
#test <- Simulate_CICA()
#test <- Simulate_CICA(N = 40, V = 500, D = 2, Q = 2, R = 4, E = 0.5)
#test <- Simulate_CICA(N = 40, V = 500, D = 64, Q = 20, R = 2, E = 0.5)







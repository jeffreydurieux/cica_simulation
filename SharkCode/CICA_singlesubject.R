# example data

singlesub <- function(ncomp=5 ,nrow=1000, ncol=8, dist="b", error=0, sd = 1){

  if(ncomp > ncol){
    stop("Too many components. Set comp <= ncol")
  }

  S <- replicate(n = ncomp, icasamp(dname = dist, query = "rnd", nsamp = nrow))
  A <- matrix( rnorm( ncol*ncomp,sd = sd), ncol, ncomp)
  X <- tcrossprod(S,A)
  Xe <- addError(X, error)

  out <- list()
  out$Xe <- Xe
  out$X <- X
  out$S  <- S
  out$A  <- A

  return(out)
}


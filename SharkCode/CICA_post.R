## Functions for post-CICA analysis

# Mon Apr 24 13:13:17 2017 ------------------------------

# add zeros

addZeros <- function(idx, data, dim){
  zeros <- matrix(data = 0, nrow = dim[4], ncol = dim[1]*dim[2]*dim[3]) # t x v
  zeros[, idx] <- data
  zeros <- to4d(zeros, dim = dim)
  nif <- as.nifti(zeros)
  return(nif)
}

# z threshold


Thres_img <- function(img, rm_zeros, abs){

  if(rm_zeros == TRUE){
    img <- (img - mean(img)) / sd(img)
    img[img < -2 | img > 2] = NA
    img[img > 0.01 & img < 0.01] = NA
    return(img)
  }else if(abs == TRUE){
    img <- (abs(img) - mean(abs(img))) / sd(abs(img))
    img[img > 2] = NA
    return(img)
  }else{
    img <- (img - mean(img)) / sd(img)
    img[img < -2 | img > 2] = NA
  }
}



ViewSimMap <- function(S, map = 1, dim = c(64,64)){
  lattice::levelplot(matrix(S[,map],nrow=dim[1], ncol=dim[2]), col.regions=gplots::colorpanel(100, low="grey",mid = 'red', high="yellow"))
}







# to nifti


# mlR_loc.R
# Author: Matthew J. Holland
# Updated: 2016/10/12
# Description:
# A few fixed-point type routines for computing M-estimates of location, similar to Algorithm 2.

ml_loc.s <- function(l, s, psi_fn = ml_psi.lgst, thres = 1e-3, iters = 50L){
  
  n <- length(l)
  diff <- Inf
  loc_old <- numeric(1L)
  loc_new <- sum(l) / n # initialize to sample mean.
  for(i in 1:iters){
    loc_old <- loc_new
    val <- {l - loc_old} / s
    loc_new <- loc_old + s * sum(psi_fn(val)) / n # update, with scale factor.
    diff <- abs({loc_new-loc_old})
    if (diff <= thres){
      return(loc_new)
    }
  }
  
  return(loc_new)
}

ml_loc.nos <- function(l, s, psi_fn = ml_psi.lgst, thres = 1e-3, iters = 50L){
  
  n <- length(l)
  diff <- Inf
  loc_old <- numeric(1L)
  loc_new <- sum(l) / n # initialize to sample mean.
  for(i in 1:iters){
    loc_old <- loc_new
    val <- {l - loc_old} / s
    loc_new <- loc_old + sum(psi_fn(val)) / n # update, without scale factor.
    diff <- abs({loc_new-loc_old})
    if (diff <= thres){
      return(loc_new)
    }
  }
  
  return(loc_new)
}




# mlR_scale.R
# Author: Matthew J. Holland
# Updated: 2016/10/12
# Description:
# A few fixed-point type routines for computing M-estimates of scale, similar to Algorithm 2.

ml_s.A <- function(l, chi_fn, ..., thres = 1e-3, iters = 50L){
  
  n <- length(l)
  diff <- Inf
  s_old <- numeric(1L)
  s_new <- sqrt(sum({l - sum(l)/n}^2)/{n-1}) # initialize to sample standard deviation.
  pivot <- sum(l)/n
  num <- l-pivot
  chi_0 <- chi_fn(0)
  
  for(i in 1:iters){
    s_old <- s_new
    val <- num/s_old
    s_new <- s_old * sqrt({1 - sum(chi_fn(val))/{chi_0*n}}) # update, with square root.
    diff <- abs({s_new-s_old})
    if (diff <= thres){
      return(s_new)
    }
  }
  
  return(s_new)
}

ml_s.B <- function(l, chi_fn, ..., thres = 1e-3, iters = 50L){
  
  n <- length(l)
  diff <- Inf
  s_old <- numeric(1L)
  s_new <- sqrt(sum({l - sum(l)/n}^2)/{n-1}) # initialize to sample standard deviation.
  pivot <- sum(l)/n
  num <- l-pivot
  chi_0 <- chi_fn(0)
  
  for(i in 1:iters){
    s_old <- s_new
    val <- num/s_old
    s_new <- s_old * {1 - sum(chi_fn(val))/{chi_0*n}} # update, without square root.
    diff <- abs({s_new-s_old})
    if (diff <= thres){
      return(s_new)
    }
  }
  
  return(s_new)
}

ml_s.C <- function(l, chi_fn, ..., thres = 1e-3, iters = 50L){
  
  n <- length(l)
  diff <- Inf
  s_old <- numeric(1L)
  s_new <- sqrt(sum({l - sum(l)/n}^2)/{n-1}) # initialize to sample standard deviation.
  pivot <- sum(l)/n
  num <- l-pivot
  chi_0 <- chi_fn(0)
  
  for(i in 1:iters){
    s_old <- s_new
    val <- num/s_old
    s_new <- s_old + sum(chi_fn(val))/n # update, of addition form.
    diff <- abs({s_new-s_old})
    if (diff <= thres){
      return(s_new)
    }
  }
  
  return(s_new)
}



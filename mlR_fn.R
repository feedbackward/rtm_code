
# mlR_fn.R
# Author: Matthew J. Holland
# Updated: 2016/10/12
# Description:
# A collection of functions useful for testing routines of the proposed form in R. 


#### Functions related to computing rho, psi, eta. Both valid (via Defn 1) and invalid examples included. ####

# NOTE: ml_psi.* functions here are ready to use with location routines ml_loc.*.

## Helper function for integrating psi functions. ##

ml_psi.toint <- function(u, psi, ...){
  
  # Just a helper function. Pass psi fn in to this function, and integrate up to the desired value (u).
  
  if (u >= 0){
    int.start <- 0
    int.end <- u
  } else{
    int.start <- u
    int.end <- 0
  }
  
  abs(integrate(f = psi, lower = int.start, upper = int.end, ...)$value)
}


## Helper functions for the M-est using truncator of Catoni (2010, 2012), Audibert and Catoni (2011). ##

ml_catlow <- function(x){
  -log1p(-{x-x^2/2})
}

ml_catlow_d1 <- function(x){
  {1-x}/{1-x+x^2/2}
}

ml_catup <- function(x){
  log1p({x+x^2/2})
}

ml_catup_d1 <- function(x){
  {1+x}/{1+x+x^2/2}
}


## The widest truncator of Catoni (2012) ##

ml_psi.catwide <- function(x){
  
  out <- x
  idx_nonneg <- which(out >= 0)
  idx_neg <- which(out < 0)
  out[idx_nonneg] <- ml_catup(out[idx_nonneg])
  out[idx_neg] <- ml_catlow(out[idx_neg])
  
  return(out)
}

ml_rho.catwide <- function(x){
  unlist(lapply(X = x, FUN = ml_psi.toint, psi = ml_psi.catwide))
}

ml_wt.catwide <- function(x){
  noz <- abs(x) + 1e-50
  ml_psi.catwide(noz)/noz
}

ml_eta.catwide <- function(x){
  out <- x
  idx_nonneg <- which(out >= 0)
  idx_neg <- which(out < 0)
  out[idx_nonneg] <- ml_catup_d1(out[idx_nonneg])
  out[idx_neg] <- ml_catlow_d1(out[idx_neg])
  return(out)
}


## The narrowest truncator of Catoni (2012). ##

ml_psi.catnar <- function(x){
  
  out <- x
  idx_nonneg <- which(out >= 0)
  idx_neg <- which(out < 0)
  idx_small <- which(out <= -1)
  idx_large <- which(out >= 1)
  out[idx_nonneg] <- ml_catlow(out[idx_nonneg])
  out[idx_neg] <- ml_catup(out[idx_neg])
  out[idx_small] <- -log(2)
  out[idx_large] <- log(2)
  
  return(out)
}

ml_rho.catnar <- function(x){
  unlist(lapply(X = x, FUN = ml_psi.toint, psi = ml_psi.catnar))
}

ml_wt.catnar <- function(x){
  noz <- abs(x) + 1e-50
  ml_psi.catnar(noz)/noz
}

ml_eta.catnar <- function(x){
  out <- x
  idx_nonneg <- which(out >= 0)
  idx_neg <- which(out < 0)
  idx_small <- which(out <= -1)
  idx_large <- which(out >= 1)
  out[idx_nonneg] <- ml_catlow_d1(out[idx_nonneg])
  out[idx_neg] <- ml_catup_d1(out[idx_neg])
  out[idx_small] <- 0
  out[idx_large] <- 0
  return(out)
}


## The usual quadratic function. ##

ml_rho.l2 <- function(x){
  x^2 / 2
}

ml_psi.l2 <- function(x){
  x
}

ml_eta.l2 <- function(x){
  rep(1, length(x))
}

ml_wt.l2 <- function(x){
  rep(1, length(x))
}


## The usual absolute value function. ##

ml_rho.l1 <- function(x){
  abs(x)
}

ml_psi.l1 <- function(x){
  sign(x)
}

ml_eta.l1 <- function(x){
  rep(0, length(x))
}

ml_wt.l1 <- function(x){
  1/abs(x)
}


## The log(cosh) function ##

ml_rho.lcosh <- function(x){
  log(cosh(x))
}

ml_psi.lcosh <- function(x){
  tanh(x)
}

ml_eta.lcosh <- function(x){
  1 / {cosh(x)^2}
}

ml_wt.lcosh <- function(x){
  noz <- abs(x) + 1e-50
  ml_psi.lcosh(noz)/noz
}


## The "square root algebraic" function. ##

ml_rho.algsq <- function(x){
  2*sqrt({1+x^2/2}) - 2
}

ml_psi.algsq <- function(x){
  x / sqrt({1+x^2/2})
}

ml_eta.algsq <- function(x){
  u <- 1 + x^2/2
  out <- {1 - x^2/{2*u}} / sqrt(u)
  return(out)
}

ml_wt.algsq <- function(x){
  1 / sqrt({1+x^2/2})
}


## The "fair" function, cited from Rey (1983, 6.4.5). ##
# Origin of "fair" is unknown, but was in well-known ROSEPACK library; see Holland and Welsch (1977).

ml_rho.fair <- function(x, c = 1.3998){
  u <- abs(x)/c
  c^2 * {u - log({1+u})}
}

ml_psi.fair <- function(x, c = 1.3998){
  x / {1+abs(x)/c}
}

ml_eta.fair <- function(x, c = 1.3998){
  u <- abs(x)/c
  out <- {1-u}/{1+u}
  return(out)
}

ml_wt.fair <- function(x, c = 1.3998){
  1 / {1+abs(x)/c}
}


## The Huber function, originally proposed in Huber (1964) ##

ml_rho.huber <- function(x, c = 1.345){
  out <- x
  u <- abs(x)/c
  idx_big <- which(u >= 1)
  idx_small <- which(u < 1)
  out[idx_big] <- c^2 * {u[idx_big] - 1/2}
  out[idx_small] <- x[idx_small]^2 / 2
  out
}

ml_psi.huber <- function(x, c = 1.345){
  out <- x
  u <- abs(x)/c
  idx_big <- which(u >= 1)
  idx_small <- which(u < 1)
  out[idx_big] <- c * sign(x[idx_big])
  out[idx_small] <- x[idx_small]
  out
}

ml_eta.huber <- function(x, c = 1.345){
  out <- x
  u <- abs(x)/c
  idx_big <- which(u >= 1)
  idx_small <- which(u < 1)
  out[idx_big] <- 0
  out[idx_small] <- 1
  out
}

ml_wt.huber <- function(x, c = 1.345){
  out <- x
  u <- abs(x)/c
  idx_big <- which(u >= 1)
  idx_small <- which(u < 1)
  out[idx_big] <- c / abs(x[idx_big])
  out[idx_small] <- 1
  out
}


## Modified Huber function, from Rey (1983, 6.4.4). Has continuous second derivative. ##

ml_rho.hmod <- function(x, c = 1.2107){
  out <- x
  u <- abs(x)/c
  idx_big <- which(u >= pi/2)
  idx_small <- which(u < pi/2)
  out[idx_big] <- u[idx_big] + 1 - pi/2
  out[idx_small] <- 1-cos(x[idx_small]/c)
  c^2 * out
}

ml_psi.hmod <- function(x, c = 1.2107){
  out <- x
  u <- abs(x)/c
  idx_big <- which(u >= pi/2)
  idx_small <- which(u < pi/2)
  out[idx_big] <- sign(x) / c
  out[idx_small] <- sin(x/c) / c
  c^2 * out
}

ml_eta.hmod <- function(x, c = 1.2107){
  out <- x
  u <- abs(x)/c
  idx_big <- which(u >= pi/2)
  idx_small <- which(u < pi/2)
  out[idx_big] <- 0
  out[idx_small] <- 2*cos(x/c)
  out
}

ml_wt.hmod <- function(x, c = 1.2107){
  noz <- abs(x) + 1e-50
  ml_psi.hmod(x = noz, c = c)/noz
}


## Gudermannian function, appears in Abramowitz and Stegun. ##

ml_psi.gud <- function(x){
  2 * atan(exp(x = x)) - pi/2
}

ml_rho.gud <- function(x){
  unlist(lapply(X = x, FUN = ml_psi.toint, psi = ml_psi.gud))
}

ml_wt.gud <- function(x){
  noz <- abs(x) + 1e-10
  ml_psi.gud(noz)/noz
}

ml_eta.gud <- function(x){
  2*exp(x)/{exp(2*x)+1}
}


## Logistic function ##

ml_psi.lgst <- function(x, c1 = 4, c2 = 1){
  c1 / {1+exp(-c2*x)} - c1/2
}

ml_rho.lgst <- function(x, c1 = 4, c2 = 1){
  unlist(lapply(X = x, FUN = ml_psi.toint, psi = ml_psi.lgst, c1 = c1, c2 = c2))
}

ml_wt.lgst <- function(x, c1 = 4, c2 = 1){
  noz <- abs(x) + 1e-10
  ml_psi.lgst(x = noz, c1 = c1, c2 = c2)/noz
}

ml_eta.lgst <- function(x, c1 = 4, c2 = 1){
  c1*c2*exp(-c2*x) / {1+exp(-c2*x)}^2
}


## The rho function whose output is inverse tangent (atan) ##

ml_rho.atan <- function(x){
  x*atan(x) - log({1+x^2})/2
}

ml_psi.atan <- function(x){
  atan(x)
}

ml_wt.atan <- function(x){
  noz <- abs(x) + 1e-50
  ml_psi.atan(noz)/noz
}

ml_eta.atan <- function(x){
  1 / {1+x^2}
}


#### Functions for carrying out the scaling procedure denoted chi in paper. ####

# NOTE: functions here of form ml_chi.* are ready to be used with scale estimation routines ml_s.*.

ml_chi.tukey <- function(x, beta, c = 1.547){
  out <- x
  idx_big <- which(abs(x) >= c)
  idx_small <- which(abs(x) < c)
  out[idx_big] <- c^2/6 - beta
  out[idx_small] <- x^6/{6*c} - x^4/{2*c} + x^2/{2*c} - beta
}

# From the example in paper.
ml_chi.lp <- function(x, beta, c = 0.5){
  c*x^{1+c} - beta
}

# Quadratic psi, originally given for Huber function in Huber (1964).
ml_chi.quad <- function(x, psi = ml_psi.l2, beta, ...){
  psi(x, ...)^2 - beta
}

# Log-quadratic psi.
ml_chi.lquad <- function(x, psi = ml_psi.l2, beta, ...){
  log({1+psi(x, ...)^2}) - beta
}

# Standard form, as given in Appendix.
ml_chi.huber <- function(x, rho = ml_rho.l2, psi = ml_psi.l2, beta, ...){
  psi(x, ...)*x - rho(x, ...) - beta
}

# Median absolute deviations from zero.
ml_chi.madzero <- function(x){
  sign({abs(x = x)-1})
}





# mlR_rtn.R
# Author: Matthew J. Holland
# Updated: 2016/10/12
# Description:
# Code for the computational routine described in the paper. All parameters are preset as described.


#### Code to be executed prior to running the main routine. ####

# Required packages.
require(stats)
require(Rcpp)
require(RcppArmadillo)

# Required functions.
load(file = "./mlR_all.RData") # fill in as needed.

# Compile the associated C++ source and load R wrappers into the environment.
sourceCpp(file = "./mlC_loc.cpp") # fill in as needed.
sourceCpp(file = "./mlC_scale.cpp" ) # fill in as needed.

# Function for carrying out calibration on an arbitrary rho fn (given eta = second derivative).
k_q <- function(q, eta, ...){
  toroot <- function(x, myfn = eta, ...){
    myfn(x, ...) - q
  }
  int <- c(0, 100) # if required, make the interval wider.
  uniroot(f = toroot, interval = int, ...)$root
}
ref_val <- 1.345
q_val <- 0.75

# OLS Initializer.
init.ols <- function(y, X){
  coef <- unname(lm(formula = y ~ X - 1)$coefficients) # Do not add an intercept term.
  resid <- y - drop(X %*% coef)
  list(coef = coef, resid = resid)
}

init.sd <- function(resid, alt.loss){
  sd(x = alt.loss(resid))
}


#### Main routine to pass carrier matrix and response vector. ####

# Procedure for robust target minimizing routine.
main_routine <- function(y, X, thres = 1e-3, alt.psiwls = ml_wt.lgst,
                         altC.est = est_lgst, altC.scale = scale_madmean,
                         alt.loss = function(x){x^2},
                         alt.init = init.ols,
                         adjval = k_q(q = q_val, eta = ml_eta.lgst),
                         iters = 50L){
  
  n <- length(y)
  d <- ncol(X)
  
  # Initialize weights.
  fit <- alt.init(y = y, X = X)
  w_init <- fit$coef
  
  w_old <- w_init
  y_est <- drop(X %*% w_old)
  res <- y - y_est
  s <- altC.scale(x = alt.loss(res), n = n)
  objval_old <- altC.est(x = alt.loss(res), n = n, s = s/adjval, thres = thres, iters = iters)
  diff <- Inf
  num_iter <- 0
  
  while (diff > thres){
    num_iter <- num_iter + 1
    if (num_iter > iters){
      return(w_old) # so it doesn't run too long.
    }
    pivot <- sum(alt.loss(res))/n # gamma in the paper.
    vals <- {alt.loss(res)-pivot}/s # values for doing the re-weighting.
    weights <- alt.psiwls(adjval*vals) # new weights using psi function.
    fit <- lm.wfit(x = X, y = y, w = weights) # run IRLS using computed weights.
    w_new <- fit$coefficients # new candidate parameters for which we check theta-hat criterion.
    y_est <- drop(X %*% w_new) # new estimate.
    res <- y - y_est # update the residuals using new estimate.
    s <- altC.scale(x = alt.loss(res), n = n)
    objval_new <- altC.est(x = alt.loss(res), n = n, s = s/adjval, thres = thres, iters = iters)
    diff <- objval_old - objval_new
    if (diff <= 0){
      return(w_old) # if does not improve the theta-hat objective, stop.
    }
    w_old <- w_new
    objval_old <- objval_new
  }
  
  return(w_old) # If do one or more updates, return the updated parameter estimates.
}





















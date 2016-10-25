
# mlR_competitors.R
# Author: Matthew J. Holland
# Updated: 2016/10/12
# Description:
# Code for the competitive benchmarks described in the paper.

# Required packages
require(MASS)
require(quantreg)
require(Rcpp)
require(RcppArmadillo)
require(robustbase)
require(stats)


## Ordinary least squares (ERM-L2) ##

linreg_ols <- function(y, X){
  ml_ols(y = y, X = X) # return the coefficients, that's it.
}

ml_ols <- function(y, X){
  
  # Do not fit an intercept term.
  unname(lm(formula = y ~ X - 1)$coefficients)
}


## Least absolute deviations (ERM-L1) ##

linreg_lad <- function(y, X){
  
  # The special case of quantile regression for the 0.5 quantile, also known
  # as median regression, is precisely least absolute deviations regression.
  
  fit <- rq(formula = y ~ X - 1, tau = 0.5, method = "br") # no intercept, as model (known) has no intercept.
  
  unname(fit$coefficients)
}



## Geometric median routine ##

linreg_geomed <- function(y, X){
  
  fname <- "linreg_geomed"
  
  n <- dim(X)[1]
  d <- dim(X)[2]
  
  # Set the number of subsets to partition into, num_segments.
  if (floor(n/2) <= d){
    cat("Error (", fname, "): not enough data to do a partition.", "\n", sep = "")
    return(NULL)
  } else {
    num_segments <- max(2, floor(n/{2*d}))
  }
  
  w_Candidates <- matrix(0, nrow = num_segments, ncol = d)
  
  # split observations into NUM_SEGMENTS segments, (via Sec 5 of Minsker (2015)), and make segment-wise estimates.
  init_seglen <- floor(x = {n / num_segments})
  remainder <- n - num_segments*init_seglen
  
  if (remainder == 0){
    
    for (j in 1:num_segments){
      idx_start <- {j-1}*init_seglen + 1
      idx_stop <- j*init_seglen
      idx <- idx_start:idx_stop
      y_sub <- y[idx]
      X_sub <- X[idx, , drop = FALSE]
      
      w_Candidates[j, ] <- ml_ols(y = y_sub, X = X_sub)
      
      #cat("Length:", length(idx), "Start:", idx_start, "Stop:", idx_stop, "\n", sep = " ")
    }
    
  } else {
    
    for (j in 1:remainder){
      idx_start <- {j-1}*{init_seglen+1} + 1
      idx_stop <- j*{init_seglen+1}
      idx <- idx_start:idx_stop
      y_sub <- y[idx]
      X_sub <- X[idx, , drop = FALSE]
      
      w_Candidates[j, ] <- ml_ols(y = y_sub, X = X_sub)
      
      #cat("Length:", length(idx), "Start:", idx_start, "Stop:", idx_stop, "\n", sep = " ")
    }
    for (j in {remainder+1}:num_segments){
      idx_start <- idx_stop + 1
      idx_stop <- idx_stop + init_seglen
      idx <- idx_start:idx_stop
      y_sub <- y[idx]
      X_sub <- X[idx, , drop = FALSE]
      
      w_Candidates[j, ] <- ml_ols(y = y_sub, X = X_sub)
      
      #cat("Length:", length(idx), "Start:", idx_start, "Stop:", idx_stop, "\n", sep = " ")
    }
    
  }
  
  # Output the geometric median of the segment-wise estimates
  ml_geomedian(X = w_Candidates)
}


# Subroutine for computing the geomedian.

ml_geomedian <- function(X, thres = 1e-3, max_iter = 100){
  
  fname <- "ml_geomedian"
  
  # Assumes X is a n by d matrix of d-dimensional real vectors.
  # Here we seek the "Geometric Median." Use algorithm of Vardi and Zhang (2000). More refs in Minsker (2015), sec 5.
  # The eta weights (in Vardi and Zhang, 2000) can all be set to 1 as a convenient special case,
  # matching Minsker's formulation.
  
  t_X <- t(X)
  
  # first, let's check that not all the observations are identical. If so, return to sender that vector.
  check <- sum(abs(X[1, ] - t_X)) == 0
  
  if (check){
    return(X[1, ])
  }
  
  # Some helper functions.
  
  weiszfeld_term <- function(u, X){
    
    diffs <- 1 / sqrt(colMeans({{t(X) - u}^2}))
    
    out <- colSums({X * diffs}) / sum(diffs)
    
    return(out)
  }
  
  r_val <- function(u, t_X){
    
    diffs <- 1 / sqrt(colMeans({{t_X - u}^2}))
    big_R <- colSums({t({t_X - u}) * diffs})
    out <- sqrt(sum(big_R^2)/length(big_R))
    
    return(out)
  }
  
  # The main routine.
  
  old_u <- colMeans(x = X)
  
  diff <- Inf
  
  iter <- 1
  
  while (diff > thres && iter <= max_iter){
    
    check <- which(colSums(abs(old_u - t_X)) == 0)
    if (length(check) > 0){
      
      # In this case, the current iteration hit one of the samples.
      r_inv <- 1/r_val(u = old_u, t_X = t_X)
      new_u <- max(0, {1-r_inv}) * weiszfeld_term(u = old_u, X = X) + min(r_inv, 1) * old_u
      
    } else {
      
      
      # In this case, current iteration different from all samples.
      new_u <- weiszfeld_term(u = old_u, X = X)
      
    }
    
    diff <- sqrt(sum({new_u - old_u}^2)/length(new_u))
    old_u <- new_u
    #cat("Iteration #", iter, "\n", sep = "")
    #cat("Diff val: ", diff, "\n", sep = "")
    iter <- iter + 1
  }
  
  return(new_u)
}


## Hsu and Sabato median-of-means generalization for linear regression. ##

linreg_hs <- function(y, X, lambda = 0, estMethod = 2, Sigma = NULL){
  
  fname <- "linreg_hs"
  
  # We follow the code provided by Sivan Sabato on her website (filename median_regression.m) precisely.
  
  n <- dim(X)[1]
  d <- dim(X)[2]
  
  # Set the number of subsets to partition into, k.
  if (floor(n/2) <= d){
    cat("Error (", fname, "): not enough data to do a partition.", "\n", sep = "")
    return(NULL)
  } else {
    k <- max(2, floor(n/{2*d}))
  }
  
  # Shuffle the order of observations.
  
  idx <- sample(x = 1:n, size = n, replace = FALSE)
  y <- y[idx]
  X <- X[idx,]
  
  # Size of each segment of the partition (some samples may be unused)
  nsmall <- floor(x = {n/k})
  
  W <- matrix(0, nrow = k, ncol = d)
  
  for (i in 1:k){
    
    Xsmall <- X[{nsmall*{i-1}+1}:{nsmall*i}, ]
    ysmall <- y[{nsmall*{i-1}+1}:{nsmall*i}]
    
    z <- linreg_hs_myregression(y = ysmall, X = Xsmall, lambda = lambda)$w # note this fn returns a list.
    
    W[i,] <- z
    
  }
  
  if (estMethod == 1){
    Sigma <- t(X) %*% X # d x d
    
    covariances <- W %*% {Sigma + lambda*diag(x = 1, nrow = d, ncol = d)} %*% t(W) # k x k
    
    distances <- diag(covariances) %*% t(rep(1, k)) + t({diag(covariances) %*% t(rep(1, k))}) - 2*covariances # k x k
    
  } # exit estMethod 1.
  
  if (estMethod == 3){
    
    # only difference from case 1 is that we assume Sigma is provided as an argument.
    
    covariances <- W %*% {Sigma + lambda*diag(x = 1, nrow = d, ncol = d)} %*% t(W) # k x k
    
    distances <- diag(covariances) %*% t(rep(1, k)) + t({diag(covariances) %*% t(rep(1, k))}) - 2*covariances # k x k
    
  } # exit estMethod 3.
  
  if (estMethod == 2){
    
    # This is the "full" random procedure proposed in their work.
    
    numsplits <- k-1
    # the extra re-shuffling that they heuristically mention is not carried out here.
    nsmall <- floor(n/numsplits)
    distances <- matrix(0, nrow = k, ncol = k)
    
    covlist <- list()
    
    for (i in 1:numsplits){
      Xsmall <- X[{nsmall*{i-1}+1}:{nsmall*i}, ]
      covlist[[i]] <- t(Xsmall) %*% Xsmall # d x d
    }
    
    for (i in 1:k){
      
      for (j in 1:k){
        
        if (i==j){
          distances[i,j] <- 0
        } else {
          
          if (j < i){
            numcov <- j
          } else {
            numcov <- j-1
          }
          
          distances[i,j] <- drop({W[i,]-W[j,]} %*% {{covlist[[numcov]] + lambda*diag(1, nrow = d, ncol = d)} %*% {W[i,]-W[j,]}})
          
        }
        
      }
      
    }
    
  } # exit estMethod 2.
  
  # now we process the distances to get our final output.
  index <- linreg_hs_generalized_median(distances)
  w <- W[index, ]
  
  return(w)
}

linreg_hs_myregression <- function(y, X, lambda){
  
  n <- dim(X)[1]
  d <- dim(X)[2]
  
  est_cov <- {t(X) %*% X}/n + lambda*diag(1, nrow = d, ncol = d)
  est_cor <- drop(t(X) %*% y) / n
  
  w <- drop(ginv(X = est_cov) %*% est_cor) 
  emploss <- sum({y - drop(X %*% w)}^2)/n
  
  list(w = w, emploss = emploss)
}

linreg_hs_generalized_median <- function(distances){
  
  # distances is a k x t (for whatever t) matrix.
  # Take the median from each column, and then the minimum of each.
  # What this functions returns, though, is just the INDEX (which, of 1,2,...,t cols had the smallest median val)
  
  medvec <- apply(X = distances, MARGIN = 2, FUN = median)
  
  which(medvec == min(medvec))[1] # return the index.
}



## linreg using the Brownlees et al. (2015) method, with estimated variance bound. ##

linreg_bjl <- function(y, X, thres = 1e-3, altC.est = est_catwide, 
                       alt.scale = init.sd,
                       alt.loss = function(x){x^2},
                       alt.init = init.ols,
                       method = "Nelder-Mead", iters = 50L){
  
  # Compute these here and pass down to all further computations, so only compute once.
  n <- length(y)
  d <- ncol(X)
  
  # Initialize weights.
  fit <- alt.init(y = y, X = X)
  w_init <- fit$coef
  resid_init <- fit$resid
  
  # Generate a pre-fixed scale estimate to be used throughout, of the proper order.
  scale_est <- alt.scale(resid = resid_init, alt.loss = alt.loss)
  var_est <- scale_est^2
  alpha <- sqrt(2/{n*var_est}) # confidence-free setting; BJL p.4.
  
  # Optimize the weights using a pre-built nonlinear solver.
  optim_out <- optim(par = w_init, fn = linreg_bjl_obj, method = method,
                     y = y, X = X, n = n,
                     thres = thres, altC.est = altC.est, s = {1/alpha},
                     alt.loss = alt.loss, iters = iters)
  
  w_est <- optim_out$par # output from the pre-built.
  
  return(w_est)
}


linreg_bjl_obj <- function(w, y, X, n, thres, altC.est, s, alt.loss, iters){
  
  # The objective function, a fn of the weights, to be optimized.
  # It returns an estimate of the quasi-risk parameter, theta-star.
  
  resids <- y - drop(X %*% w)
  
  # Compute the new objective fn value, the solution to M-estimate of location, i.e., solves the psi-condition.
  altC.est(x = alt.loss(resids), n = n, s = s, thres = thres, iters = iters)
}












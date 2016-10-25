
# mlR_data.R
# Author: Matthew J. Holland
# Updated: 2016/10/12
# Description:
# Functions related to data generation.


#### Regarding use of seeds for reproducing data. ####

# In our tests, for each experimental condition, have a unique seed (specified in exp_df_*.rda).
# For the jth condition, we take myseed <- exp_df$seed[j].
# Then we run set.seed(seed = myseed) at the following points in the data generation process:
## ONCE prior to generating all trials (under condition j) worth of carrier matrices (X).
## ONCE prior to generating all trials (under condition j) worth of additive noise.
## ONCE for EACH NOISE LEVEL, prior to generating all trials worth of true weights.
# Thus, for a single (n,d) setting and single noise condition, if there are 15 noise levels,
# then set.seed(seed = myseed) will be executed precisely 1+1+15 times.


#### Regarding the experiment parameters (exp_df_*.rda files) ####

# Each row of exp_df corresponds to one setting of n, d, and the noise distribution family.
# n: total number of samples to be generated
# num_trials: number of trials per condition
# num_train: number of samples to be allocated for training model parameters.
# Xdist, Xcorr: always set to "norm" and "zero" respectively. Uncorrelated Normal carrier matrices.
# Xdim: the dimension of the inputs, i.e., the number of columns in the carrier matrix.
# noisedist: the name of the noise distribution used, works with data_noise function.
# sparsity: irrelevant here, always set to 1.
# seed: the seed used for the experiments captured by a given row, set.seed used as described above.


#### A collection of random-data generating functions, also other quantities of interest. ####


## Arcsine ##

ml_pasin <- function(x, shift = 0, scale = 1){
  
  # Probability (distribution) function for Arcsine distro.
  
  # Note that the support is (shift, shift+scale).
  
  2 * asin(sqrt({x- shift}/scale)) / pi
  
}

ml_dasin <- function(x, shift = 0, scale = 1){
  
  # Density function for Arcsine distro.
  
  # Note that the support is (shift, shift+scale).
  
  1 / {pi * sqrt({x - shift}*{shift + scale - x})}
  
}

ml_qasin <- function(p, shift = 0, scale = 1){
  
  # Quantile function for Arcsine distro.
  
  shift + scale * sin(pi*p/2)^2
  
}

ml_rasin <- function(n, shift = 0, scale = 1){
  
  # Sample n iid observations from Arcsine distro.
  
  ml_qasin(p = runif(n = n, min = 0, max = 1), shift = shift, scale = scale)
  
}

ml_masin <- function(shift, scale){
  
  # Expected value of Arcsine distro.
  
  shift + scale/2
  
}

ml_vasin <- function(shift, scale){
  
  # Variance of Arcsine distro.
  
  scale^2 / 8
  
}

ml_sasin <- function(shift, scale){
  
  # Skewness of Arcsine distro.
  
  0
  
}


## Beta Prime ##

ml_pbpri <- function(x, shape1, shape2){
  
  # Probability (distribution) function for Beta Prime distro.
  
  require(stats)
  
  pbeta(q = {x/{x+1}}, shape1 = shape1, shape2 = shape2)
  
}

ml_dbpri <- function(x, shape1, shape2){
  
  # Density function for Beta Prime distro.
  
  x^{shape1 - 1} / {beta(a = shape1, b = shape2) * {1+x}^{shape1 + shape2}}
  
}

ml_qbpri <- function(p, shape1, shape2){
  
  # Quantile function for Beta Prime distro.
  
  require(stats)
  
  qbeta(p = p, shape1 = shape1, shape2 = shape2) / {1 - qbeta(p = p, shape1 = shape1, shape2 = shape2)}
  
}

ml_rbpri <- function(n, shape1, shape2){
  
  # Sample n iid observations from Beta Prime distro.
  
  ml_qbpri(p = runif(n = n, min = 0, max = 1), shape1 = shape1, shape2 = shape2)
  
}

ml_mbpri <- function(shape1, shape2){
  
  # Expected value of Beta Prime distro.
  
  fname <- "ml_mbpri"
  
  idx <- which(shape2 <= 1)
  
  # Remember that mean is Inf when shape2 <= 1. 
  if (length(idx) > 0){
    cat(paste("Note (", fname, "): there are shape values less/eq 1. This means infinite mean.", "\n", sep = ""))
    
    out <- shape1/{shape2-1}
    out[idx] <- Inf
    return(out)
  } else {
    return(shape1/{shape2-1})
  }
  
}

ml_vbpri <- function(shape1, shape2){
  
  # Variance of Beta Prime distro.
  
  fname <- "ml_vbpri"
  
  # Remember that variance is Inf when shape2 <= 2 (more precisely, if shape2 <= 1 then it doesn't exist)
  
  idx2 <- which(shape2 <= 2)
  idx1 <- which(shape2 <= 1)
  
  if (length(idx2) > 0){
    cat(paste("Note (", fname, "): there are shape2 values <= 2. This means infinite variance.", "\n", sep = ""))
    
    out <- shape1 * {shape1 + shape2 - 1} / {{shape2 - 1}^2 * {shape2 - 2}}
    out[idx2] <- Inf
    if (length(idx1) > 0){
      out[idx1] <- NaN
    }
    return(out)
  } else {
    return({shape1 * {shape1 + shape2 - 1} / {{shape2 - 1}^2 * {shape2 - 2}}})
  }
  
}

ml_sbpri <- function(shape1, shape2){
  
  # Skewness of Beta Prime distro.
  
  fname <- "ml_sbpri"
  
  idx <- which(shape2 <= 3)
  
  # Remember that skewness is Inf when shape2 <= 3.
  if (length(idx) > 0){
    cat(paste("Note (", fname, "): there are shape2 values <= 3. This means infinite skewness.", "\n", sep = ""))
    
    out <- 2 * {2*shape1 + shape2 - 1} * sqrt({shape2 - 2}/{shape1*{shape1 + shape2 - 1}}) / {shape2 - 3}
    out[idx] <- Inf
    return(out)
  } else {
    return({2 * {2*shape1 + shape2 - 1} * sqrt({shape2 - 2}/{shape1*{shape1 + shape2 - 1}}) / {shape2 - 3}})
  }
  
}


## Chi-squared ##

ml_mchisq <- function(df){
  
  # Expected value of Chi-squared distro.
  
  df
  
}

ml_vchisq <- function(df){
  
  # Variance of Chi-squared distro.
  
  2*df
  
}

ml_schisq <- function(df){
  
  # Skewness of Chi-squared distro.
  
  2*sqrt(2/df)
  
}


## Exponential ##

ml_mexp <- function(r){
  
  # Expected value of exponential distro.
  
  1/r
  
}

ml_vexp <- function(r){
  
  # Variance of exponential distro.
  
  1/r^2
  
}

ml_sexp <- function(r){
  
  # Skewness of exponential distro.
  
  2 # it's constant, regardless of rate.
}


## Exponential-logarithmic, shape in (0,1), scale > 0 ##

ml_pexplog <- function(x, shape, scale = 1){
  
  # Probability (distribution) function of exponential-logarithmic distro.
  
  1 - log({1 - {1-shape}*exp(-x/scale)}) / log(shape)
  
}

ml_dexplog <- function(x, shape, scale = 1){
  
  # Density function of exponential-logarithmic distro.
  
  {shape - 1} * exp(-x/scale) / {scale*log(shape)*{1 - {1-shape}*exp(-x/scale)}}
  
}

ml_qexplog <- function(p, shape, scale = 1){
  
  # Quantile function (inverse distribution fn) of exponential-logarithmic distro.
  
  scale * {log({1-shape}) - log({1-shape^{1-p}})}
  
}

ml_rexplog <- function(n, shape, scale = 1){
  
  # Sampling independent random values from exponential-logarithmic distro.
  
  require(stats)
  
  ml_qexplog(p = runif(n, min = 0, max = 1), shape = shape, scale = scale)
}

ml_mexplog <- function(shape, scale = 1){
  
  # Expected value of exponential-logarithmic distro.
  
  require(pracma)
  
  -scale * polylog(z = {1-shape}, n = 2) / log(shape)
  
}

ml_vexplog <- function(shape, scale = 1){
  
  # Variance of exponential-logarithmic distro.
  
  require(pracma)
  
  scale^2 * {-2*polylog(z = {1-shape}, n = 3)/log(shape) - {polylog(z = {1-shape}, n = 2)/log(shape)}^2}
  
}


## F ##

ml_mf <- function(df1, df2){
  
  # Expected value for F distro.
  
  fname <- "ml_mf"
  
  # Remember that expectation is Inf when 0 < df2 <= 2
  
  idx <- which(df2 <= 2)
  
  if (length(idx) > 0){
    cat(paste("Note (", fname, "): there are df2 values <= 2. This means infinite mean.", "\n", sep = ""))
    
    out <- df2 / {df2 - 2}
    out[idx] <- Inf
    
    return(out)
  } else {
    return({df2 / {df2 - 2}})
  }
  
}


ml_vf <- function(df1, df2){
  
  # Variance for F distro.
  
  fname <- "ml_vf"
  
  # Remember that variance is Inf when 2 < df2 <= 4, and undefined when df2 <= 2.
  
  idx4 <- which(df2 <= 4)
  idx2 <- which(df2 <= 2)
  
  if (length(idx4) > 0){
    cat(paste("Note (", fname, "): there are df2 values <= 4. This means Inf or NaN variance.", "\n", sep = ""))
    
    out <- 2 * {df2/{df2 - 2}}^2 * {df1 + df2 - 2} / {df1*{df2-4}}
    out[idx4] <- Inf
    
    if (length(idx2) > 0){
      out[idx2] <- NaN
    }
    
    return(out)
  } else {
    return({2 * {df2/{df2 - 2}}^2 * {df1 + df2 - 2} / {df1*{df2-4}}})
  }
  
}


ml_sf <- function(df1, df2){
  
  # Skewness for F distro.
  
  fname <- "ml_sf"
  
  # Remember that skewness isn't defined for df2 <= 6.
  
  idx <- which(df2 <= 6)
  
  if (length(idx) > 0){
    cat(paste("Note (", fname, "): there are df2 values <= 6. This means NaN skewness.", "\n", sep = ""))
    
    out <- {2*df1 + df2 - 2} * sqrt(8*{df2-4}) / {{df2-6}*sqrt(df1*{df1+df2-2})}
    out[idx] <- NaN
    
    return(out)
  } else {
    return({{2*df1 + df2 - 2} * sqrt(8*{df2-4}) / {{df2-6}*sqrt(df1*{df1+df2-2})}})
  }
  
}


## Folded Normal ##

ml_pfnorm <- function(x, shift = 0, scale = 1){
  
  # Probability (distribution) function of Folded Normal distro.
  
  # Domain is [0, Inf)
  pnorm(q = {x-shift}/scale) - pnorm(q = {-x-shift}/scale)
  
}

ml_dfnorm <- function(x, shift = 0, scale = 1){
  
  # Density function of Folded Normal distro.
  
  # Domain is [0, Inf)
  {exp(-{{x+shift}/scale}^2/2) + exp(-{{x-shift}/scale}^2/2)} / {scale*sqrt(2*pi)}
  
}

ml_rfnorm <- function(n, shift = 0, scale = 1){
  
  # Sample n iid observations from the Folded Normal distro.
  
  abs(rnorm(n = n, mean = shift, sd = scale))
  
}

ml_mfnorm <- function(shift = 0, scale = 1){
  
  # Expected value of Folded Normal distro.
  
  shift*{1 - 2*pnorm(q = -{shift/scale})} + scale*sqrt(2/pi)*exp(-{shift/scale}^2/2)
  
}

ml_vfnorm <- function(shift = 0, scale = 1){
  
  # Variance of Folded Normal distro.
  
  shift^2 + scale^2 - ml_mfnorm(shift = shift, scale = scale)^2
  
}


## Frechet ##

ml_pfrec <- function(x, shift = 0, scale = 1, shape){
  
  # Probability (distribution) function of Frechet distro.
  
  exp(-{{x-shift}/scale}^{-shape})
  
}

ml_dfrec <- function(x, shift = 0, scale = 1, shape){
  
  # Density function of Frechet distro.
  
  shape * {{x-shift}/scale}^{-1-shape} * exp(-{{x-shift}/scale}^{-shape}) / scale
  
}

ml_qfrec <- function(p, shift = 0, scale = 1, shape){
  
  # Quantile function of Frechet distro.
  
  shift + scale / {-log(p)}^{1/shape}
  
}

ml_rfrec <- function(n, shift = 0, scale = 1, shape){
  
  # Sample n iid observations from the Frechet distro.
  
  ml_qfrec(p = runif(n = n, min = 0, max = 1), shift = shift, scale = scale, shape = shape)
  
}

ml_mfrec <- function(shift = 0, scale = 1, shape){
  
  # Expected value of Frechet distro.
  
  fname <- "ml_mfrec"
  
  idx <- which(shape <= 1)
  
  # Remember that mean is Inf when shape <= 1. 
  if (length(idx) > 0){
    cat(paste("Note (", fname, "): there are shape values less/eq 1. This means infinite mean.", "\n", sep = ""))
    
    out <- numeric(length(shape))
    
    out[idx] <- Inf
    if (length(idx) < length(shape)){
      out[-idx] <- shift[-idx] + scale[-idx] * gamma({1-1/shape[-idx]})
    }
    return(out)
  } else {
    return({shift + scale * gamma({1-1/shape})})
  }
  
}

ml_vfrec <- function(shift = 0, scale = 1, shape){
  
  # Variance of Frechet distro.
  
  fname <- "ml_vfrec"
  
  idx <- which(shape <= 2)
  
  # Remember that var is Inf when shape <= 2. 
  if (length(idx) > 0){
    cat(paste("Note (", fname, "): there are shape values less/eq 2. This means infinite var.", "\n", sep = ""))
    
    out <- numeric(length(shape))
    
    out[idx] <- Inf
    if (length(idx) < length(shape)){
      out[-idx] <- scale[-idx]^2 * {gamma({1-2/shape[-idx]}) - gamma({1-1/shape[-idx]})^2}
    }
    return(out)
  } else {
    return({scale^2 * {gamma({1-2/shape}) - gamma({1-1/shape})^2}})
  }
  
}

ml_sfrec <- function(shift = 0, scale = 1, shape){
  
  # Skewness of Frechet distro.
  
  fname <- "ml_sfrec"
  
  idx <- which(shape <= 3)
  
  # Remember that Skewness is Inf when shape <= 3. 
  if (length(idx) > 0){
    cat(paste("Note (", fname, "): there are shape values less/eq 3. This means infinite skewness.", "\n", sep = ""))
    
    out <- numeric(length(shape))
    
    out[idx] <- Inf
    if (length(idx) < length(shape)){
      out[-idx] <- {gamma({1-3/shape[-idx]}) - 3*gamma({1-2/shape[-idx]})*gamma({1-1/shape[-idx]}) + 2*gamma({1-1/shape[-idx]})^3} / sqrt({gamma({1-2/shape[-idx]}) - gamma({1-1/shape[-idx]})^2}^3)
    }
    return(out)
  } else {
    return({{gamma({1-3/shape}) - 3*gamma({1-2/shape})*gamma({1-1/shape}) + 2*gamma({1-1/shape})^3} / sqrt({gamma({1-2/shape}) - gamma({1-1/shape})^2}^3)})
  }
  
}


## Gamma ##

ml_mgamma <- function(scale, shape){
  
  # Mean of Gamma distribution.
  
  out <- scale * shape
  
  return(out)
}

ml_vgamma <- function(scale, shape){
  
  # Variance of Gamma distribution.
  
  out <- scale^2 * shape
  
  return(out)
}

ml_sgamma <- function(scale, shape){
  
  # Skewness of Gamma distribution (i.e., third moment of centered, StdDev-scaled RV)
  
  out <- 2 / sqrt(shape)
  
  return(out)
}


## Gaussian mixture model (k components, each dimension 1) ##

ml_hgmm_comp_idx <- function(x, j){
  
  fname <- "ml_hgmm_comp_idx" # the "h" is for "helper."
  
  # Return binary vector, turning on indices which match j.
  
  out <- numeric(length(x))
  on_idx <- which(x == j)
  out[on_idx] <- 1
  
  return(out)
}

ml_rgmm <- function(n, k, comp_weight, comp_mean, comp_sd){
  
  fname <- "ml_rgmm"
  
  # Randomly generate n values from a Gaussian mixture model with k components.
  
  component_check <- {length(comp_weight) == k} + {length(comp_mean) == k} + {length(comp_sd) == k}
  if (component_check != 3){
    cat("Error (", fname, "): the size of something isn't right.", "\n", sep = "")
    return(NULL)
  }
  
  comps <- sample(x = 1:k, replace = TRUE, size = n, prob = comp_weight)
  
  out <- numeric(n)
  
  for (i in 1:k){
    out <- out + ml_hgmm_comp_idx(x = comps, j = i) * rnorm(n, mean = comp_mean[i], sd = comp_sd[i])
  }
  
  return(out)
}


ml_mgmm <- function(k, comp_weight, comp_mean, comp_sd){
  
  fname <- "ml_mgmm"
  
  # Mean of one-dim Gaussian mixture model.
  
  component_check <- {length(comp_weight) == k} + {length(comp_mean) == k} + {length(comp_sd) == k}
  if (component_check != 3){
    cat("Error (", fname, "): not the size of something isn't right.", "\n", sep = "")
    return(NULL)
  }
  
  out <- as.numeric({comp_weight %*% comp_mean})
  
  return(out)
}


ml_vgmm <- function(k, comp_weight, comp_mean, comp_sd){
  
  fname <- "ml_vgmm"
  
  # Variance of one-dim Gaussian mixture model.
  
  component_check <- {length(comp_weight) == k} + {length(comp_mean) == k} + {length(comp_sd) == k}
  if (component_check != 3){
    cat("Error (", fname, "): not the size of something isn't right.", "\n", sep = "")
    return(NULL)
  }
  
  second_moments <- comp_sd^2 + comp_mean^2
  true_mean <- ml_mgmm(k = k, comp_weight = comp_weight, comp_mean = comp_mean, comp_sd = comp_sd)
  out <- as.numeric({comp_weight %*% second_moments}) - true_mean^2
  
  return(out)
}


ml_sgmm <- function(k, comp_weight, comp_mean, comp_sd){
  
  fname <- "ml_sgmm"
  
  # Skewness of one-dim Gaussian mixture model.
  
  component_check <- {length(comp_weight) == k} + {length(comp_mean) == k} + {length(comp_sd) == k}
  if (component_check != 3){
    cat("Error (", fname, "): not the size of something isn't right.", "\n", sep = "")
    return(NULL)
  }
  
  # for now, computations via REF file (in code/robest).
  k1 <- ml_mgmm(k = k, comp_weight = comp_weight, comp_mean = comp_mean, comp_sd = comp_sd)
  k2 <- ml_vgmm(k = k, comp_weight = comp_weight, comp_mean = comp_mean, comp_sd = comp_sd)
  
  mean2 <- comp_mean^2
  mean3 <- comp_mean^3
  sd2 <- comp_sd^2
  mean_sd2 <- comp_mean * sd2
  
  k3 <- as.numeric({comp_weight %*% {mean3 + 3*comp_mean*sd2}}) + 2*k1^3 - 3 * k1 * as.numeric({comp_weight %*% {mean2 + sd2}})
  
  out <- k3 / {k2^{3/2}}
  
  return(out)
}


ml_dgmm <- function(x, k, comp_weight, comp_mean, comp_sd){
  
  fname <- "ml_dgmm"
  
  component_check <- {length(comp_weight) == k} + {length(comp_mean) == k} + {length(comp_sd) == k}
  if (component_check != 3){
    cat("Error (", fname, "): not the size of something isn't right.", "\n", sep = "")
    return(NULL)
  }
  
  out <- numeric(length(x))
  
  for (i in 1:k){
    out <- out + comp_weight[i] * dnorm(x = x, mean = comp_mean[i], sd = comp_sd[i])
  }
  
  return(out)
}


## Gompertz ##

ml_pgomp <- function(x, shape, scale = 1){
  
  # Probability (distribution) function for Gompertz distro.
  
  1 - exp(-shape*{exp(x/scale) - 1})
  
}

ml_dgomp <- function(x, shape, scale = 1){
  
  # Density function for Gompertz distro.
  
  shape * exp(x/scale) * exp(-shape*{exp(x/scale) - 1}) / scale
  
}

ml_qgomp <- function(p, shape, scale = 1){
  
  # Quantile function for Gompertz distro.
  
  scale * log({1 - log({1-p})/shape})
  
}

ml_rgomp <- function(n, shape, scale = 1){
  
  # Sample n iid observations from Gompertz distro.
  
  require(stats)
  
  ml_qgomp(p = runif(n = n, min = 0, max = 1), shape = shape, scale = scale)
}

ml_mgomp <- function(shape, scale = 1){
  
  # Expected value of the Gompertz distro.
  
  require(pracma)
  
  scale * exp(shape) * expint_E1(shape)
  
}

ml_vgomp <- function(shape, scale = 1){
  
  # Variance of the Gompertz distro.
  
  require(hypergeo)
  
  # note that the hypergeometric radius of convergence on the complex plane is 1, so shape
  # should be less than 1 for things to converge.
  
  scale^2 * exp(shape) * {-2*shape*genhypergeo(U = c(1,1,1), L = c(2,2,2), z = {-shape}) + digamma(1)^2 + 
      pi^2/6 + 2*{-digamma(1)}*log(shape) + log(shape)^2} - ml_mgomp(shape = shape, scale = scale)^2
  
}


## Gumbel (type-1 EV distro for maximums) ##

ml_pgum <- function(x, shift = 0, scale = 1){
  
  # Probability (distribution) function for Gumbel (type-1 EV, max) distro.
  
  exp(-exp(-{x-shift}/scale))
  
}

ml_dgum <- function(x, shift = 0, scale = 1){
  
  # Density function for Gumbel (type-1 EV, max) distro.
  
  exp(-{x-shift}/scale) * exp(-exp(-{x-shift}/scale)) / scale
  
}

ml_qgum <- function(p, shift = 0, scale = 1){
  
  # Quantile function for Gumbel (type-1 EV, max) distro.
  
  shift - scale * log(-log(p))
  
}

ml_rgum <- function(n, shift = 0, scale = 1){
  
  # Sample n iid values from Gumbel (type-1 EV, max) distro.
  
  require(stats)
  
  ml_qgum(p = runif(n = n, min = 0, max = 1), shift = shift, scale = scale)
  
}

ml_mgum <- function(shift, scale){
  
  # Expected value of Gumbel (type-1 EV, max) distro.
  
  shift + scale * {-digamma(1)}
  
}

ml_vgum <- function(shift, scale){
  
  # Variance of Gumbel (type-1 EV, max) distro.
  
  scale^2 * pi / 6
  
}

ml_sgum <- function(shift, scale){
  
  # Skewness of Gumbel (type-1 EV, max) distro.
  
  require(pracma)
  # note zeta() evaluated at 3 is called "Apery's constant"
  
  12 * sqrt(6) * zeta(3) / pi^3 # constant, regardless of parameters.
  
}


## Hyperbolic Secant ##

ml_phsec <- function(x, shift, scale = 1){
  
  # Probability (distribution) function of hyperbolic secant distro.
  
  2 * atan(exp(pi*{x-shift}/{2*scale})) / pi
  
}

ml_dhsec <- function(x, shift, scale = 1){
  
  # Density function of hyperbolic secant distro.
  
  require(pracma)
  
  sech(pi*{x-shift}/{2*scale}) / {2*scale}
  
}

ml_qhsec <- function(p, shift, scale = 1){
  
  # Quantile function of hyperbolic secant distro.
  
  shift + scale * 2 * log(tan(p*pi/2)) / pi
  
}

ml_rhsec <- function(n, shift, scale = 1){
  
  # Sample n iid observations from hyperbolic secant distro.
  
  require(stats)
  
  ml_qhsec(p = runif(n = n, min = 0, max = 1), shift = shift, scale = scale)
  
}

ml_mhsec <- function(shift, scale){
  
  # Expected value of hyperbolic secant distro.
  
  shift
  
}

ml_vhsec <- function(shift, scale){
  
  # Variance of hyperbolic secant distro.
  
  scale^2
  
}

ml_shsec <- function(shift, scale){
  
  # Skewness of hyperbolic secant distro.
  
  0
  
}


## Laplace ##

ml_plap <- function(x, shift = 0, scale = 1){
  
  # Probability (distribution) function of Lapace distro.
  
  as.numeric({x <= shift}) * exp({x-shift}/scale)/2 + as.numeric({x > shift})*{1 - exp(-{x-shift}/scale)/2}
  
}

ml_dlap <- function(x, shift = 0, scale = 1){
  
  # Density function of Lapace distro.
  
  exp(-abs({x-shift})/scale) / {2*scale}
  
}

ml_qlap <- function(p, shift = 0, scale = 1){
  
  # Quantile function of Lapace distro.
  
  as.numeric({p <= 1/2}) * {shift + scale*log(2*p)} + as.numeric({p > 1/2})*{shift - scale*log(2*{1-p})}
  
}

ml_rlap <- function(n, shift = 0, scale = 1){
  
  # Sample n iid observations from the Laplace distro.
  
  ml_qlap(p = runif(n = n, min = 0, max = 1), shift = shift, scale = scale)
  
}

ml_mlap <- function(shift, scale){
  
  # Expected value of Laplace distro.
  
  shift
  
}

ml_vlap <- function(shift, scale){
  
  # Variance of Laplace distro.
  
  2*scale^2
  
}

ml_slap <- function(shift, scale){
  
  # Skewness of Laplace distro.
  
  0
  
}


## Levy ##

ml_plevy <- function(x, shift = 0, scale = 1){
  
  # Probability (distribution) function for Levy distro.
  
  require(stats)
  
  # Note that {x > shift} is the support.
  2 * {1 - pnorm(sqrt(scale/{x-shift}))}
  
}

ml_dlevy <- function(x, shift = 0, scale = 1){
  
  # Density function for Levy distro.
  
  # Note that {x > shift} is the support.
  sqrt(scale/{2*pi}) * exp(-scale/{2*{x - shift}}) / {x - shift}^{3/2}
  
}

ml_qlevy <- function(p, shift = 0, scale = 1){
  
  # Quantile function for Levy distro.
  
  require(stats)
  
  shift + scale / qnorm({1-p/2})^2
  
}

ml_rlevy <- function(n, shift = 0, scale = 1){
  
  # Sample n iid observations from Levy distro.
  
  require(stats)
  
  ml_qlevy(p = runif(n = n, min = 0, max = 1), shift = shift, scale = scale)
  
}


## log-Logistic ##

ml_pllog <- function(x, shape, scale = 1){
  
  # Probability (distribution) function for log-Logistic distro.
  
  x^shape / {scale^shape + x^shape}
  
}

ml_dllog <- function(x, shape, scale = 1){
  
  # Density function for log-Logistic distro.
  
  scale^shape * shape * x^{shape - 1} / {scale^shape + x^shape}^2
  
}

ml_qllog <- function(p, shape, scale = 1){
  
  # Quantile function for log-Logistic distro.
  
  scale * {p / {1-p}}^{1/shape}
  
}

ml_rllog <- function(n, shape, scale = 1){
  
  # Sample n iid observations from log-Logistic distro.
  
  ml_qllog(p = runif(n = n, min = 0, max = 1), shape = shape, scale = scale)
  
}

ml_mllog <- function(shape, scale){
  
  fname <- "ml_mllog"
  
  # Expected value of the log-Logistic distro.
  
  idx <- which(shape <= 1)
  
  # Remember that mean is Inf when shape <= 1. 
  if (length(idx) > 0){
    cat(paste("Note (", fname, "): there are shape values less/eq 1. This means infinite mean.", "\n", sep = ""))
    
    out <- scale / ml_sinc(1/shape)
    out[idx] <- Inf
    return(out)
  } else {
    return({scale / ml_sinc(1/shape)})
  }
  
}

ml_vllog <- function(shape, scale){
  
  # Variance of the log-Logistic distro.
  
  idx2 <- which(shape <= 2)
  idx1 <- which(shape <= 1)
  
  # Remember that variance is Inf when shape <= 2 (more precisely, if shape <= 1 then it doesn't exist)
  if (length(idx2) > 0){
    cat(paste("Note (", fname, "): there are shape values less/eq 2. This means infinite variance.", "\n", sep = ""))
    
    out <- scale^2 * {1/ml_sinc(2/shape) - 1/ml_sinc(1/shape)^2}
    out[idx2] <- Inf
    if (length(idx1) > 0){
      out[idx1] <- NaN
    }
    return(out)
  } else {
    return({scale^2 * {1/ml_sinc(2/shape) - 1/ml_sinc(1/shape)^2}})
  }
  
}


## log-Normal ##

ml_vlnorm <- function(shift, scale){
  
  fname <- "ml_vlnorm"
  
  # Variance of log-Normal distribution.
  
  out <- {exp(scale^2) - 1} * exp(x = {2*shift + scale^2})
  
  return(out)
}


ml_mlnorm <- function(shift, scale){
  
  fname <- "ml_mlnorm"
  
  # Mean of log-Normal distribution.
  
  out <- exp(x = {shift + scale^2 / 2})
  
  return(out)
}


ml_slnorm <- function(shift, scale){
  
  fname <- "ml_slnorm"
  
  # Skewness of log-Normal distribution (i.e., third moment of centered, StdDev-scaled RV)
  
  expval <- exp(scale^2)
  
  out <- {expval + 2} * sqrt({expval - 1})
  
  return(out)
}


## Logistic ##

ml_plgst <- function(x, shift = 0, scale = 1){
  
  # Probability (distribution) function for Logistic distro.
  
  exp({x - shift}/scale) / {1+exp({x - shift}/scale)}
  
}

ml_dlgst <- function(x, shift = 0, scale = 1){
  
  # Density function for Logistic distro.
  
  exp({x - shift}/scale) / {{1+exp({x - shift}/scale)}^2 * scale}
  
}

ml_qlgst <- function(p, shift = 0, scale = 1){
  
  # Quantile function for Logistic distro.
  
  shift + scale * log(p/{1-p})
  
}

ml_rlgst <- function(n, shift = 0, scale = 1){
  
  # Sample n iid observations from Logistic distro.
  
  require(stats)
  
  ml_qlgst(p = runif(n = n, min = 0, max = 1), shift = shift, scale = scale)
  
}

ml_mlgst <- function(shift, scale){
  
  # Expected value of Logistic distro.
  
  shift
  
}

ml_vlgst <- function(shift, scale){
  
  # Variance of Logistic distro.
  
  scale^2 * pi^2 / 3
  
}

ml_slgst <- function(shift, scale){
  
  # Skewness of Logistic distro.
  
  0
  
}


## Maxwell ##

ml_pmaxw <- function(x, scale = 1){
  
  # Probability (distribution) function for Maxwell distro.
  
  require(stats)
  
  2 * pnorm(x/scale) - sqrt(2/pi) * x * exp(-{x/scale}^2/2) / scale - 1
  
}

ml_dmaxw <- function(x, scale = 1){
  
  # Density function for Maxwell distro.
  
  sqrt(2/pi) * x^2 * exp(-{x/scale}^2/2) / scale^3
  
}

ml_rmaxw <- function(n, scale = 1){
  
  # Sample n iid observations from Maxwell distro.
  
  require(stats)
  
  sqrt({rnorm(n = n, mean = 0, sd = scale)^2 + rnorm(n = n, mean = 0, sd = scale)^2 + rnorm(n = n, mean = 0, sd = scale)^2})
  
}

ml_mmaxw <- function(scale){
  
  # Expected value of Maxwell distro.
  
  scale * 2 * sqrt(2/pi)
  
}

ml_vmaxw <- function(scale){
  
  # Variance of Maxwell distro.
  
  scale^2 * {3 - 8/pi}
  
}

ml_smaxw <- function(scale){
  
  # Skewness of Maxwell distro.
  
  2 * sqrt(2) * {16 - 5*pi} / {3*pi - 8}^{3/2}
  
}


## Normal (multivariate) ##

ml_rmvnorm <- function(n, mu, Sigma){
  
  fname <- "ml_rmvnorm"
  
  require(MASS)
  
  out <- mvrnorm(n = n, mu = mu, Sigma = Sigma) # returns a dim-col matrix with n rows.
}


## Pareto ##

ml_dpareto <- function(x, a, b = 1){
  
  fname <- "ml_dpareto"
  
  # Pareto density function
  
  out <- a * {b/x}^a / x
  
  return(out)
}

ml_ppareto <- function(x, a, b = 1){
  
  fname <- "ml_ppareto"
  
  # Pareto prob distro fn
  
  out <- 1 - {b/x}^a
  
  return(out)
}

ml_qpareto <- function(p, a, b = 1){
  
  fname <- "ml_qpareto"
  
  # Pareto quantile function
  
  out <- b / {1-p}^{1/a}
  
  return(out)
}

ml_rpareto <- function(n, a, b = 1){
  
  fname <- "ml_rpareto"
  
  # Pareto random variable sampler
  
  require(stats)
  
  p <- runif(n)
  
  out <- b / p^{1/a} # note: not using qpareto since 1-p and p are the same if unif-sampled. Saves an operation.
  
  return(out)
}

ml_mpareto <- function(a, b = 1){
  
  fname <- "ml_mpareto"
  
  # Expectation of a Pareto rv
  bad_idx <- which(a <= 1)
  
  out <- b * a / {a-1}
  
  if (length(bad_idx) > 0){
    cat("Warning (", fname, "): parameter a is too small, some outputs undefined.", "\n", sep="")
    out[bad_idx] <- NaN
  }
  
  return(out)
}

ml_vpareto <- function(a, b = 1){
  
  fname <- "ml_vpareto"
  
  # Variance of a Pareto rv
  bad_idx <- which(a <= 2)
  
  out <- b^2 * a / {{a-1}^2 * {a-2}}
  
  if (length(bad_idx) > 0){
    cat("Warning (", fname, "): parameter a is too small, some outputs undefined.", "\n", sep="")
    out[bad_idx] <- NaN
  }
  
  return(out)
}

ml_spareto <- function(a, b = 1){
  
  fname <- "ml_spareto"
  
  # Skewness of a Pareto rv  (i.e., third moment of centered, StdDev-scaled RV)
  
  bad_idx <- which(a <= 3)
  
  out <- 2*{1+a} * sqrt({{a-2}/a}) / {a-3}
  
  if (length(bad_idx) > 0){
    cat("Warning (", fname, "): parameter a is too small, some outputs undefined.", "\n", sep="")
    out[bad_idx] <- NaN
  }
  
  return(out)
}


## Rayleigh ##

ml_prayl <- function(x, scale = 1){
  
  # Probability (distribution) function for Rayleigh distro.
  
  1 - exp(-{x/scale}^2/2)
  
}

ml_drayl <- function(x, scale = 1){
  
  # Density function for Rayleigh distro.
  
  x * exp(-{x/scale}^2/2) / scale^2
  
}

ml_qrayl <- function(p, scale = 1){
  
  # Quantile function for Rayleigh distro.
  
  scale * sqrt(-2*log({1-p}))
  
}

ml_rrayl <- function(n, scale = 1){
  
  # Sample n iid observations from Rayleigh distro.
  
  ml_qrayl(p = runif(n = n, min = 0, max = 1), scale = scale)
  
}

ml_mrayl <- function(scale){
  
  # Expected value of Rayleigh distro.
  
  scale * sqrt(pi/2)
  
}

ml_vrayl <- function(scale){
  
  # Variance of Rayleigh distro.
  
  scale^2 * {2-pi/2}
  
}

ml_srayl <- function(scale){
  
  # Skewness of Rayleigh distro.
  
  2 * sqrt(pi) * {pi - 3} / {4-pi}^{3/2}
  
}


## Semicircle ##

ml_pscir <- function(x, center = 0, rad = 1){
  
  # Probability (distribution) function for Semicircle distro.
  
  # The support is [a-r, a+r], thus the term "radius" for the scale parameter.
  
  1/2 + {x - center} * sqrt({rad^2 - {x - center}^2}) / {pi * rad^2} + asin({x-center}/rad)/pi
  
}

ml_dscir <- function(x, center = 0, rad = 1){
  
  # Density function for Semicircle distro.
  
  # The support is [a-r, a+r], thus the term "radius" for the scale parameter.
  
  2 * sqrt({rad^2 - {x - center}^2}) / {pi * rad^2}
  
}

ml_rscir <- function(n, center = 0, rad = 1){
  
  # Sample n iid observations from the Semicircle distro.
  
  # No quantile fn we can use, so a few tricks are required.
  
  # First, we want to sample "uniformly" from the unit circle.
  # Note it may be readily verified that (r*cos(theta),r*sin(theta)) as
  # defined below has constant (= uniform) pdf of 1/pi on the
  # unit circle.
  
  tmp_Mtx <- matrix(nrow = n, ncol = 2)
  tmp_Mtx[, 1] <- runif(n, min = 0, max = 1)
  tmp_Mtx[, 2] <- runif(n, min = 0, max = 1)
  
  r <- apply(X = tmp_Mtx, MARGIN = 1, FUN = max)
  
  theta <- 2*pi*runif(n = n, min = 0, max = 1)
  
  tmp_Mtx[, 1] <- r*cos(theta)
  tmp_Mtx[, 2] <- r*sin(theta)
  
  # So each row of this mtx is a 2-dim value uniformly (in PDF sense) distributed
  # on the UNIT circle. Basic calculus shows that both coordinates are respectively
  # following the standard semicircle distribution. Either one is fine, but we take the
  # first coordinate for use below.
  
  center + rad * tmp_Mtx[, 1] # shift and scale.
}

ml_mscir <- function(center, rad){
  
  # Expected value of Semicircle distro.
  
  center
  
}

ml_vscir <- function(center, rad){
  
  # Variance of Semicircle distro.
  
  rad^2 / 4
  
}

ml_sscir <- function(center, rad){
  
  # Skewness of Semicircle distro.
  
  0
  
}


## t (univariate Student-t) ##

ml_mt <- function(df){
  
  fname <- "ml_mt"
  
  # Expected value of Student-t distribution.
  
  0
}

ml_vt <- function(df){
  
  fname <- "ml_vt"
  
  # Variance of Student-t distribution.
  
  out <- df / {df - 2}
  
  return(out)
}

ml_st <- function(df){
  
  fname <- "ml_st"
  
  # Skewness of Student-t distribution.
  
  0
}


## t (multivariate Student-t) ##

ml_rmvt <- function(n, mu, Scale, df){
  
  fname <- "ml_rmvt"
  
  require(mnormt) 
  
  out <- rmt(n = n, df = df, mean = mu, S = Scale)
}


## Triangle ##

ml_ptri <- function(x, vert, shift = 0, scale = 1){
  
  # Probability (distribution) function for Triangle distro.
  
  # Note that there are a number of conditions, due to the "kink" in pdf.
  # Also, note support is [shift, shift+scale]
  
  if (vert == 0){
    return({1 - {{shift+scale-x}/scale}^2})
  }
  
  if (vert == 1){
    
    return({{x - shift}/scale}^2)
    
  }
  
  if (vert < 1 && vert > 0){
    
    idx <- which(x <= {shift+vert*scale})
    
    out <- numeric(length(x))
    
    if (length(idx) > 0){
      
      out[idx] <- {x[idx]-shift}^2 / {vert*scale^2}
      
      if (length(idx) < length(x)){
        
        out[-idx] <- 1 - {shift + scale - x[-idx]}^2 / {{1-vert}*scale^2}
        
      }
      
      return(out)
      
    } else {
      
      return({1 - {shift + scale - x}^2 / {{1-vert}*scale^2}})
      
    }
    
  }
  
}

ml_dtri <- function(x, vert, shift = 0, scale = 1){
  
  # Density function for Triangle distro.
  
  # Note that there are a number of conditions, due to the "kink".
  # Also, note support is [shift, shift+scale]
  
  if (vert == 0){
    return({2*{shift+scale-x}/scale^2})
  }
  
  if (vert == 1){
    
    return({2*{x-shift}/scale^2})
    
  }
  
  if (vert < 1 && vert > 0){
    
    idx <- which(x <= {shift+vert*scale})
    
    out <- numeric(length(x))
    
    if (length(idx) > 0){
      
      out[idx] <- 2*{x[idx]-shift} / {vert*scale^2}
      
      if (length(idx) < length(x)){
        
        out[-idx] <- 2*{shift + scale - x[-idx]} / {{1-vert}*scale^2}
        
      }
      
      return(out)
      
    } else {
      
      return({2*{shift + scale - x[-idx]} / {{1-vert}*scale^2}})
      
    }
    
  }
  
}

ml_qtri <- function(p, vert, shift = 0, scale = 1){
  
  # Quantile function for Triangle distro.
  
  # Note the conditions are embedded cleverly here.
  
  shift + scale*sqrt(p*vert)*as.numeric({p <= vert}) + scale*{1 - sqrt({1-p}*{1-vert})}*as.numeric({p > vert})
  
}

ml_rtri <- function(n, vert, shift = 0, scale = 1){
  
  # Sample n iid observations from the Triangle distro.
  
  ml_qtri(p = runif(n = n, min = 0, max = 1), vert = vert, shift = shift, scale = scale)
  
}

ml_mtri <- function(vert, shift = 0, scale = 1){
  
  # Expected value of the Triangle distro.
  
  shift + scale * {1 + vert} / 3
  
}

ml_vtri <- function(vert, shift = 0, scale = 1){
  
  # Variance of the Triangle distro.
  
  scale^2 * {1 - vert*{1-vert}} / 18
  
}

ml_stri <- function(vert, shift = 0, scale = 1){
  
  # Skewness of the Triangle distro.
  
  sqrt(2) * {1 - 2*vert} * {1+vert} * {2-vert} / {5 * {1 - vert*{1-vert}}^{3/2}}
  
}


## Uniform (on interval) ##

ml_punif <- function(x, a = 0, b = 1){
  
  # Probability function of interval Uniform distribution.
  
  punif(q = x, min = a, max = b)
  
}

ml_dunif <- function(x, a = 0, b = 1){
  
  # Density function of interval Uniform distribution.
  
  dunif(x = x, min = a, max = b)
  
}

ml_qunif <- function(p, a = 0, b = 1){
  
  # Quantile function of interval Uniform distribution.
  
  qunif(p = p, min = a, max = b)
  
}

ml_munif <- function(a = 0, b = 1){
  
  # Expected value of interval Uniform distribution.
  
  {a+b}/2
  
}

ml_vunif <- function(a = 0, b = 1){
  
  # Variance of interval Uniform distribution.
  
  {b-a}^2 / 12
  
}

ml_sunif <- function(a = 0, b = 1){
  
  # Skewness of interval Uniform distribution.
  
  {b-a} / {2*sqrt(3)}
  
}


## U-Power ##

ml_pupwr <- function(x, shape, shift = 0, scale = 1){
  
  # Probability (distribution) function of U-Power distro.
  
  if (shape != ceiling(shape)){
    stop("We can't set SHAPE parameter to be non-integer.")
  }
  
  # Support is [shift-scale, shift+scale]. Shape is integer.
  
  {1 + {{x - shift}/scale}^{2*shape+1}} / 2
  
}

ml_dupwr <- function(x, shape, shift = 0, scale = 1){
  
  # Density function of U-Power distro.
  
  if (shape != ceiling(shape)){
    stop("We can't set SHAPE parameter to be non-integer.")
  }
  
  # Support is [shift-scale, shift+scale]. Shape is integer.
  
  {2*shape + 1} * {{x-shift}/scale}^{2*shape} / {2*scale}
  
}

ml_qupwr <- function(p, shape, shift = 0, scale = 1){
  
  # Quantile function of U-Power distro.
  
  if (shape != ceiling(shape)){
    stop("We can't set SHAPE parameter to be non-integer.")
  }
  
  shift + scale * sign({2*p-1}) * abs({2*p-1})^{1/{2*shape+1}}
  
}

ml_rupwr <- function(n, shape, shift = 0, scale = 1){
  
  # Sample n iid observations from U-Power distro.
  
  if (shape != ceiling(shape)){
    stop("We can't set SHAPE parameter to be non-integer.")
  }
  
  ml_qupwr(p = runif(n = n, min = 0, max = 1), shape = shape, shift = shift, scale = scale)
  
}

ml_mupwr <- function(shape, shift = 0, scale = 1){
  
  # Expected value of U-Power distro.
  
  if(sum({shape != ceiling(shape)}) > 0){
    stop("We can't set SHAPE parameter to be non-integer.")
  }
  
  shift
  
}

ml_vupwr <- function(shape, shift = 0, scale = 1){
  
  # Variance of U-Power distro.
  
  if(sum({shape != ceiling(shape)}) > 0){
    stop("We can't set SHAPE parameter to be non-integer.")
  }
  
  scale^2 * {2*shape + 1} / {2*shape + 3}
  
}

ml_supwr <- function(shape, shift = 0, scale = 1){
  
  # Skewness of U-Power distro.
  
  if(sum({shape != ceiling(shape)}) > 0){
    stop("We can't set SHAPE parameter to be non-integer.")
  }
  
  0
  
}


## Wald (aka "inverse Gaussian") ##

ml_pwald <- function(x, mean = 1, shape = 1){
  
  # Probability (distribution) function of Wald distro.
  
  require(stats)
  
  pnorm(sqrt(shape/x)*{x/mean-1}) + exp(2*shape/mean)*pnorm(-sqrt(shape/x)*{x/mean+1})
  
}

ml_dwald <- function(x, mean, shape = 1){
  
  # Density function of Wald distro.
  
  sqrt(shape/{2*pi*x^3}) * exp(-shape*{x-mean}^2/{2*mean^2*x})
  
}

ml_mwald <- function(mean, shape = 1){
  
  # Expected value of Wald distro.
  
  mean # constant regardless of shape.
  
}

ml_vwald <- function(mean, shape = 1){
  
  # Variance of Wald distro.
  
  mean^3 / shape
  
}

ml_swald <- function(mean, shape = 1){
  
  # Skewness of Wald distro.
  
  3 * sqrt(mean/shape)
  
}

ml_rwald <- function(n, mean, shape = 1){
  
  # Randomly generate n values from Wald distro.
  # We follow technique of Michael et al. (1976).
  
  require(stats)
  
  v <- rnorm(n = n)^2 # squared standard Normals are chi-square with df=1.
  
  # the smaller root
  x1 <- mean + mean^2*v/{2*shape} - mean*sqrt({4*mean*shape*v + {mean*v}^2})/{2*shape}
  
  # the larger root
  x2 <- mean^2 / x1
  
  thres <- mean / {mean + x1} # the probability assigned to choosing root x1
  
  rv <- runif(n = n)
  idx <- as.numeric({rv <= thres}) # the obs which will take the smaller root
  
  idx*x1 + {1-idx}*x2
}


## Weibull ##

ml_vweibull <- function(scale, shape){
  
  fname <- "ml_vweibull"
  
  # Variance of Weibull distribution.
  
  out <- scale^2 * {gamma(x = {1+2/shape}) - {gamma(x = {1+1/shape})}^2}
  
  return(out)
}


ml_mweibull <- function(scale, shape){
  
  fname <- "ml_mweibull"
  
  # Mean of Weibull distribution.
  
  out <- scale * gamma(x = {1+1/shape})
  
  return(out)
}


ml_sweibull <- function(scale, shape){
  
  fname <- "ml_sweibull"
  
  # Skewness of Weibull distribution (i.e., third moment of centered, StdDev-scaled RV)
  
  meanval <- ml_mweibull(scale = scale, shape = shape)
  sd <- sqrt(ml_vweibull(scale = scale, shape = shape))
  
  out <- {gamma({1+3/shape}) * scale^3 - 3*meanval*sd^2 - meanval^3} / {sd^3}
  
  return(out)
}


#### A routine for generating noise data easily using same parameters as in the paper. ####

data_noise <- function(n, num_trials = 1, dist, seed = NULL){
  
  fname <- "data_noise"
  
  set.seed(seed = seed) # For reproducing results, set seed.
  
  require(mnormt)
  require(MASS)
  require(stats)
  
  # This function generates the additive, zero-mean noise from the specified distribution.
  # dist in {"arcsine", "betaprime", "chisq", "exp", "explog", "f", "frechet", "gamma", "gmm", "gompertz",
  #          "gumbel", "hsec", "laplace", "llog", "lnorm", "logistic", "maxwell", "norm", "pareto",
  #          "rayleigh", "semicircle", "t", "tri_a", "tri_s", "upower", "wald", "weibull"}
  
  # Arcsine distro, control scale para.
  if (dist == "arcsine"){
    
    paras <- seq(0.3*sqrt(8), 20*sqrt(8), {20*sqrt(8)-0.3*sqrt(8)}/14)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_masin(shift = 0, scale = paras[j])
      noise_true_var <- ml_vasin(shift = 0, scale = paras[j])
      noise_true_skew <- ml_sasin(shift = 0, scale = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rasin(n = n, shift = 0, scale = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Beta-Prime distro, shape2 para.
  if (dist == "betaprime"){
    
    paras <- numeric(15)
    tmp <- seq(0.3, 20, {20-0.3}/14)
    s1 <- 1.5
    for(i in 1:length(paras)){
      paras[i] <- Re(polyroot(z = c({-2-{s1/tmp[i]}^2+s1/{tmp[i]^2}}, {5-{s1/{tmp[i]^2}}}, {-4}, {1}))[3]) # third is real one in cubic case.
    }
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mbpri(shape1 = 1.5, shape2 = paras[j])
      noise_true_var <- ml_vbpri(shape1 = 1.5, shape2 = paras[j])
      noise_true_skew <- ml_sbpri(shape1 = 1.5, shape2 = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rbpri(n = n, shape1 = 1.5, shape2 = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Chi-squared distro, df para.
  if (dist == "chisq"){
    
    paras <- seq(0.3^2/2, 20^2/2, {20^2/2-0.3^2/2}/14)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mchisq(df = paras[j])
      noise_true_var <- ml_vchisq(df = paras[j])
      noise_true_skew <- ml_schisq(df = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        rchisq(n = n, df = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Exponential distro, rate para.
  if (dist == "exp"){
    
    paras <- 1/sqrt(seq(0.3^2, 20^2, {20^2 - 0.3^2}/14))
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mexp(r = paras[j])
      noise_true_var <- ml_vexp(r = paras[j])
      noise_true_skew <- ml_sexp(r = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        rexp(n = n, rate = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Exponential-logarithmic distro, rate para.
  if (dist == "explog"){
    
    paras <- seq(0.31, 20.01, {20.01-0.31}/14)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mexplog(shape = 0.95, scale = paras[j])
      noise_true_var <- ml_vexplog(shape = 0.95, scale = paras[j])
      noise_true_skew <- NaN
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rexplog(n = n, shape = 0.95, scale = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # F distro, df2 para.
  if (dist == "f"){
    
    paras <- numeric(15)
    tmp <- seq(0.4, 20, {20-0.4}/14)
    d1 <- 15
    for(i in 1:length(paras)){
      paras[i] <- Re(polyroot(z = c({-16*tmp[i]^2*d1}, {tmp[i]^2*d1*20}, {d1*{tmp[i]^2*{-4}*2-2}+4}, {tmp[i]^2*d1-2}))[3]) # third is real one in cubic case.
    }
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mf(df1 = 15, df2 = paras[j])
      noise_true_var <- ml_vf(df1 = 15, df2 = paras[j])
      noise_true_skew <- NaN
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        rf(n = n, df1 = 15, df2 = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Frechet distro, shape para.
  if (dist == "frechet"){
    
    paras <- c(5.5, 4, 3, 2.5, 2.3, 2.25, 2.2, 2.15, 2.1, 2.05, 2.025, 2.0125, 2.01, 2.0075, 2.005)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mfrec(shift = 0, scale = 1, shape = paras[j])
      noise_true_var <- ml_vfrec(shift = 0, scale = 1, shape = paras[j])
      noise_true_skew <- ml_sfrec(shift = 0, scale = 1, shape = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rfrec(n = n, shift = 0, scale = 1, shape = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Gamma distro, shape para.
  if (dist == "gamma"){
    
    paras <- {seq(0.3, 20, {20-0.3}/14) / 5}^2
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mgamma(scale = 5, shape = paras[j])
      noise_true_var <- ml_vgamma(scale = 5, shape = paras[j])
      noise_true_skew <- ml_sgamma(scale = 5, shape = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        rgamma(n = n, scale = 5, shape = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # GMM distro, weight para.
  if (dist == "gmm"){
    
    comp_sd <- c(0.3, 43)
    comp_mean <- c(-15, 15)
    paras <- c(0.00001, 0.0025, seq(0.005, 0.1625, 0.0125))
    weight_1 <- 1 - paras
    weight_2 <- 0 + paras
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mgmm(k = length(comp_mean), comp_weight = c(weight_1[j], weight_2[j]), comp_mean = comp_mean, comp_sd = comp_sd)
      noise_true_var <- ml_vgmm(k = length(comp_mean), comp_weight = c(weight_1[j], weight_2[j]), comp_mean = comp_mean, comp_sd = comp_sd)
      noise_true_skew <- ml_sgmm(k = length(comp_mean), comp_weight = c(weight_1[j], weight_2[j]), comp_mean = comp_mean, comp_sd = comp_sd)
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rgmm(n = n, k = length(comp_mean), comp_weight = c(weight_1[j], weight_2[j]), comp_mean = comp_mean, comp_sd = comp_sd)
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Gompertz distro, shape para.
  if (dist == "gompertz"){
    
    paras <- seq(1, 0.0001, -{1-0.0001}/14)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mgomp(shape = paras[j], scale = 15)
      noise_true_var <- ml_vgomp(shape = paras[j], scale = 15)
      noise_true_skew <- NaN
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rgomp(n = n, shape = paras[j], scale = 15)
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Gumbel distro, scale para.
  if (dist == "gumbel"){
    
    paras <- seq(0.3, 20, {20-0.3}/14) * sqrt(6 / pi)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mgum(shift = 0, scale = paras[j])
      noise_true_var <- ml_vgum(shift = 0, scale = paras[j])
      noise_true_skew <- ml_sgum(shift = 0, scale = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rgum(n = n, shift = 0, scale = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Hyperbolic Secant distro, scale para.
  if (dist == "hsec"){
    
    paras <- seq(0.3, 20, {20 - 0.3}/14)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mhsec(shift = 0, scale = paras[j])
      noise_true_var <- ml_vhsec(shift = 0, scale = paras[j])
      noise_true_skew <- ml_shsec(shift = 0, scale = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rhsec(n = n, shift = 0, scale = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Laplace distro, scale para.
  if (dist == "laplace"){
    
    paras <- seq(0.3/sqrt(2), 20/sqrt(2), {20/sqrt(2) - 0.3/sqrt(2)}/14)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mlap(shift = 0, scale = paras[j])
      noise_true_var <- ml_vlap(shift = 0, scale = paras[j])
      noise_true_skew <- ml_slap(shift = 0, scale = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rlap(n = n, shift = 0, scale = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Log-Logistic distro, shape para.
  if (dist == "llog"){
    
    paras <- c(5.5, 4, 3, 2.5, 2.3, 2.25, 2.2, 2.15, 2.1, 2.05, 2.025, 2.0125, 2.01, 2.0075, 2.005)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mllog(shape = paras[j], scale = 1)
      noise_true_var <- ml_vllog(shape = paras[j], scale = 1)
      noise_true_skew <- NaN
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rllog(n = n, shape = paras[j], scale = 1)
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Log-Normal distro, scale para.
  if (dist == "lnorm"){
    
    paras <- c(0.3, 0.5, 0.75, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.55, 1.6, 1.625, 1.65, 1.7, 1.73)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mlnorm(shift = 0, scale = paras[j])
      noise_true_var <- ml_vlnorm(shift = 0, scale = paras[j])
      noise_true_skew <- ml_slnorm(shift = 0, scale = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        rlnorm(n = n, meanlog = 0, sdlog = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Logistic distro, scale para.
  if (dist == "logistic"){
    
    paras <- seq(0.3/sqrt(pi^2/3), 20/sqrt(pi^2/3), {20/sqrt(pi^2/3) - 0.3/sqrt(pi^2/3)}/14)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mlgst(shift = 0, scale = paras[j])
      noise_true_var <- ml_vlgst(shift = 0, scale = paras[j])
      noise_true_skew <- ml_slgst(shift = 0, scale = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rlgst(n = n, shift = 0, scale = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Maxwell distro, scale para.
  if (dist == "maxwell"){
    
    paras <- seq(0.3/sqrt({3 - 8/pi}), 20/sqrt({3 - 8/pi}), {20/sqrt({3 - 8/pi}) - 0.3/sqrt({3 - 8/pi})}/14)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mmaxw(scale = paras[j])
      noise_true_var <- ml_vmaxw(scale = paras[j])
      noise_true_skew <- ml_smaxw(scale = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rmaxw(n = n, scale = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Normal distro, scale para.
  if (dist == "norm"){
    
    paras <- seq(0.3, 20, {20-0.3}/14)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- 0
      noise_true_var <- paras[j]
      noise_true_skew <- 0
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        rnorm(n = n, mean = 0, sd = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Pareto distro, shape para.
  if (dist == "pareto"){
    
    paras <- numeric(15)
    tmp <- seq(0.3, 20, {20-0.3}/14)
    b <- 1
    for(i in 1:length(paras)){
      paras[i] <- Re(polyroot(z = c({-2}, {5-{b/tmp[i]}^2}, {-4}, {1}))[3]) # third is real one in cubic case.
    }
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mpareto(a = paras[j], b = 1)
      noise_true_var <- ml_vpareto(a = paras[j], b = 1)
      noise_true_skew <- ml_spareto(a = paras[j], b = 1)
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rpareto(n = n, a = paras[j], b = 1)
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Rayleigh distro, scale para.
  if (dist == "rayleigh"){
    
    paras <- seq(0.3/sqrt({2-pi/2}), 20/sqrt({2-pi/2}), {20/sqrt({2-pi/2}) - 0.3/sqrt({2-pi/2})}/14)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mrayl(scale = paras[j])
      noise_true_var <- ml_vrayl(scale = paras[j])
      noise_true_skew <- ml_srayl(scale = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rrayl(n = n, scale = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Semicircle distro, scale (radius) para.
  if (dist == "semicircle"){
    
    paras <- seq(0.3*2, 20*2, {20*2 - 0.3*2}/14)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mscir(center = 0, rad = paras[j])
      noise_true_var <- ml_vscir(center = 0, rad = paras[j])
      noise_true_skew <- ml_sscir(center = 0, rad = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rscir(n = n, center = 0, rad = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Student-t distro, df para.
  if (dist == "t"){
    
    tmp <- seq(1.1, 20, {20-1.1}/14)^2
    paras <- -2*tmp/{1-tmp}
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mt(df = paras[j])
      noise_true_var <- ml_vt(df = paras[j])
      noise_true_skew <- ml_st(df = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        rt(n = n, df = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Triangle (asymmetric) distro, scale para.
  if (dist == "tri_a"){
    
    paras <- seq(0.3/sqrt({1 - 0.9*{1-0.9}} / 18), 20/sqrt({1 - 0.9*{1-0.9}} / 18),
                 {20/sqrt({1 - 0.9*{1-0.9}} / 18) - 0.3/sqrt({1 - 0.9*{1-0.9}} / 18)}/14)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mtri(vert = 0.9, shift = 0, scale = paras[j])
      noise_true_var <- ml_vtri(vert = 0.9, shift = 0, scale = paras[j])
      noise_true_skew <- ml_stri(vert = 0.9, shift = 0, scale = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rtri(n = n, vert = 0.9, shift = 0, scale = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Triangle (symmetric) distro, scale para.
  if (dist == "tri_s"){
    
    paras <- seq(0.3/sqrt({1 - 0.5*{1-0.5}} / 18), 20/sqrt({1 - 0.5*{1-0.5}} / 18),
                 {20/sqrt({1 - 0.5*{1-0.5}} / 18) - 0.3/sqrt({1 - 0.5*{1-0.5}} / 18)}/14)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mtri(vert = 0.5, shift = 0, scale = paras[j])
      noise_true_var <- ml_vtri(vert = 0.5, shift = 0, scale = paras[j])
      noise_true_skew <- ml_stri(vert = 0.5, shift = 0, scale = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rtri(n = n, vert = 0.5, shift = 0, scale = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # U-power distro, scale para.
  if (dist == "upower"){
    
    paras <- seq(0.3/sqrt({2*2 + 1} / {2*2 + 3}), 20/sqrt({2*2 + 1} / {2*2 + 3}),
                 {20/sqrt({2*2 + 1} / {2*2 + 3}) - 0.3/sqrt({2*2 + 1} / {2*2 + 3})}/14)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mupwr(shape = 2, shift = 0, scale = paras[j])
      noise_true_var <- ml_vupwr(shape = 2, shift = 0, scale = paras[j])
      noise_true_skew <- ml_supwr(shape = 2, shift = 0, scale = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rupwr(n = n, shape = 2, shift = 0, scale = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Wald distro, shape para.
  if (dist == "wald"){
    
    paras <- 1/seq(0.3, 20, {20-0.3}/14)^2
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mwald(mean = 1, shape = paras[j])
      noise_true_var <- ml_vwald(mean = 1, shape = paras[j])
      noise_true_skew <- ml_swald(mean = 1, shape = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        ml_rwald(n = n, mean = 1, shape = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # Weibull distro, shape para.
  if (dist == "weibull"){
    
    paras <- c(0.99, 0.9, 0.8, 0.7, 0.65, 0.6, 0.55, 0.525, 0.5, 0.475, 0.45, 0.425, 0.4, 0.375, 0.35)
    
    noise <- list()
    noise_oracle_mean <- list()
    noise_oracle_var <- list()
    noise_oracle_skew <- list()
    
    for (j in 1:length(paras)){
      
      noise_true_mean <- ml_mweibull(scale = 1, shape = paras[j])
      noise_true_var <- ml_vweibull(scale = 1, shape = paras[j])
      noise_true_skew <- ml_sweibull(scale = 1, shape = paras[j])
      
      # generate a matrix (n by num_trials) for each "level".
      noise[[j]] <- mapply(FUN = function(t){
        rweibull(n = n, scale = 1, shape = paras[j])
      }, 1:num_trials) - noise_true_mean
      noise_oracle_mean[[j]] <- noise_true_mean
      noise_oracle_var[[j]] <- noise_true_var
      noise_oracle_skew[[j]] <- noise_true_skew
    }
    
    noise_list <- list(noise, noise_oracle_mean, noise_oracle_var, noise_oracle_skew)
    
    return(noise_list)
  }
  
  # If we come this far, then there must have been an issue with the dist argument.
  cat("Error (", fname, "): the distribution wasn't properly specified.", "\n", sep = "")
  return(NULL)
}
















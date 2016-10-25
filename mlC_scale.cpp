
// C++ routines related to scale estimation, for robust target minimizer.

#define _USE_MATH_DEFINES // to get access to M_PI etc from math.h.

//[[Rcpp::depends(RcppArmadillo)]]
#include<RcppArmadillo.h> // With this, don't need armadillo header.
#include<stdlib.h> //  Use fabs() when want to return real vals.
#include<math.h> // for log(), log1p(), pow(), tanh(), fabs(), sqrt(), atan(), exp().

using namespace Rcpp;
using namespace std; // for string, and i/o related names.

//// CONTENTS ////

// - Helper functions, not exported to R.
// - Scale estimation routines, exported to R.


//// MAIN BODY ////

// Helper functions, not exported to R.

//// A routine for turning the values of an array into their abs vals.
arma::vec absArr(arma::vec x, int n){

  arma::vec u(n);
  for (int i=0; i<n; i++){
    u(i) = fabs(x(i));
  }

  return u;
}

//// Huber's proposal 2.
double chi_huber2(double x, double c, double low){

  if (fabs(x) <= c){
    return pow(x,2) - low;
  }
  else {
    return pow(c,2) - low;
  }

}

//// Sum divided by sample size, Huber's proposal 2.
double chi_huber2_n(arma::vec x, int n, double c, double low){

  double out = 0;
  for (int i=0; i<n; i++){
    out += chi_huber2(x(i),c,low);
  }

  return out/n;
}


// Scale computation routines to be exported to R.

//// MAD about the mean.
// [[Rcpp::export]]
double scale_madmean(arma::vec x, int n){
  double meanval = mean(x);
  return median(absArr(x-meanval,n))/0.8263717;
}

//// Iterative scale estimate using Huber proposal 2.
// [[Rcpp::export]]
double scale_huber2(arma::vec x, int n, double thres, int iters, double c, double low){

  double s_new = stddev(x); // Initialize to sd.
  double s_old;
  double diff = 1; // noting thres << 1.
  double pivot = mean(x); // Pivot set as sample mean.
  arma::vec num = x - pivot;
  arma::vec val;

  for (int t=0; t<iters; t++){
    s_old = s_new;
    val = num/s_old;
    s_new = s_old * sqrt(1 + chi_huber2_n(val,n,c,low)/low);
    diff = fabs(s_new - s_old);
    if (diff <= thres){
      break;
    }
  }

  return s_new;
}




// C++ routines related to location estimation, for robust target minimizer.

#define _USE_MATH_DEFINES // to get access to M_PI etc from math.h.

//[[Rcpp::depends(RcppArmadillo)]]
#include<RcppArmadillo.h> // With this, don't need armadillo header.
#include<stdlib.h> //  Use fabs() when want to return real vals.
#include<math.h> // for log(), log1p(), pow(), tanh(), fabs(), sqrt(), atan(), exp().

using namespace Rcpp;
using namespace std; // for string, and i/o related names.

//// CONTENTS ////

// - Helper functions, not exported to R.
// - Location estimation routines, exported to R.


//// MAIN BODY ////

// Helper functions, not exported to R.

//// Lower/upper sides of the Catoni truncator class.
double catlow(double x){
  return -log1p(-x+pow(x,2)/2);
}
double catup(double x){
  return log1p(x+pow(x,2)/2);
}

//// The narrowest member of the Catoni truncator class.
double psi_catnar(double x){

  double out;

  if(x >= 0){
    out = catlow(x);
  }
  else {
    out = catup(x);
  }
  if (fabs(x) >= 1){
    out = copysign(log(2), x);
  }

  return out;
}

//// The widest member of the Catoni truncator class.
double psi_catwide(double x){

  double out;

  if(x >= 0){
    out = catup(x);
  }
  else {
    out = catlow(x);
  }

  return out;
}

//// Convex comb. of wide/narrow Catoni class limits.
//// Weighting is [s]ymmetric, ie same on both sides of the origin.
double psi_catsmix(double x, double wt){

  double out = (1-wt)*psi_catnar(x) + wt*psi_catwide(x);

  return out;
}
//// In contrast to catsmix, this fn allows for [a]symmetric weighting.
double psi_catamix(double x, double wt1, double wt2){

  double out;
  if (x>=0){
    out = psi_catsmix(x, wt1);
  }
  else {
    out = psi_catsmix(x, wt2);
  }

  return out;
}

//// The "standard" logistic function, odd, rate 1, slope 1 at origin.
double psi_lgst(double x){
  return 4*(1/(1+exp(-x)) - 0.5);
}

//// A general logistic function, odd, slope 1 at origin, choose rate.
double psi_lgst_gen(double x, double rate){
  double c1 = 4/rate;
  return c1*(1/(1+exp(-x*rate)) - 0.5);
}

//// The Gudermannian function, fixed to have slope 1 at origin.
double psi_gud(double x){
  return 2*atan(exp(x))-M_PI/2;
}

//// The value of the "psi-condition" term, for Catoni mixture.
double psi_catamix_n(arma::vec x, int n, double wt1, double wt2){

  double out = 0;
  for (int i=0; i < n; i++){
    out += psi_catamix(x(i), wt1, wt2);
  }

  return out/n;
}

//// The value of the "psi-condition" term, for widest Catoni truncator.
double psi_catwide_n(arma::vec x, int n){

  double out = 0;
  for (int i=0; i < n; i++){
    out += psi_catwide(x(i));
  }

  return out/n;
}

//// The value of the psi-condition term, for narrowest Catoni truncator.
double psi_catnar_n(arma::vec x, int n){

  double out = 0;
  for (int i=0; i < n; i++){
    out += psi_catnar(x(i));
  }

  return out/n;
}

//// The value of the psi-condition term, for rho=log(cosh).
double psi_lcosh_n(arma::vec x, int n){

  double out = 0;
  for (int i=0; i<n; i++){
    out += tanh(x(i));
  }

  return out/n;
}

//// The value of the psi-condition term, for std logistic.
double psi_lgst_n(arma::vec x, int n){

  double out = 0;
  for (int i=0; i<n; i++){
    out += psi_lgst(x(i));
  }

  return out/n;
}

//// The value of the psi-condition term, for general logistic.
double psi_lgst_gen_n(arma::vec x, int n, double rate){

  double out = 0;
  for (int i=0; i<n; i++){
    out += psi_lgst_gen(x(i),rate);
  }

  return out/n;
}

//// The value of the psi-condition term, for Gudermannian.
double psi_gud_n(arma::vec x, int n){

  double out = 0;
  for (int i=0; i<n; i++){
    out += psi_gud(x(i));
  }

  return out/n;
}


// Location estimation routines to be exported to R.

//// Solution to psi-condition term under Catoni mixture.
// [[Rcpp::export]]
double est_catmix(arma::vec x, int n, double s, double thres, int iters, double wt1, double wt2){

  double new_theta = mean(x); // Initialize at sample mean.
  double old_theta, diff=1; // noting thres << 1.

  for (int t=0; t<iters; t++){
    old_theta = new_theta;
    new_theta = psi_catamix_n((x-old_theta)/s,n,wt1,wt2)*s + old_theta;
    diff = fabs(new_theta - old_theta);
    if (diff <= thres){
      break;
    }
  }

  return new_theta;
}

//// Solution to psi-condition term under narrow Catoni.
// [[Rcpp::export]]
double est_catnar(arma::vec x, int n, double s, double thres, int iters){

  double new_theta = mean(x); // Initialize at sample mean.
  double old_theta, diff=1; // noting thres << 1.

  for (int t=0; t<iters; t++){
    old_theta = new_theta;
    new_theta = psi_catnar_n((x-old_theta)/s,n)*s + old_theta;
    diff = fabs(new_theta - old_theta);
    if (diff <= thres){
      break;
    }
  }

  return new_theta;
}

//// Solution to psi-condition term under wide Catoni.
// [[Rcpp::export]]
double est_catwide(arma::vec x, int n, double s, double thres, int iters){

  double new_theta = mean(x); // Initialize at sample mean.
  double old_theta, diff=1; // noting thres << 1.

  for (int t=0; t<iters; t++){
    old_theta = new_theta;
    new_theta = psi_catwide_n((x-old_theta)/s,n)*s + old_theta;
    diff = fabs(new_theta - old_theta);
    if (diff <= thres){
      break;
    }
  }

  return new_theta;
}

//// Solution to psi-condition term under Gudermannian.
// [[Rcpp::export]]
double est_gud(arma::vec x, int n, double s, double thres, int iters){

  double new_theta=mean(x); // Initialize at sample mean.
  double old_theta, diff=1; // noting thres << 1.

  // Solve the psi-condition.
  for (int t=0; t<iters; t++){
    old_theta = new_theta;
    new_theta = psi_gud_n((x-old_theta)/s,n)*s + old_theta;
    diff = fabs(new_theta - old_theta);
    if (diff <= thres){
      break;
    }
  }

  return new_theta;
}

//// Solution to psi-condition term under logistic function.
// [[Rcpp::export]]
double est_lgst(arma::vec x, int n, double s, double thres, int iters){

  double new_theta=mean(x); // Initialize at sample mean.
  double old_theta, diff=1; // noting thres << 1.

  // Solve the psi-condition.
  for (int t=0; t<iters; t++){
    old_theta = new_theta;
    new_theta = psi_lgst_n((x-old_theta)/s,n)*s + old_theta;
    diff = fabs(new_theta - old_theta);
    if (diff <= thres){
      break;
    }
  }

  return new_theta;
}

//// Solution to psi-condition term under generalized logistic function.
// [[Rcpp::export]]
double est_lgst_gen(arma::vec x, int n, double s, double thres, int iters, double rate){

  double new_theta=mean(x); // Initialize at sample mean.
  double old_theta, diff=1; // noting thres << 1.

  // Solve the psi-condition.
  for (int t=0; t<iters; t++){
    old_theta = new_theta;
    new_theta = psi_lgst_gen_n((x-old_theta)/s,n,rate)*s + old_theta;
    diff = fabs(new_theta - old_theta);
    if (diff <= thres){
      break;
    }
  }

  return new_theta;
}


//// Solution to psi-condition term under tanh.
// [[Rcpp::export]]
double est_lcosh(arma::vec x, int n, double s, double thres, int iters){

  double new_theta=mean(x); // Initialize at sample mean.
  double old_theta, diff=1; // noting thres << 1.

  // Solve the psi-condition.
  for (int t=0; t<iters; t++){
    old_theta = new_theta;
    new_theta = psi_lcosh_n((x-old_theta)/s,n)*s + old_theta;
    diff = fabs(new_theta - old_theta);
    if (diff <= thres){
      break;
    }
  }

  return new_theta;
}




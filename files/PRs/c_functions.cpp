#include <Rcpp.h>
//#include <RcppGSL.h>
#include <random>
#include <iostream>
#include <vector>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_eigen.h>
using namespace Rcpp;


// Function below generates random AR(1) sequences
// [[Rcpp::export]]
NumericVector arp_C(int n, double c, double phi, double eps) {
  NumericVector x(n);
  double value;
  // Generate a random number from the normal distribution
  value = eps/sqrt(1-phi*phi) * R::rnorm(0,1);
  x[0] = value + c/(1-phi);
  // Loop from p to n
  for(int i = 1; i < n; i++) {
    // Generate a random number from the normal distribution
    value = eps * R::rnorm(0,1);
    x[i] = value + c + phi*x[i-1];
  }
  return x;
}

// Function below generates bivariate correlated AR(1) sequences
// [[Rcpp::export]]
NumericMatrix ar_C_2d(int n, double c_x, double c_y, double phi_x, double phi_y, double eps_x, double eps_y, double cor) {
  NumericMatrix z(n,2);
  // Initialize
  double noise_x = R::rnorm(0,1);
  double noise_y = cor * noise_x + sqrt(1-cor*cor)*R::rnorm(0,1);
  double value_x = eps_x/sqrt(1-phi_x*phi_x) * noise_x;
  double value_y = eps_y/sqrt(1-phi_y*phi_y) * noise_y;
  z(0,0) = value_x + c_x/(1-phi_x);
  z(0,1) = value_y + c_y/(1-phi_y);
  // Loop from p to n
  for(int i = 1; i < n; i++) {
    noise_x = R::rnorm(0,1);
    noise_y = cor * noise_x + sqrt(1-cor*cor)*R::rnorm(0,1);    
    value_x = eps_x * noise_x;
    value_y = eps_y * noise_y;
    z(i,0) = value_x + c_x + phi_x*z(i-1,0);
    z(i,1) = value_y + c_y + phi_y*z(i-1,1);
  }
  return z;
}



// Function below computes the mean of a vector
// [[Rcpp::export]]
double mean_C(NumericVector x) {
  int p = x.size();
  double m = 0;
  for(int j = 0; j < p; j++){
    m += x[j];
  }
  m = m/p;
  return m;
}

// Function below computes the covariance between two vectors
// [[Rcpp::export]]
double cov_C(NumericVector x, NumericVector y) {
  int p = x.size();
  double m = 0;
  double mx = mean_C(x);
  double my = mean_C(y);
  for(int j = 0; j < p; j++){
    m = m + (x[j]-mx)*(y[j]-my);
  }
  m = m/(p-1);
  return m;
}

// Function below computes the covariance between two vectors with bias
// [[Rcpp::export]]
double cov_C_bias(NumericVector x, NumericVector y) {
  int p = x.size();
  double m = 0;
  double mx = mean_C(x);
  double my = mean_C(y);
  for(int j = 0; j < p; j++){
    m += (x[j]-mx)*(y[j]-my);
  }
  m = m/p;
  return m;
}

// Function below computes the variance of a vector
// [[Rcpp::export]]
double var_C(NumericVector x) {
  int p = x.size();
  double m = 0;
  double mx = mean_C(x);
  for(int j = 0; j < p; j++){
    m += (x[j]-mx)*(x[j]-mx);
  }
  m = m/(p-1);
  return m;
}

// Function below computes the variance of a vector with bias
// [[Rcpp::export]]
double var_C_bias(NumericVector x) {
  int p = x.size();
  double m = 0;
  double mx = mean_C(x);
  for(int j = 0; j < p; j++){
    m += (x[j]-mx)*(x[j]-mx);
  }
  m = m/p;
  return m;
}

// Function below computes the regression coefficients
// [[Rcpp::export]]
NumericVector reg_C(NumericVector y, NumericVector x){
  NumericVector z(2);
  z[1] = cov_C(x,y)/var_C(x);
  z[0] = mean_C(y) - z[1]*mean_C(x);
  return z;
}

// Function below computes MSE in predictive regressions
// [[Rcpp::export]]
double quad_err_C(double r_x, double r_y, double s_x, double s_y, int TT, int n_sim){
  double mse = 0;
  double x_new;
  double y_new;
  double err;
  NumericVector x(TT);
  NumericVector y(TT);
  NumericVector coef(2);
  for(int i = 0; i<n_sim; i++){
    x = arp_C(TT, 0, r_x, s_x ); // * sqrt(1-r_x*r_x)
    y = arp_C(TT, 0, r_y, s_y ); // * sqrt(1-r_y*r_y)
    x_new = x[TT-1] * r_x + s_x * R::rnorm(0,1); // * sqrt(1-r_x*r_x)
    y_new = y[TT-1] * r_y + s_y * R::rnorm(0,1); // * sqrt(1-r_y*r_y)
    coef = reg_C(y, x);
    err = y_new - coef[0] - coef[1] * x_new;
    mse += err*err;
  }
  mse = mse/n_sim;
  return mse;
}

// [[Rcpp::export]]
double quad_err_C_2d(double a_x, double a_y, double r_x, double r_y, double s_x, double s_y, double cor, int TT, int n_sim){
  double mse = 0;
  double x_new;
  double y_new;
  double err;
  NumericMatrix z(TT, 2); // Where the bivariate process is stored
  NumericVector x(TT);    // Values for x (first col of z)
  NumericVector y(TT);    // Values for y (second col of z)
  NumericVector coef(2);  // Values for reg estimate (a and b)
  for(int i = 0; i<n_sim; i++){
    z = ar_C_2d(TT, a_x, a_y, r_x, r_y, s_x, s_y, cor);
    x = z(_,0); 
    y = z(_,1); 
    double noise_x = R::rnorm(0,1);
    double noise_y = cor * noise_x + sqrt(1-cor*cor)*R::rnorm(0,1);
    double value_x = s_x * noise_x;
    double value_y = s_y * noise_y;
    x_new = value_x + a_x + r_x * x(TT-1);
    y_new = value_y + a_y + r_y * y(TT-1);
    coef = reg_C(y, x);
    err = y_new - coef[0] - coef[1] * x_new;
    mse += err*err;
  }
  mse = mse/n_sim;
  return mse;
}


// [[Rcpp::export]]
NumericVector slicer(NumericVector v, int begin, int end){ // This one is ugly!
  NumericVector z(end-begin+1);
  for(int i = begin; i < end+1; i++){
    z[i-begin] = v[i];
  }
  return z;
}

// [[Rcpp::export]]
double quad_err_C_auto(double r_x, double s_x, int TT, int n_sim){
  double mse = 0;
  double x_new = 0;
  double err = 0;
  NumericVector x(TT);
  NumericVector coef(2);
  for(int i = 0; i<n_sim; i++){
    x = arp_C(TT, 0, r_x, s_x );  // * sqrt(1-r_x*r_x)
    x_new = x[TT-1] * r_x + s_x * R::rnorm(0,1); // * sqrt(1-r_x*r_x)
    NumericVector x_0 = slicer(x,0,TT-2);
    NumericVector x_1 = slicer(x,1,TT-1);
    coef = reg_C(x_0, x_1);
    err = x_new - coef[0] - coef[1] * x[TT-1];
    mse += err*err;
  }
  mse = mse/n_sim;
  return mse;
}

// [[Rcpp::export]]
NumericMatrix mmult(const NumericMatrix& m1, const NumericMatrix& m2){
  if (m1.ncol() != m2.nrow()) stop ("Incompatible matrix dimensions");
  NumericMatrix out(m1.nrow(),m2.ncol());
  NumericVector rm1, cm2;
  for (size_t i = 0; i < m1.nrow(); ++i) {
    rm1 = m1(i,_);
    for (size_t j = 0; j < m2.ncol(); ++j) {
      cm2 = m2(_,j);
      out(i,j) = std::inner_product(rm1.begin(), rm1.end(), cm2.begin(), 0.);              
    }
  }
  return out;
}


// [[Rcpp::export]]
NumericMatrix call_mvt(int n, NumericVector mu, NumericMatrix Sigma, Function f) {
  NumericMatrix res = f(n, mu, Sigma);
  return res;
}

// [[Rcpp::export]]
NumericVector subind(NumericVector x, int beg, int end) {
  IntegerVector idx = seq(beg, end);  
  return x[idx];
}

// [[Rcpp::export]]
NumericVector quad_err_C_2d2(double a_x, double a_y, double r_x, double r_y, double s_x, double s_y, double cor, int TT, int k, int n_sim){
  double mse = 0;
  double err;
  NumericMatrix z(TT+2*k, 2); // Where the bivariate process is stored
  NumericVector x_train(TT);  // Values for x (first col of z)
  NumericVector y_train(TT);  // Values for y (second col of z)
  NumericVector x(TT+2*k);    // Values for x (first col of z)
  NumericVector y(TT+2*k);    // Values for y (second col of z)
  NumericVector coef(2);      // Values for reg estimate (a and b)
  NumericVector e(n_sim);   // Storing individual errors
  NumericVector ret(2);

  for(int i = 0; i<n_sim; i++){
    z = ar_C_2d(TT+2*k, a_x, a_y, r_x, r_y, s_x, s_y, cor);

    x = z(_,0); 
    y = z(_,1); 
    
    x_train = subind(x, 0, TT-1);
    y_train = subind(y, k, TT-1+k);

    coef = reg_C(y_train, x_train);
    err = z(TT-1+2*k,1) - coef[0] - coef[1] * z(TT-1+k,0);
    
    e(i) = err*err;
    mse += err*err;
  }
  mse = mse/n_sim;
  ret(0) = mse;
  ret(1) = var_C(e);
  return ret;
}

// [[Rcpp::export]]
NumericVector quad_err_C_2d_vec(double a_x, double a_y, double r_x, double r_y, double s_x, double s_y, double cor, int TT, int n_sim){
  NumericVector mse(n_sim);
  double x_new;
  double y_new;
  double err;
  NumericMatrix z(TT, 2); // Where the bivariate process is stored
  NumericVector x(TT);    // Values for x (first col of z)
  NumericVector y(TT);    // Values for y (second col of z)
  NumericVector coef(2);  // Values for reg estimate (a and b)
  for(int i = 0; i<n_sim; i++){
    z = ar_C_2d(TT, a_x, a_y, r_x, r_y, s_x, s_y, cor);
    x = z(_,0); 
    y = z(_,1); 
    double noise_x = R::rnorm(0,1);
    double noise_y = cor * noise_x + sqrt(1-cor*cor)*R::rnorm(0,1);
    double value_x = s_x * noise_x;
    double value_y = s_y * noise_y;
    x_new = value_x + a_x + r_x * x(TT-1);
    y_new = value_y + a_y + r_y * y(TT-1);
    coef = reg_C(y, x);
    err = y_new - coef[0] - coef[1] * x_new;
    mse(i) = err;
  }
  return mse;
}

// [[Rcpp::export]]
double quad_err_C_uni(double a_y, double r_y, double s_y, int TT, int k, int n_sim){
  double mse = 0;
  double err;
  NumericVector y_train0(TT);  // Values for y (past)
  NumericVector y_train1(TT);  // Values for y  (future)
  NumericVector y(TT+2*k);    // Values for y (second col of z)
  NumericVector coef(2);      // Values for reg estimate (a and b)
  
  for(int i = 0; i<n_sim; i++){
    y = arp_C(TT+2*k, a_y, r_y, s_y); 
    y_train0 = subind(y, 0, TT-1);
    y_train1 = subind(y, k, TT+k-1);
    coef = reg_C(y_train1, y_train0);
    err = y(TT+2*k-1) - coef[0] - coef[1] * y(TT+k-1);
    
    mse += err*err;
  }
  mse = mse/n_sim;
  return mse;
}
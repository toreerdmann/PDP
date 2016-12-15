/*
 * Implementations of the function for evaluation of the posterior predictive
 * density and needed utility functions (mahalanobis distance, t-pdf).
 */

// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

const double logpi = std::log(M_PI); 
const double log2pi = std::log(2.0 * M_PI);

using namespace Rcpp;
  
// [[Rcpp::export]]
arma::vec Mahalanobis(arma::mat& x, arma::rowvec& center, arma::mat& cov) {
    int n = x.n_rows;
    arma::mat x_cen;
    x_cen.copy_size(x);
    for (int i=0; i < n; i++) {
        x_cen.row(i) = x.row(i) - center;
    }
    return sum((x_cen * cov.i()) % x_cen, 1);
}

// [[Rcpp::export]]
arma::vec dmvt_arma(arma::mat x, arma::rowvec mean, arma::mat sigma, double nu) {
  // double n = (double) x.n_rows;
  double d = (double) x.n_cols;
  arma::vec distval = Mahalanobis(x, mean, sigma);
  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  return lgamma((nu + d) / 2) -
    lgamma(nu / 2) -
    d/2 * (log(nu) + logpi) +
    (-.5 * logdet) -
    ((nu + d) / 2) * log(1 + (1/nu) * distval);
}


// [[Rcpp::export]]
arma::vec Mahalanobis2(const arma::mat& x, const arma::rowvec& center, 
                       const arma::mat& chol_S) {
    int n = x.n_rows;
    arma::mat x_cen;
    x_cen.copy_size(x);
    for (int i=0; i < n; i++) {
        x_cen.row(i) = x.row(i) - center;
    }
    return sum((x_cen * (chol_S.t() * chol_S).i()) % x_cen, 1);
}

// [[Rcpp::export]]
arma::vec Mahalanobis_chol(const arma::mat& x, 
                           const arma::rowvec& center, 
                           const arma::mat& chol_S) {
    int n = x.n_rows;
    int d = x.n_cols;
    arma::mat cholDec = arma::trimatl(chol_S.t());
    arma::vec D = cholDec.diag();
    
    arma::vec out(n);
    arma::vec tmp(d); 
    double acc;
    uint32_t icol, irow, ii;  
    // For each of the "n" random vectors, forwardsolve the corresponding linear system.
    // Forwardsolve because I'm using the lower triangle Cholesky.
    for(icol = 0; icol < n; icol++) {
      for(irow = 0; irow < d; irow++) {
        acc = 0.0;
        for(ii = 0; ii < irow; ii++) acc += tmp.at(ii) * cholDec.at(irow, ii);
        tmp.at(irow) = ( x.at(icol, irow) - center.at(irow) - acc ) / D.at(irow);
      }
      out.at(icol) = sum(square(tmp)); 
    }
    return out;
}

// [[Rcpp::export]]
arma::vec dmvt_chol(const arma::mat& x, const arma::rowvec& mean, 
                    const arma::mat& chol_S, const double& nu) {
  // double n = (double) x.n_rows;
  double d = (double) x.n_cols;
  arma::vec distval = Mahalanobis2(x, mean, chol_S);
  double logdet = sum(2 * arma::log(chol_S.diag()));
  return lgamma((nu + d) / 2) -
    lgamma(nu / 2) -
    d/2 * (log(nu) + logpi) +
    (-.5 * logdet) -
    ((nu + d) / 2) * log(1 + (1/nu) * distval);
}
/*** R

testthat::test_that("mahalanobis and dmvt_arma work in 1D", {
  ## 1D example
  logpi = log(pi)
  d = 1; nu = 2
  m = 2; sigma = matrix(1)
  x = matrix(1:5)
  logdet = determinant(sigma)$modulus
  distval = mahalanobis(x, m, sigma)
  res_man = lgamma((nu + d) / 2) - lgamma(nu / 2) - d/2 * log(nu) - d/2 * logpi +
    1/2 * logdet - (nu + d) / 2 * log(1 + (1/nu) * matrix(distval))
  res_cpp = dmvt_arma(x, m, sigma, nu)
  res_cpp2 = dmvt_chol(x, m, chol(sigma), nu)
  
  ## check against manual implementation and mvtnorm::dmvt
  testthat::expect_equal(res_cpp, res_man)
  testthat::expect_equal(res_cpp, res_cpp2)
  testthat::expect_equal(res_cpp, matrix(mvtnorm::dmvt(x, m, sigma, nu)))
})

testthat::test_that("mahalanobis and dmvt_arma work in 3D", {
  ## 3D example
  logpi = log(pi)
  d = 3; nu = 4
  m = c(2,2,2); sigma = matrix(c(1,0,0,0,1,0,0,0,1), 3, 3)
  x = matrix(rnorm(30), 10, 3)
  logdet = determinant(sigma)$modulus
  distval = mahalanobis(x, m, sigma)
  res_man = lgamma((nu + d) / 2) - lgamma(nu / 2) - d/2 * log(nu) - d/2 * logpi +
    1/2 * logdet - (nu + d) / 2 * log(1 + (1/nu) * matrix(distval))
  res_cpp = dmvt_arma(x, m, sigma, nu)
  testthat::expect_equal(res_cpp, res_man)
  testthat::expect_equal(res_cpp, matrix(mvtnorm::dmvt(x, m, sigma, nu)))
})
*/


// [[Rcpp::export]]
double log_pred(arma::mat x,
                        arma::rowvec mean,
                        double kappa, double nu,
                        arma::mat sigma) {
  using namespace Rcpp;
  double n = x.n_rows;
  double d = x.n_cols;
  return  sum(dmvt_arma(x, mean,
                        (kappa + n) / (kappa * (nu - d + n)) * sigma,
                        nu - d + n));
}

/*** R
testthat::test_that("log pred works", {
  set.seed(1)
  x = matrix(rnorm(30), 10, 3)
  n <- nrow(x); d = ncol(x)
  m = c(2,2,2); sigma = matrix(c(1,0,0,0,1,0,0,0,1), 3, 3)
  kappa = 2; nu = 2
  ## transformed sigma
  sigi = ((kappa + n) / (kappa * (nu - d + n)) * sigma)
  testthat::expect_equal(dmvt_arma(x, m, sigma, 2),
                         matrix(mvtnorm::dmvt(x, m, sigma, 2)))
  testthat::expect_equal(dmvt_arma(x, m, sigi, 2),
                         matrix(mvtnorm::dmvt(x, m, sigi, 2)))
  testthat::expect_equal(sum(mvtnorm::dmvt(x, m, (kappa + n) / (kappa * (nu - d + n)) * sigma, nu - d + n, log = TRUE)),
                         log_pred(x, m, kappa, nu, sigma))
})
## and some benchmarking
library(microbenchmark)
x = mvtnorm::rmvnorm(9999, m, sigma)
res = microbenchmark(times = 1000,
  mahalanobis(x, m, sigma),
  Mahalanobis(x, m, sigma),
  Mahalanobis2(x, m, sigma),
  Mahalanobis_chol(x, m, chol_S)
)
res
# m = c(2,2,2); sigma = MCMCpack::rwish(3, matrix(c(1,0,0,0,1,0,0,0,1), 3, 3))
# kappa = 2; nu = 2; chol_S = chol(sigma)
# x = mvtnorm::rmvnorm(99999, m, sigma)
# res = microbenchmark(mvtnorm::dmvt(x, m, sigma, nu),
#                      mvnfast::dmvt(x, m, sigma, nu),
#                      dmvt_arma(x, m, sigma, nu),
#                      dmvt_chol(x, m, chol_S, nu))
# res
*/


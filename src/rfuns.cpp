// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
mat rmvnorm(int n, rowvec mu, mat sigma, bool is_cholesky) {
  int d = sigma.n_cols; 
  mat x(n, d);
  x.randn(n, d);
  if (is_cholesky)
    return repmat(mu, n, 1) + (sigma.t() * x.t()).t();
  else
    stop("Need cholesky");
}

// [[Rcpp::export]]
mat rwish(int nu, mat sigma, bool is_cholesky) {
  int d = sigma.n_cols;
  rowvec mu(d); mu.fill(0.0);
  mat x = rmvnorm(nu, mu, sigma, is_cholesky);
  mat rval = mat(d, d); rval.fill(0.0);
  for (int l=0; l<nu; l++) 
    rval += x.row(l).t() * x.row(l);
  return rval;
}

/*** R
testthat::test_that("my rmvnorm works", {
  sigma = matrix(c(1, .6, .6, 1), 2, 2)
  chol_S = chol(sigma)
  plot(rmvnorm(100, c(1, 1), chol_S, TRUE))
  points(mvtnorm::rmvnorm(100, c(1, 1), sigma), col = 2)
})

testthat::test_that("my rinvwish works", {
  nu = 4
  sigma = matrix(c(1, .6, .6, 1), 2, 2)
  chol_S = chol(sigma)
  rinvwish(nu, sigma, TRUE)
  MCMCpack::rwish(nu, sigma)
})
*/

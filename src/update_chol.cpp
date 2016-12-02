// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Function to perform the rank-1 Cholesky update, needed for updating the
// covariance matrix of a normal-inverse wishart prior
//[[Rcpp::export]]
void update_chol(arma::mat& L, const arma::rowvec& v1, bool downdate)
{
  arma::vec v = v1.t();
  double sign = 1.0;
  if (downdate) {
    // Perform the downdate instead
    sign = -1.0;
  }
  for (int k=0; k<L.n_rows; k++) {
    double r = sqrt( L(k,k) * L(k,k) + sign * v(k) * v(k) );
    double c = r / L(k,k);
    double s = v(k) / L(k,k);
    L(k,k) = r;
    if (k < L.n_rows-1) {
      L(k,arma::span(k+1,L.n_rows-1)) = (L(k,arma::span(k+1,L.n_rows-1)) +
        sign * s * v(arma::span(k+1,v.n_elem-1)).t()) / c;
      v(arma::span(k+1,v.n_elem-1)) = c * v(arma::span(k+1,v.n_elem-1)) -
        s * L(k,arma::span(k+1,L.n_rows-1)).t();
    }
  }
  // if diagonal elements are negative, flip the sign of that row
  for (int k=0; k<L.n_rows; k++) {
    if (L(k,k) < 0)
      L(k,arma::span(0,L.n_cols)) = L(k,arma::span(0,L.n_cols)) * -1;
  }
}


/*** R
library(SamplerCompare)
# set.seed(1)
## sample a random covariance matrix
sigma = matrix(c(1, .3, .3, 1), 2, 2)

testthat::test_that("chud and my implementation agree", {
  S0 = MCMCpack::rwish(2, sigma)
  ## compute the cholesky decomposition
  c0 = chol(S0)
  c00 = chol(S0)
  testthat::expect_equal(t(c0) %*% c0, S0)
  ## now update with a data vector
  x0 = x = rnorm(2)
  c10 = c1 = chud(c0, x)
  update_chol(c0, x, FALSE)
  ## data is not altered
  testthat::expect_equal(x0, x)
  testthat::expect_equal(c10, c1)
  # testthat::expect_equal(c1, c0) # these are not equal, since chud does not correct when the diagonal was negative
  S1 =  S0 + x %*% t(x)
  testthat::expect_equal(t(c0) %*% c0, S1)
  testthat::expect_equal(t(c1) %*% c1, S1)
})

testthat::test_that("up- and downdate return original", {
  x = rnorm(2)
  S0 = MCMCpack::rwish(2, sigma)
  c0 = chol(S0)
  c00 = chol(S0)
  testthat::expect_equal(c00, c0)
  update_chol(c0, x, FALSE)
  testthat::expect_false(all(c00 == c0))
  update_chol(c0, x, TRUE)
  testthat::expect_equal(c00, c0)
})

testthat::test_that("chud works", {
  sigma = matrix(c(1, .3, .3, 1), 2, 2)
  S0 = MCMCpack::rwish(2, sigma)
  ## compute the cholesky decomposition
  c0 = chol(S0)
  testthat::expect_equal(t(c0) %*% c0, S0)
  ## now update with a data vector
  x = rnorm(2)
  c1 = chud(c0, x)
  S1 =  S0 + x %*% t(x)
  testthat::expect_equal(t(c1) %*% c1, S1)
})

testthat::test_that("updating with a factor works", {
  fac = .5
  S0 = MCMCpack::rwish(2, sigma)
  x0 = x = rnorm(2)
  c0 = chol(S0)
  testthat::expect_equal(t(c0) %*% c0, S0)
  update_chol(c0, sqrt(fac) * x, FALSE)
  S1 =  S0 + fac * x %*% t(x)
  testthat::expect_equal(t(c0) %*% c0, S1)
})
*/

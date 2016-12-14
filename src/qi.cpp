#include <RcppArmadillo.h>
#include <log_pred.cpp>

using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double qik_cpp(arma::mat x, int k,
                  IntegerVector z,
                  List components,
                  NumericVector Pi,
                  double lambda) {
  LogicalVector nk = z == k;
  if (sum(nk) == 0)
    stop("Cannot compute probability for empty component");
  IntegerVector neighbor_classes = z[Pi == 1];
  double nc_sum = sum(neighbor_classes == k);

  // extract relevant component
  List ll = components[k-1];
  arma::rowvec prior_mean   = as<arma::rowvec>(ll[1]);
  double prior_kappa        = ll[2];
  double prior_nu           = ll[3];
  arma::mat prior_sigma     = as<arma::mat>(ll[4]);

  return log_pred(x, prior_mean, prior_kappa,
                  prior_nu, prior_sigma) +
                    log(sum(nk)) + nc_sum * lambda;
}

/*** R
B = Matrix::bandSparse(10, k = c(-1, 1, -1)); I = diag(10)
P = kronecker(B, I) + kronecker(I, B)
lambda = 2
components = list(list(m = 0, kappa = 2, nu = 4, sigma = matrix(1)),
                  list(m = 2, kappa = 2, nu = 4, sigma = matrix(1)))
x = matrix(rnorm(100), nrow = 50)
z = rep(0, 100)
k = 1


testthat::test_that("postpred works", {
  source("~/projects/Masterthesis/code/PDP/R/nw_component.R")
  xi = matrix(x[1,1,drop=F])
  n <- nrow(x); d = ncol(x)
  prior = components[[1]]
  nv1 = normal_component(prior$m, prior$kappa, prior$nu, prior$sigma)
  z[1] = 1
  i = 1
  expect_equal(qik_cpp(xi, k = 1, z = z, components[[k]], P[i,], lambda),
               as.numeric(lambda * sum(z[P[i, ] == 1] == k) +
                            log(table(z[z != 0])) +
                            nv1$log_post_pred(xi)))
  xi = matrix(x[1:10,1,drop=F])
  n = nrow(x); d = ncol(x)
  expect_equal(qik_cpp(xi, k = 1, z = z, components[[k]], P[i,], lambda),
               as.numeric(lambda * sum(z[P[i, ] == 1] == k) +
                            log(table(z[z != 0])) +
                            nv1$log_post_pred(xi)))
})
*/

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double qi0_cpp(arma::mat x, List prior, double alpha) {
  return alpha +
    log_pred(x, prior["m"], prior["kappa"], prior["nu"], prior["sigma"]);
}


/*** R
## compare log prior pred
testthat::test_that("prior pred works", {
  set.seed(1)
  x = matrix(rnorm(100))
  d = 1; n = 1
  prior = components[[1]]
  expect_equal(mvtnorm::dmvt(x[1,,drop=F], prior$m,
                             (prior$kappa + n) / (prior$kappa * (prior$nu - d + n)) * prior$sigma,
                             prior$nu - d + n, log = TRUE),
               qi0_cpp(x[1,,drop=F], components[[1]], 1))
})

## multiple dimensions
testthat::test_that("prior pred works with multiple obs", {
  set.seed(1)
  x = matrix(rnorm(100), nrow = 50)
  d = ncol(x); n = 3
  prior =  list(m = c(0, 0), kappa = 2, nu = 2, sigma = matrix(c(1,0,0,1),2,2))
  expect_equal(
    sum(mvtnorm::dmvt(x[1:3,,drop=F], prior$m,
                      (prior$kappa + n) / (prior$kappa * (prior$nu - d + n)) * prior$sigma,
                      prior$nu - d + n, log = TRUE)),
    qi0_cpp(x[1:3,,drop=F], prior, 1))
})
*/


// [[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::export]]
arma::vec get_assignment_prob(arma::mat x,
                             IntegerVector z,
                             List components,
                             List prior,
                             NumericVector Pi,
                             double lambda, double alpha,
                             int debug) {
  int K = components.size();
  arma::vec q(K+1);
  if (debug == 1)
    print(wrap(0));
  if (K > 20) // set a hard threshold here
    q[0] = -1 * arma::datum::inf;
  else
    q[0] = qi0_cpp(x, prior, alpha);
  for (int j=1; j <= K; j++) {
    if (debug == 1)
      print(wrap(j));
    q[j] = qik_cpp(x, j, z, components, Pi, lambda);
  }
  return q;
}

/*** R
alpha = 1; lambda = 1
xi = x[1,,drop=F]
i = 1
z = rep(1:2, each = 50)
sigma = matrix(c(1,0,0,1),2,2)
components = list(
  comp1 = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma,
               prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)),
  comp2 = list(m = c(2, 2), kappa = 2, nu = 4, sigma = sigma,
               prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)))

prior = components[[1]][[5]]  
get_assigment_prob(xi, z, components, prior, P[i,], lambda, alpha)
 # # get_assigment_prob(x[4,,drop=F], z, components, P[i,], lambda, alpha)
 #
 # ## get assigments via R
 # prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)
 # nv1 = normal_component(c(0, 0), 2, 4, sigma)
 # nv2 = normal_component(c(2, 2), 2, 4, sigma)
 # ## qi0 from R
 # qR = unname(c(log(alpha) + nv1$log_prior_pred(xi),
 #               qik(x, i, z, list(nv1, nv2), P, lambda)))
 #
 # expect_equal(get_assigment_prob(xi, z, components, P[i,], lambda, alpha)[,1],
 #              qR)
 # xi = x[1:3,,drop=F]
 # qR = unname(c(log(alpha) + nv1$log_prior_pred(xi),
               # qik(x, i, z, list(nv1, nv2), P, lambda)))
 # expect_equal(get_assigment_prob(xi, z, components, P[i,], lambda, alpha)[,1],
              # qR)
*/

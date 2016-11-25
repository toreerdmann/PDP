#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void update(List components, int k, const arma::mat x, int debug) {

  // extract kth component
  List ll = components[k-1];
  arma::rowvec prior_mean   = as<rowvec>(ll["m"]);
  double prior_kappa        = ll["kappa"];
  double prior_nu           = ll["nu"];
  arma::mat prior_sigma     = as<mat>(ll["sigma"]);

  // data stats
  double n = x.n_rows;
  double d = x.n_cols;

  // compute colmeans
  arma::rowvec x_mean(d);
  for (int j=0; j<d; j++)
    for (int i=0; i<n; i++)
      x_mean[j] += x(i, j) / n;

  // centered data crossprod / empirical covariance matrix
  arma::mat xc =
    trans(x - repmat(x_mean, n, 1)) * (x - repmat(x_mean, n, 1));

  if (debug == 1) {
    print(wrap(x_mean));
    print(wrap(prior_mean));
    print(wrap("xc:"));
    print(wrap(xc));
    print(wrap("x_mean - prior_mean * ..."));
    print(wrap(trans(x_mean - prior_mean) * (x_mean - prior_mean)));
    print(wrap("fac * x_mean - prior_mean * ..."));
    print(wrap(prior_kappa / (prior_kappa + n) *
      trans(x_mean - prior_mean) * (x_mean - prior_mean)));
  }

  // update stats
  ll["m"] =
    prior_kappa /  (prior_kappa + n) * prior_mean +
    n / (prior_kappa + n) * x_mean;
  ll["kappa"] = prior_kappa + n;
  ll["nu"]    = prior_nu    + n;
  ll["sigma"] =
    prior_sigma + xc +
    prior_kappa / (prior_kappa + n) *
    trans(x_mean - prior_mean) * (x_mean - prior_mean);

  if (debug == 1) {
    print(wrap("After updating:"));
    print(wrap(ll));
  }

  // replace with updated comp
  components[k-1] = ll;
}



/*** R
set.seed(1)
x = matrix(rnorm(100), nrow=50)
xi = x[1:10,,drop=F]
sigma = matrix(c(1,0,0,1),2,2)
components = list(
  comp1 = list(m = c(1, 1), kappa = 2, nu = 4, sigma = sigma,
               prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)),
  comp2 = list(m = c(2, 2), kappa = 2, nu = 4, sigma = sigma,
               prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)),
  comp3 = list(m = c(3, 3), kappa = 2, nu = 4, sigma = sigma,
               prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)))

# x_mean = colMeans(xi)
# crossprod(xi-x_mean, xi-x_mean) # equals t(x) %*% x

# update(components, k = 1, xi, 1)
# components[[1]]
*/



// [[Rcpp::export]]
List reinit_component(List components, int k, const arma::mat x, int debug) {

  // extract kth component
  List ll = components[k-1];
  ll = ll["prior"];
  ll["prior"] = ll;
  arma::rowvec prior_mean   = as<rowvec>(ll["m"]);
  double prior_kappa        = ll["kappa"];
  double prior_nu           = ll["nu"];
  arma::mat prior_sigma     = as<mat>(ll["sigma"]);

  // data stats
  double n = x.n_rows;
  double d = x.n_cols;

  // compute colmeans
  arma::rowvec x_mean(d);
  for (int j=0; j<d; j++)
    for (int i=0; i<n; i++)
      x_mean[j] += x(i, j) / n;

  // centered data crossprod / empirical covariance matrix
  arma::mat xc =
    trans(x - repmat(x_mean, n, 1)) * (x - repmat(x_mean, n, 1));

  if (debug == 1) {
    print(wrap(x_mean));
    print(wrap(prior_mean));
    print(wrap("xc:"));
    print(wrap(xc));
    print(wrap("x_mean - prior_mean * ..."));
    print(wrap(trans(x_mean - prior_mean) * (x_mean - prior_mean)));
    print(wrap("fac * x_mean - prior_mean * ..."));
    print(wrap(prior_kappa / (prior_kappa + n) *
      trans(x_mean - prior_mean) * (x_mean - prior_mean)));
  }

  // update stats
  ll["m"] =
    prior_kappa /  (prior_kappa + n) * prior_mean +
    n / (prior_kappa + n) * x_mean;
  ll["kappa"] = prior_kappa + n;
  ll["nu"]    = prior_nu    + n;
  ll["sigma"] =
    prior_sigma + xc +
    prior_kappa / (prior_kappa + n) *
    trans(x_mean - prior_mean) * (x_mean - prior_mean);

  if (debug == 1) {
    print(wrap("After updating:"));
    print(wrap(ll));
  }

  // return new prior updated with data as new component
  return ll;
}

/*** R
set.seed(1)
x = matrix(rnorm(100), nrow=50)
xi = x[1:10,,drop=F]
sigma = matrix(c(1,0,0,1),2,2)
components = list(
  comp1 = list(m = c(1, 1), kappa = 2, nu = 4, sigma = sigma,
               prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)),
  comp2 = list(m = c(2, 2), kappa = 2, nu = 4, sigma = sigma,
               prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)),
  comp3 = list(m = c(3, 3), kappa = 2, nu = 4, sigma = sigma,
               prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)))

print(reinit_component(components, k = 1, xi, 1))
*/


// [[Rcpp::export]]
List sample_new_component(List prior, const arma::mat x, Function samplefun2, int debug) {

  // extract parameters 
  arma::rowvec prior_mean   = as<rowvec>(prior["m"]);
  double prior_kappa        = prior["kappa"];
  double prior_nu           = prior["nu"];
  arma::mat prior_sigma     = as<mat>(prior["sigma"]);

  // call sampling function from R
  List ll = samplefun2(prior_mean, prior_kappa, prior_nu, prior_sigma);

  // data stats
  double n = x.n_rows;
  double d = x.n_cols;

  // compute colmeans
  arma::rowvec x_mean(d);
  for (int j=0; j<d; j++)
    for (int i=0; i<n; i++)
      x_mean[j] += x(i, j) / n;

  // centered data crossprod / empirical covariance matrix
  arma::mat xc =
    trans(x - repmat(x_mean, n, 1)) * (x - repmat(x_mean, n, 1));

  if (debug == 1) {
    print(wrap(x_mean));
    print(wrap(prior_mean));
    print(wrap("xc:"));
    print(wrap(xc));
    print(wrap("x_mean - prior_mean * ..."));
    print(wrap(trans(x_mean - prior_mean) * (x_mean - prior_mean)));
    print(wrap("fac * x_mean - prior_mean * ..."));
    print(wrap(prior_kappa / (prior_kappa + n) *
      trans(x_mean - prior_mean) * (x_mean - prior_mean)));
  }

  // update stats
  ll["m"] =
    prior_kappa /  (prior_kappa + n) * prior_mean +
    n / (prior_kappa + n) * x_mean;
  ll["kappa"] = prior_kappa + n;
  ll["nu"]    = prior_nu    + n;
  ll["sigma"] =
    prior_sigma + xc +
    prior_kappa / (prior_kappa + n) *
    trans(x_mean - prior_mean) * (x_mean - prior_mean);

  if (debug == 1) {
    print(wrap("After updating:"));
    print(wrap(ll));
  }

  // return new prior updated with data as new component
  return ll;
}

/*** R
x = matrix(rnorm(100), nrow=50)
B = Matrix::bandSparse(nrow(x), k = c(-1, 1, -1)); I = diag(nrow(x))
P = kronecker(B, I) + kronecker(I, B)
alpha = 1; lambda = 1
xi = x[1,,drop=F]
i = 1
z = rep(1:2, each = 50)
sigma = matrix(c(1,0,0,1),2,2)
components = list(
  comp1 = list(m = c(1, 1), kappa = 2, nu = 4, sigma = sigma,
               prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)),
  comp2 = list(m = c(2, 2), kappa = 2, nu = 4, sigma = sigma,
               prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)),
  comp3 = list(m = c(3, 3), kappa = 2, nu = 4, sigma = sigma,
               prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)))
prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)
z = rep(0, 100)
z[1] = 1
z[2:4] = 2
z[5] = 3
inds = 1:5
xi = x[inds,,drop=F]
samplefun2 = function(prior_mean, prior_kappa, prior_nu, prior_sigma) {
        S0 = MCMCpack::rwish(prior_nu, prior_sigma)
	m0 = mvtnorm::rmvnorm(1, prior_mean, 1 / prior_kappa * S0)
        list(m = m0, kappa = prior_kappa, nu = prior_nu, sigma = S0)
}
sample_new_component(prior, xi, samplefun2, 1)
*/

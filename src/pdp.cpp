// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <one_pass.cpp>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List two_pass() {

  List ret;
  rowvec m = Rcpp::rnorm(2);
  mat sigma(2,2);
  std::fill(sigma.begin(), sigma.end(), R::rnorm(1, 1));
  ret["m"] = m;
  ret["sigma"] = sigma;
  return ret;
}

// [[Rcpp::export]]
List pdp_test(int nrep, NumericVector z, mat y) {
  int d = y.n_cols;
  NumericMatrix z_path(nrep, z.length());
  NumericMatrix m_path(nrep, d);
  NumericMatrix sigma_path(nrep, d*d);

  List comp_i = two_pass();

  for (int rep = 0; rep < nrep; rep++) {
    comp_i = two_pass();
    z_path(rep,_) = z;
    m_path(rep,_) = as<NumericVector>(comp_i["m"]);
    sigma_path(rep,_) = as<NumericVector>(comp_i["sigma"]);
  }

  List sampling_paths;
  sampling_paths["z"] = z_path;
  sampling_paths["m"] = m_path;
  sampling_paths["sigma"] = sigma_path;
  return sampling_paths;
}


/*** R
# z = rep(0, 100)
# nrep = 2
# d = 2
# y = mvtnorm::rmvnorm(10, c(2, 2))
#
# sampling_paths = list(z = matrix(NA, nrep, length(z)),
#                       m = matrix(NA, nrep, d),
#                       sigma = matrix(NA, nrep, d^2))
# two_pass()
# ret = pdp_test(nrep, z, y)
# str(ret$m)
# str(ret$sigma)
*/


// =====================================================
// Now the real deal
// =====================================================

// [[Rcpp::export]]
List pdp(int nrep, IntegerVector z, mat y,
         List components, List prior,
         NumericMatrix P, double alpha, double lambda,
         int debug, Function samplefun, Function samplefun2) {

  int n = z.length();
  int d = y.n_cols;
  mat z_path(nrep, z.length());
  NumericMatrix m_path(nrep, d);
  NumericMatrix sigma_path(nrep, d*d);
  List comp_i;

  if (sum(z == 0) == n) {
    print(wrap("init with first pass"));
    z[0] = 1;
    arma::mat ysub = y.rows(1, 1);
    // print(wrap(ysub));
    // update(components, 1, ysub, 0);
    comp_i = initialize_components(z[z != 0], y, prior, samplefun2);
    // first pass
    IntegerVector ind0 = Rcpp::seq_len(n-2) + 1;
    print(wrap(comp_i));
    comp_i = one_pass_cpp(ind0, z, y, comp_i, prior,
                          P, alpha, lambda, 0,
                          samplefun, samplefun2);
    print(wrap("first pass done"));
  } else if (sum(z == 0) == 0) {
    print(wrap("initialize components according to z"));
    comp_i = initialize_components(z, y, prior, samplefun2);
    print(wrap("init done"));
  }

  IntegerVector ind = Rcpp::seq_len(n);
  for (int rep = 0; rep < nrep; rep++) {
    std::random_shuffle(ind.begin(), ind.end(), randWrapper);
    comp_i = one_pass_cpp(ind, z, y, comp_i, prior,
                          P, alpha, lambda, debug,
                          samplefun, samplefun2);
    // print(wrap("almost done iter:"));
    for (int i=0; i<n; i++)
      z_path(rep, i) = z[i];
    // m_path.row(rep) = as<NumericVector>(comp_i["m"]);
    // sigma_path(rep+1,_) = as<NumericVector>(comp_i["sigma"]);
    // print(wrap("done iter:"));
    // print(wrap(rep+1));
  }

  List sampling_paths;
  sampling_paths["z"] = z_path;
  // sampling_paths["m"] = m_path;
  // sampling_paths["sigma"] = sigma_path;
  return sampling_paths;
}


/*** R
# set.seed(1)
# d = 2
# sigma = matrix(c(1,0,0,1),2,2)
# x = rbind(
#     mvtnorm::rmvnorm(50, c(-2, -2), sigma),
#     mvtnorm::rmvnorm(50, c(2, 2), sigma))
# B = Matrix::bandSparse(10, k = c(-1, 1, -1)); I = diag(10)
# P = as.matrix(kronecker(B, I) + kronecker(I, B))
# alpha = 1; lambda = 1
#
# ## new init
# prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma,
#              prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma))
# components = list(comp1 = list(m = c(1, 1), kappa = 2, nu = 4, sigma = sigma,
#                                prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)))
# z = rep(0, 100)
# source("~/projects/Masterthesis/code/PDP/R/samplefuns.R")
# nrep = 2
# ret = pdp(nrep, z, x, components, prior,
#           P, alpha, lambda, 0, samplefun, samplefun2)
# ret
*/


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <update_chol.cpp>
#include <projects/Masterthesis/code/PDP/src/log_pred.cpp>

using namespace Rcpp;

struct nw_component {
  // current values
  arma::rowvec m;
  double nu;
  double kappa;
  arma::mat S;
  arma::mat chol_S;
  
  // prior
  arma::rowvec prior_m;
  double prior_nu;
  double prior_kappa;
  arma::mat prior_S;
   
  void initialize(arma::rowvec mean, double n_mean, double n_sigma,
                  arma::mat sigma) {
    m = mean;
    kappa = n_mean;
    nu = n_sigma;
    S = sigma;
    chol_S = arma::chol(sigma);
    
    prior_m = mean;
    prior_kappa = n_mean;
    prior_nu = n_sigma;
    prior_S = sigma;
  };
  void initialize(List ll) {
    m =       as<arma::rowvec>(ll[1]);
    kappa =   ll[2];
    nu =      ll[3];
    S =       as<arma::mat>(ll[4]);
    chol_S = arma::chol(S);
    
    prior_m = as<arma::rowvec>(ll[1]);
    prior_kappa = ll[2];
    prior_nu = ll[3];
    prior_S = as<arma::mat>(ll[4]);
  };
  void init_sample(List prior, Function samplefun2) {
    List ll = samplefun2(prior);
    initialize(ll["m"], ll["kappa"], ll["nu"], ll["sigma"]);
  };
  void update(const arma::mat& x, bool downdate, int debug) {
    // data stats
    double n = x.n_rows;
    int d = x.n_cols;
    
    // compute colmeans
    arma::rowvec x_mean = arma::mean(x, 0); 
    
    // centered data crossprod / empirical covariance matrix
    arma::mat xc =
      trans(x - repmat(x_mean, n, 1)) * (x - repmat(x_mean, n, 1));
    
    if (debug == 1) {
      print(wrap(x_mean));
      print(wrap(m));
      print(wrap("xc:"));
      print(wrap(xc));
      print(wrap("x_mean - prior_mean * ..."));
      print(wrap(trans(x_mean - m) * (x_mean - m)));
      print(wrap("fac * x_mean - prior_mean * ..."));
      print(wrap(kappa / (kappa + n) *
        trans(x_mean - m) * (x_mean - m)));
    }
    
    // update stats
    if (downdate) {
      // needs this order
      kappa = kappa - n;
      nu    = nu    - n;
      m =  
        (kappa + n) / kappa * m - 
        (n * x_mean / kappa);
      double fac1 = kappa / (kappa + n);
      arma::rowvec xi_m = sqrt(fac1) * x_mean - m;
      S = S - xc - xi_m.t() * xi_m;
      // update chol
      update_chol(chol_S, xi_m, downdate);
    } else {
      double fac1 = kappa / (kappa + n);
      arma::rowvec xi_m = sqrt(fac1) * x_mean - m;
      S = S + xc + xi_m.t() * xi_m;
      m = (kappa * m + n * x_mean) / (kappa + n);
      kappa = kappa + n;
      nu    = nu    + n;
      // update chol
      update_chol(chol_S, xi_m, downdate);
    }
  };
  double log_pred(const arma::mat& x) const {
    double n = x.n_rows;
    double d = x.n_cols;
    return  sum(dmvt_chol(x, m,
                          (kappa + n) / (kappa * (nu - d + n)) * chol_S,
                          nu - d + n));
  };
  double log_pred2(const arma::mat& x) const {
    double n = x.n_rows;
    double d = x.n_cols;
    return  ::log_pred(x, m, kappa, nu, S);
  };
  void print_params() const {
    print(wrap(m));
    print(wrap(nu));
    print(wrap(kappa));
    print(wrap(S));
    print(wrap("chol and chol.t * chol:"));
    print(wrap(chol_S));
    print(wrap(chol_S.t() * chol_S));
  };
  void print_mean() const {
    print(wrap(m));
  };
  void print_prior_params() const {
    print(wrap("prior params:"));
    print(wrap(prior_m));
    print(wrap(prior_nu));
    print(wrap(prior_kappa));
    print(wrap(prior_S));
  };
  List ret_params() const{
    List rval; 
    rval["m"] = m;
    rval["kappa"] = kappa;
    rval["nu"] = nu;
    rval["S"] = S;
    return rval;
  };
  
};

// [[Rcpp::export]]
void test(List& ll, const arma::mat& x, Function samplefun2) {
  
  std::vector<nw_component> mixture(1);
  nw_component c1;
  c1.initialize(ll["m"], ll["kappa"], ll["nu"], ll["S"]); 
  mixture[0] = c1;
  
  nw_component c2;
  c2.init_sample(ll, samplefun2);
  mixture.push_back(c2);
  
  int K = mixture.size();
  for (int j = 0; j < K; j++) 
    mixture[j].print_mean();
  
  // replace first element with last and pop
  mixture[0] = mixture[K - 1];
  mixture.pop_back();
  K--;
  
  for (int j = 0; j < K; j++) 
    mixture[j].print_mean();
  
  for (int j = 0; j < 4; j++) {
    nw_component nw_new;
    nw_new.init_sample(ll, samplefun2);
    mixture.push_back(nw_new);
    K += 1;
  }
  print(wrap(mixture.size()));
  print(wrap(K));
  
  for (int j = 0; j < K; j++) 
    mixture[j].print_mean();
  return ;
}

// [[Rcpp::export]]
NumericVector test_logpred(List& ll, arma::mat& x) {
  nw_component c1;
  c1.initialize(ll[0], ll[1], ll[2], ll[3]);
  // c1.update(x.row(1), false,  0); // update
  int n = x.n_rows;
  NumericVector res(n);
  for (int i=0; i<n; i++) {
    res[i] = c1.log_pred(x.row(i));
  }
  return res;
}

// [[Rcpp::export]]
NumericVector test_logpred2(List& ll, arma::mat& x) {
  nw_component c1;
  c1.initialize(ll[0], ll[1], ll[2], ll[3]);
  // c1.update(x.row(1), false,  0); // update
  int n = x.n_rows;
  NumericVector res(n);
  for (int i=0; i<n; i++) { 
    res[i] = c1.log_pred2(x.row(i));
  }
  return res;
}


/*** R
source("~/projects/Masterthesis/code/PDP/R/samplefuns.R")
x = mvtnorm::rmvnorm(100, mean = c(2, 2), sigma = matrix(c(1, .3, .3, 1), 2, 2))
prior = list(m = c(0,0), kappa = 2, nu = 4, S = matrix(c(1, 0, 0, 1), 2, 2))
xi = x[1, , drop=FALSE]
test(prior, xi, samplefun2)

sapply(1:10, function(i)
  log_pred(x[i,,drop=F], prior$m, prior$kappa, prior$nu, prior$S) )
test_logpred2(prior, x[1:10, ,drop=F])

dmvt_chol(x[1,,drop=F], prior$m, prior$S, prior$nu)
test_logpred(prior, x[1:10, ,drop=F])
# xi = x[1:10, , drop=FALSE]
# test(prior, xi)

*/

/*
 * Testing 
 * 
 */

// [[Rcpp::export]]
List test_update(List& ll, arma::mat& x) {
  nw_component c1;
  c1.initialize(ll[0], ll[1], ll[2], ll[3]);
  c1.update(x, false,  0); // update
  List rval = c1.ret_params();
  return rval;
}

// [[Rcpp::export]]
List test_up_down_date(List& ll, arma::mat& x) {
  nw_component c1;
  c1.initialize(ll[0], ll[1], ll[2], ll[3]);
  c1.update(x, false,  0); // update
  c1.update(x, true,  0);  // downdate again
  List rval = c1.ret_params();
  return rval;
}

/*** R
testthat::test_that("updating works", {
  x = mvtnorm::rmvnorm(100, mean = c(2, 2), 
                       sigma = matrix(c(1, .3, .3, 1), 2, 2))
  prior = list(mean = matrix(c(0,0), 1), kappa = 2, nu = 4, S = matrix(c(1, 0, 0, 1), 2, 2))
  xi = x[1, , drop = F]; n = nrow(xi) 
  res = test_update(prior, xi)
  m2 = (prior$kappa*prior$m + n * colMeans(xi)) / (prior$kappa + n)
  xi_m = sqrt(prior$kappa / (prior$kappa + n)) * colMeans(xi) - prior$m
  S2 = prior$S + crossprod(xi_m, xi_m)
  testthat::expect_equal(res$m, m2)
  testthat::expect_equal(res$S, S2)
  
  ## multiple obs
  xi = x[1:10, , drop = F]; n = nrow(xi) 
  res = test_update(prior, xi)
  m2 = (prior$kappa*prior$m + n * colMeans(xi)) / (prior$kappa + n)
  xi_m = sqrt(prior$kappa / (prior$kappa + n)) * colMeans(xi) - prior$m
  S2 = prior$S + crossprod(xi_m, xi_m)
  testthat::expect_equal(res$m, m2)
  testthat::expect_equal(res$S, S2)
  
  res = test_up_down_date(prior, x[1:10, , drop = FALSE])
  testthat::expect_equal(res$m, prior$m)
  testthat::expect_equal(res$S, prior$S)
})
*/
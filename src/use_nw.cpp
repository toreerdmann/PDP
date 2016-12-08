// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <nw_component.cpp>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
void bla (IntegerVector z, mat x, List prior) {
  
  int n = x.n_rows();
  std::vector<nw_component> components;
  
  int K = Rcpp::unique(z).length();
    nw_compnent nw_new;
    nw_new.init_sample(prior, samplefun2);
    
  List ll(K);
  arma::colvec z2 = as<arma::colvec>(z);
  for (int j=0; j<K; j++) {
    ll[j] = sample_new_component(prior);
    mat ysub = y.rows(find(z2 == (j+1)));
    update(ll, j+1, ysub, 0);  // use k with indexing from 1
  }
    
  // for (int i=0; i<n; i++) {
  // }
}



/*** R
set.seed(1)
sigma = matrix(c(1,0,0,1),2,2)
x = rbind(
    mvtnorm::rmvnorm(50, c(-2, -2), sigma),
    mvtnorm::rmvnorm(50, c(2, 2), sigma))
B = Matrix::bandSparse(10, k = c(-1, 1, -1)); I = diag(10)
P = kronecker(B, I) + kronecker(I, B)
alpha = 1; lambda = 1
plot(x)

## new init
prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma,
             prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma))
z = rep(0, 100)

bla(z, x)

*/

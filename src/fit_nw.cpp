/*
 * Here, we fit a mixture of normals, given a partition of the data.
 */

// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <nw_component.cpp>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
void fit_nw(std::vector<nw_component> components, colvec & z2, arma::mat & x, List & prior) {
  int n = x.n_rows;
  // int K = components.size();
  // for (int j=0; j<K; j++) {
  //     nw_component nw_new;
  //     nw_new.init_sample(prior);
  //     print(wrap(find(z2 == j+1)));
  //     mat xsub = x.rows(find(z2 == j+1));
  //     for (int i=0; i<xsub.n_rows; i++)
  //       nw_new.update(xsub.row(i), false, 0);
  //     components[j] = nw_new;
  // }
}

// [[Rcpp::export]]
List test(mat & x, IntegerVector & z, int K, List & prior) {
  // std::vector<nw_component> components(K);
  // colvec z2 = as<colvec>(z);
  std::cout << "bla";
  // fit_nw(components, z2, x, prior);
  List ll;
  // ll["c1"] = components[0].ret_params();
  // ll["c2"] = components[1].ret_params();
  return ll;
}
/*** R
# set.seed(1)
sigma = matrix(c(1,0,0,1),2,2)
x = rbind(
    mvtnorm::rmvnorm(50, c(-2, -2), sigma),
    mvtnorm::rmvnorm(50, c(2, 2), sigma))
B = Matrix::bandSparse(10, k = c(-1, 1, -1)); I = diag(10)
P = kronecker(B, I) + kronecker(I, B)
alpha = 1; lambda = 1

## new init
prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma, chol_S = sigma,
             prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma))
z = rep(1:2, each = 50)

res = test(z, x, 2, prior)
# xseq = seq(-10,10, len = 10)
# xgrid = expand.grid(xseq, xseq)
# plot(x, xlim = c(-10, 10), ylim = c(-10, 10), col = rep(1:2, each = 50))
# contour(xseq, xseq, matrix(mvtnorm::dmvnorm(xgrid, res[[1]]$m, res[[1]]$S, log = FALSE), 10, 10), add = TRUE)
# contour(xseq, xseq, matrix(mvtnorm::dmvnorm(xgrid, res[[2]]$m, res[[2]]$S, log = FALSE), 10, 10), col = 2, add = TRUE)
        

*/
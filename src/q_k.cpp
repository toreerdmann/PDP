// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <nw_component.cpp>
#include <graph.cpp>
using namespace Rcpp;
using namespace arma;

/*
* Note to self: you cannot use a std::vector<myType> in an exported function.
* Below is an example of how to deal with this.
*/

void print_k(int k, std::vector<nw_component>  components) {
  components.at(k).print_mean();
}

// [[Rcpp::export]]
void use_print_k(int k, List prior) {
  nw_component nw1;
  nw1.init_sample(prior);
  std::vector<nw_component> components(3);
  components.at(0) = nw1;
  print_k(k-1, components);
}
/*** r
prior = list(m = c(0, 0), nu = 2, kappa = 2, sigma = diag(2), chol_S = diag(2))
use_print_k(1, prior)
*/

/* 
* Compute the assignment probability for class k
*/

double qi_k1(int i, int k, mat &  x, IntegerVector & z, Graph & g, std::vector<nw_component>  components) {
  
  int n = x.n_rows;
  int d = x.n_cols;
  
  return 
    // nj = number of obs in class k
    // log(sum(z == (k+1))) + 
    sum(z == (k+1)) +
      // nc_sum = number of neighbors of class k
      sum(as<IntegerVector>(z[g.getEdges(i)]) == (k+1));
      // nc_sum * lambda; +
      // log posterior predictive
      // sum(dmvt_chol(x,
      //               components.at(k).m,
      //               (components.at(k).kappa + n) /
      //                 (components.at(k).kappa * (components.at(k).nu - d + n)) * components.at(k).chol_S,
      //                 components.at(k).nu - d + n));
}

// [[Rcpp::export]]
double test_qi_k1(int K, int i, int k, mat & x, IntegerVector & z, List & adjlist, List & prior) {
  
  int n = x.n_rows;
  
  // init neighborhood
  Graph g(n);
  g.init_w_adjlist(adjlist, n);
  
  // init mixture
  std::vector<nw_component> components(K);
  for (int j=0; j<K; j++) {
    nw_component nw1;
    nw1.init_sample(prior);
    components.at(j) = nw1;
  }
  
  // print assignment prob
  return qi_k1(i-1, k-1, x, z, g, components);
}

/*** R
testthat::test_that("computing nj and nc in qi_k works", {
  set.seed(2)
  x = mvtnorm::rmvnorm(12, c(0, 0), diag(2)); z = sample(1:3, 12, T)
  prior = list(m = c(0, 0), nu = 2, kappa = 2, sigma = diag(2), chol_S = diag(2))
  B = Matrix::bandSparse(3, k = c(-1, 1, -1)); I = diag(4); P = kronecker(B, I) + kronecker(I, B)
  adjlist = sapply(1:nrow(P), function(i) which(P[i,] == 1))
  matrix(1:12, 4, 3) ## labelmat
  matrix(z, 4, 3)    ## z
  expect_equal(test_qi_k1(3, i = 1, k = 1, x, z, adjlist, prior), sum(z == 1))
  expect_equal(test_qi_k1(3, i = 1, k = 2, x, z, adjlist, prior), sum(z == 2))
  expect_equal(test_qi_k1(3, i = 1, k = 3, x, z, adjlist, prior), sum(z == 3) + sum(z[c(2,5)] == 3))
})

*/

arma::vec qi_k(int i, arma::mat &  x, IntegerVector & z, double lambda, Graph & g, std::vector<nw_component>  components) {
  
  int n = x.n_rows;
  int d = x.n_cols;
  int K = components.size();
  
  arma::vec qout(K);
  for (int k=0; k<K; k++)
    qout[k] = 
      // nj = number of obs in class k
      log(sum(z == (k+1))) +
      // nc_sum = number of neighbors of class k
      sum(as<IntegerVector>(z[g.getEdges(i)]) == (k+1)) * lambda +
      // log posterior predictive
      sum(dmvt_chol(x,
                    components.at(k).m,
                    (components.at(k).kappa + n) /
                    (components.at(k).kappa * (components.at(k).nu - d + n)) * components.at(k).chol_S,
                    components.at(k).nu - d + n));
  return qout;
}

// [[Rcpp::export]]
vec test_qi_k2(int K, int i, mat & x, IntegerVector & z, 
                  double lambda, List & adjlist, List & prior) {
  
  int n = x.n_rows;
  
  // init neighborhood
  Graph g(n);
  g.init_w_adjlist(adjlist, n);
  
  // init mixture
  std::vector<nw_component> components(K);
  for (int j=0; j<K; j++) {
    nw_component nw1;
    nw1.init_sample(prior);
    components.at(j) = nw1;
  }
  
  // return assignment prob for k=1,...,K
  return qi_k(i, x, z, lambda, g, components);
}
  
  
/*** R
testthat::test_that("computing nj and nc in qi_k works", {
  set.seed(2)
  x = mvtnorm::rmvnorm(12, c(0, 0), diag(2)); z = sample(1:3, 12, T)
  prior = list(m = c(0, 0), nu = 2, kappa = 2, sigma = diag(2), chol_S = diag(2))
  B = Matrix::bandSparse(3, k = c(-1, 1, -1)); I = diag(4); P = kronecker(B, I) + kronecker(I, B)
  adjlist = sapply(1:nrow(P), function(i) which(P[i,] == 1))
  matrix(1:12, 4, 3) ## labelmat
  matrix(z, 4, 3)    ## z
  test_qi_k2(3, i = 1, x, z, 2, adjlist, prior)
  
  # expect_equal(test_qi_k(3, i = 1, k = 1, x, z, adjlist, prior), sum(z == 1))
  # expect_equal(test_qi_k(3, i = 1, k = 2, x, z, adjlist, prior), sum(z == 2))
  # expect_equal(test_qi_k(3, i = 1, k = 3, x, z, adjlist, prior), sum(z == 3) + sum(z[c(2,5)] == 3))
})

*/
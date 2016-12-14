// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <nw_component.cpp>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
void bla (IntegerVector & z, const mat & x, List prior) {
  
  int n = x.n_rows;
  
  IntegerVector classes = unique(z).sort(false);
  int K = classes.length();
  std::vector<nw_component> components(K);
  print(wrap(classes));
  for (int j=0; j<K; j++) {
    if (classes[j] != 0) {
      print(wrap(classes[j]));
      nw_component nw_new;
      nw_new.init_sample(prior);
      arma::colvec z2 = as<arma::colvec>(z);
      print(wrap(find(z2 == classes[j])));
      mat xsub = x.rows(find(z2 == classes[j]));
      nw_new.update(xsub, false, 0);
      components[j] = nw_new;
    }
  }
  print(wrap(components.size()));
  std::cout << "init done.";
}

/* 
 * Function for fitting a given partition, setting 
 * the parameters of a vector of components.
 */
void fit_partition(std::vector<nw_component> & components, 
                   colvec & z2, mat & x, List & prior) {
  int K = components.size();
  for (int j=0; j<K; j++) {
      nw_component nw_new;
      nw_new.init_sample(prior);
      mat xsub = x.rows(find(z2 == j+1));
      // std::cout << "Fit data to comp" << j+1 << std::endl;
      print(wrap(find(z2 == j+1) + 1));
      for (int i=0; i<xsub.n_rows; i++)
        nw_new.update(xsub.row(i), false, 0);
      components[j] = nw_new;
  }
}

// [[Rcpp::export]]
List fit_nw(IntegerVector & partition, arma::mat & x, List & prior) {
  
  // get classes
  IntegerVector z = partition[partition != 0];
  IntegerVector classes = unique(z).sort(false);
  int K = classes.length();
  std::cout << "K:" << K <<  std::endl;
  colvec z2 = as<colvec>(partition);
  std::vector<nw_component> components(K);
  fit_partition(components, z2, x, prior); // assign data to components
  
  List ll;
  ll["c1"] = components[0].ret_params();
  ll["c2"] = components[1].ret_params();
  return ll;
}


/*** R
set.seed(2)
sigma = matrix(c(1,0,0,1),2,2) * 1/5
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
# bla(z, x, prior)
z[1] = 0
res = fit_nw(z, x, prior)
xseq = seq(-10,10, len = 10)
xgrid = expand.grid(xseq, xseq)
plot(x, xlim = c(-10, 10), ylim = c(-10, 10), col = rep(1:2, each = 50))
contour(xseq, xseq, matrix(mvtnorm::dmvnorm(xgrid, res[[1]]$m, res[[1]]$S, log = FALSE), 10, 10), add = TRUE)
contour(xseq, xseq, matrix(mvtnorm::dmvnorm(xgrid, res[[2]]$m, res[[2]]$S, log = FALSE), 10, 10), col = 2, add = TRUE)
        
*/

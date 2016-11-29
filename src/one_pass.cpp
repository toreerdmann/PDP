// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <qi.cpp>
#include <update_list.cpp>
#include <logsumexp.cpp>
#include <resize.cpp>

using namespace Rcpp;


// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }


// [[Rcpp::export]]
List one_pass_cpp (IntegerVector ind, IntegerVector z,
                   arma::mat y, List components, List prior,
                   NumericMatrix P, double alpha, double lambda,
                   int debug, Function samplefun, Function samplefun2) {

  int n = ind.length();
  // int d = y.n_cols;
  int K = components.length();

  // // create vector with shuffled indeces
  IntegerVector s = Rcpp::seq_len(n);
  // NumericVector ind (s.begin(), s.end());
  // std::random_shuffle(ind.begin(), ind.end(), randWrapper);
  // Rcpp::print(wrap(ind));
  int i;

  for (int i1=0; i1 < n; i1++) {

    i = ind[i1] - 1;

    if (debug == 1) {
      Rcpp::print(wrap("================================="));
      Rcpp::print(wrap("=========== new iter ============"));
      Rcpp::print(wrap("================================="));
      Rcpp::print(wrap("i ="));
      Rcpp::print(wrap(i));
      Rcpp::print(wrap("z[i] ="));
      Rcpp::print(wrap(z[i]));
    }

    // use sub labels
    int old_class = z[i];
    // set current labels to zero
    z[i] = 0;

    // cache the old_component, so it can possibly be restored
    // init here, so that it has global scope?
    List old_component;

    if (debug == 1) {
      Rcpp::print(wrap("z after setting z[i] = 0: "));
      Rcpp::print(wrap(z));
    }

    // remove data from current spin cluster from its component
    if (old_class != 0) {

      if (debug == 1)
        Rcpp::print(wrap("old_class != 0"));

      if (sum(z == old_class) == 0) {
        // if only observation in that cluster, remove the component

        if (debug == 1)
          Rcpp::print(wrap("sum(z == old_class) == 0"));

        // cache the old_component, so it can possibly be restored
        List old_component = (components[old_class-1]);

        if (debug == 1) {
          Rcpp::print(wrap("only obs in cluster: remove component (old_class = )"));
          Rcpp::print(wrap(old_class));
          // CharacterVector nam = components.names();
          Rcpp::print(wrap("components before: "));
          for(int j=0;j<K;j++) {
            List compj = components[j];
            rowvec m0 = compj["m"];
            print(wrap(m0));
          }
          // Rcpp::print(wrap(nam));
          Rcpp::print(wrap("Move Kth component"));
          // Rcpp::print(wrap(nam[K-1]));
          Rcpp::print(wrap("to spot of old_class - 1 "));
          // Rcpp::print(wrap(nam[old_class-1]));
        }

        // move last component forward to empty spot
        components[old_class-1] = components[K-1];
        // Rcpp::print(wrap(z[z == K]));
        z[z == K] = old_class;


        if (debug == 1) {
          Rcpp::print(wrap("moved indeces from last to old_class"));
          Rcpp::print(wrap("new indeces:"));
          Rcpp::print(wrap(z));
        }

        // remove last component from list
        IntegerVector idx_int = seq_len(components.length());
        components = components[idx_int != K];
        K = K - 1;

        if (debug == 1) {
          Rcpp::print(wrap("component list after removal:"));
          Rcpp::print(wrap(components.names()));
          Rcpp::print(wrap(components));
          Rcpp::print(wrap("K after removal:"));
          Rcpp::print(wrap(K));
        }
      } else {
        // remove x[i] from old component
        // by new initialization

        if (debug == 1)
          print(wrap("remove x[i] from component[old_class-1]"));

        // update with all data except for current
        arma::colvec z2 = as<arma::colvec>(z);
        mat ysub = y.rows(find(z2 == old_class));

        // old:
        // List new_comp = clone(prior);
        // components[old_class - 1] = new_comp;
        // even older: changed this to above
        // List new_comp = reinit_component(components, old_class-1, ysub, 0);
        // components[old_class - 1] = new_comp;

        // new version
        set_to_prior(components, old_class); // use k=old_class
        update(components, old_class, ysub, debug);  // use k=old_class

        if (debug == 1) {
          Rcpp::print(wrap("set parameters of this component to prior values:"));
          // Rcpp::print(wrap(components));
          Rcpp::print(wrap("used this data to update component"));
          Rcpp::print(wrap(ysub));
        }
      }
    }

    if (debug == 1)
      Rcpp::print(wrap("compute q"));
    // compute assignment probabilities
    // arma::colvec z2 = as<arma::colvec>(z);
    // arma::mat ysub = y.rows(arma::find(z2 == 0));
    arma::mat ysub = y.rows(i, i);
    // Rcpp::print(wrap("ysub:"));
    // Rcpp::print(wrap(ysub));
    arma::vec q = get_assignment_prob(ysub, z, components, prior,
                                     P(i,_), lambda, alpha, debug);
    // normalize
    q = exp(q - logsumexp(q));

    // draw class
    // double ret = R::runif(0, 1);
    z[i] = as<int>(samplefun(K, q));
    // z[i] = as<int>(R::sample.int(K+1, 1, prob = q) - 1);
    // z[i] = 0;

    if (debug == 1) {
      print(wrap("compute assignment prob using:"));
      print(wrap(ysub));
      print(wrap("q:"));
      print(wrap(q));
      print(wrap("z[i] drawn:"));
      print(wrap(z[i]));
    }

    // if equal to 0, create new cluster
    if (z[i] == 0) {
      if (sum(z == old_class) == 0) {
        // just use the old singleton cluster if it had one

        // alternative: (causing segfault?)
        // components.push_back(components[K-1]);
        // components[old_class] = old_component;

        // components.push_back(old_component);
        components = resize(components, K+1);
        components[K] = old_component;
        K++;
        z[i] = K;
        // z[i] = old_class;
      } else {
        // else, make a new one and update with the
        // single observation

        // sample new prior values
        List new_comp = sample_new_component(prior, samplefun2);
        // add
        components = resize(components, K+1);
        components[K] = new_comp;
        // components.push_back(new_comp);
        K++;
        z[i] = K;
      }
      if (debug == 1) {
        print(wrap("push_back"));
        print(wrap(components));
        print(wrap("K:"));
        print(wrap(K));
        print(wrap("old_class:"));
        print(wrap(old_class));
        print(wrap("z:"));
        print(wrap(z));
      }
    }

    // update stats of cluster z[i]
    // arma::colvec z2 = as<arma::colvec>(z);
    // arma::mat ysub = y.rows(arma::find(z2 == old_class) - 1);
    if (debug == 1) {
      print(wrap("update component[z[i]]:"));
      print(wrap(components));
    }
    update(components, z[i], ysub, 0);
  }
  // list(z = z, components = components)
  return components;
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
# xi = x[1,,drop=F]
# i = 1
# z = rep(1:2, each = 50)
# components = list(
#   comp1 = list(m = c(1, 1), kappa = 2, nu = 4, sigma = sigma,
#                prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)),
#   comp2 = list(m = c(2, 2), kappa = 2, nu = 4, sigma = sigma,
#                prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)),
#   comp3 = list(m = c(3, 3), kappa = 2, nu = 4, sigma = sigma,
#                prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)))
# prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma,
#              prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma))
# z = rep(0, 100)
# z[1] = 1
# z[2:4] = 2
# z[5] = 3
# inds = 1:5
# xi = x[inds,,drop=F]


## new init
prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma,
             prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma))
components = list( comp1 = list(m = c(1, 1), kappa = 2, nu = 4, sigma = sigma,
                                prior = list(m = c(0, 0), kappa = 2, nu = 4, sigma = sigma)))
z = rep(0, 100)
z[1] = 1
# update(components, 1, x[1,,drop=F], 1)
samplefun = function(K, q){
  sample.int(K+1, 1, prob = q) - 1
}
samplefun2 = function(prior) {
  prior_mean = prior[[1]]
  prior_kappa = prior[[2]]
  prior_nu = prior[[3]]
  prior_sigma = prior[[4]]
  S0 = MCMCpack::rwish(prior_nu, prior_sigma)
  m0 = mvtnorm::rmvnorm(1, prior_mean, 1 / prior_kappa * S0)
  list(m = m0, kappa = prior_kappa, nu = prior_nu, sigma = S0)
}

# inds = 1:10
# one_pass_cpp(z[inds], x[inds,,drop=F], components, prior,
#              as.matrix(P)[inds,inds,drop=F], alpha, lambda, 1,
# 	     samplefun, samplefun2)
# z
#
# alpha = 2
# ## ohne subset
# ind = 2:100
# components = one_pass_cpp(ind, z, x, components, prior,
#                           as.matrix(P), alpha, lambda, 0,
#                           samplefun, samplefun2)
# components
# z
# ind = sample(1:100)
# components = one_pass_cpp(ind, z, x, components, prior,
#                           as.matrix(P), alpha, lambda, 0,
#                           samplefun, samplefun2)
# z
# # components

*/

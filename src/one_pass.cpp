// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <qi.cpp>
#include <update_list.cpp>
#include <logsumexp.cpp>

using namespace Rcpp;


// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }


// [[Rcpp::export]]
void one_pass_cpp (NumericVector z, arma::mat y, List components, List prior,
                         NumericMatrix P, double alpha, double lambda,
                         int debug) {
  int n = y.n_rows;
  int d = y.n_cols;
  int K = components.length();
  // Rcpp::print(wrap(K));

  // // create vector with shuffled indeces
  Range s = Rcpp::seq(1, n);
  NumericVector ind (s.begin(), s.end());
  // std::random_shuffle(ind.begin(), ind.end(), randWrapper);
  // Rcpp::print(wrap(ind));
  //

  int i = 1;
  // for (int i=0; i < n; i++) {
    // i = ind[i];

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

        if (debug == 1) {
          Rcpp::print(wrap("only obs in cluster: remove component (old_class = )"));
          Rcpp::print(wrap(old_class));
          CharacterVector nam = components.names();
          Rcpp::print(wrap("components before: "));
          Rcpp::print(wrap(nam));
          Rcpp::print(wrap("Move Kth component"));
          Rcpp::print(wrap(nam[K-1]));
          Rcpp::print(wrap("to spot of old_class - 1 "));
          Rcpp::print(wrap(nam[old_class-1]));
        }

        // cache the old_component, so it can possibly be restored
        List old_component = components[old_class-1];
        // move last component forward to empty spot
        components[old_class - 1] = components[K-1];
        Rcpp::print(wrap(z[z == K]));
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

        // update with all data except for current
        arma::colvec z2 = as<arma::colvec>(z);
        arma::mat ysub = y.rows(arma::find(z2 == old_class) - 1);
        List new_comp = new_component(components, old_class-1, ysub, 0);
        components[old_class - 1] = new_comp;

        if (debug == 1) {
          Rcpp::print(wrap("new component"));
          Rcpp::print(wrap(new_comp));
          Rcpp::print(wrap("used this data to update component"));
          Rcpp::print(wrap(ysub));
        }
      }
    }

    // compute assignment probabilities
    arma::colvec z2 = as<arma::colvec>(z);
    arma::mat ysub = y.rows(arma::find(z2 == 0));
    arma::vec q = get_assigment_prob(ysub, z, components, prior,
                                     P(i,_), lambda, alpha);
    // normalize
    q = exp(q - logsumexp(q));

    if (debug == 1) {
      print(wrap("compute assignment prob using:"));
      print(wrap(ysub));
      print(wrap("q:"));
      print(wrap(q));
    }

    // draw class
    double ret = R::runif(0, 1);
    // int ret = sample(seq_len(K), 1, false, q);
    print(wrap(ret));
    print(wrap(ret < q));


  // }
        // ## draw assignment
        // z[i] = sample(0:K, 1, prob = q)
        // ## if equal to 0, create new cluster
        // if (z[i] == 0) {
        //     if (sum(z == old_class) == 0) {
        //         ## just use the old singleton cluster if it had one
        //         components =
        //             purrr::splice(components,
        //                           list(old_component))
        //         ## components = append(components, components[[old_class]])
        //         ## components[[old_class]] = old_component
        //         K = K + 1
        //         z[i] = K
        //     } else {
        //         ## else, make a new one and update with the
        //         ## single observation
        //
        //         ## sample new prior values
        //         S0 = MCMCpack::rwish(prior_nu, prior_precision)
        //         m0 = mvtnorm::rmvnorm(1, prior_mean, 1 / prior_kappa * S0)
        //         components =
        //             purrr::splice(components,
        //                           list(normal_component(m0, prior_kappa, prior_nu, S0)))
        //         K = K + 1
        //         z[i] = K
        //
        //
        //     }
        // }
        // ## update stats of cluster z[i]
        // components[[z[i][1]]]$update(y[i, , drop = FALSE])
    // }
    // list(z = z, components = components)
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
one_pass_cpp(z[inds], x[inds,,drop=F], components, prior,
             as.matrix(P)[inds,inds,drop=F], alpha, lambda, 1)
*/
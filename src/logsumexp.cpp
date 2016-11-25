// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

using namespace Rcpp;

#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#  define EXPL expl
#else
#  define LDOUBLE double
#  define EXPL exp
#endif

// [[Rcpp::export]]
double logsumexp(const arma::vec& x) {
  unsigned int maxi = x.index_max();
  LDOUBLE maxv = x(maxi);
  LDOUBLE cumsum = 0.0;

  if (!(maxv > -arma::datum::inf)) {
    return -arma::datum::inf;
  }
  for (unsigned int i = 0; i < x.n_elem; i++) {
    if ((i != maxi) & (x(i) > -arma::datum::inf)) {
      cumsum += EXPL(x(i) - maxv);
    }
  }

  return maxv + log1p(cumsum);
}

/*** R
test_that("normalization works", {
  q = c(-4, -4.8, -0.9)
  p = exp(q - logsumexp(q))
  expect_equal(sum(p), 1)
})
*/

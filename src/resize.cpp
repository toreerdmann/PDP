#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List resize(List x, int new_size) {
  int K = x.length();
  List tmp(new_size);
  for (int i=0; i<K; i++) {
    tmp[i] = x[i];
  }
  return tmp;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
L = list(1, 2, 3)
L = resize(L, 5)
L
*/

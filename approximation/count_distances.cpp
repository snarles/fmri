#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector countDist(NumericMatrix Xm, NumericMatrix Ym, NumericVector rSq) {
  int n = Xm.nrow(), 
    m = Ym.nrow(),
    k = Xm.ncol();
  NumericVector counts(m);
  for (int i1 = 0; i1 < m; i1++) {
    for (int j1 = 0; j1 < n; j1++) {
      //compute distance between Y[i1, ] and X[j1, ]
      double sum = 0;
      for (int d1 = 0; d1 < k; d1++) {
        sum += pow(Xm(j1, d1) - Ym(i1, d1), 2.0);
      }
      if (sum < rSq[i1]) {
        counts[i1]++;
      }
    }
  }
  return counts; 
}


#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector countDistEx(NumericMatrix Xm, NumericMatrix Ym, NumericVector rSq) {
  int n = Xm.nrow(), 
    k = Xm.ncol();
  NumericVector counts(n);
  for (int i1 = 0; i1 < n; i1++) {
    for (int j1 = 0; j1 < n; j1++) {
      if (j1 != i1) {
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
  }
  return counts; 
}


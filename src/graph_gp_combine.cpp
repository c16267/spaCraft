// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericVector graph_gp_combine_cpp(NumericVector z,
                                   IntegerMatrix nn_index,   // 1-based, n x k
                                   NumericMatrix nn_dist,    //          n x k
                                   double rho,
                                   bool include_self = true) {
  const int n = nn_index.nrow();
  const int k = nn_index.ncol();
  
  if ((int) z.size() != n)
    stop("length(z) must equal nrow(nn_index).");
  if (nn_dist.nrow() != n || nn_dist.ncol() != k)
    stop("nn_dist and nn_index must have identical dimensions.");
  
  NumericVector eta(n);
  if (!R_finite(rho) || rho <= 0.0) rho = 10.0;
  
  for (int i = 0; i < n; ++i) {
    double wsum = 0.0;
    double acc  = 0.0;
    
    if (include_self) {        // self-weight = exp(0) = 1
      wsum += 1.0;
      acc  += z[i];
    }
    
    for (int j = 0; j < k; ++j) {
      const int idx = nn_index(i, j) - 1;       // to 0-based
      if (idx < 0 || idx >= n) continue;
      const double w = std::exp(-nn_dist(i, j) / rho);
      wsum += w;
      acc  += w * z[idx];
    }
    
    eta[i] = (wsum > 0.0 && R_finite(wsum)) ? (acc / wsum) : z[i];
  }
  
  return eta;
}
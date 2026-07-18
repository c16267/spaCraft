// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <limits>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector greedy_assign_cpp(NumericMatrix query_pts, 
                                NumericMatrix ref_pts, 
                                IntegerMatrix nn_idx) {
  int n_query = query_pts.nrow();
  int n_ref = ref_pts.nrow();
  int k_search = nn_idx.ncol();
  
  std::vector<bool> used(n_ref, false);
  IntegerVector assigned(n_query);
  
  for(int i = 0; i < n_query; ++i) {
    bool found = false;
    
    for(int j = 0; j < k_search; ++j) {
      int cand = nn_idx(i, j) - 1; // R(1-based) to C++(0-based)
      if(!used[cand]) {
        assigned[i] = cand + 1;    // C++(0-based) to R(1-based)
        used[cand] = true;
        found = true;
        break;
      }
    }
    
    if(!found) {
      double min_d2 = std::numeric_limits<double>::max();
      int best_j = -1;
      double qx = query_pts(i, 0);
      double qy = query_pts(i, 1);
      
      for(int j = 0; j < n_ref; ++j) {
        if(!used[j]) {
          double dx = ref_pts(j, 0) - qx;
          double dy = ref_pts(j, 1) - qy;
          double d2 = dx * dx + dy * dy;
          if(d2 < min_d2) {
            min_d2 = d2;
            best_j = j;
          }
        }
      }
      assigned[i] = best_j + 1;
      used[best_j] = true;
    }
  }
  
  return assigned;
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

// 1. Conditionally include omp.h for Mac (Apple Clang) compatibility
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List compute_expr_stats_cpp(const List& logcounts_list, 
                            const List& dom_idx_list, 
                            int n_genes, 
                            int n_domains, 
                            int n_cores) {
  
  int K = logcounts_list.size();
  
  // 1. Extract inputs into fast C++ std::vectors (SAFE: Outside parallel region)
  std::vector<arma::mat> Y(K);
  std::vector<arma::uvec> D_assign(K);
  std::vector<int> N_k(K);
  
  for(int k = 0; k < K; ++k) {
    Y[k] = as<arma::mat>(logcounts_list[k]);
    D_assign[k] = as<arma::uvec>(dom_idx_list[k]);
    N_k[k] = Y[k].n_cols;
  }
  
  // 2. Prepare Output Containers (Pure C++ Armadillo objects)
  arma::mat mu_grand(n_genes, n_domains, fill::zeros);
  arma::vec sigma_bio(n_genes, fill::zeros);
  
  // Create std::vector of matrices to avoid Rcpp List API inside OpenMP loop
  std::vector<arma::mat> res_vec(K);
  for(int k = 0; k < K; ++k) {
    res_vec[k] = arma::mat(n_genes, N_k[k], fill::zeros);
  }
  
  // 3. OpenMP Parallelization across GENES
  // Conditionally apply pragma for Mac compatibility
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_cores)
#endif
  for(int g = 0; g < n_genes; ++g) {
    
    arma::mat mu_k(K, n_domains, fill::zeros);
    arma::mat count_k(K, n_domains, fill::zeros);
    
    // Step A: Calculate Mean Expression per Group per Domain (mu_k_mat)
    for(int k = 0; k < K; ++k) {
      for(int i = 0; i < N_k[k]; ++i) {
        int d = D_assign[k][i];
        mu_k(k, d) += Y[k](g, i);
        count_k(k, d) += 1.0;
      }
      for(int d = 0; d < n_domains; ++d) {
        if(count_k(k, d) > 0) {
          mu_k(k, d) /= count_k(k, d);
        } else {
          mu_k(k, d) = datum::nan; // Mark NA if domain missing in sample
        }
      }
    }
    
    // Step B: Calculate mu_grand and sigma_bio
    double var_sum = 0;
    int valid_d_var = 0;
    
    for(int d = 0; d < n_domains; ++d) {
      double sum_mu = 0;
      int valid_k = 0;
      
      for(int k = 0; k < K; ++k) {
        if(!std::isnan(mu_k(k, d))) {
          sum_mu += mu_k(k, d);
          valid_k++;
        }
      }
      
      if(valid_k > 0) {
        double mg = sum_mu / valid_k;
        mu_grand(g, d) = mg;
        
        if(valid_k > 1) { // Need >=2 samples for variance
          double v = 0;
          for(int k = 0; k < K; ++k) {
            if(!std::isnan(mu_k(k, d))) {
              v += std::pow(mu_k(k, d) - mg, 2);
            }
          }
          v /= (valid_k - 1);
          var_sum += v;
          valid_d_var++;
        }
      } else {
        mu_grand(g, d) = 0.0; 
      }
    }
    
    sigma_bio(g) = valid_d_var > 0 ? var_sum / valid_d_var : 0.0;
    
    // Step C: Pre-compute Residuals (Pure C++ memory access, thread-safe)
    for(int k = 0; k < K; ++k) {
      for(int i = 0; i < N_k[k]; ++i) {
        int d = D_assign[k][i];
        res_vec[k](g, i) = Y[k](g, i) - mu_grand(g, d);
      }
    }
  }
  
  // 4. Repackage into Rcpp List safely OUTSIDE the parallel region
  List residuals_list(K);
  for(int k = 0; k < K; ++k) {
    residuals_list[k] = res_vec[k];
  }
  
  return List::create(Named("mu_grand") = mu_grand,
                      Named("sigma_bio") = sigma_bio,
                      Named("residuals") = residuals_list);
}
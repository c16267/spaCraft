// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// Helper: Log normalizing constant for vMF (exponentially scaled for stability)
double log_Cd_cpp(double tau, int d) {
  if (tau <= 0.0) {
    return -std::log(2.0) - (d / 2.0) * std::log(M_PI) + R::lgammafn(d / 2.0);
  }
  double nu = d / 2.0 - 1.0;
  // R::bessel_i with opt=2 returns I_nu(x) * exp(-x)
  double b = R::bessel_i(tau, nu, 2); 
  b = std::max(b, 1e-300);
  return nu * std::log(tau) - std::log(b) - (d / 2.0) * std::log(2.0 * M_PI);
}

// [[Rcpp::export]]
List FG_EM_cpp(const arma::mat& X, List initial, int iter_max, double tol, int n_cores = 1) {
  

#ifdef _OPENMP
  omp_set_num_threads(n_cores);
#endif
  
  int n = X.n_rows;
  int d = X.n_cols;
  int M = as<arma::vec>(initial["pi_hat"]).n_elem;
  
  arma::vec pi_hat = as<arma::vec>(initial["pi_hat"]);
  arma::mat c_hat = as<arma::mat>(initial["c_hat"]);
  arma::vec r_hat = as<arma::vec>(initial["r_hat"]);
  arma::mat phi_hat = as<arma::mat>(initial["phi_hat"]);
  arma::vec tau_hat = as<arma::vec>(initial["tau_hat"]);
  double sigma_sq = initial["sigma_sq"];
  
  std::vector<double> ll_history;
  ll_history.reserve(iter_max);
  
  arma::mat W(n, M, fill::zeros);
  std::vector<arma::mat> y_expect(M);
  for(int m = 0; m < M; ++m) y_expect[m] = arma::mat(n, d, fill::zeros);
  
  int iter = 0;
  for (; iter < iter_max; ++iter) {
    
    // ==========================================
    // E-step
    // ==========================================
    arma::mat logW(n, M, fill::zeros);
    
    for (int k = 0; k < M; ++k) {
      arma::rowvec ck = c_hat.row(k);
      arma::rowvec phik = phi_hat.row(k);
      double rk = r_hat[k];
      double tauk = tau_hat[k];
      
      arma::mat Z = X.each_row() - ck;
      arma::vec z_norm2 = arma::sum(Z % Z, 1);
      
      arma::mat V = Z * (rk / sigma_sq);
      V.each_row() += tauk * phik;
      arma::vec v_norm = arma::sqrt(arma::sum(V % V, 1));
      
      double logCd_tau = log_Cd_cpp(tauk, d);
      double const_term = logCd_tau - (d * std::log(2.0 * M_PI * sigma_sq) / 2.0) - (rk * rk) / (2.0 * sigma_sq) - tauk;
      
      
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (int i = 0; i < n; ++i) {
        double vn = v_norm[i];
        double logCd_v = log_Cd_cpp(vn, d);
        
        // Log density for W
        logW(i, k) = std::log(std::max(pi_hat[k], 1e-300)) + const_term - logCd_v - (z_norm2[i] / (2.0 * sigma_sq)) + vn;
        
        // y_expect
        double vn_safe = std::max(vn, 1e-10);
        double b1 = std::max(R::bessel_i(vn_safe, d/2.0, 2), 1e-300);
        double b2 = std::max(R::bessel_i(vn_safe, d/2.0 - 1.0, 2), 1e-300);
        double Ad_k = b1 / b2;
        Ad_k = std::min(std::max(Ad_k, 1e-10), 1.0 - 1e-10);
        
        y_expect[k].row(i) = Ad_k * (V.row(i) / vn_safe);
      }
    }
    
    // Log-Sum-Exp Trick
    arma::vec log_max = arma::max(logW, 1);
    arma::mat centered_logW = logW.each_col() - log_max;
    arma::vec sum_exp = arma::sum(arma::exp(centered_logW), 1);
    
    double ll = arma::sum(log_max + arma::log(sum_exp));
    ll_history.push_back(ll);
    
    // [Fix] Armadillo Delayed Evaluation Fix
    arma::mat exp_centered = arma::exp(centered_logW);
    W = exp_centered.each_col() / sum_exp;
    
    // Check Convergence
    if (iter > 1 && std::abs(ll - ll_history[iter - 1]) < tol) {
      iter++;
      break;
    }
    
    // ==========================================
    // M-step
    // ==========================================
    arma::rowvec W_sums = arma::sum(W, 0);
    pi_hat = W_sums.t() / n;
    
    arma::mat new_c_hat(M, d, fill::zeros);
    arma::vec new_r_hat(M, fill::zeros);
    arma::mat new_phi_hat(M, d, fill::zeros);
    arma::vec new_tau_hat(M, fill::zeros);
    double acc_sigma = 0.0;
    
    for (int k = 0; k < M; ++k) {
      arma::vec Wk = W.col(k);
      arma::mat yk = y_expect[k];
      double sum_wk = std::max(W_sums[k], 1e-12);
      
      arma::rowvec x_bar_w = arma::sum(X.each_col() % Wk, 0);
      arma::rowvec y_bar_w = arma::sum(yk.each_col() % Wk, 0);
      
      new_c_hat.row(k) = (x_bar_w - r_hat[k] * y_bar_w) / sum_wk;
      
      arma::mat res = X.each_row() - new_c_hat.row(k);
      arma::vec proj = arma::sum(res % yk, 1);
      new_r_hat[k] = arma::sum(Wk % proj) / sum_wk;
      
      arma::rowvec phi_k_w = arma::sum(yk.each_col() % Wk, 0);
      double phi_norm = arma::norm(phi_k_w, 2);
      if (phi_norm < 1e-8) {
        new_phi_hat.row(k) = arma::rowvec(d, fill::zeros);
        new_phi_hat(k, 0) = 1.0;
      } else {
        new_phi_hat.row(k) = phi_k_w / phi_norm;
      }
      
      double R_bar = phi_norm / sum_wk;
      double denom = std::max(1e-6, 1.0 - R_bar * R_bar);
      new_tau_hat[k] = std::max((R_bar * d - R_bar * R_bar * R_bar) / denom, 1e-6);
      
      arma::vec norm_sq = arma::sum(res % res, 1);
      acc_sigma += arma::sum(Wk % (norm_sq - 2.0 * new_r_hat[k] * proj + new_r_hat[k] * new_r_hat[k]));
    }
    
    c_hat = new_c_hat;
    r_hat = new_r_hat;
    phi_hat = new_phi_hat;
    tau_hat = new_tau_hat;
    sigma_sq = acc_sigma / (n * d);
  }
  
  return List::create(
    Named("pi") = pi_hat,
    Named("c") = c_hat,
    Named("r") = r_hat,
    Named("phi") = phi_hat,
    Named("tau") = tau_hat,
    Named("sigma_sq") = sigma_sq,
    Named("log_likelihood") = ll_history,
    Named("n_iter") = iter,
    Named("W") = W
  );
}
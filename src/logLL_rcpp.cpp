#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double logLL_rcpp(NumericVector par, List locus_info, NumericVector log_freq) {
  
  // 1. Parameter validation
  if (is_true(any(par < 0)) || sum(par) >= 1) {
    return -1.0e10; // Return a large negative number
  }
  
  // 2. Complete parameter vector and compute log ratios
  int n_par = par.size();
  NumericVector full_par(n_par + 1);
  for(int i = 0; i < n_par; ++i) {
    full_par[i] = par[i];
  }
  full_par[n_par] = 1.0 - sum(par);
  
  NumericVector log_par = log(full_par);
  NumericVector log_ratios = log_par - log_freq;
  
  // 3. Loop over loci
  int nloci = locus_info.size();
  double total_log_like = 0.0; // Initialize total log-likelihood
  
  for(int i = 0; i < nloci; ++i) {
    // Extract locus information
    List info = locus_info[i];
    int n_models = as<int>(info["n_models"]);
    
    if (n_models == 0) continue; // Skip if no models for this locus
    
    NumericVector log_probs = as<NumericVector>(info["log_probs"]);
    List model_cats = as<List>(info["model_cats"]);
    NumericVector model_log_likes(n_models);
    
    // 4. Inner loop over models (this is very fast in C++)
    for(int m = 0; m < n_models; ++m) {
      IntegerVector cats = model_cats[m];
      double sum_log_ratios = 0.0;
      
      for(int k = 0; k < cats.size(); ++k) {
        // IMPORTANT: R indices start at 1, C++ indices start at 0.
        // We must subtract 1 from the R index to get the correct C++ index.
        sum_log_ratios += log_ratios[cats[k] - 1]; 
      }
      model_log_likes[m] = log_probs[m] + sum_log_ratios;
    }
    
    // 5. Log-sum-exp for numerical stability
    double max_val = Rcpp::max(model_log_likes);
    if (R_finite(max_val)) {
      total_log_like += max_val + log(sum(exp(model_log_likes - max_val)));
    } else {
      // If max_val is not finite, add a large penalty
      total_log_like += -1.0e10; 
    }
  }
  
  return total_log_like;
}
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Last revised: 6/11/2023

// [[Rcpp::export]]
double LogLikeC(mat DQP, vec all_quantiles, vec y, vec mu_x, vec sigma_x, vec group_ind) 
{
  int n = y.n_elem;
  double nrows = all_quantiles.n_elem;
  double log_like = 0;
  
  for(int i = 0; i < n; ++i) {
    double y_unif_i = R::pnorm(y(i), mu_x(i), sigma_x(i), true, false);

    vec DQP_col(nrows);
    double t = 0;
    for (int j = 0; j < nrows; j++){
    
    // DQP_col calculation //
      DQP_col(j) = DQP(j, group_ind(i));
      
      // t calculation //
      if(y_unif_i > DQP_col(j)){
        t++;
      }
    } 
    
    if((t == 0) | (t == nrows)){ 
      // if t value is odd due to extreme a y_unif value //
      log_like = sqrt(-1); // NaN
      break;
    }

    double num = all_quantiles(t) - all_quantiles(t-1);
    double denom = DQP_col(t) - DQP_col(t-1);
    double log_phi_y = R::dnorm(y(i), mu_x(i), sigma_x(i), true); // ** Fixed as of 3/5/2023
    double log_f = log(num) - log(denom) + log_phi_y;

    log_like += log_f;
  }
  
  return log_like;
}




// // Last revised: 5/31/2023

// // [[Rcpp::export]]
// double logf(double yi, double mu_xi, double sigma_xi, int group_indi, mat DQP, vec all_quantiles)
// {
//   double y_unif_i = R::pnorm(yi, mu_xi, sigma_xi, true, false);
//   double nrows = all_quantiles.n_elem;
//   double t = 0;
//   vec DQP_col(nrows);
   
//   // #pragma omp parallel for reduction(+:t) // Let's try this later
//   for (int j = 0; j < nrows; j++){
    
//     // copying DQP_col //
//     DQP_col(j) = DQP(j, group_indi); // ** check here
    
//     // t calculation //
//     if(y_unif_i > DQP_col(j)){
//       t++;
//     }
//   } 

//   double num = all_quantiles(t) - all_quantiles(t-1);
//   double denom = DQP_col(t) - DQP_col(t-1);
//   double log_phi_y = R::dnorm(yi, mu_xi, sigma_xi, true); // ** Fixed as of 3/5/2023
//   double log_f = log(num) - log(denom) + log_phi_y;

//   return log_f;
// } 

// // [[Rcpp::export]]
// double LogLikeC(mat DQP, vec all_quantiles, vec y, vec mu_x, vec sigma_x, vec group_ind) 
// {
//   // DQP = log version (for all x)
//   // mu_x = long version (for all x)
//   // sigma_x = long version (for all x)
  
//   int n = y.n_elem;
//   double log_like = 0.0;

//   for(int i = 0; i < n; ++i) {
//     double y_i = y(i);
//     double mu_x_i = mu_x(i);
//     double sigma_x_i = sigma_x(i);
//     double grp_i = group_ind(i);
//     double log_f = logf(y_i, mu_x_i, sigma_x_i, grp_i, DQP, all_quantiles);
//     log_like += log_f;
//   }
//   return log_like;
// }



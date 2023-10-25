#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Last revised: 5/29/2023

// [[Rcpp::export]]
double LogLikeC(mat DQP, vec all_quantiles, vec y, vec mu_x, vec sigma_x, vec group_ind) 
{
  int n = y.n_elem;
  double y_unif_i;

  double nrows = all_quantiles.n_elem;
  double log_like = 0;
  
  for(int i = 0; i < n; ++i) {
    y_unif_i = R::pnorm(y(i), mu_x(group_ind(i)), sigma_x(group_ind(i)), true, false);

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
    
    double num = all_quantiles(t) - all_quantiles(t-1);
    double denom = DQP_col(t) - DQP_col(t-1);
    double log_phi_y = R::dnorm(y(i), mu_x(group_ind(i)), sigma_x(group_ind(i)), true); // ** Fixed as of 3/5/2023
    double log_f = log(num) - log(denom) + log_phi_y;

    log_like += log_f;
  }
  return log_like;
}

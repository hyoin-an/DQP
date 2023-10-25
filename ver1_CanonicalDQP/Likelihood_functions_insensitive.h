#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Last revised: 3/5/2023

// [[Rcpp::export]]
double LogLikeC(mat DQP, vec all_quantiles, vec y, vec mu_x, vec sigma_x, vec group_ind) 
{
  int n = y.n_elem;
  vec y_unif(n);

  double nrows = all_quantiles.n_elem;
  double log_like = 0;
  
  for(int i = 0; i < n; ++i) {
    y_unif(i) = R::pnorm(y(i), mu_x(group_ind(i)), sigma_x(group_ind(i)), true, false);

    vec DQP_col(nrows);
    double t = 0;
    for (int j = 0; j < nrows; j++){
    
    // DQP_col calculation //
      DQP_col(j) = DQP(j, group_ind(i));
      
      // t calculation //
      if(y_unif(i) > DQP_col(j)){
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
    double log_phi_y = R::dnorm(y(i), mu_x(group_ind(i)), sigma_x(group_ind(i)), true); // ** Fixed as of 3/5/2023
    double log_f = log(num) - log(denom) + log_phi_y;

    log_like += log_f;
  }
  
  return log_like;
}


// // Last revised: 3/2/2023
// // need to change all the other vectors to arma::vec & check the output values

// // [[Rcpp::export]]
// double LogLikeC(mat DQP, vec group_ind, vec y_unif, vec tau) 
// {
//   int n = y_unif.size();
//   double nrows = tau.size();
//   double log_like = 0;
  
//   for(int i = 0; i < n; ++i) {
  
//     vec DQP_col(nrows);
//     double t = 0;
      
//     for (int j = 0; j < nrows; j++){
    
//     // DQP_col calculation //
//       DQP_col(j) = DQP(j, group_ind(i));
      
//       // t calculation //
//       if(y_unif(i) > DQP_col(j)){
//         t++;
//       }
//     } 
    
//     double num = tau(t) - tau(t-1);
//     double denom = DQP_col(t) - DQP_col(t-1);
//     double log_phi_y = R::dnorm(R::qnorm(y_unif(i), 0, 1, true, false), 0, 1, true); // ** Fixed as of 2/17/2023
//     double log_f = log(num) - log(denom) + log_phi_y;

//     log_like += log_f;
//   }
//   return log_like;
// }

// // [[Rcpp::export]]
// double LogLikelihood_Discrete(mat DQP, vec all_quantiles, vec y, vec mu_x, vec sigma_x, vec group_ind)
// {
//     // vec group_ind_cpp;
//     // group_ind_cpp = group_ind - 1; // for C++ indexing

//     vec y_unif(y.n_elem);
//     for(int i=0; i<y.n_elem; i++){
//         y_unif(i) = R::pnorm(y(i), mu_x(group_ind(i)), sigma_x(group_ind(i)), true, false);
//     }
//     vec tau = all_quantiles;

//     double log_like;
//     log_like = LogLikeC(DQP, group_ind, y_unif, tau);

//     return log_like;
// }

// // [[Rcpp::export]]
// double LogLikeC(NumericMatrix DQP, vec group_ind_cpp,
//                 vec y_unif, vec tau) 
// {
//   int n = y_unif.size();
//   double nrows = tau.size();
//   double log_like = 0;
  
//   for(int i = 0; i < n; ++i) {
  
//     vec DQP_col(nrows);
//     double t = 0;
      
//     for (int j = 0; j < nrows; j++){
    
//     // DQP_col calculation //
//       DQP_col(j) = DQP(j, group_ind_cpp(i));
      
//       // t calculation //
//       if(y_unif(i) > DQP_col(j)){
//         t++;
//       }
//     } 
    
//     double num = tau(t) - tau(t-1);
//     double denom = DQP_col(t) - DQP_col(t-1);
//     double log_phi_y = R::dnorm(R::qnorm(y_unif(i), 0, 1, false), 0, 1, true); // ** Fixed as of 2/17/2023
//     double log_f = log(num) - log(denom) + log_phi_y;

//     log_like += log_f;
//   }
//   return log_like;
// }

// // [[Rcpp::export]]
// vec as_numeric(CharacterVector v)
// {
//     int n = v.size();
//     vec out(n);
//     std::string s;
//     for(int i=0; i<n; i++){
//         s = std::string(v[i]);
//         out[i] = std::stod(s);
//     }
//     return out;
// }

// // [[Rcpp::export]]
// double LogLikelihood_Discrete(NumericMatrix DQP, vec y, vec mu_x, vec sigma_x, vec group_ind)
// {
//   vec group_ind_cpp;
//   group_ind_cpp = group_ind - 1; // for C++ indexing

//   vec y_unif = exp(Rcpp::pnorm(y, mu_x(group_ind_cpp), sigma_x(group_ind_cpp)));
//   vec tau = as_numeric(rownames(DQP));

//   double log_like;
//   log_like = LogLikeC(DQP, group_ind_cpp, y_unif, tau);
  
//   return log_like;
// }



// // [[Rcpp::export]]
// double LogLikeC(NumericMatrix DQP, NumericVector group_ind_cpp,
//                 NumericVector y_unif, NumericVector tau) 
// {
//   int n = y_unif.size();
//   double nrows = tau.size();
//   double log_like = 0;
  
//   for(int i = 0; i < n; ++i) {
  
//     NumericVector DQP_col(nrows);
//     double t = 0;
      
//     for (int j = 0; j < nrows; j++){
    
//     // DQP_col calculation //
//       DQP_col[j] = DQP(j, group_ind_cpp[i]);
      
//       // t calculation //
//       if(y_unif[i] > DQP_col[j]){
//         t++;
//       }
//     } 
    
//     double num = tau[t] - tau[t-1];
//     double denom = DQP_col[t] - DQP_col[t-1];
//     double log_phi_y = R::dnorm(R::qnorm(y_unif[i], 0, 1, false), 0, 1, true); // ** Fixed as of 2/17/2023
//     double log_f = log(num) - log(denom) + log_phi_y;

//     log_like += log_f;
//   }
//   return log_like;
// }


// // [[Rcpp::export]]
// NumericVector as_numeric(CharacterVector v)
// {
//     int n = v.size();
//     NumericVector out(n);
//     std::string s;
//     for(int i=0; i<n; i++){
//         s = std::string(v[i]);
//         out[i] = std::stod(s);
//     }
//     return out;
// }

// // [[Rcpp::export]]
// double LogLikelihood_Discrete(NumericMatrix DQP, NumericVector y, NumericVector mu_x, NumericVector sigma_x, NumericVector group_ind)
// {
//   NumericVector group_ind_cpp;
//   group_ind_cpp = group_ind - 1; // for C++ indexing

//   NumericVector y_unif = exp(Rcpp::pnorm(y, mu_x[group_ind_cpp], sigma_x[group_ind_cpp]));
//   NumericVector tau = as_numeric(rownames(DQP));

//   double log_like;
//   log_like = LogLikeC(DQP, group_ind_cpp, y_unif, tau);
  
//   return log_like;
// }

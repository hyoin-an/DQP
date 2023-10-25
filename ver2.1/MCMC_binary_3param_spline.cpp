#include <RcppArmadillo.h>
#include "utilities_spline.h"
#include "DQPpriorSampling_functions.h" // not used in mcmc, but used to generate an initial value //
#include "DQPdensity_functions.h"
#include "Likelihood_functions.h"
#include "UpdateQ_binary_spline.h"
#include "UpdateBeta_binary_spline.h"
#include "UpdateGamma_binary_linear_new.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Last revised: 7/16/2023 
// need to check the output values

// [[Rcpp::export]]
mat unique_rows(const mat& x) {
  
  // Initialize matrix for unique rows
  mat unique_x = x.rows(0, 0);
  
  for (uword i = 1; i < x.n_rows; ++i) {
    bool is_unique = true;
    for (uword j = 0; j < unique_x.n_rows; ++j) {
      if (all(x.row(i) == unique_x.row(j))) {
        is_unique = false;
        break;
      }
    }
    if (is_unique) {
      unique_x = join_vert(unique_x, x.rows(i, i));
    }
  }
  
  return unique_x;
}

// [[Rcpp::export]]
vec unique_rows_identifier(mat x) {
  
  // Create a mapping of unique rows to their identifiers
  std::unordered_map<std::string, uword> unique_row_map;
  uword next_id = 0;
  
  // Create a vector of identifiers for each row
  vec id_vec(x.n_rows);
  for (uword i = 0; i < x.n_rows; ++i) {
    std::string row_str;
    for (uword j = 0; j < x.n_cols; ++j) {
      if (x(i,j) == -0){ // fixing the issue of 0 != -0
        x(i,j) = 0;
      }
      row_str += std::to_string(x(i,j));
    }

    auto iter = unique_row_map.find(row_str);
    if (iter == unique_row_map.end()) {
      // This is a new unique row
      unique_row_map[row_str] = next_id++;
    }
    id_vec(i) = unique_row_map[row_str];
  }
  
  return id_vec; 
}




// [[Rcpp::export]]
List MCMC_param3(int niter, vec y, mat x, mat x_QP, vec x_SP, List quantile_levels, 
                 vec mu, mat Sigma, vec mu_0, mat Sigma_0, vec mu_1, mat Sigma_1, double alpha_scale, 
                 double delta, mat Sigma_beta, mat Sigma_gamma, int n_chunk_Q, int n_chunk_beta, int n_chunk_gamma,
                 NumericMatrix DQP_Q_initial_rcpp, mat DQP_Z_initial, vec beta_initial, vec gamma_initial,
                 int d=0 // d = order of spline functions (knots = x_QP, 1-dimensional only)
                 ){
    // ** Data related input **
    // y = n x 1, original response variable vector
    // x = n x 1, original covariate matrix (replicates O, all locations O), 1-dimensional (due to spline mu & linear interpolation)
    //   -> has to be sorted for interpolation and should be univariate (or up to polynomial, e.g [x x^2])
    // x_QP = sorted, selected covariate matrix with QP locations (replicates X, all locations X)

    // sort x and y
    uvec sortedIndices = sort_index(x.col(0)); // obtain the order based on the first column of x
    x = x.rows(sortedIndices);  // sort mat x based on its first column
    y = y(sortedIndices);       // sort vec y based on x's first column

    // * maybe it's better to allow replicates for any inputs and remove group_ind...?
    // * NO, we have to have a group ind because we want the same parameter values for the same location
    
    // Getting 'all_quantiles' vector
    vec all_quantiles = as_numeric(rownames(DQP_Q_initial_rcpp));
    mat DQP_Q_initial = as<mat>(DQP_Q_initial_rcpp);

    // Handling design matrix //
    mat X = join_horiz(ones(x.n_rows), x);          // Design matrix with all data points
    mat X_prime = unique_rows(X);                   // Deisgn matrix without replicates
    vec group_ind = unique_rows_identifier(X);      // group indicator for QP locations
    mat X_QP = join_horiz(ones(x_QP.n_rows), x_QP); // Locations where QPs are constructed without replicates
    
    // building a special design matrix for spline //

    mat X_sp = SplineDesignMatrix(X, x_SP, d);
    mat X_sp_QP = SplineDesignMatrix(X_QP, x_SP, d);

    // mat X_sp = ones(x.n_rows);  
    // mat X_sp_QP = ones(x_QP.n_rows); 

    // if(d >= 1){ // d == 0; no spline
    //   // base of spline design matrix
    //   for(int l=0; l<d; l++){
    //     X_sp = join_horiz(X_sp, pow(X.col(1),l+1));
    //     X_sp_QP = join_horiz(X_sp_QP, pow(X_QP.col(1),l+1));
    //   }
    //   // additional part of design matrix
    //   for(int l=0; l<x_QP.n_rows; l++){
    //     X_sp = join_horiz(X_sp, pow(zero_plus(X.col(1)-X_QP(l,1)*ones(X_sp.n_rows)),d));
    //     X_sp_QP = join_horiz(X_sp_QP, pow(zero_plus(X_QP.col(1)-X_QP(l,1)*ones(X_sp_QP.n_rows)),d));
    //   }
    // }
    
    // For efficient computation //
    mat Sigma_inv = inv(Sigma);
    double log_Sigma_det = log(det(Sigma)); 
    mat Sigma_0_inv = inv(Sigma_0);
    double log_Sigma_0_det = log(det(Sigma_0));
    mat Sigma_1_inv = inv(Sigma_1);
    double log_Sigma_1_det = log(det(Sigma_1));

    // initial values //
    mat DQP_Q_post = DQP_Q_initial;
    mat DQP_Z_post = DQP_Z_initial;
    vec beta_post = beta_initial;
    vec mu_x = X_sp_QP * beta_post; // **new** 
    vec mu_x_long = X_sp * beta_post; // **new** 
    vec gamma_post = gamma_initial;

    List Q_posteriors(niter);
    List beta_posteriors(niter);
    List gamma_posteriors(niter);

    // iteration //
    for(int i=0; i<niter; i++){
        // Rcpp::Rcout << "Interation: " << i+1 << " begins..." << std::endl;
        // Rcpp::Rcout << "Updating Q..." << std::endl;
        List DQP_post = UpdateQ(DQP_Q_post, DQP_Z_post, mu_x, mu_x_long, gamma_post, 
                                quantile_levels, all_quantiles, y, X, X_prime, group_ind, X_QP, 
                                mu, Sigma, Sigma_inv, log_Sigma_det, delta, alpha_scale, n_chunk_Q);
        
        DQP_Q_post = as<mat>(DQP_post["Q"]);
        DQP_Z_post = as<mat>(DQP_post["Z"]);
        // Rcpp::Rcout << DQP_Q_post << std::endl;

        // Rcpp::Rcout << "Updating beta..." << std::endl;
        beta_post = UpdateBeta(DQP_Q_post, beta_post, gamma_post, 
                               quantile_levels, all_quantiles, y, X, X_prime, group_ind, X_QP, X_sp, X_sp_QP,
                               mu, Sigma, Sigma_inv, log_Sigma_det, alpha_scale, Sigma_beta, 
                               mu_0, Sigma_0_inv, log_Sigma_0_det, n_chunk_beta);
        // Rcpp::Rcout << beta_post << std::endl;
        mu_x = X_sp_QP * beta_post; // **new** 
        mu_x_long = X_sp * beta_post; // **new** 

        // Rcpp::Rcout << "Updating gamma..." << std::endl;
        gamma_post = UpdateGamma(DQP_Q_post, mu_x, mu_x_long, gamma_post, 
                                 quantile_levels, all_quantiles, y, X, X_prime, group_ind, X_QP,  
                                 mu, Sigma, Sigma_inv, log_Sigma_det, alpha_scale, Sigma_gamma,
                                 mu_1, Sigma_1_inv, log_Sigma_1_det, n_chunk_gamma);
        // Rcpp::Rcout << gamma_post << std::endl;

        Q_posteriors[i] = DQP_Q_post; // we need to copy, not using reference
        beta_posteriors[i] = beta_post;
        gamma_posteriors[i] = gamma_post;
    }

    List posterior_all = List::create(Named("Q") = Q_posteriors, Named("beta") = beta_posteriors, Named("gamma") = gamma_posteriors);

    return posterior_all;
}


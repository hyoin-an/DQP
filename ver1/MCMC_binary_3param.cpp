#include <RcppArmadillo.h>
#include "utilities.h"
#include "DQPpriorSampling_functions.h" // not used in mcmc, but used to generate an initial value //
#include "DQPdensity_functions.h"
#include "Likelihood_functions_insensitive.h"
#include "UpdateQ_binary.h"
#include "UpdateBeta_binary.h"
#include "UpdateSigmaX_binary.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Last revised: 9/13/2023 
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
vec unique_rows_identifier(const mat& x) {
  
  // Create a mapping of unique rows to their identifiers
  std::unordered_map<std::string, uword> unique_row_map;
  uword next_id = 1;
  
  // Create a vector of identifiers for each row
  vec id_vec(x.n_rows);
  for (uword i = 0; i < x.n_rows; ++i) {
    std::string row_str;
    for (uword j = 0; j < x.n_cols; ++j) {
      row_str += std::to_string(x(i,j));
    }
    auto iter = unique_row_map.find(row_str);
    if (iter == unique_row_map.end()) {
      // This is a new unique row
      unique_row_map[row_str] = next_id++;
    }
    id_vec(i) = unique_row_map[row_str];
  }
  
  return id_vec-1; // to make the first index = 0, instead of 1
}

// [[Rcpp::export]]
List MCMC_param3(int niter, vec y, mat x, List quantile_levels, 
                 vec mu, mat Sigma, double delta, double alpha_scale, 
                 vec mu_0, mat Sigma_0, mat Sigma_beta, vec mu_1, mat Sigma_1, mat Sigma_eta, 
                 int n_chunk_Q, int n_chunk_beta, int n_chunk_eta,
                 NumericMatrix DQP_Q_initial_rcpp, mat DQP_Z_initial, vec beta_initial, vec sigma_x_initial
                 ){
    
    // Getting 'all_quantiles' vector
    vec all_quantiles = as_numeric(rownames(DQP_Q_initial_rcpp));
    mat DQP_Q_initial = as<mat>(DQP_Q_initial_rcpp);

    // Handling design matrix //
    mat X = join_horiz(ones(y.n_elem), x); // Design matrix with replicates
    mat X_prime = unique_rows(X); // Design matrix without replicates
    vec group_ind = unique_rows_identifier(X); 

    // For efficient computation //
    mat Sigma_inv = inv(Sigma);
    double log_Sigma_det = log(det(Sigma)); // this det() function needs to be checked
    mat Sigma_0_inv = inv(Sigma_0);
    double log_Sigma_0_det = log(det(Sigma_0));
    mat Sigma_1_inv = inv(Sigma_1);
    double log_Sigma_1_det = log(det(Sigma_1));

    // initial values //
    mat DQP_Q_post = DQP_Q_initial;
    mat DQP_Z_post = DQP_Z_initial;
    vec beta_post = beta_initial;
    vec mu_x = X_prime * beta_post; // # need this for each group, not individually
    vec sigma_x_post = sigma_x_initial;

    List Q_posts(niter);
    List beta_posts(niter);
    List sigma_x_posts(niter);

    // iteration //
    for(int i=0; i<niter; i++){
        // A few major changes
        // 1. Default for X is X_prime (without duplicates) for everywhere
        //    If the full X is needed, we do that by X_prime[group_ind]
        // 2. Likelihood function does not take beta as an input

        // Rcpp::Rcout << "Iteration: " << i+1 << " begins..." << std::endl;
        List DQP_post = UpdateQ(DQP_Q_post, DQP_Z_post, y, mu_x, sigma_x_post, group_ind, quantile_levels, all_quantiles,
                                mu, Sigma, Sigma_inv, log_Sigma_det, delta, alpha_scale, n_chunk_Q);
        
        DQP_Q_post = as<mat>(DQP_post["Q"]);
        DQP_Z_post = as<mat>(DQP_post["Z"]);
        Q_posts[i] = DQP_Q_post; // fixed as of 8/27
        // Rcpp::Rcout << "Update Q: " << DQP_Q_post << std::endl;

        beta_post = UpdateBeta(DQP_Q_post, beta_post, X_prime, y, sigma_x_post, group_ind, quantile_levels, all_quantiles,
                                     mu, Sigma, Sigma_inv, log_Sigma_det, alpha_scale, Sigma_beta, mu_0, Sigma_0_inv, log_Sigma_0_det, n_chunk_beta);
        mu_x = X_prime * beta_post; // # need this for each group, not individually
        beta_posts[i] = beta_post;
        // Rcpp::Rcout << "Update beta: " << beta_post << std::endl;

        sigma_x_post = UpdateSigmaX(DQP_Q_post, y, mu_x, sigma_x_post, group_ind, quantile_levels, all_quantiles,
                                          mu, Sigma, Sigma_inv, log_Sigma_det, alpha_scale, Sigma_eta, mu_1, Sigma_1_inv, log_Sigma_1_det, n_chunk_eta);
        sigma_x_posts[i] = sigma_x_post;
        // Rcpp::Rcout << "Update sigma_x: " << sigma_x_post << std::endl;
    }

    List posterior_all = List::create(Named("Q") = Q_posts, Named("beta") = beta_posts, Named("sigma_x") = sigma_x_posts);

    return posterior_all;
}



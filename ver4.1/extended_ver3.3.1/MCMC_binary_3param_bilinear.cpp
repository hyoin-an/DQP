#include <RcppArmadillo.h>
#include "utilities_bilinear.h"
#include "DQPpriorSampling_functions.h" // not used in mcmc, but used to generate an initial value //
#include "DQPdensity_functions.h"
#include "Likelihood_functions.h"
#include "UpdateQ_binary_bilinear.h"
#include "UpdateBeta_binary_bilinear.h"
#include "UpdateEta_binary_bilinear.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Last revised: 12/4/2023 
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
List MCMC_param3(int niter, vec y, mat x, mat x_QP, List quantile_levels, uvec locs, // locs = {loc1, loc2}
                 vec mu, mat Sigma, vec mu_0, mat Sigma_0, vec mu_1, mat Sigma_1, double alpha_scale, 
                 mat Sigma_p, double alpha_scale_p, mat Sigma_beta, mat Sigma_eta, 
                 int n_chunk_Q, int n_chunk_beta, int n_chunk_eta,
                 NumericMatrix DQP_Q_initial_rcpp, mat DQP_Z_initial, vec beta_initial, vec eta_initial
                 ){
    // ** Data related input **
    // y = n x 1, original response variable vector
    // x = n x 2, original covariate matrix (replicates O, all locations O)
    // x_QP = has to be sorted, selected covariate matrix with QP locations (replicates X, all locations X)

    // sort x and y -- I don't think this is necessary, but the logic becomes easier if we sort these too (**new**)
    //  -> applying sorting twice using Excel sorting principle for 2 variables (tested) 
    uvec sortedIndices_1 = sort_index(x.col(1));
    x = x.rows(sortedIndices_1);  
    y = y(sortedIndices_1);       
    uvec sortedIndices_0 = sort_index(x.col(0)); 
    x = x.rows(sortedIndices_0);  
    y = y(sortedIndices_0);      

    uvec alocs = AddIntercept(locs); //** new: (loc1, loc2) -> (0, loc1, loc2) ** as of 10/21/2023   

    // MAKE SURE THAT x_QP is sorted in a correct format; this is necessary (**new**) 
    //  -> applying sorting twice using Excel sorting principle for 2 variables 
    //  -> quantiles will be generated for each pair, bestowed the same order
    //  -> this occurs an issue on Unity (when there are ties, they are not ordered properly)
    //  -> for now, make sure the input already has the correct order
    // uvec sortedIndices_QP1 = sort_index(x_QP.col(1));
    // x_QP = x_QP.rows(sortedIndices_QP1);
    // uvec sortedIndices_QP0 = sort_index(x_QP.col(0));
    // x_QP = x_QP.rows(sortedIndices_QP0); 

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
    vec eta_post = eta_initial;
    vec eta_post_long = bilinearInterpolator(X_QP.cols(1,2), eta_post, X_prime.cols(1,2));
    vec sigma_x_post = exp(eta_post);
    vec sigma_x_post_long = exp(eta_post_long);

    List Q_posteriors(niter);
    List beta_posteriors(niter);
    List sigma_x_posteriors(niter);

    // iteration //
    for(int i=0; i<niter; i++){
        // Rcpp::Rcout << "Iteration: " << i << " begins..." << std::endl;
        // Rcpp::Rcout << "Updating Q..." << std::endl;
        List DQP_post = UpdateQ(DQP_Q_post, DQP_Z_post, beta_post, sigma_x_post, sigma_x_post_long, 
                                quantile_levels, all_quantiles, alpha_scale, y, X, X_prime, group_ind, X_QP, locs, alocs,
                                mu, Sigma, Sigma_inv, log_Sigma_det, Sigma_p, alpha_scale_p, n_chunk_Q);
        
        DQP_Q_post = as<mat>(DQP_post["Q"]);
        DQP_Z_post = as<mat>(DQP_post["Z"]);
        // Rcpp::Rcout << DQP_Q_post << std::endl;

        // Rcpp::Rcout << "Updating beta..." << std::endl;
        beta_post = UpdateBeta(DQP_Q_post, beta_post, sigma_x_post, sigma_x_post_long, 
                               quantile_levels, all_quantiles, y, X, X_prime, group_ind, X_QP, locs, alocs,
                               mu, Sigma, Sigma_inv, log_Sigma_det, alpha_scale, Sigma_beta, 
                               mu_0, Sigma_0_inv, log_Sigma_0_det, n_chunk_beta);
        // Rcpp::Rcout << beta_post << std::endl;

        // Rcpp::Rcout << "Updating eta..." << std::endl;
        eta_post = UpdateEta(DQP_Q_post, beta_post, eta_post, 
                             quantile_levels, all_quantiles, y, X, X_prime, group_ind, X_QP, locs, alocs,
                             mu, Sigma, Sigma_inv, log_Sigma_det, alpha_scale, Sigma_eta,
                             mu_1, Sigma_1_inv, log_Sigma_1_det, n_chunk_eta);
        eta_post_long = bilinearInterpolator(X_QP.cols(1,2), eta_post, X.cols(locs));
        sigma_x_post = exp(eta_post);
        sigma_x_post_long = exp(eta_post_long);

        Q_posteriors[i] = DQP_Q_post; // we need to copy, not using reference
        beta_posteriors[i] = beta_post;
        sigma_x_posteriors[i] = sigma_x_post;
    }

    List posterior_all = List::create(Named("Q") = Q_posteriors, Named("beta") = beta_posteriors, Named("sigma_x") = sigma_x_posteriors);

    return posterior_all;
}


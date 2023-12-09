#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Last revised: 5/31/2023 
// need to check the output values


// [[Rcpp::export]]
vec UpdateBeta(mat DQP, vec beta_c, vec gamma, List quantile_levels, vec all_quantiles, 
               vec y, mat X, mat X_prime, vec group_ind, mat X_QP, // vec group_ind_QP, 
               vec mu, mat Sigma, mat Sigma_inv, double log_Sigma_det,  double alpha_scale, mat Sigma_beta, 
               vec mu_0, mat Sigma_0_inv, double log_Sigma_0_det, int n_chunk=2)
{
    // X = design matrix with replicates including intercept column
    // X_prime = design matrix without replicates including intercept column
    // Sigma_beta = This is a covariance matrix for proposal, could be used in Adaptive metropolis later!

    uvec ind_chunks = chunk_indices(beta_c.n_elem, n_chunk);
    // update chunks (lower-dimensional update) //
    for (uword i = 0; i < ind_chunks.n_elem-1; i++){
        vec beta_p = beta_c;

        // (1) Propose a new candidate
        uvec ind = regspace<uvec>(ind_chunks(i), ind_chunks(i+1)-1);
        vec beta_p_ind = proposeNormal_chunk(ind, beta_c, Sigma_beta);
        beta_p(ind) = beta_p_ind;
        
        // (2) Calculate (log) acceptance probability: a_star
        // Log prior density of beta: Normal
        double log_beta_prior_p = LogNormalPrior(beta_p, mu_0, Sigma_0_inv, log_Sigma_0_det);
        double log_beta_prior_c = LogNormalPrior(beta_c, mu_0, Sigma_0_inv, log_Sigma_0_det);

        // Log Likelihood
        vec mu_x_p = X_QP * beta_p;
        vec mu_x_c = X_QP * beta_c;
        vec mu_x_p_long = X * beta_p;   // **NEW**
        vec mu_x_c_long = X * beta_c;   // **NEW**
        vec sigma_x = exp(X_QP * gamma);     // **NEW**
        vec sigma_x_long = exp(X * gamma);   // **NEW**

        mat DQP_u_p(DQP.n_rows, DQP.n_cols);
        mat DQP_u_c(DQP.n_rows, DQP.n_cols);
        for(int r=0; r < DQP.n_rows; r++){
            DQP_u_p.row(r) = pnorm_cpp(DQP.row(r).t(), mu_x_p, sigma_x).t();
            DQP_u_c.row(r) = pnorm_cpp(DQP.row(r).t(), mu_x_c, sigma_x).t();
        }

        // **NEW: interpolationg Q in uniform scale **
        mat DQP_u_p_long = interpolater_mat(X_QP.col(1), DQP_u_p, X_prime.col(1)); // interpolating based on the lowest order of polynomials
        mat DQP_u_c_long = interpolater_mat(X_QP.col(1), DQP_u_c, X_prime.col(1)); // interpolating based on the lowest order of polynomials
       
        double log_like_p = LogLikeC(DQP_u_p_long, all_quantiles, y, mu_x_p_long, sigma_x_long, group_ind);
        double log_like_c = LogLikeC(DQP_u_c_long, all_quantiles, y, mu_x_c_long, sigma_x_long, group_ind);

        // we need to include these lines b/c prior of Q (uniform/y-scale) does depend on beta 
        double log_Q_density_given_beta_p = DQPlogPrior_binary(DQP_u_p, quantile_levels, all_quantiles, mu_x_p, sigma_x, mu, Sigma, Sigma_inv, log_Sigma_det, alpha_scale);
        double log_Q_density_given_beta_c = DQPlogPrior_binary(DQP_u_c, quantile_levels, all_quantiles, mu_x_c, sigma_x, mu, Sigma, Sigma_inv, log_Sigma_det, alpha_scale);

        // Log posterior of beta (proposal density ratio is 1 due to symmetry)
        double log_posterior_p = log_like_p + log_beta_prior_p + log_Q_density_given_beta_p;
        double log_posterior_c = log_like_c + log_beta_prior_c + log_Q_density_given_beta_c;

        double log_a_star = log_posterior_p - log_posterior_c;

        // (3) Decide to update
        // [Decide to accept vs. reject]
        double u = R::runif(0, 1);
        if (!Rcpp::NumericVector::is_na(log_a_star)) {
            // consider only when Q_p is reasonably generated
            if (log(u) <= log_a_star) {
                beta_c = beta_p;
                }
            }
    }
    return beta_c;
}

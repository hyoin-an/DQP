#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Last revised: 11/26/2023
// need to check the output values


// [[Rcpp::export]]
vec UpdateGamma(mat DQP, vec beta, vec gamma_c, List quantile_levels, vec all_quantiles, 
                vec y, mat X, mat X_prime, vec group_ind, mat X_QP, uvec locs, uvec alocs, // ** new **
                vec mu, mat Sigma, mat Sigma_inv, double log_Sigma_det, double alpha_scale, mat Sigma_gamma, 
                vec mu_1, mat Sigma_1_inv, double log_Sigma_1_det, int n_chunk=5)
{
    // X = design matrix with replicates including intercept column
    // X_prime = design matrix without replicates including intercept column
    // Sigma_gamma = This is a covariance matrix for proposal, could be used in Adaptive metropolis later!

    uvec ind_chunks = chunk_indices(gamma_c.n_elem, n_chunk);

    // update chunks (lower-dimensional update) //
    for (uword i = 0; i < ind_chunks.n_elem-1; i++){
        vec gamma_p = gamma_c;

        // (1) Propose a new candidate  
        uvec ind = regspace<uvec>(ind_chunks(i), ind_chunks(i+1)-1);
        vec gamma_p_ind = proposeNormal_chunk(ind, gamma_c, Sigma_gamma);
        gamma_p(ind) = gamma_p_ind;

        // (2) Calculate (log) acceptance probability: a_star

        // Log Prior densities: Normal
        double log_gamma_prior_p = LogNormalPrior(gamma_p, mu_1, Sigma_1_inv, log_Sigma_1_det);
        double log_gamma_prior_c = LogNormalPrior(gamma_c, mu_1, Sigma_1_inv, log_Sigma_1_det);

        // Log Likelihood
        vec mu_x = X_QP * beta.elem(alocs); // **NEW**
        vec mu_x_long = X * beta;  // **NEW**
        vec sigma_x_p = exp(X_QP * gamma_p.elem(alocs));
        vec sigma_x_c = exp(X_QP * gamma_c.elem(alocs));
        vec sigma_x_p_long = exp(X * gamma_p); // **NEW**
        vec sigma_x_c_long = exp(X * gamma_c); // **NEW**

        // Rcpp::Rcout << "sigma_x_p: " << sigma_x_p << std::endl;
        // Rcpp::Rcout << "sigma_x_c: " << sigma_x_c << std::endl;
        
        mat DQP_u_p(DQP.n_rows, DQP.n_cols);
        mat DQP_u_c(DQP.n_rows, DQP.n_cols);
        for(int r=0; r < DQP.n_rows; r++){
            DQP_u_p.row(r) = pnorm_cpp(DQP.row(r).t(), mu_x, sigma_x_p).t();
            DQP_u_c.row(r) = pnorm_cpp(DQP.row(r).t(), mu_x, sigma_x_c).t();
        }

        // **NEW: interpolationg Q in uniform scale **
        mat DQP_u_p_long = interpolator_mat(X_QP.col(1), DQP_u_p, X_prime.cols(locs)); // interpolating based on the lowest order of polynomials
        mat DQP_u_c_long = interpolator_mat(X_QP.col(1), DQP_u_c, X_prime.cols(locs)); // interpolating based on the lowest order of polynomials

        // Rcpp::Rcout << "DQP_u_p: " << DQP_u_p << std::endl;
        // Rcpp::Rcout << "DQP_u_c: " << DQP_u_c << std::endl;

        double log_like_p = LogLikeC(DQP_u_p_long, all_quantiles, y, mu_x_long, sigma_x_p_long, group_ind);
        double log_like_c = LogLikeC(DQP_u_c_long, all_quantiles, y, mu_x_long, sigma_x_c_long, group_ind);

        // we need to include these lines b/c prior of Q (uniform/y-scale) does depend on sigma 
        double log_Q_density_given_eta_p = DQPlogPrior_binary(DQP_u_p, quantile_levels, all_quantiles, mu_x, sigma_x_p, mu, Sigma, Sigma_inv, log_Sigma_det, alpha_scale);
        double log_Q_density_given_eta_c = DQPlogPrior_binary(DQP_u_c, quantile_levels, all_quantiles, mu_x, sigma_x_c, mu, Sigma, Sigma_inv, log_Sigma_det, alpha_scale);

        // Log posterior of eta (RandomWalk Metropolis: proposal density ratio is 1 due to symmetry)
        double log_posterior_p = log_like_p + log_gamma_prior_p + log_Q_density_given_eta_p;
        double log_posterior_c = log_like_c + log_gamma_prior_c + log_Q_density_given_eta_c;

        double log_a_star = log_posterior_p - log_posterior_c;

        // (3) Decide to update
        // [Decide to accept vs. reject]
        double u = R::runif(0, 1);
        if (!Rcpp::NumericVector::is_na(log_a_star)) {
            // consider only when Q_p is reasonably generated
            if (log(u) <= log_a_star) {
                gamma_c = gamma_p;
            }
        }
    }
    
    return gamma_c;
}

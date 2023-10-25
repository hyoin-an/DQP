#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Last revised: 9/13/2023 
// need to check the output values

// [[Rcpp::export]]
vec UpdateSigmaX(mat DQP, vec y, vec mu_x, vec sigma_x_c, vec group_ind, List quantile_levels, vec all_quantiles,
                 vec mu, mat Sigma, mat Sigma_inv, double log_Sigma_det,  double alpha_scale, mat Sigma_eta,
                 vec mu_1, mat Sigma_1_inv, double log_Sigma_1_det, int n_chunk=5)
{
    // X = design matrix without replicates including intercept column
    // Lambda = This is a covariance matrix for proposal, could be used in Adaptive metropolis later!

    uvec ind_chunks = chunk_indices(sigma_x_c.n_elem, n_chunk);
    vec eta_c = log(sigma_x_c);

    // update chunks (lower-dimensional update) //
    for (uword i = 0; i < ind_chunks.n_elem-1; i++){
        
        // (1) Propose a new candidate  
        uvec ind = regspace<uvec>(ind_chunks(i), ind_chunks(i+1)-1);
        vec eta_p_ind = proposeNormal_chunk(ind, eta_c, Sigma_eta);

        vec eta_p = eta_c;
        eta_p(ind) = eta_p_ind;
        vec sigma_x_p = exp(eta_p);

        // (2) Calculate (log) acceptance probability: a_star

        // Log prior density of eta: dependent log normal
        double log_eta_prior_p = LogNormalPrior(eta_p, mu_1, Sigma_1_inv, log_Sigma_1_det);
        double log_eta_prior_c = LogNormalPrior(eta_c, mu_1, Sigma_1_inv, log_Sigma_1_det);

        mat DQP_u_p(DQP.n_rows, DQP.n_cols);
        mat DQP_u_c(DQP.n_rows, DQP.n_cols);
        for(int r=0; r < DQP.n_rows; r++){
            DQP_u_p.row(r) = pnorm_cpp(DQP.row(r).t(), mu_x, sigma_x_p).t();
            DQP_u_c.row(r) = pnorm_cpp(DQP.row(r).t(), mu_x, sigma_x_c).t();
        }

        // Log Likelihood
        double log_like_p = LogLikeC(DQP_u_p, all_quantiles, y, mu_x, sigma_x_p, group_ind);
        double log_like_c = LogLikeC(DQP_u_c, all_quantiles, y, mu_x, sigma_x_c, group_ind);

        // we need to include these lines b/c prior of Q (uniform/y-scale) does depend on eta 
        double log_Q_density_given_eta_p = DQPlogPrior_binary(DQP_u_p, quantile_levels, all_quantiles, mu_x, sigma_x_p, mu, Sigma, Sigma_inv, log_Sigma_det, alpha_scale);
        double log_Q_density_given_eta_c = DQPlogPrior_binary(DQP_u_c, quantile_levels, all_quantiles, mu_x, sigma_x_c, mu, Sigma, Sigma_inv, log_Sigma_det, alpha_scale);

        // Log posterior of eta (RandomWalk Metropolis: proposal density ratio is 1 due to symmetry)
        double log_posterior_p = log_like_p + log_eta_prior_p + log_Q_density_given_eta_p;
        double log_posterior_c = log_like_c + log_eta_prior_c + log_Q_density_given_eta_c;

        double log_a_star = log_posterior_p - log_posterior_c;

        // (3) Decide to update
        // [Decide to accept vs. reject]
        double u = R::runif(0, 1);
        if (!Rcpp::NumericVector::is_na(log_a_star)) {
            // consider only when Q_p is reasonably generated
            if (log(u) <= log_a_star) {
                eta_c = eta_p;
                sigma_x_c = sigma_x_p;
            }
        }
    }
    
    return sigma_x_c;
}

// // [[Rcpp::export]]
// vec proposeEta_chunk(uvec ind, vec eta_c, mat Lambda) // checked //
// {
//     // Subsetting Sigma_beta //
//     vec rowind = zeros<vec>(Lambda.n_rows);
//     vec colind = zeros<vec>(Lambda.n_cols);
//     rowind.elem(ind).fill(1);
//     colind.elem(ind).fill(1);
//     mat Sigma11 = Lambda.submat(ind, ind); 
//     mat Sigma12 = Lambda.submat(arma::find(rowind==1), arma::find(colind==0));
//     mat Sigma21 = Sigma12.t();
//     mat Sigma22 = Lambda.submat(arma::find(rowind==0), arma::find(colind==0));
//     mat Sigma22_inv = inv(Sigma22);

//     // conditional parameters //
//     vec cond_mean = eta_c.elem(ind);
//     mat cond_var = Sigma11 - Sigma12 * Sigma22_inv * Sigma21;

//     // sample a chunk //
//     vec eta_p_ind = sampleZ(cond_mean, cond_var);

//     return eta_p_ind;
// }

// // [[Rcpp::export]]
// double LogEtaPrior(vec eta, vec nu, double lambda_eta, mat Sigma_inv, double log_Sigma_det) // checked //
// {
//     int n = eta.size();

//     mat Lambda_inv = Sigma_inv/lambda_eta;
//     double log_lambda_det = n * log(lambda_eta) + log_Sigma_det;
    
//     double log_density = logMVNdensity(eta, nu, Lambda_inv, log_lambda_det);
//     return log_density;
// }

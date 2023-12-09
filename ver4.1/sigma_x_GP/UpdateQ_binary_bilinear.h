#include <iostream>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Last revised: 12/4/2023 
// need to check the output values

// [[Rcpp::export]]
mat CondVar(uvec ind, mat Sigma)
{
    // Subsetting Sigma //
    vec rowind = zeros<vec>(Sigma.n_rows);
    vec colind = zeros<vec>(Sigma.n_cols);
    rowind.elem(ind).fill(1);
    colind.elem(ind).fill(1);
    mat Sigma11 = Sigma.submat(ind, ind); 
    mat Sigma12 = Sigma.submat(arma::find(rowind==1), arma::find(colind==0));
    mat Sigma21 = Sigma12.t();
    mat Sigma22 = Sigma.submat(arma::find(rowind==0), arma::find(colind==0));
    mat Sigma22_inv = inv(Sigma22);

    mat cond_var = Sigma11 - Sigma12 * Sigma22_inv * Sigma21;
    return cond_var;
}

// [[Rcpp::export]]
List PyramidProposal(uvec ind, mat CondMu, mat CondSigma, List quantile_levels, vec all_quantiles, double alpha_scale = 1)
{
    // propose whole quantile pyramids (uniform scale) at a time (default is when ind is a scalar) //
    int n = CondMu.n_cols;
    vec Q_left = zeros<vec>(n);
    vec Q_right = ones<vec>(n);
    mat Q_mat = join_rows(Q_left, Q_right).t();
    mat Z_mat = join_rows(Q_left, Q_right).t(); // just a placeholder

    vec quantiles_estimated = {0, 1};

    int quantile_levels_length = quantile_levels.size();
    for(int m = 0; m < quantile_levels_length; m++) {
        vec quantile_targets = as<vec>(quantile_levels[m]);

        int quantile_targets_length = quantile_targets.size();
        for(int j = 0; j < quantile_targets_length; j++) {
            double tau_l = quantiles_estimated[j];
            double tau_r = quantiles_estimated[j+1];
            
            vec Q_left = Q_mat.row(2*j).t(); // row vector
            vec Q_right = Q_mat.row(2*j+1).t(); // row vector
            
            double tau = quantile_targets[j];
            double a1 = (tau - tau_l) * pow(m + alpha_scale, 2);
            double a2 = (tau_r - tau) * pow(m + alpha_scale, 2);

            int tau_ind = arma::as_scalar(find(all_quantiles == tau)); 
            vec mu = CondMu.row(tau_ind).t(); // mu vector for corresponding quantile (8/31/2023)

            vec Z = sampleZ(mu, CondSigma);
            vec U = sampleU(Z, mu, CondSigma);
            vec V = sampleV(U, a1, a2);
            vec Q = sampleQ(Q_left, Q_right, V); // column vector (default)

            int nr = Q_mat.n_rows;
            Q_mat = join_cols(Q_mat.rows(0,2*j), Q.t(), Q_mat.rows(2*j+1,nr-1));
            Z_mat = join_cols(Z_mat.rows(0,2*j), Z.t(), Z_mat.rows(2*j+1,nr-1));
        }
        quantiles_estimated = sort(join_cols(quantiles_estimated, quantile_targets));
    }
    // NumericMatrix Q_mat_rcpp = clone(wrap(Q_mat));
    // CharacterVector Q_names = clone(wrap(quantiles_estimated));
    // rownames(Q_mat_rcpp) = Q_names;
    // List outputList = List::create(Named("Q") = Q_mat_rcpp, Named("Z") = Z_mat);

    return List::create(Named("Q") = Q_mat, Named("Z") = Z_mat);
}


// [[Rcpp::export]]
mat update_columns(mat DQP, uvec ind, mat Q_p){

    for(uword i=0; i<ind.n_elem; i++){
        DQP.col(ind(i)) = Q_p.col(i); // replace each column (entire pyramid)
    }
    return DQP;
}


// [[Rcpp::export]]
List UpdateQ(mat DQP_c_y, mat DQP_Z_c, vec beta, vec sigma_x, vec sigma_x_long, List quantile_levels, vec all_quantiles, double alpha_scale, 
             vec y, mat X, mat X_prime, vec group_ind, mat X_QP, uvec locs, uvec alocs, // ** new **
             vec mu, mat Sigma, mat Sigma_inv, double log_Sigma_det, 
             mat Sigma_p, double alpha_scale_p, int n_chunk=5)
{
    // [input: mat DQP_c (y scale), all_quantiles] //
    // vec all_quantiles = as_numeric(rownames(DQP_c_rcpp)); // we get all_quantiles from rownames
    
    vec mu_x = X_QP * beta.elem(alocs);           // on QP locations
    vec mu_x_long = X * beta;               // on all locations

    // transform Q to uniform scale
    mat DQP_c(DQP_c_y.n_rows, DQP_c_y.n_cols);
    for (uword r = 0; r < DQP_c_y.n_rows; r++) {
        DQP_c.row(r) = pnorm_cpp(DQP_c_y.row(r).t(), mu_x, sigma_x).t();
    }

    vec quantiles_estimated = {0, 1};
    uvec ind_chunks = chunk_indices(DQP_c.n_cols, n_chunk);

    for (uword i = 0; i < ind_chunks.n_elem-1; i++){
        uvec ind = regspace<uvec>(ind_chunks(i), ind_chunks(i+1)-1);
        mat Q_cu = DQP_c.cols(ind);
        mat Z_c = DQP_Z_c.cols(ind);

        // (i) Propose a new candidate (Pyramid)
        mat condmu_p = Z_c;
        mat condSigma = CondVar(ind, Sigma_p); // common for current/proposed pyramids
        mat condSigma_inv = inv(condSigma);
        double log_condSigma_det = log(det(condSigma));

        List DQP_QZ = PyramidProposal(ind, condmu_p, condSigma, quantile_levels, all_quantiles, alpha_scale_p);
        mat Q_pu = as<arma::mat>(DQP_QZ["Q"]);      // proposed Q matrix (uniform scale)
        mat Z_p = as<arma::mat>(DQP_QZ["Z"]);       // corresponding Z matrix

        mat DQP_p = update_columns(DQP_c, ind, Q_pu);
        mat DQP_Z_p = update_columns(DQP_Z_c, ind, Z_p);
        mat condmu_c = Z_p;
        
        // (ii) Calculate the ratio of conditional proposal density
        // this is a joint density (instead of calculating 1-dim conditional density)
        vec mu_x_ind = mu_x(ind);
        vec sigma_x_ind = sigma_x(ind);
        double LogQdensity_c = LogCondDensityBinaryDQP(Q_cu, quantile_levels, all_quantiles, mu_x_ind, sigma_x_ind, condmu_c, condSigma, condSigma_inv, log_condSigma_det, alpha_scale_p);
        double LogQdensity_p = LogCondDensityBinaryDQP(Q_pu, quantile_levels, all_quantiles, mu_x_ind, sigma_x_ind, condmu_p, condSigma, condSigma_inv, log_condSigma_det, alpha_scale_p);

        // (iii) Calculate the (log)prior densities
        double logprior_c = DQPlogPrior_binary(DQP_c, quantile_levels, all_quantiles, mu_x, sigma_x, mu, Sigma, Sigma_inv, log_Sigma_det, alpha_scale);
        double logprior_p = DQPlogPrior_binary(DQP_p, quantile_levels, all_quantiles, mu_x, sigma_x, mu, Sigma, Sigma_inv, log_Sigma_det, alpha_scale);

        mat DQP_c_long = bilinearInterpolator_mat(X_QP.cols(1,2), DQP_c, X_prime.cols(locs)); // interpolating in uniform scale based on 2 columns
        mat DQP_p_long = bilinearInterpolator_mat(X_QP.cols(1,2), DQP_p, X_prime.cols(locs)); // interpolating in uniform scale based on 2 columns
        
        // (iv) Calculate the (log)likelihoods (interpolated Q in uniform scale)
        double loglike_c = LogLikeC(DQP_c_long, all_quantiles, y, mu_x_long, sigma_x_long, group_ind);
        double loglike_p = LogLikeC(DQP_p_long, all_quantiles, y, mu_x_long, sigma_x_long, group_ind);

        // (v) Calculate the (log) acceptance probability
        double log_a_star = (loglike_p + logprior_p) - (loglike_c + logprior_c) + (LogQdensity_c - LogQdensity_p);

        // [Decide to accept vs. reject]
        double u = R::runif(0, 1);
        if (!Rcpp::NumericVector::is_na(log_a_star)) {
            // consider only when Q_p is reasonably generated
            if (log(u) <= log_a_star) {
                DQP_c = DQP_p;
                DQP_Z_c = DQP_Z_p;
                // Rcpp::Rcout << "accepted..." << std::endl;
            }
        }
    }

    // transform back to y scale //
    for (uword r = 0; r < DQP_c.n_rows; r++) {
        DQP_c_y.row(r) = qnorm_cpp(DQP_c.row(r).t(), mu_x, sigma_x).t();
    }

    // [output Q =  arma::mat, Z = arma::mat]
    return List::create(Named("Q") = DQP_c_y, Named("Z") = DQP_Z_c);
}


// // [[Rcpp::export]]
// List proposeQ_chunk(uvec ind, vec Z_c, mat Sigma, double delta, double a1, double a2, vec Q_Lu, vec Q_Ru)
// {
//     // vector Z_c = normal scale vector related to current quantile vector //
//     // vector Q_Lu = left-endpoints of Q given its neighbors (uniform scale) //
//     // vector Q_Ru = right-endpoints of Q given its neighbors (uniform scale) //

//     // Subsetting Sigma //
//     vec rowind = zeros<vec>(Sigma.n_rows);
//     vec colind = zeros<vec>(Sigma.n_cols);
//     rowind.elem(ind).fill(1);
//     colind.elem(ind).fill(1);
//     mat Sigma11 = delta * Sigma.submat(ind, ind); 
//     mat Sigma12 = delta * Sigma.submat(arma::find(rowind==1), arma::find(colind==0));
//     mat Sigma21 = Sigma12.t();
//     mat Sigma22 = delta * Sigma.submat(arma::find(rowind==0), arma::find(colind==0));
//     mat Sigma22_inv = inv(Sigma22);

//     // conditional parameters //
//     vec cond_mean = Z_c.elem(ind); // + Sigma12 %*% Sigma22_inv %*% (Z_p_given - mu2) // = always zero
//     mat cond_var = Sigma11 - Sigma12 * Sigma22_inv * Sigma21;

//     // sample a chunk //
//     vec Z_p_ind = sampleZ(cond_mean, cond_var);

//     // calculate U, V, Q chunks //
//     vec U_p_ind = sampleU(Z_p_ind, cond_mean, cond_var); // element-wise
//     vec V_p_ind = sampleV(U_p_ind, a1, a2);
//     vec Q_pu_ind = (1-V_p_ind) % Q_Lu.elem(ind) + V_p_ind % Q_Ru.elem(ind);

//     return List::create(Named("Q_pu_ind") = Q_pu_ind, Named("Z_p_ind") = Z_p_ind);
// }


// // [[Rcpp::export]]
// double LogQdensity(vec Q_pu, vec Q_Lu, vec Q_Ru, vec Z_p, vec Z_c, mat Sigma, vec mu_x, vec sigma_x,
//                    double delta, double a1, double a2, mat Sigma_inv, double log_Sigma_det) {

//     vec V_p = (Q_pu - Q_Lu) / (Q_Ru - Q_Lu);
//     double logJ = LogJacobianC(Q_pu, Q_Lu, Q_Ru, mu_x, sigma_x, V_p, a1, a2, Z_p, Z_c, delta*Sigma);
//     mat Sigma_inv_p = Sigma_inv/delta;
//     double log_Sigma_det_p = Z_p.n_elem * log(delta) + log_Sigma_det;
//     double logq = logMVNdensity(Z_p, Z_c, Sigma_inv_p, log_Sigma_det_p) + logJ;

//     return logq;
// }


// // [[Rcpp::export]]
// mat update_row_part(mat DQP_p, int tau_ind, vec Q_pu, uvec ind){

//     rowvec DQP_p_row = DQP_p.row(tau_ind);  // step1. find tau row
//     DQP_p_row(ind) = Q_pu(ind);             // step2. replace the corresponding elements
//     DQP_p.row(tau_ind) = DQP_p_row;         // step3. update the entire matrix
//     return DQP_p;
// }
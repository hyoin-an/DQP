library(Rcpp)
sourceCpp('MCMC_binary_3param.cpp')


# *************************** #
#     DISCRETE/UNIVARIATE X   #
# *************************** #

### Data generation ###

n = 300
rep_seq = sort(sample(1:25, n, replace=T))
set.seed(156); x1 = 1:5;  x2 = 1:5
x = as.matrix(expand.grid(x1=x1, x2=x2))[rep_seq,]
y = x[,1] + x[,2] + rnorm(n) 


### hyper parameters ###
# mu and Sigma for Gaussian process
#   (length of mu and dimension of Sigma determine the number of quantile pyramids)

x_prime = unique_rows(x)

mu = rep(0, nrow(x_prime))
tausq = 0 # zero nugget effect for smooth covariance function
sigsq = 1 # we want the correlation matrix
r = 1 # exponential
phi = 2.5 # for discrete chose the distance parameter so the nearby points would have meaningful dependence
Sigma = sigsq * exp(-abs(as.matrix(dist(x_prime)))^r/phi) + diag(rep(tausq), nrow(x_prime))

ols_res = lm(y ~ x[,1] + x[,2])
mu_0 = ols_res$coef
Sigma_0 = diag(summary(ols_res)$coef[,2])
mu_1 = rep(0, nrow(x_prime))
Sigma_1 = diag(rep(0.1, nrow(x_prime)))


# proposal parameters
delta = 1
Sigma_beta = Sigma_beta = vcov(ols_res)
Sigma_eta = Sigma_1/3



### Quantile pyramid structure specification ###

quantile_levels = list(c(0.5), c(0.25, 0.75))

# initial values
beta_initial = mu_0
mu_x_initial = cbind(1, x_prime) %*% beta_initial
sigma_x_initial = exp(mu_1)


DQP_res = DQPbinarySampling(mu, Sigma, quantile_levels, alpha_scale=5)
DQP_initial_unit = DQP_res$Q
DQP_Z_initial = DQP_res$Z
DQP_Q_initial = DQP_initial_unit
for(r in 1:nrow(DQP_Q_initial)){
  DQP_Q_initial[r,] = qnorm(DQP_initial_unit[r,], mu_x_initial, sigma_x_initial)
}


### Run MCMC ###


### Run MCMC ###

niter = 300000
set.seed(957)
tic = Sys.time()
## [Rcpp version]
# QP_grid: locations on the covariate space where QPs are constructed
res = MCMC_param3(niter, y, as.matrix(x), quantile_levels, mu, Sigma, delta, alpha_scale=5, 
                  mu_0, Sigma_0, Sigma_beta, mu_1, Sigma_1, Sigma_eta, 
                  n_chunk_Q=nrow(x_prime), n_chunk_beta=1, n_chunk_eta=1,
                  DQP_Q_initial, DQP_Z_initial, beta_initial, sigma_x_initial)
toc = Sys.time()
toc - tic


save(res, file = 'ver1_twodim1_m2_nchunkQ25.RData')


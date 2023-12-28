# DQP

A process of dependent quantile pyramids (DQP) is a Bayesian nonparametric approach that offers a flexible modeling framework for the simultaneous estimation of multiple quantile regression curves. 

This repository contains code for implementing Markov Chain Monte Carlo (MCMC) simulations for simultaneous quantile regression using a binary DQP.
To accommodate diverse modeling needs with DQP, various code versions have been developed.

Version 1 is the code for the canonical version of DQP. More examples and results to come for the other versions. 


### ver 1. Canonical DQP
This is the simplest version of DQP. Conditional and dependent quantile pyramids are constructed at all available data points.

*	Normal tail for transformation
*	Linear regression model for mu_x
*	Log-linear regression model for sigma_x
*	Nonparametric model for sigma_x with a log-normal proposal with dependence

To use DQP, hyperparameters and initial values of the parameters are required. Consider the following example:
```
n = 300
rep_seq = sort(sample(1:25, n, replace=T))
set.seed(156); x1 = 1:5;  x2 = 1:5
x = as.matrix(expand.grid(x1=x1, x2=x2))[rep_seq,]
y = x[,1] + x[,2] + rnorm(n) 
```

First, hyperparameters for the Gaussian process need to be specified. The distance parameter `phi` decides the degree of correlation between the data points. 

```
# Hyperparameters for the Gaussian process
x_prime = unique_rows(x)
mu = rep(0, nrow(x_prime))
tausq = 0 
sigsq = 1
r = 1 
phi = 2.5 
Sigma = sigsq * exp(-abs(as.matrix(dist(x_prime)))^r/phi) + diag(rep(tausq), nrow(x_prime))
```

Also, the hyperparameters for prior distributions can be specified as below.
```
ols_res = lm(y ~ x[,1] + x[,2])
mu_0 = ols_res$coef
Sigma_0 = diag(summary(ols_res)$coef[,2])
mu_1 = rep(0, nrow(x_prime))
Sigma_1 = diag(rep(0.1, nrow(x_prime)))
```

Next, the parameters for the proposal distribution should be pre-determined.
```
Sigma_beta = Sigma_beta = vcov(ols_res)
Sigma_eta = Sigma_1/3
```

Lastly, the initial values for the parameters need to be chosen. 
```
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
```

Then, MCMC run can be performed using the following code. 
```
quantile_levels = list(c(0.5), c(0.25, 0.75))
res = MCMC_param3(niter, y, as.matrix(x), quantile_levels, mu, Sigma, delta, alpha_scale=5, 
                  mu_0, Sigma_0, Sigma_beta, mu_1, Sigma_1, Sigma_eta, 
                  n_chunk_Q=nrow(x_prime), n_chunk_beta=1, n_chunk_eta=1,
                  DQP_Q_initial, DQP_Z_initial, beta_initial, sigma_x_initial)
```

An illustrative example of detailed code usage can be found [here](ver1/Example).
Simulation studies and an illustrative application of this version of DQP can be referenced in the paper below.

*Reference*: An, H., & MacEachern, S. N. (2023). A Process of Dependent Quantile Pyramids. arXiv preprint arXiv:2306.02126. [[link]](https://doi.org/10.48550/arXiv.2306.02126)



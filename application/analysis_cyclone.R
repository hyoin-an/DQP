Rcpp::sourceCpp("ver0/MCMC_binary.cpp")


### Data set up ###
dat <- read.table("data/globalTCmax4.txt") # 2098 obs, over the globe, 1981-2006 (test data)
dat <- dat[is.na(dat$Basin),] # NA = North Atlantic

y <- dat$WmaxST
x <- dat$Year

x_centered <- x - mean(x)
x_prime <- unique(x)

##### set up #####

# mu and Sigma for Gaussian process
mu <- rep(0, length(x_prime))
tausq <- 0  # zero nugget effect for smooth covariance function
sigsq <- 1  # we want a correlation matrix
r <- 1      # exponential
phi <- 5    # for discrete chose the distance parameter so the nearby points would have meaningful dependence
Sigma <- sigsq * exp(-abs(as.matrix(dist(x_prime)))^r/phi) + diag(rep(tausq), length(x_prime))
# View(Sigma)

# beta prior & proposal parameters
mu_0 <- c(71, 0.5)
Sigma_0 <- matrix(c(15, 0, 0, 2), ncol = 2)
ols_res <- lm(y~x_centered)
Sigma_beta <- solve((t(cbind(1, x_centered))) %*% as.matrix(cbind(1, x_centered))) * summary(ols_res)$sigma^2

# sigma_x
sigma_x <- aggregate(y, by=list(x), FUN=sd)[,2] # empirical sd vector
sigma_x <- lm(sigma_x ~ x_prime)$fit



# Quantile pyramid structure specification
quantile_levels = list(c(0.5), c(0.25, 0.75), c(0.1, 0.35, 0.65, 0.9), c(0.05, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.95))

# initial values
beta_initial <- ols_res$coefficients
mu_x_initial <- cbind(1, unique(x_centered)) %*% beta_initial
DQP_res <- DQPbinarySampling(mu, Sigma, quantile_levels, alpha_scale=5)
DQP_initial_unit <- DQP_res$Q
DQP_Z_initial <- DQP_res$Z
DQP_Q_initial <- DQP_initial_unit
for(r in 1:nrow(DQP_Q_initial)){
  DQP_Q_initial[r,] <- qnorm(DQP_initial_unit[r,], mu_x_initial, sigma_x)
}


##### Run MCMC #####

niter = 210000
set.seed(803)
res <- MCMC_param2(niter, y, as.matrix(x_centered), sigma_x, delta=1, Sigma_beta, 
                   DQP_Q_initial, DQP_Z_initial, beta_initial, quantile_levels, 
                   mu, Sigma, mu_0, Sigma_0, alpha_scale=5, n_chunk_Q=13, n_chunk_beta=1)

save(res, file = "cyclone_result.RData")



# Load the data
load("ver1_1_twodim1_m2_nchunkQ25.RData")

quantile_levels = list(c(0.5), c(0.25, 0.75))
quantiles_of_interest = sort(unlist(quantile_levels))
n_quantiles = length(quantiles_of_interest)


# burning/thinning
niter = length(res$beta)
nburn = 10000
nthin = 100
ind = seq(from = nburn + 1, to = niter, by = nthin)


# acceptance ratio of proposals
accept_ratio_Q = NULL
for(t in 1:n_quantiles){
  accept_Q = NULL
  for(i in 1:(niter-1)){
    accept_Q[i] = !all(round(res$Q[[i]][t+1,],6) == round(res$Q[[i+1]][t+1,],6))
  }
  accept_ratio_Q[t] = mean(accept_Q)
}
accept_ratio_Q 


accept_beta = NULL
for(i in 1:(niter-1)){
  accept_beta[i] = !all(res$beta[[i]] == res$beta[[i+1]])
}
accept_ratio_beta = mean(accept_beta)
accept_ratio_beta 


accept_sigma_x = NULL
for(i in 1:(niter-1)){
  accept_sigma_x[i] = !all(res$sigma_x[[i]] == res$sigma_x[[i+1]])
}
accept_ratio_sigma_x = mean(accept_sigma_x)
accept_ratio_sigma_x



## Posterior samples 
Q_post = list()
for(j in 1:(length(quantiles_of_interest))){
  Q_post_each = matrix(NA, ncol=ncol(res$Q[[1]]), nrow=3*niter)
  for(i in 1:niter){
    Q_post_each[i,] = res$Q[[i]][j+1,]
  }
  Q_post[[j]] = Q_post_each
}

beta_post = matrix(NA, nrow=nrow(res$beta[[1]]), ncol=niter)
for(i in 1:niter){
  for(j in 1:nrow(res$beta[[1]])){
    beta_post[j,i] = res$beta[[i]][j,]
  }
}

sigma_x_post = matrix(NA, nrow=nrow(res$sigma_x[[1]]), ncol=niter)
for(i in 1:niter){
  for(j in 1:nrow(res$sigma_x[[1]])){
    sigma_x_post[j,i] = res$sigma_x[[i]][j,]
  }
}


## Trace plots
for(t in 1:length(Q_post)){
  par(mfrow=c(5, 5), mar=c(4, 2, 1, 1))
  for(i in 1:ncol(res$Q[[1]])){
    plot(ind, Q_post[[t]][ind,i], xlab=paste("Q(",quantiles_of_interest[t],") ",i,sep=""), cex=0.5)
  }
}

par(mfrow=c(1, 3), mar=c(4, 2, 1, 1))
for(j in 1:nrow(res$beta[[1]])){
  plot(ind, beta_post[j,ind], cex=0.5, xlab=paste("beta",j-1,sep=""))
}

par(mfrow=c(5, 5), mar=c(4, 2, 1, 1))
for(j in 1:nrow(res$sigma_x[[1]])){
  plot(ind, sigma_x_post[j,ind], cex=0.5, xlab=paste("sigma_x",j-1,sep=""))
}



### Recall: Data generation ###
n = 300
rep_seq = sort(sample(1:25, n, replace=T))
set.seed(156); x1 = 1:5;  x2 = 1:5
x = as.matrix(expand.grid(x1=x1, x2=x2))[rep_seq,]
y = x[,1] + x[,2] + rnorm(n) # Fully Linear


Rcpp::sourceCpp('MCMC_binary_3param.cpp')
x_prime = unique_rows(x)

true_q = data.frame(x = x_prime)
for(l in quantiles_of_interest){
  true_q = data.frame(true_q, q = x_prime[,1] + x_prime[,2] + qnorm(l))
}


# Fitted plots
par(mfrow=c(1,2), mar=c(4, 4, 1, 1))
for(t in 1:length(quantiles_of_interest)){
  y_quantiles = Q_post[[t]][ind,]
  plot(x[,1], y, pch=16, col='grey', cex=0.8) # 1-dim plot (x1)
  points(true_q[,1], true_q[,t+2], col='red', pch=16, cex=2)
  points(true_q[,1], colMeans(y_quantiles), col = 'black', pch=16, cex=2)
  
  y_quantiles = Q_post[[t]][ind,]
  plot(x[,2], y, pch=16, col='grey', cex=0.8) # 1-dim plot (x2)
  points(true_q[,2], true_q[,t+2], col='red', pch=16, cex=2)
  points(true_q[,2], colMeans(y_quantiles), col = 'black', pch=16, cex=2)
}


### Linearized fit
par(mfrow=c(1,2), mar=c(4, 4, 1, 1))
for(t in 1:length(quantiles_of_interest)){
  qs = colMeans(Q_post[[t]][ind,])
  qs_linear = lm(qs ~ x_prime[,1] + x_prime[,2])
  plot(x[,1], y, pch=16, col='grey', cex=0.8) # 1-dim plot (x1)
  points(true_q[,1], true_q[,t+2], col='red', pch=16, cex=2)
  points(true_q[,1], qs_linear$fitted.values, col = 'black', pch=16, cex=2)
  
  y_quantiles = Q_post[[t]][ind,]
  plot(x[,2], y, pch=16, col='grey', cex=0.8) # 1-dim plot (x2)
  points(true_q[,2], true_q[,t+2], col='red', pch=16, cex=2)
  points(true_q[,2], qs_linear$fitted.values, col = 'black', pch=16, cex=2)
}





library(ggplot2)
library(tidyverse)


load('cyclone_result.RData')

quantile_levels = list(c(0.5), c(0.25, 0.75), c(0.1, 0.35, 0.65, 0.9), c(0.05, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.95)) 
quantiles_of_interest = sort(unlist(quantile_levels))
n_quantiles = length(quantiles_of_interest)

niter = length(res$Q)

# burning/thinning
nburn = 10000
nthin = 100
ind = seq(from = nburn + 1, to = niter, by = nthin)


DQP_post = list()

for(j in 1:(length(quantiles_of_interest))){
  DQP_post_each = matrix(NA, ncol=ncol(res$Q[[1]]), nrow=niter)
  for(i in 1:niter){
    DQP_post_each[i,] = res$Q[[i]][j+1,]
  }
  DQP_post[[j]] = DQP_post_each
}

## Trace plots
quantiles_of_interest

for(j in 1:length(quantiles_of_interest)){
  par(mfrow=c(3, 9), mar=c(4, 2, 1, 1))
  for(i in 1:ncol(res$Q[[1]])){
    plot(ind, DQP_post[[j]][ind,i], xlab=i, cex=0.5)
  }
}

beta0_post = NULL
beta1_post = NULL
for(i in 1:niter){
  beta0_post[i] = res$beta[[i]][1]
  beta1_post[i] = res$beta[[i]][2]
}

par(mfrow=c(1,2), mar=c(4, 2, 1, 1))
plot(ind, beta0_post[ind], cex=0.5, xlab='beta0')
plot(ind, beta1_post[ind], cex=0.5, xlab='beta1')
mean(beta0_post[ind]); mean(beta1_post[ind])




### Recall: Data generation ###
dat = read.table('data/globalTCmax4.txt') # 2098 obs, over the globe, 1981-2006 (test data)
dat = dat[is.na(dat$Basin),] # NA = North Atlantic


### Data set up ###
y = dat$WmaxST
x = dat$Year
n = length(y)


par(mfrow=c(4,4), mar=c(4, 4, 1, 1))
for(t in 1:(length(quantiles_of_interest))){
  y_quantiles = DQP_post[[t]][ind,]
  boxplot(y_quantiles, cex = 0.5, names = unique(x), ylim=c(25, 160),
          xlab = paste(quantiles_of_interest[t],"Quantile"), ylab='WmaxST')
}

plot(x, y, pch=4, cex=0.8, col='grey', xlab="Year", ylab="WmaxST")
for(t in 1:(length(quantiles_of_interest))){
  y_quantiles = DQP_post[[t]][ind,]
  lines(unique(x), colMeans(y_quantiles), lwd=2)
}



#### Figures ####

# (1) Fitted plots

Qs = Qs_lm = matrix(NA, nrow=ncol(DQP_post[[1]]), ncol=length(quantiles_of_interest))
for(t in 1:ncol(Qs)){
  y_quantiles = colMeans(DQP_post[[t]][ind,])
  Qs[,t] = y_quantiles
  Qs_lm[,t] = lm(y_quantiles ~ unique(x))$fitted
}
colnames(Qs) = colnames(Qs_lm) = quantiles_of_interest
Qs = as.data.frame(cbind(Year=unique(x), Qs))
Qs_lm = as.data.frame(cbind(Year=unique(x), Qs_lm))


fit_DQP_long = reshape2::melt(Qs, id.vars = "Year", variable.name = "Level", value.name="Quantile")
head(fit_DQP_long)

fit_DQP_lm_long = reshape2::melt(Qs_lm, id.vars = "Year", variable.name = "Level", value.name="Quantile")
head(fit_DQP_lm_long)

gg1 = ggplot() + geom_point(aes(x, y), col="grey", size=0.5) + 
  geom_line(data = fit_DQP_long, aes(Year, Quantile, group=Level), col="grey") + 
  geom_line(data = fit_DQP_lm_long, aes(Year, Quantile, group=Level)) + 
  theme_bw() + labs(x="Year", y="WmaxST", title="(a) DQP")


# (2) 95% Credible intervals of beta1

# linearized coefficients
beta1_linearized = matrix(NA, nrow=length(ind), ncol=(length(quantiles_of_interest)))
for(t in 1:(length(quantiles_of_interest))){
  y_quantiles = DQP_post[[t]][ind,]
  for(i in 1:length(ind)){
    beta1_linearized[i,t] = coef(lm(y_quantiles[i,] ~ unique(x)))[2]
  }
}

q025 = apply(beta1_linearized, 2, function(x){quantile(x, 0.025)})
q975 = apply(beta1_linearized, 2, function(x){quantile(x, 0.975)})
beta1 = data.frame(Level=quantiles_of_interest, Slope=colMeans(beta1_linearized), q025=q025, q975=q975)

gg2 = ggplot(beta1, aes(Level, Slope)) + geom_point(size=4) + 
  geom_errorbar(aes(ymin=q025, ymax=q975)) + ylim(-0.3, 2.6) + 
  theme_bw() + geom_hline(aes(yintercept=0), linetype=2) + labs(title="95% Credible intervals")

gridExtra::grid.arrange(gg1, gg2, ncol=2)



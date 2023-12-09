# DQP

A process of dependent quantile pyramids (DQP) is a Bayesian nonparametric approach that offers a flexible modeling framework for the simultaneous estimation of multiple quantile regression curves. 

This repository contains code for implementing Markov Chain Monte Carlo (MCMC) simulations for simultaneous quantile regression using DQP.
To accommodate diverse modeling needs with DQP, various code versions have been provided and uploaded.

### ver 1. Canonical DQP
This is the simplest version of DQP. Conditional and dependent quantile pyramids are constructed at all available data points.

*	Normal tail for transformation
*	Linear regression model for mu_x
*	Log-linear regression model for sigma_x
*	Nonparametric model for sigma_x with a log-normal proposal with dependence


### ver2.0 One-dimensional DQP with linear interpolation

##### Log-linear modeling for sigma_x
* One-dimensional DQP
* Linear interpolation across X
* Normal tail for transformation
* Linear regression model for mu_x
* Log-linear regression model for sigma_x
* Pyramid proposal

##### Nonparametric modeling for sigma_x
* Two-dimensional DQP
* Linear interpolation across X
* Normal tail for transformation
* Linear regression model for mu_x
* Nonparametric model for sigma_x with a log-normal proposal with dependence
* Pyramid proposal

### ver3.0 Two-dimensional DQP with bilinear interpolation

##### Log-linear modeling for sigma_x
* Two-dimensional DQP
* Bilinear interpolation across X
* Normal tail for transformation
* Linear regression model for mu_x
* Log-linear regression model for sigma_x
* Pyramid proposal

##### Nonparametric modeling for sigma_x
* Two-dimensional DQP
* Bilinear interpolation across X
* Normal tail for transformation
* Linear regression model for mu_x
* Nonparametric model for sigma_x with a log-normal proposal with dependence
* Pyramid proposal


### ver4.0 Location/Scale family of a One-dimensional DQP

##### Log-linear modeling for sigma_x
*	Extension of One dimensional DQP (ver2.0)
*	Linear interpolation across X
*	Normal tail for transformation
*	Linear regression model for mu_x
*	Log-linear regression model for sigma_x
*	Pyramid proposal for quantiles
*	One predictor is designated for whose space QPs are constructed; Location/scale shift for all the other predictors
  
##### Nonparametric modeling for sigma_x
*	Extension of One dimensional DQP (ver2.1)
*	Linear interpolation across X
*	Normal tail for transformation
*	Linear regression model for mu_x
*	Nonparametric model for sigma_x with a log-normal proposal with dependence
*	Pyramid proposal for quantiles
*	One predictor is designated for whose space QPs are constructed; Location/scale shift for all the other predictors


### ver4.1 Location/Scale family of a Two-dimensional DQP

##### Log-linear modeling for sigma_x
* Extension of One dimensional DQP (ver3.1)
* Bilinear interpolation across X
* Normal tail for transformation
* Linear regression model for mu_x
* Log-linear regression model for sigma_x
* Pyramid proposal for quantiles
* Two predictors are designated for whose space QPs are constructed; Location/scale shift for all the other predictors


##### Nonparametric modeling for sigma_x
* Extension of One dimensional DQP (ver3.1)
* Bilinear interpolation across X
* Normal tail for transformation
* Linear regression model for mu_x
* Nonparametric for sigma_x with a log-normal proposal with dependence
* Pyramid proposal for quantiles
* Two predictors are designated for whose space QPs are constructed; Location/scale shift for all the other predictors


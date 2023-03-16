# STAT 597
# Homework 3
# Veronica A. Winter
# 03/16/2023

# clean env
rm(list = ls())
gc()
ls()

# set for reproducibility
set.seed(1234)

# Libraries
library(coda)
library(mvtnorm)
library(maptools)
library(maps)
library(mgcv)
library(glmnet)
library(MASS)
#library(nimble)

# Source fun
source("Homework3/mcmc.lm.lasso.r")

# Question 1 ----
# Let 𝑦𝑖 be the number of sick days that a person takes due to an illness, 
# Let 𝑥𝑖 be the number of months that person has been taking part in a 
# treatment program.  

# Data model : yi ~ Pois(exp{a + bxi})

# a Poisson regression model with two regression parameters 
# a and b.  
# 
# We will assume that both of these regression parameters have prior distributions 
# that are Gaussian with mean zero and variance equal to 4. 

x=c(8,14,11,7,32,8,28,21,27,15,26,13,19,22,15,
    12,15,7,9,15,26,22,16,12,6)

y=c(5,2,5,4,1,3,0,2,1,2,2,5,3,2,1,2,2,8,5,2,1,1,6,4,3)

# a. Construct a Metropolis Hasting sampler that jointly proposes (a,b) from a 
# bivariate normal distribution with the current values of both parameters as the mean.  
# Begin with a proposal distribution that has a diagonal covariance matrix with 
# 0.01 on both diagonal elements. 

# Use Shaby and Wells’ log-adaptive tuning approach to adaptively tune the 
# proposal distribution, with adaptation happening every 100 MCMC iterations.  
# Run this until you are sure that your algorithm converges to the stationary 
# distribution.

# Load in Question 1 data
x=c(8,14,11,7,32,8,28,21,27,15,26,13,19,22,15,
    12,15,7,9,15,26,22,16,12,6)

y=c(5,2,5,4,1,3,0,2,1,2,2,5,3,2,1,2,2,8,5,2,1,1,6,4,3)

## Define the log-likelihood function
loglike <- function(a, b, x, y) {
  return(sum(dpois(y, exp(a + b*x), log=TRUE)))
}

# Set up the prior distribution parameters
prior_mean <- c(0, 0)
prior_var <- diag(2)*4

# Set up the proposal distribution
proposal_var <- diag(2)*0.01

# Set up the initial values
a <- 0
b <- 0

# Set up the MCMC algorithm
niter <- 10000
accept <- 0
cov_mat <- proposal_var
tune_interval <- 100
tune_count <- 0
tune_iter <- tune_interval

# Sample storage
samples <- matrix(0, niter, 2)

# Run the MCMC algorithm
for(i in 1:niter) {
  # Generate a proposed value for a and b
  proposal <- mvrnorm(1, c(a, b), cov_mat)
  a_prop <- proposal[1]
  b_prop <- proposal[2]
  
  # Calculate the log-likelihoods for the current and proposed values
  loglike_curr <- loglike(a, b, x, y)
  loglike_prop <- loglike(a_prop, b_prop, x, y)
  
  # Calculate the log-prior densities for the current and proposed values
  logprior_curr <- dmvnorm(c(a, b), prior_mean, prior_var, log=TRUE)
  logprior_prop <- dmvnorm(c(a_prop, b_prop), prior_mean, prior_var, log=TRUE)
  
  # Calculate the log-acceptance ratio
  log_accept_ratio <- loglike_prop + logprior_prop - loglike_curr - logprior_curr
  
  # Accept or reject the proposed value
  if (log(runif(1)) < log_accept_ratio) {
    a <- a_prop
    b <- b_prop
    accept <- accept + 1
  }
  
  # Save the current values of a and b
  samples[i, 1] <- a
  samples[i, 2] <- b
  
  # Update the covariance matrix and tune the proposal distribution
  if (i == tune_iter) {
    # Calculate the acceptance rate
    acceptance_rate <- accept/tune_interval
    
    # Update the covariance matrix
    cov_mat <- cov(samples[(i-tune_interval+1):i,])
    
    # Tune the proposal distribution
    if (acceptance_rate < 0.2) {
      proposal_var <- proposal_var/2
    } else if (acceptance_rate > 0.3) {
      proposal_var <- proposal_var*2
    }
    
    # Reset the acceptance count
    accept <- 0
    
    # Increment the tune count and tune iteration
    tune_count <- tune_count
  }}  

# Print summary of results
cat("Mean of a:", mean(samples[,1]), "\n")
cat("Mean of b:", mean(samples[,2]), "\n")
cat("95% credible interval for a:", quantile(samples[,1], c(0.025, 0.975)), "\n")
cat("95% credible interval for b:", quantile(samples[,2], c(0.025, 0.975)), "\n")

### Accept ----
accept

# Plot results
par(mfrow=c(2,1))
plot(samples[,1], type="l", ylab="a", main="Traceplot of a")
plot(samples[,2], type="l", ylab="b", main="Traceplot of b")

# Question 3 ----
#############################################
##
## #3 Linear Model with Ridge Regression priors
#log(tn[i]) ~ N(mu[i], sigma^2)
#mu[i] = beta0 + beta1log(tp[i]) + beta2log(secchi[i]) + beta3ag[i] + beta4forest[i] + beta5x[i] + beta6y[i] + beta7x[i]^2 + beta8y[i]^2

#where i = 1, 2, ..., n are the lake observations, 
# and beta0, beta1, ..., beta8 are the regression coefficients to be estimated.

#To incorporate ridge regression prior on the regression coefficients, 
#we assume that:
  
#beta1, beta2, ..., beta8 ~ N(0, tau^2)
#beta0 ~ N(0, 100)
#where tau is a hyperparameter to be specified.
#############################################
# Load in data
load("Homework3/lagos.Rdata")
source("Homework3/mcmc.lm.rr.r")
#source("Homework3/mcmc.lm.lasso.r")
source("Homework3/ShabyWellsLogAdapt.r")
library(statmod)

## Create Indicator Vars
lagos$ag=0
lagos$ag[lagos$SurLand=="agricultural"] <- 1
lagos$forest=0
lagos$forest[lagos$SurLand=="forest"] <- 1

## Format the data for 
df <- data.frame(
  #int = rep(1,nrow(lagos)),
  log.tp = log(lagos$tp),
  log.sec = log(lagos$secchi),
  x = lagos$lon,
  y = lagos$lat,
  x2 = lagos$lon^2,
  y2 = lagos$lat^2,
  ag = lagos$ag,
  forest = lagos$forest
)

# Set up response variable
y.lag=lagos[,5]
# There is an NA value here that I think is causing issues:
y.lag[1]
# Setting this to be the mean value (probably not correct)
mean(y.lag)
y.lag[1] <-  951.8362

# Scale data frame
x_scale <- scale(df)

# Run RR mcmc using EH code
lagos.out = mcmc.lm.rr(y = y.lag, 
                       X = x_scale, 
                       n.mcmc = 10000,
                       k2.start = 1,
                       beta.prior.var = 100,
                       s2.prior.var = var(y.lag),
                       k2.exp.rate = 1,
                       # Tuning parameter starting values
                       tune.k2 = 0.01,
                       tune.s2=.01,
                       # adapting every i'th iteration
                       adapt.iter=100)

## trace plots
par(mfrow=c(3,4))
for(k in 2:ncol(lagos.out)){
  plot(lagos.out[,k],main=colnames(lagos.out)[k],type="l")
  abline(v=0,col="red")
}

# Removing first 5000 for burn-in
burnin = 5000
effectiveSize(lagos.out[-c(1:burnin),])

for(k in 2:ncol(lagos.out)){
  hist(lagos.out[-c(1:burnin),k],main=colnames(lagos.out)[k])
  abline(v=0,col="red")
}

# Calculate the means
apply(lagos.out,2,mean)
# Calculate the credible intervals
t(apply(lagos.out, 2, function(x) quantile(x, c(0.025, 0.975))))

# Question 4 (a)-----
# Define the data and prior distributions
y <- c(9, 15, 14, 5, 6, 4, 4, 4, 8, 0, 10, 22, 7, 2, 6, 2, 18, 
       9, 4, 2, 5, 7, 7, 5, 8, 7, 2, 15, 17, 7, 1, 4, 8, 5, 8, 
       9, 25, 6, 6, 4, 22, 3, 2, 5, 3, 4, 8, 4, 17, 14)

# Define the log-likelihood function
loglik <- function(lambda0, lambda1, p, y, n) {
  lambda <- lambda0 + lambda1 * p
  sum(dpois(y, lambda, log = TRUE))
}

# Define the log-prior function for lambda0
logprior_lambda0 <- function(lambda0, shape, rate) {
  dgamma(lambda0, shape, rate, log = TRUE)
}

# Define the log-prior function for lambda1
logprior_lambda1 <- function(lambda1, shape, rate) {
  dgamma(lambda1, shape, rate, log = TRUE)
}

# Set the initial values for the MCMC sampler
lambda0 <- rgamma(1, shape = 1, rate = 0.1)
lambda1 <- rgamma(1, shape = 1, rate = 0.1)
p <- runif(1)

# Set the parameters for the MCMC sampler
num_iter <- 5000
burn_in <- 1000
prop_sd_lambda0 <- 0.1
prop_sd_lambda1 <- 0.1
n <- length(y)
shape <- 1
rate <- 0.1
mean_gamma <- 10
var_gamma <- 100

# Initialize acceptance rates
accept_lambda0 <- 0
accept_lambda1 <- 0

# Run the MCMC sampler
for (i in 2:num_iter) {
  # Update p using Metropolis-Hastings
  p_prop <- rnorm(1, p[i-1], 0.1)
  log_alpha_p <- loglik(lambda0[i-1], lambda1[i-1], p_prop, y, n) - loglik(lambda0[i-1], lambda1[i-1], p[i-1], y, n)
  if (log(runif(1)) < log_alpha_p) {
    p[i] <- p_prop
  } else {
    p[i] <- p[i-1]
  }
  
  # Update lambda0 using Metropolis-Hastings
  lambda0_prop <- rnorm(1, lambda0[i-1], prop_sd_lambda0)
  log_alpha0 <- loglik(lambda0_prop, lambda1[i-1], p[i], y, n) - loglik(lambda0[i-1], lambda1[i-1], p[i], y, n) + 
    logprior_lambda0(lambda0_prop, shape, rate) - logprior_lambda0(lambda0[i-1], shape, rate)
  if (log(runif(1)) < log_alpha0) {
    lambda0[i] <- lambda0_prop
    accept_lambda0 <- accept_lambda0 + 1
  } else {
    lambda0[i] <- lambda0[i-1]
  }
  
  # Update lambda1 using Metropolis-Hastings
  lambda1_prop <- rnorm(1, lambda1[i-1], prop_sd_lambda1)
  log_alpha1 <- loglik(lambda0[i], lambda1_prop, p[i], y, n) - loglik(lambda0[i], lambda1[i-1], p[i], y, n) + 
    logprior_lambda1(lambda1_prop, shape, rate) - logprior_lambda1(lambda1[i-1], shape, rate)
  if (log(runif(1)) < log_alpha1) {
    lambda1[i] <- lambda1_prop
    accept_lambda1 <- accept_lambda1 + 1
  } else {
    lambda1[i] <- lambda1[i-1]
  }
}

# Find the posterior mean and 95% credible intervals for lambda_high
lambda_high <- lambda0 + lambda1
post_mean_lambda_high <- mean(lambda_high)
post_ci_lambda_high <- quantile(lambda_high, c(0.025, 0.975))

# Print the results
cat("Posterior mean of lambda_high: ", post_mean_lambda_high, "\n")
cat("95% credible interval for lambda_high: ", post_ci_lambda_high[1], "-", post_ci_lambda_high[2], "\n")


# Question 4 (b) ----
library(rjags)

# Define the data
y <- c(9, 15, 14, 5, 6, 4, 4, 4, 8, 0, 10, 22, 7, 2, 6, 2, 18, 9, 4, 2, 5, 7, 7, 5, 8, 7, 2, 15, 17, 7, 1, 4, 8, 5, 8, 9, 25, 6, 6, 4, 22, 3, 2, 5, 3, 4, 8, 4, 17, 14)
T <- length(y)

# Define the JAGS model
model_string <- "model{
  # Priors
  lambda_0 ~ dgamma(1, 0.1)
  lambda_1 ~ dgamma(1, 0.1)
  p ~ dunif(0, 1)

  # Data augmentation
  for(t in 1:T) {
    y[t] ~ dpois(nt0[t] + nt1[t])
    nt0[t] ~ dpois(lambda_0)
    nt1[t] ~ dpois(lambda_1 * x[t])
    x[t] ~ dbern(p)
  }
}"

# Prepare the data list for JAGS
data_list <- list(y = y, T = T)

# Set the JAGS parameters
params <- c("lambda_0", "lambda_1", "p", "x")

# Set the number of iterations and burn-in period
n_iter <- 5000
n_burnin <- 1000

# Initialize the MCMC chains
n_chains <- 3
jags_inits <- list(
  list(lambda_0 = 10, lambda_1 = 10, p = 0.5, x = rbinom(T, 1, 0.5)),
  list(lambda_0 = 10, lambda_1 = 10, p = 0.5, x = rbinom(T, 1, 0.5)),
  list(lambda_0 = 10, lambda_1 = 10, p = 0.5, x = rbinom(T, 1, 0.5))
)

# Run the MCMC sampler
jags_model <- jags.model(textConnection(model_string), data = data_list, inits = jags_inits, n.chains = n_chains)
jags_samples <- coda.samples(jags_model, variable.names = params, n.iter = n_iter, n.burnin = n_burnin, thin = 1)

# Calculate the posterior means and credible intervals
post_means <- apply(jags_samples[[1]], 2, mean)
post_cis <- t(apply(jags_samples[[1]], 2, function(x) quantile(x, c(0.025, 0.975))))
post_lambdas <- post_means[1] + post_means[2]

# Print the results
cat("Posterior mean of lambda_0:", post_means[1], "\n")
cat("95% credible interval of lambda_0:", post_cis[1], "\n")




# Start to Q5 (unfinished):
# # Define the model
# code <- nimbleCode({
#   for (i in 1:n) {
#     log(tn) ~ dnorm(mu[i], tau)
#     mu[i] <- beta0 + beta1*ag[i] + beta2*forest[i] + f[i] + g[i]
#     f[i] <- splines::bs(log_tp, knots = knots_tp, degree=3) %*% beta_f
#     g[i] <- splines::bs(log_secchi, knots = knots_secchi, degree=3) %*% beta_g
#   }
#   beta0 ~ dnorm(0, 1.0E-6)
#   beta1 ~ dnorm(0, 1.0E-6)
#   beta2 ~ dnorm(0, 1.0E-6)
#   beta_f ~ dnorm(0, 1.0E-6)
#   beta_g ~ dnorm(0, 1.0E-6)
#   tau ~ dgamma(1.0E-3, 1.0E-3)
# })
# 
# # Set up the data and parameters
# data <- list(
#   n = nrow(lagos),
#   log_tp = log(lagos$tp),
#   log_secchi = log(lagos$secchi),
#   ag = lagos$ag,
#   forest = lagos$forest,
#   tn = log(lagos$tn)
# )
# params <- c("beta0", "beta1", "beta2", "beta_f", "beta_g", "tau")
# 
# # Set up the priors
# inits <- list(
#   beta0 = 0,
#   beta1 = 0,
#   beta2 = 0,
#   beta_f = rnorm(length(knots_tp) + 1),
#   beta_g = rnorm(length(knots_secchi) + 1),
#   tau = 1
# )
# prior <- list(
#   beta0 = list(distribution = "dnorm", mean = 0, sd = 1.0E6),
#   beta1 = list(distribution = "dnorm", mean = 0, sd = 1.0E6),
#   beta2 = list(distribution = "dnorm", mean = 0, sd = 1.0E6),
#   beta_f = list(distribution = "dnorm", mean = 0, sd = 1.0E6),
#   beta_g = list(distribution = "dnorm", mean = 0, sd = 1.0E6),
#   tau = list(distribution = "dgamma", alpha = 0.001, beta = 0.001)
# )
# 
# # Set up the knots and degrees for the B-spline basis functions
# knots_tp <- mean(log(lagos$tp))
# knots_secchi <- mean(log(lagos$secchi))
# 
# # Compile and run the model
# model <- nimbleModel(code, data = data, inits = inits, 
#                      constants = list(knots_tp=knots_tp, knots_secchi=knots_secchi))
# 
# samples <- nimbleMCMC(model, params, prior = prior, niter = 20000, nburn = 5000, nthin = 2)
# 
# # Summarize the results
# summary(samples)
# 

# Unfinished

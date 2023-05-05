# STAT 597
# Homework 4
# Veronica A. Winter
# 04/03/2023

# clean env
rm(list = ls())
gc()
ls()

# load required packages
library(coda)
library(loo)
library(rjags)
library(invgamma)

# set for reproducibility
set.seed(123)

# Question 1 ----
## Part 1:
# Now construct a MCMC sampler BY HAND to fit this.  
# You should work out full conditional distributions for all 
# parameters in the following model:

#       y_t∼N(z_t,τ^2 )
#       z_(t+1)∼N(z_t-β(z_t-u)+g_1 h_1t+g_2 h_2t,σ^2 ),t=2,3,…,300
#       u~N(〖0,100〗^2 )
#       β~N(〖0,100〗^2 )
#       g_1~N(〖0,100〗^2 )
#       g_2~N(〖0,100〗^2 )
#       z_1~N(〖0,100〗^2 )
#       σ^2∼IG(10,100)
#       τ^2∼IG(10,100)


# Steps to implement by hand:
# 1. Initialize the values of all the parameters.
# 2. For i = 1 to the desired number of iterations:
#   Sample z_1 from N((y_1+u)/(1+beta), tau^2/(1+beta^2))

# For t = 2 to 300:
#   Sample z_t from N(mu_t, sigma^2), where
#       mu_t = (y_t + betaz_{t-1} + g_1h_1t + g_2h_2*t)/(1 + beta^2)
#   Sample u from N((beta*z_1)/(1+beta), tau^2/(1+beta^2))
#   Sample beta from N((z_1-u)/(z_2-z_1), tau^2/(z_2-z_1)^2)
#   Sample g_1 from N((z_t-z_{t-1})/(h_1t), sigma^2/(h_1^2t^2))
#   Sample g_2 from N((z_t-z_{t-1})/(h_2t), sigma^2/(h_2^2t^2))
#   Sample sigma^2 from IG(a_sigma, b_sigma), where a_sigma = 10 + 1500/2 
#       and b_sigma = 100 + 0.5sum((z_t - z_{t-1} - beta(z_{t-1} - u) - g_1h_1t - g_2h_2t)^2)
#   Sample tau^2 from IG(a_tau, b_tau), where a_tau = 10 + 150/2 and 
#       b_tau = 100 + 0.5*sum((y_t - z_t)^2)


# Note that we have used the notation N(a,b) to denote a normal distribution 
# with mean a and variance b, and IG(a,b) to denote an inverse gamma distribution 
#   with shape parameter a and scale parameter b.


## part 1: ----
# set initial values and parameters
u <- 74
beta <- 0.06
s2 <- 3
z0 <- 70
g1 <- 40
g2 <- 75
T <- 300

# simulate time series for z and y
z <- rep(z0, T+1)
for(t in 1:T){
  z[t+1] <- z[t] - beta*(z[t]-u) + rnorm(1, mean=0, sd=sqrt(s2))
  if(t == 100){
    z[t+1] <- z[t+1] + g1
  }
  if(t == 200){
    z[t+1] <- z[t+1] + g2
  }
}
y <- rnorm(length(z), mean=z, sd=4)

# thin out y
y[-seq(10, 300, by=10)] <- NA
t.obs=which(y>-Inf)
y.obs=y[t.obs]

# set initial values 
sigma2 <- rinvgamma(1, 10, 100)
tau2 <- rinvgamma(1, 10, 100)
h1 <- rep(0, T)
h1[100] <- 1
h2 <- rep(0, T)
h2[200] <- 1
z[1] <- rnorm(1, 0, 100)

# set up priors
n=length(y.obs)
# create data list
data_list <- list(
  y = y.obs,
  t.obs = t.obs,
  n = length(y.obs),
  T = 300,
  h1 = h1,
  h2 = h2
)

# define model
ssm_model <- "model {
  # observed data
  for(i in 1:n){
    y[i] ~ dnorm(z[t.obs[i]], sqrt(tau2))
  }
  
  # OU state space process
  for(t in 2:T){
    z[t] ~ dnorm(z[t - 1] - beta * (z[t - 1] - u) + g2 * h1[t - 1] + g1 * h2[t - 1], 
                  sqrt(s2))
  }
  
  # prior for t = 1
  z[1] ~ dnorm(0, 100)
  
  # priors for parameters
  s2 ~ dgamma(10,1/ 100)
  tau2 ~ dgamma(10, 1/100)
  beta ~ dnorm(0, 100)
  g1 ~ dnorm(0, 100)
  g2 ~ dnorm(0, 100)
  u ~ dnorm(0, 100)
}"

# compile model
ssm_jags <- jags.model(textConnection(ssm_model), data = data_list, n.chains = 4)

# burn-in
update(ssm_jags, 1000)

# obtain posterior samples
ssm_samples <- coda.samples(model = ssm_jags, variable.names = c("u", "beta", "g1", "g2", "s2", "tau2"), n.iter = 5000)

# extract posterior summaries
summary(ssm_samples)
# Extract the posterior samples from the MCMC chains
posterior_samples <- as.matrix(ssm_samples)

plot(posterior_samples[, "beta"], main = "Trace plot for beta")
plot(posterior_samples[, "g1"], main = "Trace plot for g1")
plot(posterior_samples[, "g2"], main = "Trace plot for g2")
plot(posterior_samples[, "s2"], main = "Trace plot for sigma2")
plot(posterior_samples[, "tau2"], main = "Trace plot for tau2")
plot(posterior_samples[, "u"], main = "Trace plot for u")

# Calculate the log-likelihood for each posterior sample                         
loglik <- function(theta, data) {
  y <- na.omit(data$y)
  h1 <- data$h1
  h2 <- data$h2
  t.obs <- data$t.obs
  n <- length(y)
  T <- 300
  u <- theta[1]
  beta <- theta[2]
  g1 <- theta[3]
  g2 <- theta[4]
  s2 <- theta[5]
  tau2 <- theta[6]
  
  # observed data
  for(i in 1:n){
    y[i] ~ dnorm(z[t.obs[i]], sqrt(tau2))
  }
  
  # OU state space process
  for(t in 2:T){
    z[t] ~ dnorm(z[t - 1] - beta * (z[t - 1] - u) + g2 * h1[t - 1] + g1 * h2[t - 1], 
                 sqrt(s2))
  }
  
  # prior for t = 1
  z[1] ~ dnorm(0, 100)
  
  # priors for parameters
  s2 ~ dgamma(10,1/ 100)
  tau2 ~ dgamma(10, 1/100)
  beta ~ dnorm(0, 100)
  g1 ~ dnorm(0, 100)
  g2 ~ dnorm(0, 100)
  u ~ dnorm(0, 100)
  
  # initial values
  z[1] <- rnorm(1, 0, 100)
  z[2:(T+1)] <- rep(0, T)
  
  # calculate log-likelihood
  target <- 0
  for(i in 1:n){
    target <- target + dnorm(y[i], z[t.obs[i]], sqrt(tau2), log = TRUE)
  }
  
  return(target)
}


ll <- apply(posterior_samples, 1, loglik, data=data_list)

# Calculate DIC
dic <- dic.samples(model = ssm_jags,
                   log_likelihood = ll, 
                   n_eff = effectiveSize(posterior_samples), 
                   y = y.obs,
                   n.iter = 10000)

dic

# Calculate WAIC
# compute posterior mean of log-likelihood
lppd <- mean(na.omit(ll))

# compute effective number of parameters
pwaic <- var(na.omit(ll))

# compute WAIC
waic <- -2 * (lppd - pwaic)
waic

# Here, we compute the log-likelihood for each observation based on the posterior 
# distribution of the parameters obtained from the MCMC samples. 
# We then use these log-likelihoods to compute the deviance and effective number 
# of parameters for the WAIC and DIC. 

# Note that for the DIC, we first need to compute the posterior mean of the 
# parameters and use this to compute the log-likelihood of the held-out dataset.

# Question 2: ----
# Now, for the same data, fit the following model:
# y_t∼N(∑_k▒〖ϕ_kt α_k 〗,σ^2 )
# α∼N(0,τ^2 K)
# σ^2∼IG(10,100)
# τ^2∼IG(10,100)

# In the above model, ϕ_kt is the k-th B-spline basis function in a semiparametric model for the mean of y_t.  Use 20 B-spline basis functions of order=4
# B=create.bspline.basis(rangeval=c(mn,mx),nbasis=20,norder=4)
# Phi=eval.basis(1:300,B)
# The matrix K is K=D’D where D is a matrix that calculates the 2nd differences of α. See BsplineEx.r for examples.
# You may fit this model in nimble, or by hand.  Either way, compute WAIC and DIC for this model, and then compare this model with the model in Q1.  Which is better for this data? 
#   


# Q 2
rm(list = ls())
gc()

library(fda)
library(mgcv)
library(invgamma)
library(splines)
library(mvtnorm)

# set initial values and parameters
beta <- 0.06
s2 <- 3
z0 <- 70
g1 <- 40
g2 <- 75
T <- 300
u <- 74
# simulate time series for z and y
z <- rep(z0, T+1)
for(t in 1:T){
  z[t+1] <- z[t] - beta*(z[t]-u) + rnorm(1, mean=0, sd=sqrt(s2))
  if(t == 100){
    z[t+1] <- z[t+1] + g1
  }
  if(t == 200){
    z[t+1] <- z[t+1] + g2
  }
}
y <- rnorm(length(z), mean=z, sd=4)

# thin out y
y[-seq(10, 300, by=10)] <- NA
t.obs=which(y>-Inf)
y.obs=y[t.obs]

# set initial values 
sigma2 <- rinvgamma(1, 10, 100)
tau2 <- rinvgamma(1, 10, 100)
h1 <- rep(0, T)
h1[100] <- 1
h2 <- rep(0, T)
h2[200] <- 1

# predictor values
pred.vals = seq(mn,mx, by=.1)
Zpred=eval.basis(pred.vals,B)

# set up priors
n = length(y.obs)
mn <- min(y.obs)
mx <-  max(y.obs)
B=create.bspline.basis(rangeval=c(mn,mx),nbasis=20,norder=4)
Phi=eval.basis(y.obs,B)
n_basis <- dim(Phi)[2]
mean_prior <- rep(0, n_basis)
shape_prior <- 0.001
rate_prior <- 0.001

# Define D
p=ncol(Phi)
D=matrix(0, nrow=p-2, ncol=p)
for(i in 1:(p-2)){
  D[i,i]=1
  D[i,i+1]=-2
  D[i,i+2]=1
}
D


## K
K=t(D)%*%D
dim(K)
K

# Bspline funciton
mcmc.pen.spline = function(y, Z, K, a, b, c, d, n.mcmc) {
  
  # set up starting values
  n=length(y.obs)
  alpha=rep(0,ncol(Z))
  s2=var(y.obs)
  tau2=var(y.obs)
  
  # set up storage for posterior samples
  alpha.save <- matrix(0, n.mcmc, p)
  sigma2.save <- rep(0, n.mcmc)
  tau2.save <- rep(0, n.mcmc)
  
  # MCMC loop
  for (iter in 1:n.mcmc) {
    
    # sample alpha
    V <- solve(t(Z) %*% Z / sigma2 + K / tau2)
    m <- V %*% t(Z) %*% y / sigma2
    a.inv = solve(V)
    alpha <- as.numeric(mvrnorm(1, a.inv%*%m, a.inv))

    # sample sigma2
    a.sigma <- a + n / 2
    b.sigma <- b + 0.5 * sum((y - Z %*% alpha)^2)
    sigma2 <- 1 / rgamma(1, a.sigma, b.sigma)
    
    # sample tau2
    a.tau <- a + p / 2
    b.tau <- b + 0.5 * t(alpha) %*% K %*% alpha
    tau2 <- 1 / rgamma(1, a.tau, b.tau)
    
    # store samples
    alpha.save[iter, ] <- alpha
    sigma2.save[iter] <- sigma2
    tau2.save[iter] <- tau2
    
    # progress indicator
    if (iter %% 100 == 0) cat(iter, " ")
    
  }
  
  # return posterior samples
  list(alpha = alpha.save, sigma2 = sigma2.save, tau2 = tau2.save)
  
}


# fit the model 
fit = mcmc.pen.spline(y = y.obs, Z = Phi,K = K, 
                      a = 10, b = 10, c = 10, d = 10, n.mcmc = 10000)

# extract posterior samples
alpha.samples <- fit$alpha
sigma2.samples <- fit$sigma2
tau2.samples <- fit$tau2

# compute log-likelihood of the model for each iteration of the MCMC chain
loglik <- rep(0, nrow(alpha.samples))
for (i in 1:nrow(alpha.samples)) {
  alpha <- alpha.samples[i, ]
  sigma2 <- sigma2.samples[i]
  tau2 <- tau2.samples[i]
  loglik[i] <- sum(dnorm(y.obs, mean = Phi %*% alpha, sd = sqrt(sigma2)))
}

# compute DIC
d_bar <- mean(loglik)
p_dbar <- -2 * (lppd - d_bar)
dic <- d_bar + p_dbar
dic

# compute WAIC
lppd <- sum(loglik)
p_waic <- sum(log(colMeans(alpha.samples != 0)) != -Inf)
waic <- -2 * (lppd - p_waic)
waic



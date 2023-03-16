# Libraries
library(coda)
library(mvtnorm)
library(maptools)
library(maps)
library(mgcv)
library(glmnet)
library(MASS)
library(nimble)
library(stan)
library(rjags)
library(statmod)

load("Homework3/lagos.Rdata")

data <- lagos


library(rjags)

# Load data
data(lagos)


# Define mu
beta <- rep(NA, p)
# Define model
model_string <- "
model {
# Define x2 and y2
  x2 <- x^2
  y2 <- y^2
  
  # Priors
  beta_0 ~ dnorm(0, 1.0E-6) # vague Gaussian prior for intercept

  tau ~ dgamma(1.0E-3, 1.0E-3) # hyperprior for ridge regression tuning parameter
  
  for (j in 1:p) {
    beta[j] ~ dnorm(0, tau) # ridge regression prior for regression parameters
  }
  
  
  # Define mu
  mu_i <- rep(0, n)
  
  # Likelihood
  for(i in 1:n) {
    mu_i[i] <- beta_0 + 
             beta[1]*tp[i] + 
             beta[2]*secchi[i] + 
             beta[3]*ag[i] + 
             beta[4]*forest[i] + 
             beta[5]*x[i] + 
             beta[6]*y[i] + 
             beta[7]*x2[i] + 
             beta[8]*y2[i]
             
  log_tn[i] ~ dnorm(mu_i, tau_obs)
  }
  
  # Hyperprior for residual variance
  tau_obs ~ dgamma(1.0E-3, 1.0E-3)
}"

# Data
data_list <- list(n = nrow(lagos),
     p = 8,
     log_tn = log(lagos$tn),
     tp = log(lagos$tp),
     secchi = log(lagos$secchi),
     ag = lagos$ag,
     forest = lagos$forest,
     x = lagos$lon,
     y = lagos$lat)

# Initial values
inits <- list(beta_0 = 0, beta = rep(0, 8), tau = 1.0, tau_obs = 1.0)

# Parameters to monitor
params <- c("beta_0", "beta", "tau", "tau_obs")

# Run model
model <- jags.model(textConnection(model_string), data = data_list,
                    inits = inits, n.chains = 3, n.adapt = 1000)
update(model, n.iter = 5000)
samples <- coda.samples(model, variable.names = params, n.iter = 5000)

# Summarize results
summary(samples)









#############################
# Model code
model_code <- "model {
  # priors for variance params
  s2 ~ dt(0, var.log.tn, 1)
  lambda ~ dexp(1) ## actually 1/lambda
  ## priors for regression parameters
  beta[1] ~ dnorm(0, 10)
  for(b in 2:M){
    beta[b] ~ ddexp(0, inverse(lambda))
  }
  ## data model
  for(i in 1:N){
    y[i] ~ dnorm(inprod(X[i,1:M],beta[1:M]), sqrt(s2))
  }
}"

N <- nrow(data)
M <- ncol(data)


# Initial values for MCMC
inits_list <- list(
  s2 = 1,
  lambda = 1,
  beta = rep(0, M)
)

X <- list(
  p = 9,
  tn = log(data$tn),
  tp = log(data$tp),
  secchi = log(data$secchi),
  ag = data$ag,
  forest = data$forest,
  x = data$lon,
  y = data$lat,
  x2 = data$lon^2,
  y2 = data$lat^2
)
# Parameters to monitor
params <- c("beta", "s2", "lambda")
# Run the model
model <- jags.model(textConnection(model_code), data = X, inits = inits_list)








mcmc.lm.rr2=function(y,X,n.mcmc,
                    beta.start=rep(0,ncol(X)),s2.start=1,k2.start=1,
                    beta.prior.var=100,s2.prior.var=var(y),k2.exp.rate=1,
                    tune.k2=.01,tune.s2=.01,
                    adapt.iter=100, tau2 = 1){
  
  beta.save=matrix(NA,n.mcmc,ncol(X))
  colnames(beta.save)=colnames(X)
  s2.save=rep(NA,n.mcmc)
  k2.save=rep(NA,n.mcmc)
  n=length(y)
  p=ncol(X)-1
  beta=beta.start
  k2=k2.start
  s2=s2.start
  accept.k2=0
  accept.s2=0
  adapt.t=0
  for(iter in 1:n.mcmc){
    if(iter%%100==0){
      cat(iter," ")
    }
    ## update s2 (RW proposal on log(s2))
    s2.star=exp(rnorm(1,log(s2),tune.s2))
    mh1=sum(dnorm(y,X%*%beta,sqrt(s2.star),log=TRUE))+
      dnorm(s2.star,0,sqrt(s2.prior.var),log=TRUE)+log(s2.star)
    mh2=sum(dnorm(y,X%*%beta,sqrt(s2),log=TRUE))+
      dnorm(s2,0,sqrt(s2.prior.var),log=TRUE)+log(s2)
    if(runif(1)<exp(mh1-mh2)){
      s2=s2.star
      accept.s2=accept.s2+1
    }
    
    ## update k2 (RW proposal on log(k2))
    k2.star=exp(rnorm(1,log(k2),tune.k2))
    mh1=sum(dnorm(beta[-1],0,sqrt(k2.star),log=TRUE))+
      dnorm(k2.star,0,k2.exp.rate,log=TRUE)+log(k2.star)
    mh2=sum(dnorm(beta[-1],0,sqrt(k2),log=TRUE))+
      dnorm(k2,0,k2.exp.rate,log=TRUE)+log(k2)
    if(runif(1)<exp(mh1-mh2)){
      k2=k2.star
      accept.k2=accept.k2+1
    }
    ## update beta
    D.inv=diag(c(1/beta.prior.var, rep(1/tau2, p)))
    A.inv=solve(t(X)%*%X/s2+D.inv)
    b=1/s2*t(X)%*%y + solve(D.inv) %*% beta
    beta=t(rmvnorm(1,A.inv%*%b,A.inv))
    
    ## save
    beta.save[iter,]=beta
    s2.save[iter]=s2
    k2.save[iter]=k2
    
    ##
    ## Adaptive tuning
    ##
    if(iter%%adapt.iter==0){
      ## move adapt counter up 1
      adapt.t=adapt.t+1
      ## new tuning parameters for k2
      adapt.vals=get.sigma(tune.k2,1,data=k2.save[(iter-adapt.iter)+1:adapt.iter],accept=accept.k2,t.adapt=adapt.t)
      tune.k2=adapt.vals$s2
      ## new tuning parameters for s2
      adapt.vals=get.sigma(tune.s2,1,data=s2.save[(iter-adapt.iter)+1:adapt.iter],accept=accept.s2,t.adapt=adapt.t)
      tune.s2=adapt.vals$s2
      ## resetting acceptances to 0
      accept.k2=0
      accept.s2=0
    }
    
    
  }
  cat("\n")
  ##output
  out=cbind(beta.save,s2.save,k2.save)
  colnames(out)=c(colnames(X),"s2","k2")
  data.frame(out)
}




#Define x2 and y2
x2 <- x^2
y2 <- y^2

# Priors
beta_0 ~ dnorm(0, 1.0E-6) # vague Gaussian prior for intercept

tau ~ dgamma(1.0E-3, 1.0E-3) # hyperprior for ridge regression tuning parameter

for (j in 1:p) {
  beta[j] ~ dnorm(0, tau) # ridge regression prior for regression parameters
}


# Define mu
mu_i <- rep(0, n)

# Likelihood
for(i in 1:n) {
  mu_i[i] <- beta_0 + 
    beta[1]*tp[i] + 
    beta[2]*secchi[i] + 
    beta[3]*ag[i] + 
    beta[4]*forest[i] + 
    beta[5]*x[i] + 
    beta[6]*y[i] + 
    beta[7]*x2[i] + 
    beta[8]*y2[i]
  
  log_tn[i] ~ dnorm(mu_i, tau_obs)
}

# Hyperprior for residual variance
tau_obs ~ dgamma(1.0E-3, 1.0E-3)



Under this new model, show that 𝑛𝑡0 and 𝑛𝑡1 have conjugate updates, and implement
an MCMC sampler to draw samples from the posterior distribution, and report posterior
means and 95% credible intervals of the mean number of calls in the “high” state which
has rate (𝜆0 + 𝜆1).

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
library(nimble)

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
samples <- matrix(0, niter, 2)
accept <- 0
cov_mat <- proposal_var
tune_interval <- 100
tune_count <- 0
tune_iter <- tune_interval

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

### Accept/Reject ----
accept
# [1] 4770

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
#xy=lagos[,2:1]
#plot(xy,pch=".",cex=2,col="blue")
map("state",add=T,lwd=2)

## Create Indicator Vars
lagos$ag=0
lagos$ag[lagos$SurLand=="agricultural"] <- 1
lagos$forest=0
lagos$forest[lagos$SurLand=="forest"] <- 1


## model
lasso.lake <- nimbleCode({
  ## priors for variance params
  s2 ~ T(dnorm(mean = 0, var = var.log.tn), min=0, max=Inf)
  lambda ~ dexp(1) ## actually 1/lambda
  ## priors for regression parameters
  beta[1] ~ dnorm(mean = 0, var=1/100)
  for(b in 2:M){
    beta[b] ~ ddexp(location = 0, scale = lambda)
  }
  ## data model
  for(i in 1:N){
    y[i] ~ dnorm(mean = inprod(X[i,1:M],beta[1:M]), sd = sqrt(s2))
  }
})


## create data frame of covariates
X <- data.frame(
  int = rep(1,nrow(lagos)),
  log.tp = log(lagos$tp),
  log.sec = log(lagos$secchi),
  x = lagos$lon,
  y = lagos$lat,
  x.2 = lagos$lon^2,
  y.2 = lagos$lat^2
)
## create vector of response variable
y <- log(lagos$tn)

## number of observations
N <- nrow(X)
## number of covariates
M <- ncol(X)

lake.data <- list(
  y = y,
  X = X
)
lake.const <- list(
  N = N,
  M = M,
  var.log.tn = var(log(y))
)
lake.inits <- list(
  s2 = 1,
  lambda = 1,
  beta = rep(0, M)
)

## run model
## compile model and run MCMC
lake.out.lasso <- nimbleMCMC(
  lasso.lake,
  constants = lake.const,
  data = lake.data,
  inits = lake.inits,
  niter = 1000,
  nburnin = 600,
  thin = 1,
  monitors = c("beta","s2","lambda")
)





#################################################################

# Predictor variables
N=nrow(X.lagos)
## dim of data
M=ncol(X.lagos)

x=cbind(1,as.matrix(lagos[,-1]))
y=lagos[,1]
class(x)

n <- nrow(lagos)
tp  <-  log(lagos$tp)
secchi <- log(lagos$secchi) 
ag <- lagos$ag
forest <- lagos$forest
tn <- log(lagos$tn)

## set up 
M = 10000
mu = rep(NA, M)
b = rep(NA, M)

# Use a ridge regression prior on the regression parameters, 
# but use a vague Gaussian prior for the intercept. 

## Reminder: alpha=0 gives ridge regression
# Define data and parameters
data_list <- list(n, p, tp, secchi, ag, forest, x, y)
params <- c("b", "sigma", "tau")

# Define the model in Nimble syntax
# define the model
model_code <- "
  model {
    # define priors
    for (i in 1:length(lagos)) {
      beta[i] ~ dnorm(0, tau_beta)
    }
    beta[9] ~ dnorm(0, 100) # vague prior for intercept
    tau_beta ~ dgamma(0.001, 0.001) # ridge regression prior on beta
  
    # define likelihood
    for (i in 1:n) {
      tn_store[i] = log(tn) ~ dnorm(mu[i], tau)
      mu[i] <- beta[1]*log(tp[i]) + beta[2]*log(secchi[i]) + beta[3]*ag[i] + 
               beta[4]*forest[i] + beta[5]*x[i] + beta[6]*y[i] +
               beta[7]*x[i]^2 + beta[8]*y[i]^2 + beta[9]
    }
  
    # define hyperparameters
    tau ~ dgamma(0.001, 0.001)
  }"

# compile the model
model <- nimbleModel(code=model_code, 
                     constants=list(n=nrow(lagos)), 
                     data = list(tn = log(lagos$tn), tp=log(lagos$tp), secchi=log(lagos$secchi),
                                 ag=lagos$ag, forest=lagos$forest, x=x, y=y), 
                     inits=list(beta=rep(0,9), tau_beta=1, tau=1))

# set MCMC options
options <- list(numChains=4, numSamples=5000, thin=2, adaptDelta=0.99)

# compile and build the MCMC object
mcmc <- buildMCMC(model, options)

# run the MCMC
runMCMC(mcmc)

# summarize the MCMC output
print(mcmc, summary=TRUE)

# extract posterior samples
post <- as.matrix(mcmc)
# compute 95% credible intervals for all parameters
cred_int <- apply(post, 2, function(x) quantile(x, c(0.025, 0.975)))

# compute 95% credible intervals for functions of x, y, x^2, and y^2
b <- post[, 5]
c <- post[, 7]
d <- post[, 8]
e <- post[, 9]

cred_int_b <- quantile(b*x + c*x^2, c(0.025, 0.975))
cred_int_d <- quantile(d*x + e*x^2, c(0.025, 0.975))

# print the results
cat("95% credible intervals for all parameters:\n")
print(cred_int)

cat("95% credible interval for function bx + cx^2:\n")
print(cred_int_b)

cat("95% credible interval for function dx + ex^2:\n")
print(cred_int_d)




# Question 4-----
set.seed(123)

# Define the data and prior distributions
y <- c(9, 15, 14, 5, 6, 4, 4, 4, 8, 0, 10, 22, 7, 2, 6, 2, 18, 
       9, 4, 2, 5, 7, 7, 5, 8, 7, 2, 15, 17, 7, 1, 4, 8, 5, 8, 9, 25, 6, 6, 4, 22, 3, 2, 5, 3, 4, 8, 4, 17, 14)
n <- length(y)
shape <- 1
rate <- 0.1
mean_gamma <- 10
var_gamma <- 100

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
shape <- 10
rate <- 1

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

# Compute the posterior mean and 95% credible intervals for lambda_high
lambda_high <- lambda0 + lambda1
post_mean_lambda_high <- mean(lambda_high)
post_ci_lambda_high <- quantile(lambda_high, c(0.025, 0.975))

# Print the results
cat("Posterior mean of lambda_high: ", post_mean_lambda_high, "\n")
# Posterior mean of lambda_high:  14.02394 
cat("95% credible interval for lambda_high: ", post_ci_lambda_high[1], "-", post_ci_lambda_high[2], "\n")
# 95% credible interval for lambda_high:  9.169818 - 18.39378 



# The code implements an MCMC sampler using the Metropolis-Hastings algorithm to
# draw samples from the posterior distribution of the mean number of calls in the
# "high" state of a Poisson model with two latent states.
# 
# The model assumes that there are two latent states,
# one with a high call rate and one with a lower call rate.
# The observed data are the number of calls to a help line each hour.
# The model assumes that the observed data follow a Poisson distribution with
# parameter λ_t, where λ_t is the sum of two parameters, λ_0 and λ_1 times a
# binary indicator variable x_t. The indicator variable x_t follows a Bernoulli
# distribution with parameter p, which is itself assumed to follow a uniform
# distribution on the interval (0, 1). The parameters λ_0 and λ_1 are assumed
# to follow gamma distributions with shape parameters equal to 1 and rate
# parameters equal to 0.1, which gives them means of 10 and variances of 100.
# 
# The Metropolis-Hastings algorithm is used to draw samples from the
# joint posterior distribution of the parameters λ_0, λ_1, and p. The algorithm
# iteratively generates proposals for the parameters and accepts or rejects them
# based on the posterior density of the proposed values relative to the current values.
# 
# The code initializes the parameters λ_0, λ_1, and p and sets the
# number of MCMC iterations. It then iteratively generates proposals for
# the parameters and evaluates their posterior density using the observed data
# and the model specification. The proposal distributions are chosen to be normal
# distributions centered at the current values of the parameters, with standard
# deviations chosen to ensure reasonable acceptance rates. The code stores the
# accepted parameter values and reports the posterior means and 95% credible
# intervals for the mean number of calls in the "high" state.

# Question 5 -----
# Define the model
code <- nimbleCode({
  for (i in 1:n) {
    log(tn) ~ dnorm(mu[i], tau)
    mu[i] <- beta0 + beta1*ag[i] + beta2*forest[i] + f[i] + g[i]
    f[i] <- splines::bs(log_tp, knots = knots_tp, degree=3) %*% beta_f
    g[i] <- splines::bs(log_secchi, knots = knots_secchi, degree=3) %*% beta_g
  }
  beta0 ~ dnorm(0, 1.0E-6)
  beta1 ~ dnorm(0, 1.0E-6)
  beta2 ~ dnorm(0, 1.0E-6)
  beta_f ~ dnorm(0, 1.0E-6)
  beta_g ~ dnorm(0, 1.0E-6)
  tau ~ dgamma(1.0E-3, 1.0E-3)
})

# Set up the data and parameters
data <- list(
  n = nrow(lagos),
  log_tp = log(lagos$tp),
  log_secchi = log(lagos$secchi),
  ag = lagos$ag,
  forest = lagos$forest,
  tn = log(lagos$tn)
)
params <- c("beta0", "beta1", "beta2", "beta_f", "beta_g", "tau")

# Set up the priors
inits <- list(
  beta0 = 0,
  beta1 = 0,
  beta2 = 0,
  beta_f = rnorm(length(knots_tp) + 1),
  beta_g = rnorm(length(knots_secchi) + 1),
  tau = 1
)
prior <- list(
  beta0 = list(distribution = "dnorm", mean = 0, sd = 1.0E6),
  beta1 = list(distribution = "dnorm", mean = 0, sd = 1.0E6),
  beta2 = list(distribution = "dnorm", mean = 0, sd = 1.0E6),
  beta_f = list(distribution = "dnorm", mean = 0, sd = 1.0E6),
  beta_g = list(distribution = "dnorm", mean = 0, sd = 1.0E6),
  tau = list(distribution = "dgamma", alpha = 0.001, beta = 0.001)
)

# Set up the knots and degrees for the B-spline basis functions
knots_tp <- mean(log(lagos$tp))
knots_secchi <- mean(log(lagos$secchi))

# Compile and run the model
model <- nimbleModel(code, data = data, inits = inits, 
                     constants = list(knots_tp=knots_tp, knots_secchi=knots_secchi))

samples <- nimbleMCMC(model, params, prior = prior, niter = 20000, nburn = 5000, nthin = 2)

# Summarize the results
summary(samples)

# Get the posterior means and credible intervals for the smooth

# Q3 again
## Create Indicator Vars
lagos$ag=0
lagos$ag[lagos$SurLand=="agricultural"] <- 1
lagos$forest=0
lagos$forest[lagos$SurLand=="forest"] <- 1

# Load data
Lagos <- lagos

# Define the model in nimbleCode
model_code <- nimbleCode({
  
  # Define priors
  for (i in 1:9) {
    beta[i] ~ dnorm(0, tau[i])
  }
  tau[1:8] ~ dgamma(0.001, 0.001)
  tau[9] <- 1e-6 # very vague prior for intercept
  
  # Define model
  for (i in 1:n) {
    log(tn[i]) ~ dnorm(mu[i], sigma)
    mu[i] <- beta[1] + beta[2]*log(tp[i]) + beta[3]*log(secchi[i]) + 
      beta[4]*ag[i] + beta[5]*forest[i] + beta[6]*x[i] + beta[7]*y[i] + 
      beta[8]*x_sq[i] + beta[9]*y_sq[i]
  }
  
  # Define variance
  sigma ~ dgamma(0.001, 0.001)
})

# Create nimble model
model <- nimbleModel(code = model_code, constants = list(n = nrow(Lagos)), 
                     data = list(tn = log(Lagos$tn), tp = log(Lagos$tp), secchi = log(Lagos$secchi), 
                                 ag = Lagos$ag, 
                                 forest = Lagos$forest,
                                 x = Lagos$lon, y = Lagos$lat,
                                 x_sq = (Lagos$lon)^2, y_sq = (Lagos$lat)^2))

# Compile the model
compiled_model <- compileNimble(model)

# Set up MCMC sampler
MCMC_specs <- list(
  numChains = 1,
  thin = 1,
  nBurnin = 1000,
  nMCMC = 20000,
  DIC = FALSE,
  monitorParams = c("beta", "sigma"),
  monitorLiks = TRUE
)

# Initialize MCMC chains
MCMC_chains <- buildMCMC(compiled_model, MCMC_specs)

# Run MCMC sampler
MCMC_run <- runMCMC(MCMC_chains, inits = list(beta = rnorm(9, 0, 10), sigma = rgamma(1, 0.001, 0.001)))

# Get posterior samples
posterior_samples <- as.matrix(MCMC_run$theta)

# Get 95% credible intervals for all parameters
credible_intervals <- t(apply(posterior_samples, 2, function(x) quantile(x, c(0.025, 0.975))))

# Get 95% credible intervals for bx + cx^2 and dx + ex^2
bx_cx2_samples <- posterior_samples[,6] + posterior_samples[,8]*Lagos$lon + posterior_samples[,9]*Lagos$lon^2
dx_ex2_samples <- posterior_samples[,7] + posterior_samples[,8]*Lagos$lat + posterior_samples[,9]*Lagos$lat^2

bx_cx2_intervals <- quantile(bx_cx2_samples, c(0.025, 0.975))
dx_ex2_intervals <- quantile(dx_ex2_samples, c(0.025, 0.975))

# Print results
cat("95% credible intervals for regression coefficients:\n")
print

`
# Start to Q5 (unfinished):
```{r}
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
# Load in data
load("../Homework3/lagos.Rdata")
library(rjags)


## Create Indicator Vars
lagos$ag=0
lagos$ag[lagos$SurLand=="agricultural"] <- 1
lagos$forest=0
lagos$forest[lagos$SurLand=="forest"] <- 1

data <- lagos
# Model code
model_code <- "model {
  # priors for variance params
  s2 ~ dt(0, var.log.tn, 1)
  lambda ~ dexp(1) ## actually 1/lambda
  ## priors for regression parameters
  beta[1] ~ dnorm(0, 10)
  for(b in 2:M){
    beta[b] ~ ddexp(0, inverse(lambda))
  }
  ## data model
  for(i in 1:N){
    y[i] ~ dnorm(inprod(X[i,1:M],beta[1:M]), sqrt(s2))
  }
}"

# Initial values for MCMC
inits_list <- list(
  s2 = 1,
  lambda = 1,
  beta = rep(0, M)
)

X <- list(
  n = n,
  p = p,
  tn = log(data$tn),
  tp = log(data$tp),
  secchi = log(data$secchi),
  ag = data$ag,
  forest = data$forest,
  x = data$lon,
  y = data$lat,
  x2 = data$lon^2,
  y2 = data$lat^2
)

# Data and parameters
N <- nrow(data)
M <- ncol(data)
var.log.tn <- var(log(y))

# Parameters to monitor
params <- c("beta", "s2", "lambda")

# Run the model
model <- jags.model(textConnection(model_code), data = X, inits = inits_list)
update(model, n.iter = 1000)
samples <- coda.samples(model, variable.names = params, n.iter = 10000)

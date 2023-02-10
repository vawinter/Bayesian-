# STAT 597
# Homework 2
# Markov Chain Monte Carlo
# Veronica A. Winter
# 02/07/2023

# clean env
rm(list = ls())
gc()
ls()

# set for reproducibility
set.seed(1544)

# install and load packages ----
#install.packages(c("mvtnorm", "invgamma", "nimble"))
library(mvtnorm)
library(invgamma)
library(nimble)
library(coda)

# Q1. ----
# Parameter values
p = 0.8
u = 3
v = 7
s2 = 2
k = rbinom(200, 1, p)
mn = u*k + v*(1 - k)
y = rnorm(length(mn), mn, sqrt(s2))
n = length(y)

# plot
hist(y, breaks = 20)

## Construct a Gibbs sampler to sample from the posterior dist'n of all parameters
## 1a. Let yi ~ N(u, o2), i = 1, 2, .., n ----
#       let u ~ N(0, 100) and o2 ~ IG(10, 100)

# starting values
mu <- 0
s2 <- 1

## set the priors
mu.mean.prior <- 0
mu.s2.prior <- 100
s2.shape.prior <- 10
s2.scale.prior <- 100

## set up 
M = 10000
mu.output = rep(NA, M)
s2.output = rep(NA, M)

### for loop ----
for(i in 1:10000){
  ## print status
  print(i)
  
  ## draw mu from full conditional
  ## this needs to calculated by hand
  ## start with our starting values and updates in each iteration
  b = sum(y)/s2
  a <- n/s2 + 1/mu.s2.prior
  a.inv <- 1/a
  
  mu = rnorm(1, mean = b*a.inv, sd = sqrt(a.inv))
  
  ## draw s2 from full conditional
  s2 = rinvgamma(1, shape = s2.shape.prior + n/2,
               rate = s2.scale.prior + sum((y - mu)^2)/2)
 
   ## save out parameters
  mu.output[i] = mu
  s2.output[i] = s2
}

## means
mean(mu.output)
# [1] 3.553679

mean(s2.output)
# [1] 4.757164

## credible intervals
quantile(mu.output, c(0.025, 0.975))
#   2.5%    97.5% 
#  3.250090 3.856807 
quantile(s2.output, c(0.025, 0.975))
#  2.5%    97.5% 
#  3.910839 5.730237 


plot(s2.output, type="l")
matplot(cbind(s2.output, mu.output), type = "l")
hist(s2.output)
hist(mu.output)


## 1b. Let y_i |k_i~N(uk_i+v(1-k_i ),σ^2 ),i=1,2,…,n -----
##      let k_i~Binom(1,0.5), u~N(0,100), v~N(0,100), and σ^2~IG(10,100).

# Clean env
#rm(lists = ls())

## starting values
mu <- 2
s2 <- 1
v <- 6
k <- rbinom(n,size=1,prob=p)

## set the priors
mu.mean.prior <- 0
mu.s2.prior <- 100
v.mean.prior <- 0
v.s2.prior <- 100
s2.shape.prior <- 10
s2.scale.prior <- 100

n <- length(y)

## number of MCMC iterations
M <- 10000

## save files
mu.ouput <- rep(NA, M)
v.output <- rep(NA, M)
s2.ouput <- rep(NA, M)
k.ouput <- matrix(NA,ncol=n,nrow=M)

### MCMC loop ----
for(i in 1:M){
  ## print iteration number
  cat(i," ")
  
  # sample k
  k = rbinom(n, 1, prob = 0.5)
  
  ## sample mu
  mu.b <- sum(k*y)/s2
  mu.a <- sum(k^2) / s2 + 1/mu.s2.prior
  mu.a.inv <- 1/mu.a
  
  mu = rnorm(1, mean = mu.b * mu.a.inv, sd = sqrt(mu.a.inv))
  
  ## sample v
  bv = (sum(y*(1 - k)))/s2
  av = sum((1 - k)^2)/s2 + 1/v.s2.prior
  v.a.inv <- 1/av
  
  v = rnorm(1,mean = bv*v.a.inv, sd = sqrt(v.a.inv))
  
  ## sample s2
  s2 = rinvgamma(1, shape = s2.shape.prior + n/2, rate = s2.scale.prior + sum((y - ((mu*k) + v*(1-k)))^2)/2)
  
  ## ouput current values
  mu.ouput[i]=mu
  v.output[i]=v
  s2.ouput[i]=s2
  k.ouput[i,]=k
  
}

## means
mean(mu.ouput)
# [1] 3.554192
mean(v.output)
# [1] 3.552307
mean(s2.ouput)
# [1] 4.751697

## credible intervals
quantile(mu.ouput, c(0.025, 0.975))
#  2.5%    97.5% 
# 3.046868 4.069802 
quantile(v.output, c(0.025, 0.975))
#  2.5%    97.5% 
# 3.042426 4.063057 
quantile(s2.ouput, c(0.025, 0.975))
#  2.5%    97.5% 
#  3.938524 5.724607 


## plot results
plot(s2.ouput, type = "l")
matplot(cbind(s2.ouput,mu.output,v.output),type="l")
matplot(cbind(mu.output,v.output),type="l")

## look at the marginal posteriors with the histograms
hist(s2.ouput)
par(mfrow=c(1,2))
hist(mu.output, prob = T)
hist(v.output, prob = T)


## 1c.Let y_i |k_i~N(uk_i+v(1-k_i ),σ^2 ),i=1,2,…,n  -----
## and let k_i~Binom(1,p), u~N(0,100), u~N(0,100), p~Unif(0,1), and σ^2~IG(10,100).
# Clean env
#rm(lists = ls())

## set the priors
mu.mean.prior <- -1 ## attempted to make the priors more informative
mu.s2.prior <- 1
v.mean.prior <- 20
v.s2.prior <- 1
s2.shape.prior <- 10
s2.scale.prior <- 100
p.lower.prior <- 0
p.upper.prior <- 1

n <- length(y)
## starting values
mu <- 5
s2 <- 1
v <- 2
p <- 0.5
k <- rbinom(
  n,
  size = 1,
  prob = p
)

M <- 10000
mu.save <- rep(NA, M)
v.save <- rep(NA, M)
s2.save <- rep(NA, M)
k.save <- matrix(NA,ncol=n,nrow=M)
p.save <- rep(NA, M)

## for loop, with 100 iterations
for(iter in 1:M){
  ## draw beta from full conditional
  ## calculation provided on additional .pdf
  ## start with our starting values and updates in each iteration
  p.k1 <- exp(-(1/2*s2) * ( (-2*mu*y) + (mu)^2 ) ) * p
  p.k2 <- exp(-(1/2*s2) * ( (-2*v*y) + v^2 ) ) * (1-p)
  prob <- p.k1 / (p.k1 + p.k2)
  k <- rbinom(
    n,
    1,
    prob = prob
  )
  
  mu.b <- sum(k*y)/s2
  mu.a <- sum(k^2) / s2 + 1/mu.s2.prior
  mu.a.inv <- 1/mu.a
  v.b <- (sum(y*(1-k)))/s2
  v.a <- sum((1-k)^2)/s2 + 1/v.s2.prior
  v.a.inv <- 1/v.a
  mu = rnorm(
    1,
    mean = mu.b*mu.a.inv,
    sd = sqrt(mu.a.inv)
  )
  v = rnorm(
    1,
    mean = v.b*v.a.inv,
    sd = sqrt(v.a.inv)
  )
  mn <- (mu*k) + v*(1-k)
  ## draw s2 from full conditional 
  s2 <- invgamma::rinvgamma(
    1,
    shape = s2.shape.prior + n/2, 
    rate = s2.scale.prior + sum((y - (mu*k + v*(1-k)))^2)/2
  )
  p <- rbeta(
    1, 
    1 + sum(k),
    1 + (n - sum(k))
  )
  
  ## save out parameters
  mu.save[iter] = mu
  s2.save[iter] = s2
  # s2.save.v[iter] = s2.v
  v.save[iter] = v
  k.save[iter, ] = k
  p.save[iter] = p
}

## plot results
plot(s2.save, type = "l")
matplot(cbind(s2.save,mu.save,v.save),type="l")
matplot(cbind(mu.save,v.save),type="l")
matplot(p.save, type = "l")

## look at the marginal posteriors with the histograms
dev.off()
par(mfrow=c(1,2))
hist(s2.save)
hist(mu.save)
hist(v.save)

# Q2. Repeat 1c in nimble ----
# The model:
## y[i] | k[i] ~ N(u*k[i] + v*(1-k[i]), s2)
## k[i] ~ Binom(1, p)
## u ~ N(0, 100)
## v ~ N(0, 100)
## s2 ~ IG(10, 100)
## p ~ Unif(0, 1) (same as ~Beta(1, 1))

# Clean env
#rm(lists = ls())

## bring in the data
set.seed(1544)
p = 0.8
u = 3
v = 7
s2 = 2
k = rbinom(200 ,1, p)
mn = u * k + v * (1 - k)
y = rnorm(length(mn), mn, sqrt(s2))
hist(y, breaks=20)


### a. write a Bayesian model statement ####
c.mod <- nimbleCode({
  for(t in 1:n){
    ## likelihood
    y[t] ~ dnorm(mean = (mu*k[t] + v*(1-k[t])), var = s2)
    k[t] ~ dbinom(size = 1, prob = p)
  }
  ## priors 
  mu ~ dnorm(mean = 0, var = 100)
  v ~ dnorm(mean = 0, var = 100)
  s2 ~ dinvgamma(shape = 10, rate = 100)
  p ~ dbeta(1, 1)
})

### b. create data, constants, initial values ####
## put any data used in the "data" list
data <- list(y = y)

## put any constants you need in the "constants" list
##  - the "n=length(x)" is necessary for the "for"-loop in the data model
##  - if you have some parameters fixed at some value, put them in "constants"
consts <- list(n = length(y))

## put a starting value for each parameter in "inits"
k.init <- rep(0, length(y))
k.init[y > mean(y)] = 1
inits <- list(mu = 0, v = 0, s2 = 1, p = 0.5, k = k.init)

### c. Run MCMC ####
#### i. Construct model 
c.model.full <- nimbleModel( ## build the model 
  c.mod,
  data = data,
  constants = consts
) 

#### ii. Configure
conf.c.model <- configureMCMC( ## configure and add monitors to watch
  c.model.full,
  monitors = c("mu","v","s2","p","k")
)

#### iii. Build MCMC
mcmc.c.model <- buildMCMC( ## build MCMC
  conf.c.model
)

#### iv. Run MCMC
c.chains <- runMCMC(
  mcmc.c.model,
  inits = inits,
  niter = 100000,
  nburnin = 60000
)

#### v. Find the effective sample size
## Note: If it is too low, increase niter above and re-run MCMC
coda::effectiveSize(c.chains)

### Results: ----
# optional - Plot MCMC chains
matplot(c.chains, type="l")
hist(c.chains)

#### Get posterior means ----
apply(c.chains,2,mean)

#### Get posterior 95% Credible Intervals ----
apply(c.chains,2,quantile,c(.025,.975))

# ### Visualization w. MCMCvis ----
# MCMCsummary(
#   c.chains,
#   exc = "k",
#   round = 2
# )
# MCMCtrace(
#   c.chains,
#   params = c("mu","v"),
#   pdf = FALSE,
#   iter = 40000
# )

# Q3a. Metropolis Hastings ----
# Let y_i be the number of sick days that a person takes due to an illness and 
# let x_i be the number of months that person has been taking part in a 
# treatment program.  
 
# Assume a model for this data as y_i∼Pois(exp⁡{a+bx_i })  - 
#   a Poisson regression model with two regression parameters a and b.  

# We will assume that both of these regression parameters have prior distributions 
# that are Gaussian with mean zero and variance equal to 4.  
# The data can be read in as follows:

x = c(8,14,11,7,32,8,28,21,27,15,26,13,19,22,15,
    12,15,7,9,15,26,22,16,12,6)

z = c(5,2,5,4,1,3,0,2,1,2,2,5,3,2,1,2,2,8,5,2,1,1,6,4,3)
    
# Using Markov Chain Monte Carlo and the Metropolis Hastings algorithm,
# draw samples from the posterior distribution of [a,b│x].  
# Use a separate MH step for each parameter.  
# Use a proposal distribution for each parameter that is normally 
# distributed and centered on the current value of the parameter.  
# Choose the variance of each proposal distribution so that your acceptance 
# rate is between 20%-50%.  Draw enough samples that your effective sample size 
# is at least 1000, and higher if possible.  
# Report the posterior mean and marginal 95% credible intervals for all parameters.

# Clean env
#rm(lists = ls())

# Notes:
## z[i] ~ Pois(a + b*x[i])
## Using Markov Chain Monte Carlo and the Metropolis Hastings algorithm, draw samples 
## from the posterior distribution of [𝑎,𝑏│𝒙].

## psi ----
## a.star ~ N(a, a.tune) 
## b.star ~ N(b, b.tune)

## Starting values
a <- 1
b <- 1

## Tuning parameters
a.tune <- 0.5
b.tune <- 0.5

## Prior distributions distributions
# a ~ Normal(0, 4)
# b ~ Normal(0, 4)

# Number of MCMC iterations
M <- 100000

# Storage file
a.output <- rep(NA, M)
b.output <- rep(NA, M)
n <- length(z)
# Number of times we accept/reject
accept.a <- 0
accept.b <- 0


## For loop w. M iterations ----
for(w in 1:M){
  
  ## 1. Metropolis Hasting sampling - a
  a.star <- rnorm(1, mean = a, sd = (a.tune))
  
  # a. Numerator
  p.mh.a.num <- sum(
    dpois( ## likelihood
      z,lambda = exp(a.star + b*x),log = TRUE)) +
    dnorm( ## prior
      a.star, mean = 0, sd = sqrt(4), log = TRUE ) +
    dnorm( ## proposal
      a, mean = a.star, sd = sqrt(a.tune), log = TRUE)
  
  # b. Denominator 
  p.mh.a.denom <- sum(
    dpois( ## likelihood
      z, lambda = exp(a + b*x), log = TRUE)) +
    dnorm( ## prior
      a, mean = 0, sd = sqrt(4),log = TRUE) +
    dnorm( ## proposal
      a.star, mean = a, sd = sqrt(a.tune), log = TRUE)
  
  ## c. Find the Pmh
  p.mh.a <- exp(p.mh.a.num - p.mh.a.denom)
 
  ## d. Accept a.star or leave a
  if(runif(1) < p.mh.a){
    accept.a <- accept.a + 1
    a <- a.star
  }
  
  ## 2. Metropolis Hasting sampling - b
  ## Repeat above for b
  
  # a. Numerator
  b.star <- rnorm(1, mean = b, sd = sqrt(b.tune))
  p.mh.b.num <- sum(
    dpois( ## likelihood
      z, lambda = exp(a + b.star*x), log = TRUE)) +
    dnorm( ## prior
      b.star, mean = 0, sd = sqrt(4), log = TRUE) +
    dnorm( ## proposal
      b,mean = b.star, sd = sqrt(b.tune), log = TRUE)
  
  # b. Denominator 
  p.mh.b.denom <- sum(
    dpois( ## likelihood
      z, lambda = exp(a + b*x), log = TRUE)) +
    dnorm( ## prior
      b, mean = 0, sd = sqrt(4), log = TRUE) +
    dnorm( ## proposal
      b.star,mean = b, sd = sqrt(b.tune), log = TRUE)
  
  # c.  Find Pmh
  p.mh.b <- exp(p.mh.b.num - p.mh.b.denom)
  
  # d. Accept/reject
  if(runif(1) < p.mh.b){
    accept.b <- accept.b + 1
    b <- b.star
  }
  
  ## 3. Save output parameters
  a.output[w] = a
  b.output[w] = b
}

## Results ----
### Plots ---- 
matplot(cbind(a.output, b.output), type = "l")

### Accept/Reject ----
accept.a
# [1] 28252
accept.b
# [1] 1437

### ESS ----
coda::effectiveSize(a.output)
# var1 
# 232.3736 
coda::effectiveSize(b.output)
# var1 
# 185.6753 
### Means ----
mean(a.output)
# [1] 2.054524
mean(b.output)
# [1] -0.06999972

### Credible intervals ----
quantile(a.output, c(0.025, 0.975))
#      2.5%    97.5% 
#   1.575839 2.623520 
quantile(b.output, c(0.025, 0.975))
#      2.5%       97.5% 
#   -0.11358840 -0.03951539 

# END: HW 2

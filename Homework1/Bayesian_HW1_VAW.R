# STAT 597
# Homework 1
# Veronica A. Winter
# 01/20/2023

# Clean env
rm(lsit = ls())
gc()

set.seed(0519)

# References: "ImportanceSampling.r" & "ImporanceSampling_pt2.r"

# Function ----
logsumexp = function(z){
  ## function to compute log(sum(exp(z)))
  c=max(z)
  c+log(sum(exp(z-c)))
}

# Question 2.
# Let x_i,i=1,…,n be the time, in minutes, between successive arrivals 
# of commercial airplanes at Heathrow airport.   

# a.Find this data in “airplane.Rdata”.  

# b.Choose a reasonable parametric model for this data.  
## Do so by considering the support (possible values) of the data, 
## and through simple exploratory data analysis.  
## Then choose a reasonable prior distribution for any parameters in the 
## parametric model you chose.  

# C.Finally, use importance sampling to draw approximate samples from the 
# posterior distribution of all parameters.  
## Make sure your effective sample size is at least 1000, 
## and higher if possible.  

# d.Report the mean and marginal 95% credible intervals for all parameters.


########### Solution: Question 2 ##############
# Let x_i,i=1,…,n be the time, in minutes, between successive arrivals 
# of commercial airplanes at Heathrow airport.   

## a. Load in the data ---- 
load("Homework1/airplane.Rdata")

## b.Choose a reasonable parametric model for this data. ----
# Let's view the data
hist(airplane)

# Histogram is all positive numbers and skewed (0, ~inf)
# An exponential could be used....
# x[i] ~ exp(u)

# With a prior of Gamma
# u ~ Gamma(shape = a, rate = b)

# Therefore:
# [u|x] ~ Gamma(a+n, b+sum(x)) # via wikipedia

# C.Importance sampling ----
# Use importance sampling to draw approximate samples from the 
# posterior distribution of all parameters.  

# Let's specify prior hyperparameters
a = 1
b = 1

# Re: above
# Data: x[i] ~ exp(u)
# Prior: u ~ Gamma(shape = a, rate = b)
# Posterior: [u|x] ~ Gamma(a+n, b+sum(x))

u.vals = seq(.01,qgamma(.99,shape=a,rate=b),by=.01)

prior.vals=dgamma(u.vals, shape=a, rate=b)

post.vals=dgamma(u.vals,shape=a+length(airplane),rate=b+sum(airplane))

plot(u.vals,post.vals,type="l",col="red",lwd=3)

points(u.vals,prior.vals,type="l",col="blue",lwd=3,lty=2)

points(u.vals,post.vals,type="l",col="red",lwd=3)

legend("topright",legend=c("Prior","Posterior"),
          lwd=c(3,3),col=c("blue","red"),lty=c(2,1))  


# Now, draw M RVs from the prior
M = 100000
u.star = rgamma(M,shape=a,rate=b)


# calculate weights
log.q=rep(NA,M)

# Loop over data
for(i in 1:M){
  log.q[i]=sum(dexp(airplane, u.star[i],
                     log=TRUE))
}

## Denom
log.denom = logsumexp(log.q)

## Weights
log.w=log.q-log.denom
w=exp(log.w)
summary(w)
sum(w)

## Effective sample size
ESS = 1/sum(w^2)
ESS
# [1] 7926.4

## Resampling
N = 10000
u = sample(u.star, size = N, replace = TRUE, prob = w)
hist(u, prob = TRUE, breaks = 30, main = paste ("ESS =", round(ESS)))

# Finally, let's plot the results.....
hist(u,prob=TRUE,breaks=30,main=paste("ESS=",round(ESS)),xlim=c(0,2))
u.vals=seq(.1,70,by=.01)
prior.vals=dgamma(u.vals,shape=a,rate=b)
post.vals=dgamma(u.vals,shape=a+length(airplane),rate=b+sum(airplane))
points(u.vals,prior.vals,type="l",col="blue",lwd=3,lty=2)
points(u.vals,post.vals,type="l",col="red",lwd=3)
legend("topright",legend=c("SIS Histogram","Prior","Posterior"),
       lwd=c(10,3,3),col=c("grey","blue","red"),lty=c(1,2,1))


# D. Report the mean and marginal 95% credible intervals for all parameters. 
a.post = a+length(airplane)
b.post = b+sum(airplane)

post.mean.analytic = a.post/b.post
# [1] 0.3064645

# ANSWERS ----
# Mean: Analytics
mean(u)
# 0.3066594

# Mean: Using importance sampling weights
sum(w*u)
# 0.306537

# CI: Analytics
qgamma(c(.025,.975),shape=a.post,rate=b.post)
# 0.2496204 0.3690524

## CI: Using samples
quantile(u,c(.025,.975))
#     2.5%      97.5% 
# 0.2501406 0.3706314 

########### End Question 2 ##############

########### Solution: Question 3 ##############
# Let y_i be the number of sick days that a person takes due to an illness, 
# and let x_i be the number of months that person has been taking part in a 
# treatment program.  

# Assume a model for this data as y_i∼Pois(exp⁡{a+bx_i }) 
# -  a Poisson regression model with two regression parameters a and b.  
# We will assume that both of these regression parameters have prior 
# distributions that are Gaussian with mean zero and variance equal to 4. 

# Data = Poisson w/ 2 regression paramaters
# y_i∼Pois(exp⁡{a+bx_i })
## with priors
## b1~N(0, 4)
## b2~N(0, 4)

## Format the data ----
x = c(8,14,11,7,32,8,28,21,27,15,26,13,19,22,15,12,15,7,9,15,26,22,16,12,6)
y = c(5,2,5,4,1,3,0,2,1,2,2,5,3,2,1,2,2,8,5,2,1,1,6,4,3)

## [1] draw M RVs from the prior (normal) ----
M = 10^6
a.star = rnorm(M, 0, sqrt(4)) #a param
b.star = rnorm(M, 0, sqrt(4)) #b param

## [2] calculate weights ----
log.q = rep(NA, M)

for(m in 1:M){
  p = exp(a.star[m] + b.star[m]*x) 
  log.q[m]=sum(dpois(y, lambda = p, log=TRUE))
}

log.denom = logsumexp(log.q)

log.w = log.q-log.denom
w = exp(log.w)
summary(w)
sum(w)

## [2] Check effective sample size (ESS) ----
ESS = 1/sum(w^2)
ESS
# [1] 639.9817

## resampling
N = 500
## resample the INDEX (draw from the prior for all three parameters)
idx = sample(1:M,size = N, replace = TRUE, prob = w)

a = a.star[idx]
b = b.star[idx]

hist(a)
hist(b)

## [3a.] Posterior inference ----
post.mat = cbind(a, b)
post.mean = apply(post.mat,2,mean)
post.CI = apply(post.mat,2,quantile,c(.025,.975))

# ANSWER: Find posterior mean ----
post.mean
#       a          b 
# 2.06594295 -0.07113701 

# ANSWER: Marginal 95% CI ----
post.CI
#              a           b
# 2.5%  1.454949 -0.10781549
# 97.5% 2.587724 -0.03287357


# Using the posterior predictive distribution, provide a 95% prediction interval 
# on the number of sick days that might be taken in a year for a person who has 
# been on the treatment program for 24 months.  
mn.ppd_1 = exp(a + b*24)
ci.24 <- quantile(mn.ppd_1, c(0.025, 0.0975))
#       2.5%    9.75% 
#   0.8833687 1.0671182 
mean.24 <- mean(mn.ppd_1)
# [1] 1.472093

hist(mn.ppd_1,main="Posterior Predictive Mean: Number of sick days, 24 mo.")
summary(mn.ppd_1) 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.4831  1.2110  1.4510  1.4721  1.7039  2.8938  

# Provide a similar 95% prediction interval for a person 
# who has been on the treatment program for 36 months
mn.ppd_2 = exp(a + b*36)

ci.36 <- quantile(mn.ppd_2, c(0.025, 0.0975))
#       2.5%     9.75% 
#   0.2570704 0.3476915  
mean.36 <- mean(mn.ppd_2)
# [1] 0.6715566
hist(mn.ppd_2, main="Posterior Predictive Mean: Number of sick days, 36 mo.")
summary(mn.ppd_2)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.08509 0.45122 0.61119 0.67156 0.82567 1.92369  


### END HW 1 ###

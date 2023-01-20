# STAT 597
# Homework 1
# Veronica A. Winter
# 01/20/2023

# Clean env
rm(lsit = ls())
gc()

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
ESS=1/sum(w^2)
ESS

## Resampling
N=10000
u=sample(u.star,size=N,replace=TRUE,prob=w)
hist(u,prob=TRUE,breaks=30,main=paste("ESS=",round(ESS)))

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

# Mean: Analytics
mean(u)
# 0.3065115

# Mean: Using importance sampling weights
sum(w*u)
# 0.3062683

# CI: Analytics
qgamma(c(.025,.975),shape=a.post,rate=b.post)
# 0.2496204 0.3690524

## CI: Using samples
quantile(u,c(.025,.975))
#     2.5%      97.5% 
# 0.2491419  0.3700360

########### End Question 2 ##############

########### Solution: Question 3 ##############
# Let y_i be the number of sick days that a person takes due to an illness, 
# and let x_i be the number of months that person has been taking part in a 
# treatment program.  

# Assume a model for this data as y_i∼Pois(exp⁡{a+bx_i }) 
# -  a Poisson regression model with two regression parameters a and b.  
# We will assume that both of these regression parameters have prior 
# distributions that are Gaussian with mean zero and variance equal to 4. 

# Data = Poisson
# Prior = gaussian
m = 0
v = 4

#The data can be read in as follows:
  
x = c(8,14,11,7,32,8,28,21,27,15,26,13,19,22,15,12,15,7,9,15,26,22,16,12,6)
y = c(5,2,5,4,1,3,0,2,1,2,2,5,3,2,1,2,2,8,5,2,1,1,6,4,3)

# Using importance sampling, draw samples from the posterior distribution of 
# [a,b│x].  
# 

# Draw enough samples that your effective sample size is at least 500, 
# and higher if possible.  Report the posterior mean and marginal 95% credible 
# intervals for all parameters.

S = 500

# Using the posterior predictive distribution, provide a 95% prediction interval 
# on the number of sick days that might be taken in a year for a person who has 
# been on the treatment program for 24 months.  

# Provide a similar 95% prediction interval for a person 
# who has been on the treatment program for 36 months





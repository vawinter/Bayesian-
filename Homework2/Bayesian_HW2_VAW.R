# Homework 2
# Markov Chain Monte Carlos

# clean env
rm(list = ls())
gc()
ls()

# set for reroducibility
set.seed(1544)

# install and load packages ----
#install.packages(c("mvtnorm", "invgamma"))
library(mvtnorm)
library(invgamma)

# Question 1 ----
# Parameter values
p = 0.8
u = 3
v = 7
s2.true = 2
k = rbinom(200, 1, p)
mn = u*k + v*(1 - k)
y = rnorm(length(mn), mn, sqrt(s2.true))

# plot
hist(y, breaks = 20)

## Construct a Gibbs sampler to sample from the posterior dist'n of all parameters
## a. Let yi ~ N(u, o2), i = 1, 2, .., n ----
#       let u ~ N(0, 100) and o2 ~ IG(10, 100)

## Gibbs Sampler ----
# starting values
s2 = 1
mu = 0

## set up 
M = 1000
mu.output = matrix(NA, M, length(mu))
s2.output = rep(NA, M)

## for loop ----
for(i in 1:1000){
  ## print status
  print(i)
  
  ## draw mu from full conditional
  ## this needs to calculated by hand
  ## start with our starting values and updates in each iteration
  b = (sum(y) * (1/s2))
  A = (length(y)/s2)+ (1/100)
  
  mu = as.numeric(rnorm(1, mean = b/(A), sd = sqrt(1/A)))
  
  ## draw s2 from full conditional
  s2 = rinvgamma(1, shape = 10 + (length(y)/2),
               scale = 100 + 1/2 * sum((y - mu)^2))
 
   ## save out parameters
  mu.output[i] = mu
  s2.output[i] = s2
}

plot(s2.output, type="l")
matplot(cbind(s2.output, mu.output), type = "l")
hist(s2.output)
hist(mu.output)


## b. Let y_i |k_i~N(uk_i+v(1-k_i ),σ^2 ),i=1,2,…,n -----
##      let k_i~Binom(1,0.5), u~N(0,100), v~N(0,100), and σ^2~IG(10,100).
















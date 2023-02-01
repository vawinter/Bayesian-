########################################################
##
## Template Code for MCMC
##
########################################################


library(nimble)
library(coda)

####################################
## write a Bayesian model statement
####################################

model.code=nimbleCode({
    for(i in 1:n){
        ## Put data model here
    }
    ## put prior model for all parameters here
})


## put any data used in the "data" list
data=list(x=x)
## put any constants you need in the "constants" list
##  - the "n=length(x)" is necessary for the "for"-loop in the data model
##  - if you have some parameters fixed at some value, put them in "constants"
constants=list(n=length(x))
## put a starting value for each parameter in "inits"
inits=list(lambda=1)

###################################
## Run MCMC
###################################
## The code below can be run without changes, except for:
##   monitors=c(...) - replace ... with the names of all parameters,
##                     separated by commas, like monitors=c("mu","sigma")
##   niter=1000 is the number of MCMC iterations.  Increase if needed

model=nimbleModel(model.code,data=data,constants=constants)
## model=compileNimble(model) ## uncomment if you can compile C++ code
conf <- configureMCMC(model, monitors = c("lambda"))
mcmc <- buildMCMC(conf)
## mcmc=compileNimble(mcmc,project=model) ## uncomment if you can compile C++ code
chain=runMCMC(mcmc,inits=inits,niter=1000)

## find the effective sample size.
## If it is too low, increase niter above and re-run MCMC
effectiveSize(chain)

##############################
## optional - Plot MCMC chains
##############################
matplot(chain,type="l")
hist(chain)

## find the effective sample size
library(coda)
effectiveSize(chain)

##############################
## Get posterior means and CIs
##############################

## find posterior means
apply(chain,2,mean)

## find posterior 95% Credible Intervals
apply(chain,2,quantile,c(.025,.975))






########################################################
##
## A Poisson data example
##
########################################################

library(nimble)
library(coda)



##
## simulate data from y ~ Pois(5)
##

lambda.true=5

set.seed(123)
x=rpois(100,lambda=lambda.true)
x



##
## write a Bayesian model statement
##

model.code=nimbleCode({
    ## likelihood (in a for-loop for each data point)
    for(i in 1:n){
        x[i] ~ dpois(lambda)
    }
    ## prior distribution for all parameters
    lambda ~ dgamma(shape=1,rate=2)
})


data=list(x=x) ## put data in here
constants=list(n=length(x)) ## put any constants you need - typically need the data length
inits=list(lambda=1) ## put starting values for all parameters here

model=nimbleModel(model.code,data=data,constants=constants)
model=compileNimble(model)
conf <- configureMCMC(model, monitors = c("lambda"))
mcmc <- buildMCMC(conf)
mcmc=compileNimble(mcmc,project=model)
chain=runMCMC(mcmc,inits=inits,niter=1000)


matplot(chain,type="l")
hist(chain)

## find the effective sample size
library(coda)
effectiveSize(chain)

## find posterior means
apply(chain,2,mean)

## find posterior 95% Credible Intervals
apply(chain,2,quantile,c(.025,.975))



########################################################
##
## 2 parameter example (Normal data) with Nimble
##
########################################################

library(nimble)

##
## simulate data from y~N(5,16)
##

set.seed(1234)
x=rnorm(100,5,sqrt(16))
x

####################################
## write a Bayesian model statement
####################################

model.code=nimbleCode({
    for(i in 1:n){
        ## Put data model here
        x[i]~dnorm(mean=M,sd=sqrt(V))
    }
    ## put prior model for all parameters here
    ## prior for mean
    M~dnorm(0,1)
    V~dexp(1/10)
})


## put any data used in the "data" list
data=list(x=x)
## put any constants you need in the "constants" list
##  - the "n=length(x)" is necessary for the "for"-loop in the data model
##  - if you have some parameters fixed at some value, put them in "constants"
constants=list(n=length(x))
## put a starting value for each parameter in "inits"
inits=list(M=1,V=1)

###################################
## Run MCMC
###################################
## The code below can be run without changes, except for:
##   monitors=c(...) - replace ... with the names of all parameters,
##                     separated by commas, like monitors=c("mu","sigma")
##   niter=1000 is the number of MCMC iterations.  Increase if needed

model=nimbleModel(model.code,data=data,constants=constants)
model=compileNimble(model) ## uncomment if you can compile C++ code
conf <- configureMCMC(model, monitors = c("M","V"))
mcmc <- buildMCMC(conf)
mcmc=compileNimble(mcmc,project=model) ## uncomment if you can compile C++ code
chain=runMCMC(mcmc,inits=inits,niter=10000)

## find the effective sample size.
## If it is too low, increase niter above and re-run MCMC
effectiveSize(chain)

##############################
## optional - Plot MCMC chains
##############################
matplot(chain,type="l")

##############################
## Get posterior means and CIs
##############################

## find posterior means
apply(chain,2,mean)

## find posterior 95% Credible Intervals
apply(chain,2,quantile,c(.025,.975))



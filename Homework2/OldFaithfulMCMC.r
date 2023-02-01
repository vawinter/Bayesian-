##
##
## Old Faithful data
##
##

data(faithful)
faithful
attach(faithful)
hist(waiting,breaks=25,col="yellow")
x=waiting

##
## model:
##

n.mcmc

# MCMC prep
lambdaH.save = rep(NA, n.mcmc)
lambdaL.save = rep(NA, n.mcmc)
p.save = rep(NA, n.mcmc)
k.save = rep(NA, ncol = n, nrow = n.mcmc)

####################################
## write a Bayesian model statement
####################################
library(nimble)
library(coda)

model.code=nimbleCode({
    for(i in 1:n){
        ## Put data model here
        x[i]~dnorm(mean=mu1*k[i]+mu2*(1-k[i]),var=v1*k[i]+v2*(1-k[i]))
        k[i]~dbinom(size=1,prob=p)
    }
    ## put prior model for all parameters here
    p~dbeta(1,1)
    mu1~dnorm(mean=0,var=100)
    mu2~dnorm(mean=0,var=100)
    v1~dexp(1/10)
    v2~dexp(1/10)
})


## put any data used in the "data" list
data=list(x=x)
## put any constants you need in the "constants" list
##  - the "n=length(x)" is necessary for the "for"-loop in the data model
##  - if you have some parameters fixed at some value, put them in "constants"
constants=list(n=length(x))
## put a starting value for each parameter in "inits"
inits=list(mu1=0,mu2=0,p=.5,v1=1,v2=2,k=k.init)

k.init=rep(0,length(x))
k.init[x>mean(x)]=1

###################################
## Run MCMC
###################################
## The code below can be run without changes, except for:
##   monitors=c(...) - replace ... with the names of all parameters,
##                     separated by commas, like monitors=c("mu","sigma")
##   niter=1000 is the number of MCMC iterations.  Increase if needed

model=nimbleModel(model.code,data=data,constants=constants)
 model=compileNimble(model) ## uncomment if you can compile C++ code
conf <- configureMCMC(model, monitors = c("p","mu1","mu2","v1","v2"))
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


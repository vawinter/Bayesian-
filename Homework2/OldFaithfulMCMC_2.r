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



##
## MCMC prep
##

## starting values
lambdaH=100
lambdaL=50
p=0.5
n=length(x)
k=rbinom(n,size=1,prob=p)

## number of MCMC iterations
n.mcmc=1000

## save files
lambdaH.save=rep(NA,n.mcmc)
lambdaL.save=rep(NA,n.mcmc)
p.save=rep(NA,n.mcmc)
k.save=matrix(NA,ncol=n,nrow=n.mcmc)

## MCMC loop
for(iter in 1:n.mcmc){
    ## print iteration number
    cat(iter," ")

    ## sample k
    ## prob.k1=lambdaH^x*exp(-lambdaH)*p/(lambdaH^x*exp(-lambdaH)*p+lambdaL^x*exp(-lambdaL)*(1-p))
    ## k=rbinom(n,size=1,prob=prob.k1)

    for(i in 1:n){
      
        prob.k1=lambdaH^x[i]*exp(-lambdaH)*p/(lambdaH^x[i]*exp(-lambdaH)*p+lambdaL^x[i]*exp(-lambdaL)*(1-p))
        k[i]=rbinom(1,size=1,prob=prob.k1)
        
    }
  
    ## sample p
    p=rbeta(1,shape1=1+sum(k),shape2=1+n-sum(k))

    ## sample lambdaL
    lambdaL=rgamma(1,shape=500+sum(x*(1-k)),rate=10+n-sum(k))

    ## sample lambdaH
    lambdaH=rgamma(1,shape=1000+sum(x*(k)),rate=10+sum(k))

    ## save current values
    lambdaH.save[iter]=lambdaH
    lambdaL.save[iter]=lambdaL
    p.save[iter]=p
    k.save[iter,]=k

}

par(mfrow=c(2,2))
plot(p.save,main="p",type="l")
plot(lambdaH.save,main="lambdaH",type="l")
plot(lambdaL.save,main="lambdaL",type="l")
matplot((k.save[,1]),main="k[1] ",type="l")
matplot((k.save[,2]),main="k[2] ",type="l")
matplot((k.save[,58]),main="k[58]" ,type="l")

hist(p.save)
hist(lambdaH.save)
hist(lambdaL.save)


####################################
##now the same model in nimble
####################################
library(nimble)
library(coda)

model.code=nimbleCode({
    for(i in 1:n){
        ## Put data model here
        x[i]~dpois(lambda=mu1*k[i]+mu2*(1-k[i]))
        k[i]~dbinom(size=1,prob=p)
    }
    ## put prior model for all parameters here
    p~dbeta(1,1)
    mu1~dgamma(1000,10)
    mu2~dgamma(500,10)
})


## put any data used in the "data" list
data=list(x=x)
## put any constants you need in the "constants" list
##  - the "n=length(x)" is necessary for the "for"-loop in the data model
##  - if you have some parameters fixed at some value, put them in "constants"
constants=list(n=length(x))
## put a starting value for each parameter in "inits"
k.init=rep(0,length(x))
k.init[x>mean(x)]=1
inits=list(mu1=100,mu2=50,p=.5,k=k.init)


###################################
## Run MCMC
###################################
## The code below can be run without changes, except for:
##   monitors=c(...) - replace ... with the names of all parameters,
##                     separated by commas, like monitors=c("mu","sigma")
##   niter=1000 is the number of MCMC iterations.  Increase if needed

model=nimbleModel(model.code,data=data,constants=constants)
 model=compileNimble(model) ## uncomment if you can compile C++ code
conf <- configureMCMC(model, monitors = c("p","mu1","mu2"))
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


## compare to by-hand code:
mean(lambdaH.save)
quantile(lambdaH.save,c(.025,.975))

mean(lambdaL.save)
quantile(lambdaL.save,c(.025,.975))

mean(p.save)
quantile(p.save,c(.025,.975))


## plotting histograms of true data and data simulated from posterior mean
par(mfrow=c(1,2))
hist(x)

p.hat=mean(p.save)
lambdaH.hat=mean(lambdaH.save)
lambdaL.hat=mean(lambdaL.save)
k.sim=rbinom(n,size=1,prob=p)
lambda.sim=rep(lambdaL.hat,n)
lambda.sim[k.sim==1]=lambdaH.hat
x.sim=rpois(n,lambda=lambda.sim)
hist(x.sim,main="Simulated From posterior")

####
#### LM.examples.r
####

##
## 0. Read in datasets
## 1. Standard Bayesian Linear Model
## 2. Bayesian Ridge Regression
## 3. Bayesian LASSO Regression
## 4. Bayesian Robust Regression (post mode = minimize absolute error)

#############################################
##
## [0] Datasets
##
## diabetes data
##
## source:      Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle
##              Regression" (with discussion) _Annals of Statistics_
##
###############################################

library(coda)
library(mvtnorm)

load("diab.Rdata")
str(diab)
## note how correlated some of the predictors are
pairs(diab)
cor(diab)

X.diab=cbind(1,as.matrix(diab[,-1]))
y.diab=diab[,1]



#############################################
##
## #1 Linear Model with Ridge Regression priors
##
## y~N(X*beta,s2*I)
## beta~N(0,beta.prior.var*I)
## s2~IG(a,b)
##
#############################################

source("mcmc.lm.r")

out=mcmc.lm(y.diab,X.diab,n.mcmc=10000)
summary(out)
matplot(out[,-12],type="l")

## diagnostics
effectiveSize(out)
library(batchmeans)
bmmat(out)




#######################################################
##
## (4) Bayesian Ridge Regresssion
##
#######################################################




#############################################
##
## #2 Linear Model with Ridge Regression priors
##
## y~N(u+X*beta,s2*I)
## u~N(0,beta.prior.var)
## beta~N(0,k2*I)
## k2~Exp(exp.rate)
## s2~IG(a,b)
##
#############################################

source("mcmc.lm.rr.r")
source("ShabyWellsLogAdapt.r")
out.rr=mcmc.lm.rr(y.diab,X.diab,n.mcmc=10000,
                  beta.prior.var=100,s2.prior.var=var(y.diab),k2.exp.rate=100,
                  tune.k2=.01,tune.s2=.01,
                  adapt.iter=100)


## diagnostics

## trace plots
par(mfrow=c(3,4))
for(k in 2:ncol(out.rr)){
  plot(out.rr[,k],main=colnames(out.rr)[k],type="l")
  abline(v=0,col="red")
}

burnin=5000
effectiveSize(out.rr[-c(1:burnin),])

for(k in 2:ncol(bhm.out)){
  hist(out.rr[-c(1:burnin),k],main=colnames(out.rr)[k])
  abline(v=0,col="red")
}



##
## Same thing in nimble!
##

library(nimble)


N=nrow(X.diab)
## dim of data
M=ncol(X.diab)

bhm.data=list(y=y.diab,X=X.diab)
bhm.const=list(N=N,M=M)
bhm.inits=list(s2=1,k2=1,beta=rep(0,M))

## model
rr.code <- nimbleCode({
  ## priors for variance params
  s2 ~ T(dnorm(0,sd=77),min=0,max=Inf)
  k2~dexp(.1)
  ## k2 ~ T(dnorm(0,sd=1),min=0,max=Inf)
  ## priors for regression parameters
  beta[1]~ dnorm(0,var=100)
  for(i in 2:M){
    beta[i]~ dnorm(0,var=k2)
  }
  ## data model
  for(i in 1:N){
    y[i]~dnorm(inprod(X[i,1:M],beta[1:M]),var=s2)
  }
})



## run model
## compile model and run MCMC
bhm.out <- nimbleMCMC(rr.code,
                      constants=bhm.const,
                      data=bhm.data,
                      inits=bhm.inits,
                      niter=10000,
                      thin=1,
                      monitors=c("beta","s2","k2"))

summary(bhm.out)
apply(bhm.out,2,mean)
apply(bhm.out,2,mean)-apply(out,2,mean)

burnin=5000
for(k in 2:ncol(bhm.out)){
  hist(bhm.out[-c(1:burnin),k],main=colnames(bhm.out)[k])
  abline(v=0,col="red")
}

for(k in 2:ncol(bhm.out)){
  plot(bhm.out[,k],main=colnames(bhm.out)[k],type="l")
  abline(v=0,col="red")
}



#############################################
##
## #3 Linear Model with LASSO priors
##
## y~N(u+X*beta,s2*I)
## u~N(0,beta.prior.var)
## beta[k]~N(0,s2*tau2[k])
## tau2[k]~exp(0.5*lambda)
## lambda~Exp(exp.rate)
## s2~HalfNorm(0,var(data))
##
#############################################


library(statmod) ## needed for inverse-gaussian distribution
library(nimble)
source("mcmc.lm.lasso.r")
source("ShabyWellsLogAdapt.r")
out.lasso=mcmc.lm.lasso(y.diab,X.diab,n.mcmc=10000,
                  beta.prior.var=100,s2.prior.var=var(y.diab),lam2.exp.rate=.10,
                  tune.s2=.01,
                  adapt.iter=100)


## diagnostics

## trace plots
par(mfrow=c(3,4))
for(k in 2:ncol(out.lasso)){
  plot(out.lasso[,k],main=colnames(out.lasso)[k],type="l")
  abline(v=0,col="red")
}

burnin=2000
effectiveSize(out.lasso[-c(1:burnin),])

for(k in 2:ncol(out.lasso)){
  hist(out.lasso[-c(1:burnin),k],main=colnames(out.lasso)[k])
  abline(v=0,col="red")
}


##
## same thing in nimble!
##

N=nrow(X.diab)
## dim of data
M=ncol(X.diab)

bhm.data=list(y=y.diab,X=X.diab)
bhm.const=list(N=N,M=M)
bhm.inits=list(s2=1,lambda=1,beta=rep(0,M),tau2=rep(1,M-1))

## model
lasso.code <- nimbleCode({
  ## priors for variance params
  s2 ~ T(dnorm(0,sd=77),min=0,max=Inf)
  lambda ~ dexp(1)
  ## priors for regression parameters
  beta[1]~ dnorm(0,var=100)
  ## data augmentation version
  for(i in 2:M){
    beta[i]~ dnorm(0,var=tau2[i-1])
    tau2[i-1] ~ dexp(0.5*lambda)
  }
  ## ## not augmented version
  ## for(i in 2:M){
  ##   beta[i]~ ddexp(location=0,location=lambda)
  ## }
  ## data model
  for(i in 1:N){
    y[i]~dnorm(inprod(X[i,1:M],beta[1:M]),var=s2)
  }
})



## run model
## compile model and run MCMC
bhm.out <- nimbleMCMC(lasso.code,
                      constants=bhm.const,
                      data=bhm.data,
                      inits=bhm.inits,
                      niter=10000,
                      thin=1,
                      monitors=c("beta","s2","lambda"))

apply(bhm.out,2,mean)

burnin=5000
for(k in 2:ncol(bhm.out)){
  hist(bhm.out[-c(1:burnin),k],main=colnames(bhm.out)[k])
  abline(v=0,col="red")
}










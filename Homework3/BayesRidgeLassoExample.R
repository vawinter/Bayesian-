####
#### RidgeLassoExample.r
####
#### Outline:
####  1. Simulating data and splitting into train / test sets
####  2. AIC stepwise selection
####  3. Ridge regression
####  4. Lasso regression
####  5. Ridge / Lasso for GLMs (logistic regression example)


#######################################################
##
## (1) Simulating 50 predictors and one response
##     v1-v10 are important
##     v11-v50 are not
##
#######################################################

## generate a dummy dataset with 30 predictors (10 useful & 20 useless)
## note that each predictor has mean zero and sd=1
N=500
P=50

X=matrix(NA,nrow=N,ncol=P)


##
## do this for correlated predictors
##
covmat=matrix(rnorm(P^2,sd=2),nrow=P)
covmat=covmat+t(covmat)
U=eigen(covmat)$vectors
D=diag(rexp(P,rate=10))
covmat=U%*%D%*%t(U)

##
## Do this for uncorrelated predictors
##
covmat=diag(P)



set.seed(1234)
library(mvtnorm)
for(i in 1:N){
    X[i,]=rmvnorm(1,mean=rep(0,P),sigma=covmat)
}
X=data.frame(X)
head(X)

## true betas
betas.true=c(1,2,3,4,5,-1,-2,-3,-4,-5,rep(0,P-10))
betas.true

## simulating "y"
sigma=15.7
mu=1
X=as.matrix(X)
y=mu+X%*%betas.true+rnorm(N,mean=0,sd=sigma)


## splitting into "train" and "test" data

alldata=data.frame(cbind(y,X))
names(alldata)[1] <- "y"
head(alldata)
train=alldata[1:350,]
test=alldata[351:500,]

## y  <- rnorm(N)
## x1 <- matrix(rnorm(N*20),N,20)
## x2 <- matrix(y+rnorm(N*10),N,10)
## x  <- cbind(x1,x2)


## linear regression
fit=lm(y~.,data=train)
summary(fit)

betas.lm=coef(fit)

yhat.lm=predict(fit,newdata=test)
mspe.lm=mean((test$y-yhat.lm)^2)
mspe.lm

## AIC stepwise selection
fit.AIC=step(fit)
summary(fit.AIC)

yhat.AIC=predict(fit.AIC,newdata=test)
mspe.AIC=mean((test$y-yhat.AIC)^2)
mspe.AIC

#######################################################
##
## ridge regression
##
#######################################################

library(glmnet)
## alpha=0 gives ridge regression
## alpha=1 gives lasso regression

## fit ridge (trying 100 different lambda values)
rr=glmnet(x=as.matrix(train[,-1]),y=as.numeric(train[,1]),alpha=0,nlambda=100)
plot(rr,xvar="lambda",main="Ridge Regression Betas for Different Values of the Tuning Parameter")

## use 10-fold crossvalidation to find the best lambda
cv.rr=cv.glmnet(x=as.matrix(train[,-1]),y=as.numeric(train[,1]),alpha=0,nfolds=10,nlambda=100)

## get lambda and best rr fit
lambda.rr=cv.rr$lambda.min
lambda.rr

## some plots
par(mfrow=c(1,2))
plot(cv.rr)
abline(v=log(lambda.rr))
plot(rr,xvar="lambda",main="Ridge Regression Betas for Different Values of the Tuning Parameter")
abline(v=log(lambda.rr))

## beta estimates for best lambda
betas.rr=coef(cv.rr,s="lambda.min")

plot(betas.rr,betas.lm,xlim=c(-6,6),ylim=c(-6,6))
abline(0,1)
plot(betas.rr,c(0,betas.true),xlim=c(-6,6),ylim=c(-6,6))
abline(0,1)

yhat.rr=predict(cv.rr,s="lambda.min",newx=as.matrix(test[,-1]))
mspe.rr=mean((test$y-yhat.rr)^2)
mspe.rr



#######################################################
##
## (3) lasso regression
##
#######################################################

## alpha=0 gives ridge regression
## alpha=1 gives lasso regression

## fit lasso (trying 100 different lambda values)
lasso=glmnet(x=as.matrix(train[,-1]),y=as.numeric(train[,1]),alpha=1,nlambda=100)
plot(lasso,xvar="lambda",main="Lasso Regression Betas for Different Values of the Tuning Parameter")
plot(rr,xvar="lambda",main="Ridge Regression Betas for Different Values of the Tuning Parameter")

## use 10-fold crossvalidation to find the best lambda
cv.lasso=cv.glmnet(x=as.matrix(train[,-1]),y=as.numeric(train[,1]),alpha=1,nfolds=10)

## get lambda and best lasso fit
lambda.lasso=cv.lasso$lambda.min
lambda.lasso

## some plots
par(mfrow=c(1,2))
plot(cv.lasso)
abline(v=log(lambda.lasso))
plot(lasso,xvar="lambda")
abline(v=log(lambda.lasso))

## beta estimates for best lambda
betas.lasso=coef(cv.lasso,s="lambda.min")
betas.lasso

plot(betas.lasso,betas.lm,xlim=c(-6,6),ylim=c(-6,6))
abline(0,1)
plot(betas.lasso,c(0,betas.true),xlim=c(-6,6),ylim=c(-6,6))
abline(0,1)

yhat.lasso=predict(cv.lasso,newx=as.matrix(test[,-1]),s="lambda.min")
mspe.lasso=mean((test$y-yhat.lasso)^2)
mspe.lasso



##
## Comparison of MSPE on test set
##

mspe.lm
mspe.AIC
mspe.rr
mspe.lasso






#############################################
##
## Bayesian Linear Model
##
## y~N(u+X*beta,s2*I)
## u~N(0,beta.prior.var)
## beta~N(0,beta.prior.var*I)
## s2~HalfNorm(0,var(y))
##
#############################################

y.train=as.numeric(train[,1])
X.train=cbind(1,as.matrix(train[,-1]))


source("mcmc.lm.r")

out=mcmc.lm(y.train,X.train,n.mcmc=10000)
summary(out)
matplot(out,type="l")

## diagnostics
effectiveSize(out)
library(batchmeans)
bmmat(out)

burnin=100

colnames(out)
betas.post=out[-c(1:burnin),2:51]
mean.post=apply(betas.post,2,mean)
CI.post=apply(betas.post,2,quantile,c(.025,.975))
par(mfrow=c(1,1))
plot(0,0,type="n",xlim=c(1,50),ylim=c(-7,7),main="Posterior 95% CIs (iid Gaussian Priors)")
abline(h=0,col="green")
for(i in 1:length(mean.post)){
  points(c(i,i),CI.post[,i],type="l",)
  points(i,mean.post[i],pch=3)
  points(i,betas.true[i],col="red")
}


#############################################
##
## #2 Linear Model with Ridge Regression priors
##
## y~N(u+X*beta,s2*I)
## u~N(0,beta.prior.var)
## beta~N(0,k2*I)
## k2~Exp(exp.rate)
## s2~HalfNorm(0,var(y))
##
#############################################

source("mcmc.lm.rr.r")
source("ShabyWellsLogAdapt.r")
out.rr=mcmc.lm.rr(y.train,X.train,n.mcmc=10000,
                  beta.prior.var=100,s2.prior.var=var(y.train),k2.exp.rate=1,
                  tune.k2=.01,tune.s2=.01,
                  adapt.iter=100)


## diagnostics

## trace plots
par(mfrow=c(3,4))
for(k in 2:ncol(out.rr)){
  plot(out.rr[,k],main=colnames(out.rr)[k],type="l")
  abline(v=0,col="red")
}

par(mfrow=c(3,4))
for(k in 2:ncol(out.rr)){
  hist(out.rr[,k],main=colnames(out.rr)[k])
  abline(v=0,col="red")
}

burnin=2000
effectiveSize(out.rr[-c(1:burnin),])


colnames(out.rr)
betas.post=out.rr[-c(1:burnin),2:51]
mean.post=apply(betas.post,2,mean)
CI.post=apply(betas.post,2,quantile,c(.025,.975))
par(mfrow=c(1,1))
plot(0,0,type="n",xlim=c(1,50),ylim=c(-7,7),main="Posterior 95% CIs (Ridge Priors)")
abline(h=0,col="green")
for(i in 1:length(mean.post)){
  points(c(i,i),CI.post[,i],type="l",)
  points(i,mean.post[i],pch=3)
  points(i,betas.true[i],col="red")
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
out.lasso=mcmc.lm.lasso(y.train,X.train,n.mcmc=10000,
                  beta.prior.var=100,s2.prior.var=var(y.train),lam2.exp.rate=1,
                  tune.s2=.01,
                  adapt.iter=100)


## diagnostics

## ## trace plots
## par(mfrow=c(3,4))
## for(k in 2:ncol(out.lasso)){
##   plot(out.lasso[,k],main=colnames(out.lasso)[k],type="l")
##   abline(v=0,col="red")
## }

burnin=2000
effectiveSize(out.lasso[-c(1:burnin),])

for(k in 2:ncol(out.lasso)){
  hist(out.lasso[-c(1:burnin),k],main=colnames(out.lasso)[k])
  abline(v=0,col="red")
}





betas.post=out.lasso[-c(1:burnin),2:51]
mean.post=apply(betas.post,2,mean)
CI.post=apply(betas.post,2,quantile,c(.025,.975))
par(mfrow=c(1,1))
plot(0,0,type="n",xlim=c(1,50),ylim=c(-7,7),main="Posterior 95% CIs (LASSO Priors)")
abline(h=0,col="green")
for(i in 1:length(mean.post)){
  points(c(i,i),CI.post[,i],type="l",)
  points(i,mean.post[i],pch=3)
  points(i,betas.true[i],col="red")
}






library(nimble)


N=nrow(X.train)
## dim of data
M=ncol(X.train)

bhm.data=list(y=y.train,X=X.train)
bhm.const=list(N=N,M=M)
bhm.inits=list(s2=1,lambda=1,beta=rep(0,M),tau2=rep(1,M-1))

## model
lasso.code <- nimbleCode({
  ## priors for variance params
  s2 ~ T(dnorm(0,sd=sqrt(344)),min=0,max=Inf)
  lambda ~ dexp(.1)
  ## priors for regression parameters
  beta[1]~ dnorm(0,var=100)
  for(i in 2:M){
    beta[i]~ dnorm(0,var=tau2[i-1])
    tau2[i-1] ~ dexp(0.5*lambda)
  }
  ## data model
  for(i in 1:N){
    y[i]~dnorm(inprod(X[i,1:M],beta[1:M]),var=s2)
  }
})



## run model
## compile model and run MCMC
bhm.out.lasso <- nimbleMCMC(lasso.code,
                      constants=bhm.const,
                      data=bhm.data,
                      inits=bhm.inits,
                      niter=10000,
                      thin=1,
                      monitors=c("beta","s2","lambda"))

# apply(bhm.out,2,mean)
#
# burnin=5000
# for(k in 2:ncol(bhm.out)){
#   hist(bhm.out[-c(1:burnin),k],main=colnames(bhm.out)[k])
#   abline(v=0,col="red")
# }

betas.post=bhm.out.lasso[-c(1:burnin),2:51]
mean.post.lasso=apply(betas.post,2,mean)
CI.post=apply(betas.post,2,quantile,c(.025,.975))
par(mfrow=c(1,1))
plot(0,0,type="n",xlim=c(1,50),ylim=c(-7,7),main="Posterior 95% CIs (LASSO Priors)")
abline(h=0,col="green")
for(i in 1:length(mean.post)){
  points(c(i,i),CI.post[,i],type="l",)
  points(i,mean.post[i],pch=3)
  points(i,betas.true[i],col="red")
}














library(nimble)

X=cbind(1,as.matrix(train[,-1]))
colnames(X)[1]="int"
y=as.numeric(train[,1])

N=nrow(X)
## dim of data
M=ncol(X)

bhm.data=list(y=y,X=X)
bhm.const=list(N=N,M=M)
bhm.inits=list(s2=1,beta=rep(0,M))

## model
rr.code <- nimbleCode({
  ## data model
  for(i in 1:N){
    y[i]~dnorm(inprod(X[i,1:M],beta[1:M]),var=s2)
  }
  ## priors for regression parameters
  for(i in 1:M){
    beta[i]~ dnorm(0,var=100)
  }
  ## priors for variance params
  s2 ~ T(dnorm(0,sd=77),min=0,max=Inf)
})



## run model
## compile model and run MCMC
bhm.out <- nimbleMCMC(rr.code,
                      constants=bhm.const,
                      data=bhm.data,
                      inits=bhm.inits,
                      niter=10000,
                      thin=1,
                      monitors=c("beta","s2"))

summary(bhm.out)
apply(bhm.out,2,mean)

burnin=5000
for(k in 2:ncol(bhm.out)){
  hist(bhm.out[-c(1:burnin),k],main=colnames(bhm.out)[k])
  abline(v=0,col="red")
}
for(k in 2:ncol(bhm.out)){
  plot(bhm.out[,k],main=colnames(bhm.out)[k],type="l")
  abline(v=0,col="red")
}


colnames(bhm.out)
betas.post=bhm.out[,2:51]
mean.post=apply(betas.post,2,mean)
CI.post=apply(betas.post,2,quantile,c(.025,.975))
par(mfrow=c(1,1))
plot(0,0,type="n",xlim=c(1,50),ylim=c(-7,7),main="Posterior 95% CIs (Standard Bayesian Regression)")
abline(h=0,col="green")
for(i in 1:length(mean.post)){
  points(c(i,i),CI.post[,i],type="l",)
  points(i,mean.post[i],pch=3)
  points(i,betas.true[i],col="red")
}




#######################################################
##
## (4) Bayesian Ridge Regresssion
##
#######################################################




#############################################
##
## Linear Model with Ridge Regression priors
##
## y~N(u+X*beta,s2*I)
## u~N(0,beta.prior.var)
## beta~N(0,k2*I)
## k2~Exp(exp.rate)
## s2~IG(a,b)
##
#############################################


library(nimble)

X=cbind(1,as.matrix(train[,-1]))
colnames(X)[1]="int"
y=as.numeric(train[,1])

N=nrow(X)
## dim of data
M=ncol(X)

bhm.data=list(y=y,X=X)
bhm.const=list(N=N,M=M)
bhm.inits=list(s2=1,k2=1,beta=rep(0,M))

## model
rr.code <- nimbleCode({
  ## data model
  for(i in 1:N){
    y[i]~dnorm(inprod(X[i,1:M],beta[1:M]),var=s2)
  }
  ## priors for regression parameters
  beta[1]~ dnorm(0,var=100)
  for(i in 2:M){
    beta[i]~ dnorm(0,var=k2)
  }
  k2 ~ T(dnorm(0,sd=1),min=0,max=Inf)
  ## priors for variance params
  s2 ~ T(dnorm(0,sd=77),min=0,max=Inf)
})



## run model
## compile model and run MCMC
bhm.out.rr <- nimbleMCMC(rr.code,
                      constants=bhm.const,
                      data=bhm.data,
                      inits=bhm.inits,
                      niter=10000,
                      thin=1,
                      monitors=c("beta","s2","k2"))

# summary(bhm.out)
# apply(bhm.out,2,mean)
#
# burnin=5000
# for(k in 2:ncol(bhm.out)){
#   hist(bhm.out[-c(1:burnin),k],main=colnames(bhm.out)[k])
#   abline(v=0,col="red")
# }
# for(k in 2:ncol(bhm.out)){
#   plot(bhm.out[,k],main=colnames(bhm.out)[k],type="l")
#   abline(v=0,col="red")
# }


colnames(bhm.out.rr)
betas.post=bhm.out.rr[,2:51]
mean.post=apply(betas.post,2,mean)
CI.post=apply(betas.post,2,quantile,c(.025,.975))
par(mfrow=c(1,1))
plot(0,0,type="n",xlim=c(1,50),ylim=c(-7,7),main="Posterior 95% CIs (Ridge Priors)")
abline(h=0,col="green")
for(i in 1:length(mean.post)){
  points(c(i,i),CI.post[,i],type="l",)
  points(i,mean.post[i],pch=3)
  points(i,betas.true[i],col="red")
}



#############################################
##
## #3 Linear Model with LASSO priors
##
## y~N(u+X*beta,s2*I)
## u~N(0,beta.prior.var)
## beta[k]~N(0,tau2[k])
## tau2[k]~exp(0.5*lambda)
## lambda~Exp(exp.rate)
## s2~HalfNorm(0,var(data))
##
#############################################


library(nimble)


N=nrow(X)
## dim of data
M=ncol(X)

bhm.data=list(y=y,X=X)
bhm.const=list(N=N,M=M)
bhm.inits=list(s2=1,lambda=1,beta=rep(0,M),tau2=rep(1,M-1))

## model
lasso.code <- nimbleCode({
  ## priors for variance params
  s2 ~ T(dnorm(0,sd=77),min=0,max=Inf)
  lambda ~ dexp(1)
  ## priors for regression parameters
  beta[1]~ dnorm(0,var=100)
  for(i in 2:M){
    beta[i]~ dnorm(0,var=tau2[i-1])
    tau2[i-1] ~ dexp(0.5*lambda)
  }
  ## data model
  for(i in 1:N){
    y[i]~dnorm(inprod(X[i,1:M],beta[1:M]),var=s2)
  }
})



## run model
## compile model and run MCMC
bhm.out.lasso <- nimbleMCMC(lasso.code,
                      constants=bhm.const,
                      data=bhm.data,
                      inits=bhm.inits,
                      niter=10000,
                      thin=1,
                      monitors=c("beta","s2","lambda"))

# apply(bhm.out,2,mean)
#
# burnin=5000
# for(k in 2:ncol(bhm.out)){
#   hist(bhm.out[-c(1:burnin),k],main=colnames(bhm.out)[k])
#   abline(v=0,col="red")
# }

betas.post=bhm.out.lasso[,2:51]
mean.post=apply(betas.post,2,mean)
CI.post=apply(betas.post,2,quantile,c(.025,.975))
par(mfrow=c(1,1))
plot(0,0,type="n",xlim=c(1,50),ylim=c(-7,7),main="Posterior 95% CIs (LASSO Priors)")
abline(h=0,col="green")
for(i in 1:length(mean.post)){
  points(c(i,i),CI.post[,i],type="l",)
  points(i,mean.post[i],pch=3)
  points(i,betas.true[i],col="red")
}

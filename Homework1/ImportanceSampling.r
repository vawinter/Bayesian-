## this is a function to compute log(sum(exp)), which is numerically unstable
## unless you use this trick!

logsumexp=function(z){
  ## function to compute log(sum(exp(z)))
  c=max(z)
  c+log(sum(exp(z-c)))
}





##
## Importance sampling
##

## assume the following data are x_i ~ Pois(u)

x=c(2,0,3,3,2,2,1,6,2,2,2,8,0,6,1,3,1,3,3,2,2,1,0,1,1)

## if our Prior on u is u~Gamma(shape=a,rate=b), let's use importance sampling
##   to get samples from the posterior

## specify prior hyperparameters
a=1
b=1

## since we can get the posterior analytically, let's compare the prior and posterior
##  (try changing the prior hyperparameters a and b, and see what that does to the posterior)
##
## if prior is u~Gamma(shape=a,rate=b)
## and x[i]~Pois(u)
## then u|x~Gamma(a+sum(x),b+n)

u.vals=seq(.01,qgamma(.99,shape=a,rate=b),by=.01)
prior.vals=dgamma(u.vals,shape=a,rate=b)
post.vals=dgamma(u.vals,shape=a+sum(x),rate=b+length(x))
plot(u.vals,post.vals,type="l",col="red",lwd=3)
points(u.vals,prior.vals,type="l",col="blue",lwd=3,lty=2)
points(u.vals,post.vals,type="l",col="red",lwd=3)
legend("topright",legend=c("Prior","Posterior"),
       lwd=c(3,3),col=c("blue","red"),lty=c(2,1))

## [1] draw M RVs from the prior
M=100000
u.star=rgamma(M,shape=a,rate=b)

## (2) calculate weights
log.q=rep(NA,M)
for(i in 1:M){
  log.q[i]=sum(dpois(x,u.star[i],log=TRUE))
}

log.denom=logsumexp(log.q)

log.w=log.q-log.denom
w=exp(log.w)
summary(w)
sum(w)

ESS=1/sum(w^2)
ESS

## resampling
N=10000
u=sample(u.star,size=N,replace=TRUE,prob=w)
hist(u,prob=TRUE,breaks=30,main=paste("ESS=",round(ESS)))

## compare to analytic posterior
## if prior is u~Gamma(shape=a,rate=b)
## and x[i]~Pois(u)
## then u|x~Gamma(a+sum(x),b+n)

hist(u,prob=TRUE,breaks=30,main=paste("ESS=",round(ESS)),xlim=c(0,20))
u.vals=seq(.1,70,by=.01)
prior.vals=dgamma(u.vals,shape=a,rate=b)
post.vals=dgamma(u.vals,shape=a+sum(x),rate=b+length(x))
points(u.vals,prior.vals,type="l",col="blue",lwd=3,lty=2)
points(u.vals,post.vals,type="l",col="red",lwd=3)
legend("topright",legend=c("SIS Histogram","Prior","Posterior"),
       lwd=c(10,3,3),col=c("grey","blue","red"),lty=c(1,2,1))




##
## Inference examples
##

##
## (1) Posterior mean of u?
##

## (1a) Analytic
a.post=a+sum(x)
b.post=b+length(x)
post.mean.analytic=a.post/b.post
post.mean.analytic
## (1b) Using samples
mean(u)
## (1c) Using importance sampling weights
sum(w*u)

##
## (2) 95% Credible Interval of u (95% posterior probability interval)
##

## (2a) Analytic
qgamma(c(.025,.975),shape=a.post,rate=b.post)
## 2(b) Using samples
quantile(u,c(.025,.975))

##
## (3) Inference on a transformation
##       Let v=log(u)

## (3a) Analytic calculation is challenging (but do-able - left as exercise for student)

## (3b) Using samples
## transform samples to v
v=log(u)
## get posterior mean
post.mean.v=mean(v)
post.mean.v
## get posterior CI
CI.v=quantile(v,c(.025,.975))
CI.v

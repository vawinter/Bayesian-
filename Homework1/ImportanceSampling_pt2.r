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
M=1000
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
N=100000
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
sum(w*u.star)

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

####################################################
##
## 2-parameter example

## read in data

df <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
str(df)
df$rank=as.factor(df$rank)
summary(df)
attach(df)

pairs(df)

## consider a simple binary regression model:
## admit[i] ~ bern(p[i])
## p[i] = 1/1+exp(-( b0+b1*gre[i]+b2*gpa[i]))
## with priors
## b0~N(0,10)
## b1~N(0,10)
## b2~N(0,10)

## Goal: estimate probability of admission for someone with GPA=3.5 and GRE=760


## [1] draw M RVs from the prior
M=100000
b0.star=rnorm(M,0,sqrt(10))
b1.star=rnorm(M,0,sqrt(.001))
b2.star=rnorm(M,0,sqrt(1))

## (2) calculate weights
log.q=rep(NA,M)
for(m in 1:M){
  p=1/(1+exp(-(b0.star[m]+b1.star[m]*gre+b2.star[m]*gpa)))
  log.q[m]=sum(dbinom(admit,size=1,prob=p,log=TRUE))
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
## resample the INDEX (draw from the prior for all three parameters)
idx=sample(1:M,size=N,replace=TRUE,prob=w)
b0=b0.star[idx]
b1=b1.star[idx]
b2=b2.star[idx]
hist(b0)
hist(b1)
hist(b2)

## posterior inference
post.mat=cbind(b0,b1,b2)
post.mean=apply(post.mat,2,mean)
post.CI=apply(post.mat,2,quantile,c(.025,.975))
post.mean
post.CI

##
## Now get posterior predictive distribution of admission of someone
##   with GRE=760 and GPA=3.5

## calculate "p" for this person for each posterior draw
## result is vector of (p_1,p_2,...p_N) with each p_i~[p|y]
p=1/(1+exp(-(b0+b1*760+b2*3.5)))
## histogram of this posterior (posterior predictive dist'n of [p|y])
hist(p)

## now get posterior predictive dist'n of [admit|y]
admit.ppd=rbinom(length(p),size=1,prob=p)
hist(admit.ppd)
summary(admit.ppd)

##################################################################
##
## (1) Teacher Pay
##    "PAY" = Average public school teacher annual salary ($)
##    "SPEND" = Spending on public schools per pupil
##    "AREA" = "1" = Northeast/North Central, "2"=South, "2"=West
##
## Example goal: estimate what avg teacher salary should be if
##               state spending is raised to 10000 per student
##
##################################################################

pay=read.csv("TeacherPay.csv")
head(pay)

plot(pay$SPEND,pay$PAY,xlab="State Spending Per Student",ylab="Avg Teacher Salary")

## regression model
## pay_i ~ N(u+SPEND_i*b,s^2)
##
## priors
## u~N(10000,4000)
## b~N(0,10)
## s~Exp(mean=1000)

attach(pay)
## simulate data from the prior
M=10^5
u.star=rnorm(M,10000,sqrt(4000))
b.star=rnorm(M,0,sqrt(10))
s.star=rexp(M,1/1000)

## (2) calculate weights
log.q=rep(NA,M)
for(m in 1:M){
  log.q[m]=sum(dnorm(PAY,mean=u.star[m]+b.star[m]*SPEND,sd=s.star[m],log=TRUE))
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
## resample the INDEX (draw from the prior for all three parameters)
idx=sample(1:M,size=N,replace=TRUE,prob=w)
u=u.star[idx]
b=b.star[idx]
s=s.star[idx]
hist(u)
hist(b)
hist(s)

## posterior inference
post.mat=cbind(u,b,s)
post.mean=apply(post.mat,2,mean)
post.CI=apply(post.mat,2,quantile,c(.025,.975))
post.mean
post.CI

##
## Now posterior predictive inference on [z|y] where z~N(u+b*10000,s^2)

## get posterior predictive means
mn.ppd=u+b*10000
hist(mn.ppd,main="Posterior Predictive Mean of Avg. Teacher Pay If SPEND=10k")

## get posterior predictive Avg Teacher Pay
atp.ppd=rnorm(length(mn.ppd),mean=mn.ppd,sd=s)
hist(atp.ppd,main="Posterior Predictive Dist'n of Avg. Teacher Pay If SPEND=10k")

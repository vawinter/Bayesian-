##
## Simulate normal data
##

beta.true=c(1,-2)
s2.true=1.5

##install.packages("mvtnorm")
## install.packages("invgamma")
library(invgamma)
library(mvtnorm)
x=rnorm(50)
X=cbind(1,x)
X
y=rnorm(50,X%*%beta.true,sd=sqrt(s2.true))
plot(x,y)

##
## Gibbs Sampler
##

## starting values
beta=c(-10,-10)
s2=100

M=1000
beta.save=matrix(NA,M,length(beta))
s2.save=rep(NA,M)

## for loop
for(iter in 1:1000){
  ## draw beta from full conditional
  b=t(X)%*%y/s2
  A=diag(2)/100+t(X)%*%X/s2
  A.inv=solve(A)
  beta=as.numeric(rmvnorm(1,mean=A.inv%*%b,sigma=A.inv))
  ## draw s2 from full conditional
  s2=rinvgamma(1,shape=10+length(y)/2,
               rate=10+1/2*t(y-X%*%beta)%*%(y-X%*%beta))
  ## save out parameters
  beta.save[iter,]=beta
  s2.save[iter]=s2
}

plot(s2.save,type="l")
matplot(cbind(s2.save,beta.save),type="l")
hist(s2.save)
hist(beta.save[,2])


######################################################
##
## Hidden Markov Model Example
##


##
## Simulating data (motivation is in-nest ant trophallaxis events)
##

P.true=matrix(c(.9,.1,0,.1,.8,.1,0,.2,.8),nrow=3,byrow=TRUE)
P.true
lambda.true=c(.3,2,7)
set.seed(127)
S=1
for(t in 1:99){
    S[t+1]=sample(1:3,size=1,prob=P.true[S[t],])
}
plot(S,type="b")
y=rpois(length(S),lambda.true[S])
plot(y,type="b")
plot(cumsum(y))


##
## MCMC using nimble
##
library(nimble)
hmm.code=nimbleCode({
    ## prior for lambda
    ## note our prior for lambdas leads to identifiable lambdas with no label switching
    for(g in 1:n.states){
        lambda.tilde[g]~dexp(1)
    }
    lambda[1] <- lambda.tilde[1]
    lambda[2] <- lambda.tilde[1]+lambda.tilde[2]
    lambda[3] <- lambda.tilde[1]+lambda.tilde[2]+lambda.tilde[3]
    ## prior for P
    for(g in 1:n.states){
        P[g,1:n.states]~ddirch(alpha[1:n.states])
    }
    ## prior for initial state
    X[1]~dcat(init.prior[1:n.states])
    ## Model for state space at each time point
    for(t in 2:T){
        X[t]~dcat(P[X[t-1],1:n.states])
    }
    ## get likelihood of data, conditioned on state space
    for(t in 1:T){
        y[t]~dpois(lambda[X[t]])
    }
})


## data and constants
bhm.data=list(y=y)
bhm.const=list(T=length(y),
               n.states=nrow(P.true),
               alpha=c(1,1,1), ## prior hyperparameters for each row of P
               init.prior=c(1/3,1/3,1/3) ## prior mean state at X[1]
               )
bhm.inits=list(X=rep(1,100),lambda.tilde=c(1,1,1),P=matrix(1/3,nrow=3,ncol=3))




## run model
## compile model and run MCMC
bhm.out <- nimbleMCMC(hmm.code,
                      constants=bhm.const,
                      data=bhm.data,
                      inits=bhm.inits,
                      niter=10000,
                      thin=1,
                      monitors=c("P","lambda","X"))

colnames(bhm.out)
matplot(bhm.out[,110:112],type="l",main="trace plots for rate parameters")
abline(h=lambda.true,col=1:3,lwd=3)

matplot(bhm.out[,c(1,4,7)],type="l",main="trace plots for P[1,1:3]")
abline(h=as.numeric(P.true)[c(1,4,7)],col=1:3,lwd=3)

matplot(bhm.out[,c(2,5,8)],type="l",main="trace plots for P[2,1:3]")
abline(h=as.numeric(P.true)[c(2,5,8)],col=1:3,lwd=3)

matplot(bhm.out[,c(3,6,9)],type="l",main="trace plots for P[3,1:3]")
abline(h=as.numeric(P.true)[c(3,6,9)],col=1:3,lwd=3)

P.hat=matrix(apply(bhm.out[,1:9],2,mean),nrow=3,ncol=3)
P.hat
P.true

plot(S,type="b")
S.hat=apply(bhm.out[,10:109],2,mean)
points(S.hat,col="red",type="l")



##########################################################
##
## State Space Model Example
##

## simulate an OU process

u=74
b=.06
s2=3
z0=70
g1=40
g2=75
T=300
z=rep(z0,T+1)
for(t in 1:T){
    z[t+1]=z[t]-b*(z[t]-u)+rnorm(1,mean=0,sd=sqrt(s2))
    if(t==100){
        z[t+1]=z[t+1]+g1
    }
    if(t==200){
        z[t+1]=z[t+1]+g2
    }
}
plot(z,type="l")
abline(h=u,col="red")

y=rnorm(length(z),mean=z,sd=4)
plot(y)
## thin out y
y[-seq(10,300,by=10)]=NA
plot(y,col="red")
points(z,type="l")



##
## MCMC using nimble
##
library(nimble)
ssm.code=nimbleCode({
    ## observed data
    for(i in 1:n){
        y[i]~dnorm(z[t.obs[i]],sd=sqrt(tau2))
    }
    ## OU state space process
    for(t in 2:(T)){
        z[t]~dnorm(z[t-1]-beta*(z[t-1]-u)+g2*e200[t-1]+g1*e100[t-1],sd=sqrt(s2))
    }
    ## prior for t=1
    z[1]~dnorm(0,sd=100)
    ## prior for parameters
    s2~dinvgamma(10,100)
    tau2~dinvgamma(10,100)
    beta~dnorm(0,sd=100)
    g1~dnorm(0,sd=100)
    g2~dnorm(0,sd=100)
    u~dnorm(0,sd=100)
})

## data for nimble model
t.obs=which(y>-Inf)
y.obs=y[t.obs]
T=300
e100=rep(0,T)
e100[100]=1
e200=rep(0,T)
e200[200]=1

## data and constants
bhm.data=list(y=y.obs,
              e100=e100,
              e200=e200
)
bhm.const=list(T=T,
               t.obs=t.obs,
               n=length(y.obs)
               )
bhm.inits=list(z=rep(100,300),u=100,g1=50,g2=50,s2=5,tau2=5,beta=1)




## run model
## compile model and run MCMC
bhm.out <- nimbleMCMC(ssm.code,
                      constants=bhm.const,
                      data=bhm.data,
                      inits=bhm.inits,
                      niter=10000,
                      thin=1,
                      monitors=c("u","beta","g1","g2","tau2","s2","z"))

cn=colnames(bhm.out)

idx.beta=which(cn=="beta")
plot(bhm.out[,idx.beta],type="l",main="beta")
abline(h=b,col="red")

idx.u=which(cn=="u")
plot(bhm.out[,idx.u],type="l",main="u")
abline(h=u,col="red")

idx.z=6:305
z.out=bhm.out[,idx.z]
zmean=apply(z.out,2,mean)
z.ci.upper=apply(z.out,2,quantile,.975)
z.ci.lower=apply(z.out,2,quantile,.025)


plot(y,col="red",pch=20,cex=3)
points(z,type="l",lwd=2)
points(zmean,type="l",col="blue")
points(z.ci.upper,type="l",col="blue",lty=2)
points(z.ci.lower,type="l",col="blue",lty=2)

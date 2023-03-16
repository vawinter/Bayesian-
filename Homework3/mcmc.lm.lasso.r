mcmc.lm.lasso=function(y,X,n.mcmc,
                    beta.start=rep(0,ncol(X)),s2.start=1,lam2.start=1,
                    beta.prior.var=100,s2.prior.var=var(y),lam2.exp.rate=1,
                    tune.s2=.01,
                    adapt.iter=100){
  beta.save=matrix(NA,n.mcmc,ncol(X))
  colnames(beta.save)=colnames(X)
  s2.save=rep(NA,n.mcmc)
  lam2.save=rep(NA,n.mcmc)
  n=length(y)
  p=ncol(X)-1
  beta=beta.start
  lam2=lam2.start
  s2=s2.start
  accept.lam2=0
  accept.s2=0
  adapt.t=0
  for(iter in 1:n.mcmc){
    if(iter%%100==0){
      cat(iter," ")
    }

    ## update tau2
    #tau2=1/rinvgauss(p,mean=sqrt(lam2)*sqrt(s2)/abs(beta[-1]),shape=lam2)
    tau2=1/rinvgauss(p,mean=sqrt(lam2)/abs(beta[-1]),shape=lam2)

    ## update lam2
    lam2=rgamma(1,shape=p+1 , rate=lam2.exp.rate+1/2*sum(tau2))

    ## update s2 (RW proposal on log(s2))
    s2.star=exp(rnorm(1,log(s2),tune.s2))
    mh1=(sum(dnorm(y, X%*%beta ,sqrt(s2.star), log=TRUE))
        +sum(dexp(tau2,lam2,log=TRUE))
        ##+dnorm(s2.star,0,sqrt(s2.prior.var),log=TRUE)
        +log(s2.star)
        )
    mh2=(sum(dnorm(y,X%*%beta,sqrt(s2),log=TRUE))
        +sum(dexp(tau2,lam2,log=TRUE))
        ##+dnorm(s2,0,sqrt(s2.prior.var),log=TRUE)
        +log(s2)
        )
    if(runif(1)<exp(mh1-mh2)){
      s2=s2.star
      accept.s2=accept.s2+1
    }

    ## ## update lam2 (RW proposal on log(lam2))
    ## lam2.star=exp(rnorm(1,log(lam2),tune.lam2))
    ## mh1=sum(dnorm(beta[-1],0,sqrt(lam2.star),log=TRUE))+dnorm(lam2.star,0,lam2.exp.rate,log=TRUE)+log(lam2.star)
    ## mh2=sum(dnorm(beta[-1],0,sqrt(lam2),log=TRUE))+dnorm(lam2,0,lam2.exp.rate,log=TRUE)+log(lam2)
    ## if(runif(1)<exp(mh1-mh2)){
    ##   lam2=lam2.star
    ##   accept.lam2=accept.lam2+1
    ## }

    ## update beta
    D.inv=diag(1/c(beta.prior.var,tau2))
    A.inv=solve(t(X)%*%X/s2+D.inv)
    b=1/s2*t(X)%*%y
    beta=t(rmvnorm(1,A.inv%*%b,A.inv))


    ## save
    beta.save[iter,]=beta
    s2.save[iter]=s2
    lam2.save[iter]=lam2

    ##
    ## Adaptive tuning
    ##
    if(iter%%adapt.iter==0){
        ## move adapt counter up 1
        adapt.t=adapt.t+1
        ## ## new tuning parameters for lam2
        ## adapt.vals=get.sigma(tune.lam2,1,data=lam2.save[(iter-adapt.iter)+1:adapt.iter],accept=accept.lam2,t.adapt=adapt.t)
        ## tune.lam2=adapt.vals$s2
        ## new tuning parameters for s2
        adapt.vals=get.sigma(tune.s2,1,data=s2.save[(iter-adapt.iter)+1:adapt.iter],accept=accept.s2,t.adapt=adapt.t)
        tune.s2=adapt.vals$s2
        ## resetting acceptances to 0
        accept.lam2=0
        accept.s2=0
    }


  }
  cat("\n")
  ##output
  out=cbind(beta.save,s2.save,lam2.save)
  colnames(out)=c(colnames(X),"s2","lam2")
  data.frame(out)
}


mcmc.lm.lasso2=function(y,X,n.mcmc,
                       beta.start=rep(0,ncol(X)),s2.start=1,lam2.start=1,
                       beta.prior.var=100,s2.prior.var=var(y),lam2.exp.rate=1,
                       tune.s2=.01,
                       adapt.iter=100, lambda=NULL){
  
  beta.save=matrix(NA,n.mcmc,ncol(X))
  colnames(beta.save)=colnames(X)
  s2.save=rep(NA,n.mcmc)
  lam2.save=rep(NA,n.mcmc)
  n=length(y)
  p=ncol(X)-1
  beta=beta.start
  lam2=lam2.start
  s2=s2.start
  accept.lam2=0
  accept.s2=0
  adapt.t=0
  
  if(!is.null(lambda)){
    lam2=lambda/2
  }
  
  for(iter in 1:n.mcmc){
    if(iter%%100==0){
      cat(iter," ")
    }
    
    ## update tau2
    tau2=1/rinvgauss(p,mean=sqrt(lam2)/abs(beta[-1]),shape=lam2)
    
    ## update lam2
    lam2=rgamma(1,shape=p+1 , rate=lam2.exp.rate+1/2*sum(tau2))
    
    ## update s2 (RW proposal on log(s2))
    s2.star=exp(rnorm(1,log(s2),tune.s2))
    mh1=(sum(dnorm(y,X%*%beta,sqrt(s2.star),log=TRUE))
         +sum(dexp(tau2,lam2,log=TRUE))
         +log(s2.star)
    )
    mh2=(sum(dnorm(y,X%*%beta,sqrt(s2),log=TRUE))
         +sum(dexp(tau2,lam2,log=TRUE))
         +log(s2)
    )
    if(runif(1)<exp(mh1-mh2)){
      s2=s2.star
      accept.s2=accept.s2+1
    }
    
    ## update beta
    D.inv=diag(1/c(beta.prior.var,tau2))
    A.inv=solve(t(X)%*%X/s2+D.inv)
    b=1/s2*t(X)%*%y
    beta=t(rmvnorm(1,A.inv%*%b,A.inv))
    
    
    ## save
    beta.save[iter,]=beta
    s2.save[iter]=s2
    lam2.save[iter]=lam2
    
    ##
    ## Adaptive tuning
    ##
    if(iter%%adapt.iter==0){
      ## move adapt counter up 1
      adapt.t=adapt.t+1
      ## new tuning parameters for s2
      adapt.vals=get.sigma(tune.s2,1,data=s2.save[(iter-adapt.iter)+1:adapt.iter],accept=accept.s2,t.adapt=adapt.t)
      tune.s2=adapt.vals$s2
      ## resetting acceptances to 0
      accept.lam2=0
      accept.s2=0
    }
    
    
  }
  cat("\n")
  ##output
  out=cbind(beta.save,s2.save,lam2.save)
  colnames(out)=c(colnames(X),"s2","lam2")
  data.frame(out)
}

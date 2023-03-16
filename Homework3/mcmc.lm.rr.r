mcmc.lm.rr=function(y,X,n.mcmc,
                    beta.start=rep(0,ncol(X)),s2.start=1,k2.start=1,
                    beta.prior.var=100,s2.prior.var=var(y),k2.exp.rate=1,
                    tune.k2=.01,tune.s2=.01,
                    adapt.iter=100){
  beta.save=matrix(NA,n.mcmc,ncol(X))
  colnames(beta.save)=colnames(X)
  s2.save=rep(NA,n.mcmc)
  k2.save=rep(NA,n.mcmc)
  n=length(y)
  p=ncol(X)-1
  beta=beta.start
  k2=k2.start
  s2=s2.start
  accept.k2=0
  accept.s2=0
  adapt.t=0
  for(iter in 1:n.mcmc){
    if(iter%%100==0){
      cat(iter," ")
    }
    ## update s2 (RW proposal on log(s2))
    s2.star=exp(rnorm(1,log(s2),tune.s2))
    mh1=sum(dnorm(y,X%*%beta,sqrt(s2.star),log=TRUE))+dnorm(s2.star,0,sqrt(s2.prior.var),log=TRUE)+log(s2.star)
    mh2=sum(dnorm(y,X%*%beta,sqrt(s2),log=TRUE))+dnorm(s2,0,sqrt(s2.prior.var),log=TRUE)+log(s2)
    if(runif(1)<exp(mh1-mh2)){
      s2=s2.star
      accept.s2=accept.s2+1
    }

    ## update k2 (RW proposal on log(k2))
    k2.star=exp(rnorm(1,log(k2),tune.k2))
    mh1=sum(dnorm(beta[-1],0,sqrt(k2.star),log=TRUE))+dnorm(k2.star,0,k2.exp.rate,log=TRUE)+log(k2.star)
    mh2=sum(dnorm(beta[-1],0,sqrt(k2),log=TRUE))+dnorm(k2,0,k2.exp.rate,log=TRUE)+log(k2)
    if(runif(1)<exp(mh1-mh2)){
      k2=k2.star
      accept.k2=accept.k2+1
    }
    ## update beta
    D.inv=diag(1/c(beta.prior.var,rep(k2,p)))
    A.inv=solve(t(X)%*%X/s2+D.inv)
    b=1/s2*t(X)%*%y
    beta=t(rmvnorm(1,A.inv%*%b,A.inv))
    ## save
    beta.save[iter,]=beta
    s2.save[iter]=s2
    k2.save[iter]=k2

    ##
    ## Adaptive tuning
    ##
    if(iter%%adapt.iter==0){
        ## move adapt counter up 1
        adapt.t=adapt.t+1
        ## new tuning parameters for k2
        adapt.vals=get.sigma(tune.k2,1,data=k2.save[(iter-adapt.iter)+1:adapt.iter],accept=accept.k2,t.adapt=adapt.t)
        tune.k2=adapt.vals$s2
        ## new tuning parameters for s2
        adapt.vals=get.sigma(tune.s2,1,data=s2.save[(iter-adapt.iter)+1:adapt.iter],accept=accept.s2,t.adapt=adapt.t)
        tune.s2=adapt.vals$s2
        ## resetting acceptances to 0
        accept.k2=0
        accept.s2=0
    }


  }
  cat("\n")
  ##output
  out=cbind(beta.save,s2.save,k2.save)
  colnames(out)=c(colnames(X),"s2","k2")
  data.frame(out)
}

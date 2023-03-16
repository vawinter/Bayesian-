mcmc.lm=function(y,X,n.mcmc,beta.start=rep(0,ncol(X)),s2.start=1,beta.prior.var=100,a=10,b=100){
  beta.save=matrix(NA,n.mcmc,ncol(X))
  colnames(beta.save)=colnames(X)
  s2.save=rep(NA,n.mcmc)
  n=length(y)
  p=ncol(X)
  beta=beta.start
  for(iter in 1:n.mcmc){
    if(iter%%100==0){
      cat(iter," ")
    }
    ## update s2
    s2=1/rgamma(1,a+n/2,b+1/2*sum((y-X%*%beta)^2))
    ## update beta
    A.inv=solve(t(X)%*%X/s2+diag(p)/beta.prior.var)
    b=1/s2*t(X)%*%y
    beta=t(rmvnorm(1,A.inv%*%b,A.inv))
    ## save
    beta.save[iter,]=beta
    s2.save[iter]=s2
  }
  cat("\n")
  ##output
  out=cbind(beta.save,s2.save)
  colnames(out)=c(colnames(X),"s2")
  data.frame(out)
}

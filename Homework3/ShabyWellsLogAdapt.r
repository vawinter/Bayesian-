get.sigma=function(s2.tune,Sigma.tune,data,accept,t.adapt,c0=1,c1=0.8){
    ##
    ## code to implement log-adaptive tuning following Shaby and Wells
    ##
    ## s2.tune=current variance scaling
    ## Sigma.tune=current Covariance matrix of samples
    ## data = matrix or vector of the most recent "k" MCMC iterations
    ## accept = number of MH acceptances in the most recent "k" iterations
    ## t.adapt = the number of times, including this one, that we have
    ##           applied the adaptive tuning procedure.
    ##
    if(is.matrix(Sigma.tune)){
        ## case for multiple params
        k=nrow(data) ## number of mcmc iterations
    }else{
        k=length(data)
    }
    r.hat=accept/k ## empirical acceptance prob
    S.hat=var(data)
    gamma1=1/t.adapt^c1
    gamma2=c0*gamma1
    s2.new=exp(log(s2.tune)+gamma2*(r.hat-.234))
    Sigma.new=Sigma.tune+gamma1*(S.hat-Sigma.tune)
    list(s2=s2.new,Sigma=Sigma.new)
}

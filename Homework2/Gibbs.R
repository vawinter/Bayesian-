##
## simulate normal data
##

##load packages
library(mvtnorm)
install.packages("invgamma")
library(invgamma)
## true parameters
beta.true <- c(1,-2)
s2.true = 1.5

## draw random values 
n <- 50
x=rnorm(50)
X=cbind(1,x)
X
y=rnorm(50,X%*%beta.true, sd = sqrt(s2.true))

plot(x, y)

## 
## gibbs sampler
##

## starting values
beta = c(0,0)
s2 = 1

M=10000
beta.save=matrix(NA,M,length(beta))
s2.save=rep(NA,M)
## for loop, with 100 iterations
for(iter in 1:10000){
    ## draw beta from full ocnditional
    ## this needs to calculated by hand
    ## start with our starting values and updates in each iteration
    b <- t(X)%*%y/s2
    A <- diag(2)/100+t(X)%*%X/s2
    A.inv=solve(A)
    beta = as.numeric(
        rmvnorm(1, mean = A.inv%*%b, sigma = A.inv)
    )
    ## draw s2 from full conditional 
    s2 <- rinvgamma(
        1,
        shape=10+length(y)/2, 
        rate=10+1/2*(t(y-X%*%beta)%*%(y-X%*%beta))
    )

    ## save out parameters
    beta.save[iter, ] = beta
    s2.save[iter] = s2
}

## see what happens here
plot(s2.save, type = "l")
matplot(cbind(s2.save,beta.save),type="l")

## look at the marginal posteriors with the histograms
hist(s2.save)
hist(beta.save[, 2])

## description
## lat = latitude
## lon = longitude
## maxdepth = maximum depth of the lake
## res.time = average time it takes a particle to flow through a lake (in minutes)
## secchi = "Secchi depth" = maximum depth at which a "Secchi disk" is visible
##          See https://en.wikipedia.org/wiki/Secchi_disk
## tn = total Nitrogen in lake
## tp = total Phosphorous in lake
## state = State lake is in
## SurLand = Surrounding land classification (forest, agriculture, or other)

load("Homework3/lagos.Rdata")
library(maptools)
library(maps)
xy=lagos[,2:1]
plot(xy,pch=".",cex=2,col="blue")
map("state",add=T,lwd=2)


str(lagos)
summary(lagos)
pairs(lagos)
head(lagos)

##
## Create Indicator Vars
##

lagos$ag=0
lagos$ag[lagos$SurLand=="agricultural"] <- 1
lagos$forest=0
lagos$forest[lagos$SurLand=="forest"] <- 1


###############################################
##
## Semiparametric model fit
## (should give you some idea what to expect 
##  from a Bayesian analysis)
##    N=tn, P=tp, S=secchi, A=1_{SurLand=ag}, F=1_{SurLand=forest}
##    log(N_i) ~ b0+f(log(P_i)+g(log(S_i)+b1*A_i+b2*F_i+e_i,  e_i ~ N(0,s2) 
##
###############################################

library(mgcv)
fit=gam(log(tn)~s(log(tp))+s(log(secchi))+ag+forest,data=lagos)
summary(fit)
plot(fit)


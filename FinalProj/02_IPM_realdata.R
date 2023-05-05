##########################################################################
#
# Dead-recovery model for turkey banding data analysis of male gobblers,
#                           Pennsylvania
#
# Survival and recovery rates are modeled by age. No random effects by WMU so assumes
# Harvest rates are the same across WMUs.
#
# September 2022
#
# VAW modification and comments 4.7.2023
# Added Statewide abundance and abundance per WMU (how to  make this bayesian)

# Need: recruitment, season, sex differences (female latent state)
#########################################

# Clean env
rm(list = ls())
gc()

#Sys.setenv(JAGS_HOME="C:/Program Files (x86)/JAGS/JAGS-4.3.1")

# Load in libraries
library(jagsUI)
library(MCMCvis)
library(RODBC)
library(stringr)
library(dplyr)

# read in data
df <- read.csv("../Proposal/FromDRD/BRM/Data/Males Band Data Export Report 2000_2021_2022_8-23-22_Final.csv", header=T)


# Define function to create a matrix with information about known latent state z
known.state.mr <- function(mr){
  state <- matrix(NA, nrow = dim(mr)[1], ncol = dim(mr)[2])
  rec <- which(rowSums(mr)==2)
  for (i in 1:length(rec)){
    n1 <- min(which(mr[rec[i],]==1))
    n2 <- max(which(mr[rec[i],]==1))
    state[rec[i],n1:n2] <- 1
    state[rec[i],n1] <- NA
    state[rec[i],n2:dim(mr)[2]] <- 0
  }
  return(state)
}

# Define function to create a matrix of initial values for latent state z
mr.init.z <- function(mr){
  ch <- matrix(NA, nrow = dim(mr)[1], ncol = dim(mr)[2])
  rec <- which(rowSums(mr)==1)
  for (i in 1:length(rec)){
    n1 <- which(mr[rec[i],]==1)
    ch[rec[i],n1:dim(mr)[2]] <- 0
    ch[rec[i],n1] <- NA
  }
  return(ch)
}

########################################################################
##Read in tables from Access database for 2020 and 2021

# Note: Duane pulled directly from db but I was having RODBC issues
# Made those tables into csv's and loaded n here
db <- "C:/Users/vaw5154/OneDrive - The Pennsylvania State University/PhD/PSUTurkey/Proposal/FromDRD/BRM/"
band20 <- read.csv(paste0(db,"2020 Gobbler_Banding_Data.csv"))
band21 <- read.csv(paste0(db,"2021 Gobbler_Banding_Data.csv"))
recov20 <- read.csv(paste0(db,"2020 Band_Report_Information.csv"))
recov21 <- read.csv(paste0(db,"2021 Band_Report_Information.csv"))

### Retain needed info, fix errors, rename variables, and combine years
## recov = recovery dataframe
## band  = banding dataframe
  
# Band: 
band20 <- band20[,c(3,7,9,10,18,19)]
band21 <- band21[,c(3,7,9,10,18,19)]
band <- rbind(band20,band21)
names(band)[1] <- "WMU"
names(band)[2] <- "cdate"
names(band)[3] <- "lat"
names(band)[4] <- "long"
names(band)[5] <- "band"
names(band)[6] <- "age"
band$band <- as.numeric(band$band)
band$band <- ifelse(is.na(band$band),532,band$band)
band$age <- ifelse(band$age=="j","J",band$age)
band <- band[band$age=="J"|band$age=="A",]
band$cdate <- as.Date(band$cdate,"%m/%d/%Y") # added for consistency here

# Recovery
recov20 <- recov20[,c(1,2,14)]
recov21 <- recov21[,c(1,2,15)]
recov21$`Date of Harvest/Recovery` <- as.POSIXct(recov21$`Date of Harvest/Recovery`,format="%m/%d/%Y")
recov <- rbind(recov20,recov21)
names(recov)[1] <- "band"
names(recov)[2] <- "cause"
names(recov)[3] <- "rdate"
recov$cause <- ifelse(recov$cause=="harvest" | recov$cause=="harvested" | recov$cause=="Harvested" |
recov$cause=="harvetsed","H",recov$cause)
recov <- recov[recov$cause=="H",]
recov$rdate <- as.Date(recov$rdate,"%m/%d/%Y") # added for consistency
     
####################################################################################
## Get 2022 data imported and formatted
## Recovery data
recov22 <- read.csv("../Proposal/FromDRD/BRM/Data/2022 Winter_Spr_Band_Reports_To_PSU_08102022.csv")
recov22 <- recov22[recov22$sex=="M" & recov22$encounter_reason=="Harvest",]
names(recov22)[1] <- "band"
names(recov22)[20] <- "cause"
recov22$cause <- ifelse(recov22$cause=="Harvest","H","O")
recov22$rdate <- as.POSIXct(paste0(recov22$year,"-",recov22$month,"-",recov22$day))
#hist(recov22$killdate,"weeks")
recov22 <- recov22[,c(1,20,25)]

## Add 2022 data to 2020/21 data
recov <- rbind(recov,recov22)
## Banding data
band22 <- read.csv("../Proposal/FromDRD/BRM/Data/TurkeyBandingReport_3yrs_with Cty_Twp_20220725.csv")
band22 <- band22[,c(3,2,8,9,1,10)]
names(band22)[1] <- "WMU"
names(band22)[2] <- "cdate"
names(band22)[3] <- "lat"
names(band22)[4] <- "long"
names(band22)[5] <- "band"
names(band22)[6] <- "age"
band22$cdate <- as.Date(band22$cdate,"%m/%d/%Y")
band22$age <- ifelse(band22$age=="Adult","A",band22$age)
band22$age <- ifelse(band22$age=="Juvenile","J",band22$age)

## Add 2022 data to 2020/21 data
band <- rbind(band,band22)

###Tabulate releases by Year and age
release <- as.data.frame(table(format(band$cdate, "%Y"),band$age)) 
release

### Merge banding and recovery data
full <- left_join(band, recov, by="band") #ignores recoveries not retained in release dataframe
full$WMU <- as.factor(full$WMU)
#full$rdate <- as.Date(band22$rdate,"%m/%d/%Y")
# test <- as.numeric(year(full$rdate)) - 2019 # testing something
full$rYr <- as.numeric(format(full$rdate, "%Y"))-2019 #recovery year 99=not recovered
full$rYr <- ifelse(is.na(full$rYr), 99, full$rYr)
full$cYr <- as.numeric(format(full$cdate, "%Y"))-2019 #release year

### Tabulate recoveries by year of release/recovery and age
full %>% 
  group_by(rYr,cYr,age) %>%
  summarise(n=n()) %>%
  mutate(total=sum(n))
  
###Create mark-recapture matrix MR
n.occasions <- max(full$cYr) # number release occasions
n.release <- dim(full)[1]
MR <- matrix(NA, ncol=n.occasions+1, nrow=n.release)
# vector of occasion of marking
mark.occ <- full$cYr
recov.occ <- full$rYr

# Fill the CH matrix
for (i in 1:n.release) {
  MR[i,mark.occ[i]] <- 1 # write 1 at release occasion
  for (t in mark.occ[i]:n.occasions){
    if (t==recov.occ[i]) {MR[i,t+1] <- 1} # write 1 at recovery occasion
  } #t loop
} #i loop

MR[which(is.na(MR))] <- 0
     
### 8.2.2. Adapted model with constant parameters to allow for age-specific and time-specific rates
## Set up -----
# i. Create vector with occasion of marking 
get.first <- function(x) min(which(x!=0))
f <- apply(MR, 1, get.first)
       
# ii. Create indicator array of which occasion an animal is juvenile
I <- array(0,c(dim(MR)[1],n.occasions))
for (i in 1:dim(MR)[1]){ 
  if(full$age[i]=="J"){I[i,f[i]] <- 1}
}

# iii. Create indicator vector for which beta is the reference level
J <- c(0,rep(1,n.occasions-1))   

# iv. Create indicator vector of WMU
wmu <- as.numeric(full$WMU) #converting factor to numeric
    
## Specify model in BUGS language
sink("mr.ss.jags")
cat("
  model {
  
  # Priors and constraints
  for (i in 1:nind) { # nid = dimensions of mark recapture matrix
   for (t in f[i]:(n.occasions-1)) {
   
    # s[1,t] = survival calculations at each timestep
    #logit(s[i,t]) <- mu + alpha*I[i,t] + beta[t]*J[t] + gamma[wmu[i]]
    s[i,t] <- 1/ (1+exp(mu + alpha*I[i,t] + beta[t]*J[t] + gamma[wmu[i]]))
   
      # r[i,t] = recapture calculation
      r[i,t] <- repp[1]*I[i,t] + repp[2]*(1-I[i,t])
      } #t
    } #i

 #######  priors for time effect
    for (t in 1:(n.occasions-1)) {
      beta[t] ~ dnorm(0, 0.001)
    }

 ####### priors for WMU effect
    for (u in 1:n.wmu) {
      gamma[u] ~ dnorm(0, tau)
    }

    mu ~ dnorm(0, 0.001)                # Prior for logit of reference level survival
      sigma ~ dunif(0, 10)                # Prior for sd of logit of survival variability
      tau <- pow(sigma, -2) 

    repp[1] ~ dunif(0,1)                # Prior for juvenile reporting rate
    repp[2] ~ dunif(0,1)                # Prior for adult reporting rate
    alpha ~ dnorm(0, 0.001)             # Prior for juvenile effect on survival
    #alpha[2] ~ dnorm(0, 0.001)         # Dummy Prior for juvenile effect on survival
    
    
  ############# Derived STATEWIDE survival rates
    for (t in 1:(n.occasions-1)) {
      s.ad[t] <- 1 / (1+exp(-mu-beta[t]*J[t]))  # Back-transformed adult survival
      s.jv[t] <- 1 / (1+exp(-mu-alpha-beta[t]*J[t]))  # Back-transformed juvenile survival
    }#t
  ###### Derived MEAN STATEWIDE survival rates
      mean.s.ad <- 1 / (1+exp(-mu))  # Back-transformed adult survival
      mean.s.jv <- 1 / (1+exp(-mu-alpha))  # Back-transformed juvenile survival

  ############# Derived STATEWIDE harvest rate estimates
    for (t in 1:(n.occasions-1)) {
        h.juv[t] <- (1-s.jv[t])*repp[1]
        h.ad[t] <- (1-s.ad[t])*repp[2]
    }#t
  ########### Derived MEAN STATEWIDE harvest rate estimates
        mean.h.jv <- (1-mean.s.jv)*repp[1]
        mean.h.ad <- (1-mean.s.ad)*repp[2]


  ## WMU survival rates
    for (t in 1:(n.occasions-1)) { # over time
     for (u in 1:n.wmu) { # over WMU
      s.ad.wmu[t,u] <- 1 / (1+exp(-mu-beta[t]*J[t] - gamma[u]))  # Back-transformed adult survival
      s.juv.wmu[t,u] <- 1 / (1+exp(-mu-alpha-beta[t]*J[t]-gamma[u]))  # Back-transformed juvenile survival
     }#u
    }#t

  ### Informative priors for hunter reporting rates based on Diefenbach et al. (2012)
  rrate.j ~ dnorm(0.709240358,1/(0.071909667)^2) #  non-reward bands for juveniles
  rrate.a ~ dnorm(0.8693272,1/(0.038757301)^2) #  non-reward bands for adults
   

  ## WMU harvest rates
   for (t in 1:(n.occasions-1)) {
    for (g in 1:n.wmu) {
        h.juv.wmu[t,g] <- (1-s.juv.wmu[t,g])*repp[1]/rrate.j 
        h.ad.wmu[t,g] <- (1-s.ad.wmu[t,g])*repp[2]/rrate.a
     }#g
    }#t


# Lincoln estimator -----
## Statewide abundance
    for (t in 1:(n.occasions-1)) { # over time
      abun.juv[t] <-  (repp[1]/rrate.j)/h.juv[t]       # harvest/harvest rate
      abun.ad[t] <- (repp[2]/rrate.a)/h.ad[t]
    }#t

## WMU abundance 
    for (t in 1:(n.occasions-1)) { # over time
     for (u in 1:n.wmu) { # over WMU
      abun.juv.wmu[t,u] <-  (repp[1]/rrate.j)/h.juv.wmu[t,u]   # harvest/harvest rate
      abun.ad.wmu[t,u] <- (repp[2]/rrate.a)/h.ad.wmu[t,u]
    
     }#u
    }#t
    
# Recruitment (pg 354) -----
    
    # Pt ~ Pois(Rtft)
    # P = rate 
    # R = number of surveyed broods,
    # hp = Hens per poult
    # likelihood = Lp(P, R| hp)
    
# Latent state of F based on M
# will 'guess' at survival and abundance:
# determine number of females that can be produced per female?

    #######  priors for time effect
    for (t in 1:(n.occasions-1)) {
      hp[t] <-  mean.fec
    }#t

    mean.fec ~ dunif(0, 20)
    
    #######  Liklihood for productivity daya
    for (t in 1:(n.occasions-1)) {
      P[t] ~ dpois(rho[t])
      rho[t] <- hp[t]
    }
    
       # #######  priors for time effect
    for (t in 1:(n.occasions-1)) {
    for (u in 1:n.wmu) { # over WMU
      hp.wmu[t, u] <-  mean.fec
    }#u
    }#t

   # mean.fec ~ dunif(0, 20)

# ## WMU recruitment
    for (t in 1:(n.occasions-1)) { # over time
     for (u in 1:n.wmu) { # over WMU
      P.wmu[t, u] ~ dpois(rho.wmu[t, u])
      rho.wmu[t, u] <- hp.wmu[t, u]
     }#u
    }#t


# Likelihood 
for (i in 1:nind){
  
  # Define latent state at first capture
  z[i,f[i]] <- 1

   for (t in (f[i]+1):n.occasions){
     # State process - survival
    z[i,t] ~ dbern(mu1[i,t])
    mu1[i,t] <- s[i,t-1] * z[i,t-1]

   # Observation process
  
    y[i,t] ~ dbern(mu2[i,t])
    mu2[i,t] <- r[i,t-1] * (z[i,t-1] - z[i,t])
  } #t
} #i
}
    ", fill = TRUE)

sink()
     
# Random sample (for now) of P and R     
#R <- sample(x = 1:15, size  = 3)

# Bundle data
known.state <- known.state.mr(MR)
is.na(known.state) <- NA
jags.data <- list(y=MR, f=f, I=I, J=J, nind = dim(MR)[1], n.occasions = (dim(MR)[2]), 
                  z = known.state, n.wmu = length(unique(full$WMU)), wmu = wmu)#, hp=hp)

# Initial values
inits <- function(){list(z = mr.init.z(MR),
                         mean.fec = runif(1, 0, 1),
                         repp = runif(2, 0, 1), alpha=rnorm(1,0,0.001), 
                         beta=rnorm(n.occasions,0,0.001), 
                         mu = rnorm(1,0,0.001), sigma=runif(1,0,10), 
                         rrate.a=dnorm(1,.71,1/.072),
                         rrate.j=dnorm(1,0.87,1/0.039))}

######## Call JAGS from R (BRT 4 min)
# Parameters monitored
parameters <- c("mu","alpha","beta","gamma","rrate.a","rrate.j","sigma",
                "repp","h.ad","h.juv","h.ad.wmu","h.juv.wmu","s.jv","s.ad",
                "mean.s.ad","mean.s.jv","mean.h.ad","mean.h.jv", "abun.ad.wmu", 
                "abun.juv.wmu", "P", "P.wmu", "s.ad.wmu", "s.juv.wmu",
                "abun.ad", "abun.juv")
     
# MCMC settings
ni <- 25000
nt <- 1
nb <- 7500
nc <- 6
     
######## Call JAGS from R (BRT 4 min)
mr.ss.re.ab <- jags(jags.data, inits, parameters, "mr.ss.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                parallel=TRUE)

# See output
summary(mr.ss.re.ab)
print(mr.ss.re.ab, digits = 3)

MCMCtrace(mr.ss.re.ab, pdf=T)
     
# Print summary of results
cat("Mean of Juvenile harvest (State):", mean(mr.ss.re.ab$sims.list$h.juv), "\n")
cat("Mean of Juvenile harvest (WMU):", mean(mr.ss.re.ab$sims.list$h.juv.wmu), "\n")
cat("Mean of Adult harvest (State):", mean(mr.ss.re.ab$sims.list$h.ad), "\n")
cat("Mean of Adult harvest (WMU):", mean(mr.ss.re.ab$sims.list$h.ad.wmu), "\n")
cat("Mean of Juvenile survival (State):", mean(mr.ss.re.ab$sims.list$s.jv), "\n")
cat("Mean of Juvenile survival (WMU):", mean(mr.ss.re.ab$sims.list$s.juv.wmu), "\n")
cat("Mean of Adult survival (State):", mean(mr.ss.re.ab$sims.list$s.ad), "\n")
cat("Mean of Adult survival (WMU):", mean(mr.ss.re.ab$sims.list$s.ad.wmu), "\n")
cat("Mean of Juvenile abundance (State):", mean(mr.ss.re.ab$sims.list$abun.juv), "\n")
cat("Mean of Juvenile abundance (WMU):", mean(mr.ss.re.ab$sims.list$abun.juv.wmu), "\n")
cat("Mean of Adult abundance (State):", mean(mr.ss.re.ab$sims.list$abun.ad), "\n")
cat("Mean of Adult abundance (WMU):", mean(mr.ss.re.ab$sims.list$abun.ad.wmu), "\n")
cat("Mean of Reruitment (State):", mean(mr.ss.re.ab$sims.list$P), "\n")
cat("Mean of Reruitment (WMU):", mean(mr.ss.re.ab$sims.list$P.wmu), "\n")
cat("95% credible interval Juvenile harvest (State):", quantile(mr.ss.re.ab$sims.list$h.juv, c(0.025, 0.975)), "\n")
cat("95% credible interval Juvenile harvest (WMU):", quantile(mr.ss.re.ab$sims.list$h.juv.wmu, c(0.025, 0.975)), "\n")
cat("95% credible interval for Adult harvest (State):", quantile(mr.ss.re.ab$sims.list$h.ad, c(0.025, 0.975)), "\n")
cat("95% credible interval for Adult harvest (WMU):", quantile(mr.ss.re.ab$sims.list$h.ad.wmu, c(0.025, 0.975)), "\n")
cat("95% credible interval for Juvenile survial (State):", quantile(mr.ss.re.ab$sims.list$s.jv, c(0.025, 0.975)), "\n")
cat("95% credible interval for Juvenile survival (WMU):", quantile(mr.ss.re.ab$sims.list$s.juv.wmu, c(0.025, 0.975)), "\n")
cat("95% credible interval for Adult survival (State):", quantile(mr.ss.re.ab$sims.list$s.ad, c(0.025, 0.975)), "\n")
cat("95% credible interval for Adult survival (WMU):", quantile(mr.ss.re.ab$sims.list$s.ad.wmu, c(0.025, 0.975)), "\n")
cat("95% credible interval for Juvenile abundance (State):", quantile(mr.ss.re.ab$sims.list$abun.juv, c(0.025, 0.975)), "\n")
cat("95% credible interval for Juvenile abundance (WMU):", quantile(mr.ss.re.ab$sims.list$abun.juv.wmu, c(0.025, 0.975)), "\n")
cat("95% credible interval for Adult abundance (State):", quantile(mr.ss.re.ab$sims.list$abun.ad, c(0.025, 0.975)), "\n")
cat("95% credible interval for Adult abundance (WMU):", quantile(mr.ss.re.ab$sims.list$abun.ad.wmu, c(0.025, 0.975)), "\n")
cat("95% credible interval for Recruitment (State):", quantile(mr.ss.re.ab$sims.list$P, c(0.025, 0.975)), "\n")
cat("95% credible interval for Recruitment (WMU):", quantile(mr.ss.re.ab$sims.list$P.wmu, c(0.025, 0.975)), "\n")
  
# Get summary to check
MCMCsummary(mr.ss.re.ab, params = 's.juv.wmu', round = 2)

# test plotting
MCMCplot(mr.ss.re.ab, params = 's.ad.wmu', ci = c(50, 95), 
         horiz  = T, 
       # labels = wmu,
         
         main = "Survival of Adults per WMU", col = "lightpink",
         col2 = "red", ref_ovl = T)

wmu <- unique(wmu)


# Save outputs
saveRDS(mr.ss.re.ab, "../IPM/Data/20230503_ipm.rds")
saveRDS(mr.ss.re.ab, "Data/20230503_ipm.rds")


# Done!



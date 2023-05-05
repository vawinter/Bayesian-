library(rjags)
library(coda)

# Simulated data
set.seed(123)
N_obs <- c(500, 350, 200, 125, 75)  # observed population sizes
r <- 0.5  # recruitment rate

# IPM model in JAGS
  model_string <- "
model {
  # Priors
  N[1] ~ dunif(0, 1000)  # uniform prior for initial population size
  p ~ dbeta(1, 1)  # non-informative beta prior for survival probability

  # Transition model
  for (t in 1:4) {
    N[t+1] ~ dpois(p*N[t])  # poisson model for population dynamics
  }

  # Likelihood
  for (t in 2:6) {
    N_obs[t-1] ~ dpois(N[t-1]*r)  # Poisson likelihood for observed population sizes
  }
}
"

# Data list
data_list <- list(N=N_obs, r=r)

# Create the model object
model <- jags.model(textConnection(model_string), data=data_list)

# Specify the parameters to monitor
params <- c("N", "p", "r")

# Set the number of iterations and burn-in
n_iter <- 1000
n_burnin <- 500

# Run the model
samples <- jags.samples(model, params, n.iter=n_iter, n.burnin=n_burnin)
# Burn-in and sampling
update(model, 1000)
samples <- coda.samples(model, variable.names=params, n.iter=5000)

# Summarize the posterior distributions
summary(samples)

# Calculate the posterior means and credible intervals
post_means <- apply(samples[[1]], 2, mean)
post_cis <- t(apply(samples[[1]], 2, function(x) quantile(x, c(0.025, 0.975))))

# Print the results
cat("Posterior mean of p:", post_means[6], "\n")
cat("95% credible interval fro p:", post_cis[6], "\n")

# Create trace plots
plot(samples)
plot(samples[, "p"], main = "Trace plot for p")

      
      
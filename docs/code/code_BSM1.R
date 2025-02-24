# -----------------------------
# We need the following libraries
# -----------------------------

library(mvtnorm)
library(rstan)
library(coda)
library(ggplot2)

############################################
### HAMILTONIAN MONTE CARLO - a step   #####
############################################

# How is it working? Take a look on a single step

HMC_step_norm <-  function(epsilon, L, theta_old, A){
  
  theta = theta_old
  p = rnorm(length(theta),0,1) # independent standard normal variates
  p_old = p
  
  # Make a half step for momentum at the beginning
  p=p - (epsilon / 2) * solve(A) %*% theta 
  
  # Alternate full steps for position and momentum
  for (i in 1:L){
    theta_plot <- theta
    theta=theta + epsilon * p
    if(i != L) p = p - epsilon * solve(A) %*% theta 
    
    points(theta[1], theta[2], pch = 20)
    segments(theta_plot[1], theta_plot[2], theta[1], theta[2])
  }
  
  # Make a half step for momentum at the end.
  p = p - (epsilon / 2) * solve(A) %*% theta
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = - p
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = t(theta_old) %*% solve(A) %*% theta_old / 2
  current_K = sum(p_old^2) / 2
  proposed_U = t(theta) %*% solve(A) %*% theta / 2
  proposed_K = sum(p^2) / 2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K)){
    return (theta) # accept
  }
  else{
    return (theta_old) # reject
  }
}

#-----------------------------------------------
# plot some iterations
theta_old <- c(-1.5,-1)
A <- matrix(c(1,0.3,0.3,1), ncol = 2)

eg <- expand.grid(seq(-4, 4, length.out = 50), seq(-4, 4, length.out = 50))
dns <- apply(eg, 1, function(x) dmvnorm(x = x, sigma = A))
plot(theta_old[1], theta_old[2], xlim = c(-3,3), ylim = c(-3,3), pch = 20, xlab = "X1", ylab = "X2")
contour(seq(-4, 4, length.out = 50), seq(-4, 4, length.out = 50), matrix(dns, ncol = 50), add = T, col = 2)

set.seed(42)
for(i in 1:100){
  theta_old <- HMC_step_norm(0.2, 10, theta_old, A)
  points(theta_old[1], theta_old[2], col = 2, pch = 20)
  Sys.sleep(2)
}

############################################
### HAMILTONIAN MONTE CARLO            #####
############################################

HMC_norm <-  function(niter, burnin, theta0, mu, Sigma, epsilon, L){
  
  results <- matrix(ncol = 2, nrow = niter - burnin)
  theta_old <- theta0
  
  pb <- txtProgressBar(0, niter, style = 3)
  for(i in 1:niter){
    
    theta = theta_old
    p = rnorm(length(theta),0,1) 
    p_old = p
    
    p = p - (epsilon / 2) * solve(Sigma) %*% (theta - mu)
    for (j in 1:L){
      theta = theta + epsilon * p
      if(j != L) p = p - epsilon * solve(Sigma) %*% (theta - mu)
    }
    p = p - (epsilon / 2) * solve(Sigma) %*% (theta - mu)
    p = - p
    
    current_U = t(theta_old - mu) %*% solve(Sigma) %*% (theta_old - mu) / 2
    current_K = sum(p_old^2) / 2
    proposed_U = t(theta - mu) %*% solve(Sigma) %*% (theta - mu) / 2
    proposed_K = sum(p^2) / 2
    
    if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K)){
      theta_old <- theta
    }
    
    if(i > burnin){
      results[i - burnin, ] <- c(theta_old)
    }
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(results)
}

#-----------------------------------------------------------------
# RUN: HMC vs RW-MH

theta0 <- c(3,3)
mu <- c(2,2)
Sigma <- matrix(c(1,0.3,0.3,1), ncol = 2)
niter <- 10000
burnin <- 1000

set.seed(123)
HMC_sample <- HMC_norm(niter = niter, burnin = burnin, theta0 = theta0, 
                       mu = mu, Sigma = Sigma, epsilon = 0.2, L = 10)
plot(HMC_sample)

############################################
### STAN EXAMPLES                      #####
############################################

#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------

set.seed(42)
y <- sample(c(0,1), size = 50, replace = T)

data_STAN_BB <-list(n = 50, 
                    y = y, 
                    a0 = 1, 
                    b0 = 1)

out_STAN_BB <- stan(file = "bernoulli_beta.stan", 
                    data = data_STAN_BB,
                    chains = 1, 
                    iter = 1000, 
                    warmup = 200, 
                    seed = 42)

rstan::traceplot(out_STAN_BB, pars = c("theta"))
params_STAN_BB <- As.mcmc.list(out_STAN_BB)
summary(params_STAN_BB)
geweke.diag(params_STAN_BB)

#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------

set.seed(42)
y <- rnorm(50, 5, 2.5)

data_STAN_NNG <-list(n = 50, 
                     y = y, 
                     m0 = 0, 
                     s20 = 100,
                     a0 = 0.1, 
                     b0 = 0.1)

out_STAN_NNG <- stan(file = "normal_normal_gamma.stan", 
                     data = data_STAN_NNG,
                     chains = 2, 
                     iter = 1000, 
                     warmup = 200, 
                     seed = 42)

rstan::traceplot(out_STAN_NNG, pars = c("mu", "sigma2"))
params_STAN_NNG <- As.mcmc.list(out_STAN_NNG)
summary(params_STAN_NNG)

geweke.diag(params_STAN_NNG)
gelman.diag(params_STAN_NNG, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)
gelman.plot(params_STAN_NNG, confidence = 0.95)

#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------

set.seed(42)
y <- rpois(50, 5)

data_STAN_POI <-list(n = 50, 
                     y = y,
                     a0 = 0.1, 
                     b0 = 0.1)

out_STAN_POI <- stan(file = "poi.stan", 
                     data = data_STAN_POI,
                     chains = 2, 
                     iter = 1000, 
                     warmup = 200, 
                     seed = 42)

rstan::traceplot(out_STAN_POI, pars = c("lambda"))
params_STAN_POI <- As.mcmc.list(out_STAN_POI)
summary(params_STAN_POI)

geweke.diag(params_STAN_POI)
gelman.diag(params_STAN_POI, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)
gelman.plot(params_STAN_POI, confidence = 0.95)



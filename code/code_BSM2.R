# -----------------------------
# We need the following libraries
# -----------------------------

library(mvtnorm)
library(rstan)
library(coda)
library(ggplot2)

# -----------------------------
# The data reported in the slides and the corresponding design matrix
# -----------------------------

y <- c(-3.7, -5.0, 7.3, -4.3, -4.4, -7.3, 0.0, -5.7, 9.2, -4.0,
       -1.1, 1.5, -3.8, -7.4, 5.3, -8.8, -2.2, -5.4, 0.2, -2.6)
z <- c(1.1, 2.2, 3.6, 0.9, 1.9, 2.1, 2.7, 1.8, 4.0, 1.9,
       2.4, 3.0, 1.6, 1.0, 3.8, -0.3, 2.9, 2.0, 3.0, 2.4)
X <- cbind(rep(1, 20), z, z^2)

# -----------------------------
# In this first part, we build a Gibbs sampler to perform
# inference with linear model and the non-informative prior specification
# 
# The function is called with four arguments: 
#    - the response variable "y"
#    - the design matrix "X"
#    - total number of iterations "niter"
#    - number of burn-in iterations "nburn"
# -----------------------------

gibbs_vague_prior <- function(y, X, niter, nburn){
  
  out_sigma2 <- rep(0, niter)
  out_beta <- matrix(0, ncol = ncol(X), nrow = niter)
  
  beta_ML <- solve(t(X) %*% X) %*% t(X) %*% y
  
  for(i in 1:niter){
    an <- (nrow(X) - ncol(X)) / 2
    bn <- 0.5 * t(y - X %*% beta_ML) %*% (y - X %*% beta_ML)
    Sigman <- solve(t(X) %*% X)
      
    out_sigma2[i] <- 1 / rgamma(1, shape = an, scale = 1 / bn)
    out_beta[i,] <- rmvnorm(n = 1, mean = beta_ML, sigma = out_sigma2[i] * Sigman)
  }
  
  list(out_sigma2[-c(1:nburn)], out_beta[-c(1:nburn),])
}

# with 5 observations

out_gibbs_vague_5 <- gibbs_vague_prior(y[1:5], X[1:5,], 1200, 200)

acf(out_gibbs_vague_5[[1]])
acf(out_gibbs_vague_5[[2]][,1])
acf(out_gibbs_vague_5[[2]][,2])
acf(out_gibbs_vague_5[[2]][,3])

hist(out_gibbs_vague_5[[1]])
plot(as.data.frame(out_gibbs_vague_5[[2]]))

summary(out_gibbs_vague_5[[1]])
summary(as.data.frame(out_gibbs_vague_5[[2]]))

beta_sigma_out_5 <- as.mcmc(cbind(out_gibbs_vague_5[[2]], out_gibbs_vague_5[[1]]))
summary(beta_sigma_out_5)
effectiveSize(beta_sigma_out_5)
geweke.diag(beta_sigma_out_5)

# with the whole sample

out_gibbs_vague <- gibbs_vague_prior(y, X, 1200, 200)

acf(out_gibbs_vague[[1]])
acf(out_gibbs_vague[[2]][,1])
acf(out_gibbs_vague[[2]][,2])
acf(out_gibbs_vague[[2]][,3])

hist(out_gibbs_vague[[1]])
plot(as.data.frame(out_gibbs_vague[[2]]))

summary(out_gibbs_vague[[1]])
summary(as.data.frame(out_gibbs_vague[[2]]))

beta_sigma_out <- as.mcmc(cbind(out_gibbs_vague[[2]], out_gibbs_vague[[1]]))
summary(beta_sigma_out)
effectiveSize(beta_sigma_out)
geweke.diag(beta_sigma_out)

# -----------------------------
# In the following, we perform similar analysis with
# the model implemented in STAN and the informative prior
# 
# The STAN model is coded in the file "regression_informative_prior.stan"
# -----------------------------

# with 5 observations

data_STAN_5 <-list(n = 5, 
                   p = 3,  
                   y = y[1:5], 
                   X = X[1:5,], 
                   beta0 = rep(0, 3), 
                   Sigma0 = diag(10^2, 3), 
                   a0 = 2, 
                   b0 = 1)

out_STAN_5 <- stan(file = "regression_informative_prior.stan", 
                   data = data_STAN_5,
                   chains = 1, 
                   iter = 1000, 
                   warmup = 200, 
                   seed = 42)

rstan::traceplot(out_STAN_5, pars = c("beta", "sigma2"))
params_STAN_5 <- As.mcmc.list(out_STAN_5, pars = c("beta", "sigma2"))
summary(params_STAN_5)
effectiveSize(params_STAN_5)
geweke.diag(params_STAN_5)

df_plot <- data.frame(x = as.vector(params_STAN_5[[1]]), 
                      group = factor(rep(c("beta1", "beta2", "beta3", "sigma2"), each = 800)))
ggplot(df_plot[1:2400,]) + 
  geom_histogram(aes(x = x, fill = group), col = 1, alpha = 0.5) + 
  facet_wrap(.~group, nrow = 1) +
  theme_bw() + 
  theme(legend.position = "null")

# with 20 observations

data_STAN_20 <-list(n = 20, 
                    p = 3,  
                    y = y, 
                    X = X, 
                    beta0 = rep(0, 3), 
                    Sigma0 = diag(10^2, 3), 
                    a0 = 2, 
                    b0 = 1)

out_STAN_20 <- stan(file = "regression_informative_prior.stan", 
                    data = data_STAN_20,
                    chains = 1, 
                    iter = 1000, 
                    warmup = 200, 
                    seed = 42)

rstan::traceplot(out_STAN_20, pars = c("beta", "sigma2"))
params_STAN_20 <- As.mcmc.list(out_STAN_20, pars = c("beta", "sigma2"))
summary(params_STAN_20)
effectiveSize(params_STAN_20)
geweke.diag(params_STAN_20)

df_plot <- data.frame(x = as.vector(params_STAN_20[[1]]), 
                      group = factor(rep(c("beta1", "beta2", "beta3", "sigma2"), each = 800)))
ggplot(df_plot[1:2400,]) + 
  geom_histogram(aes(x = x, fill = group), col = 1, alpha = 0.5) + 
  facet_wrap(.~group, nrow = 1) +
  theme_bw() + 
  theme(legend.position = "null")


# -----------------------------
# Hypothesis test
# 
# we define the marginal distribution of the data, 
# then we use such a function to calculate the
# Bayes factors
# -----------------------------

marginal_distribution <- function(y, X, beta0, Sigma0, a0, b0){

  Sigman <- solve(solve(Sigma0) + t(X) %*% X)
  betan <- Sigman %*% (solve(Sigma0) %*% beta0 + t(X) %*% y)
  an <- a0 + length(y) / 2
  bn <- b0 + 0.5 * (t(y) %*% y - t(betan) %*% solve(Sigman) %*% betan + 
                      t(beta0) %*% solve(Sigma0) %*% beta0)
  
  logM <- 0.5 * log(det(Sigman)) - an * log(bn) + lgamma(an)
  return(exp(logM[1,1]))
}

# Test for a single coefficient

marginal_distribution(y[1:5], X[1:5,-2], rep(0, 2), diag(10^2, 2), 3, 2) /
  marginal_distribution(y[1:5], X[1:5,], rep(0, 3), diag(10^2, 3), 3, 2)

marginal_distribution(y, X[,-2], rep(0, 2), diag(10^2, 2), 3, 2) /
  marginal_distribution(y, X[,], rep(0, 3), diag(10^2, 3), 3, 2)

# Test vs null model

marginal_distribution(y[1:5], X[1:5,-c(2,3)], rep(0, 1), diag(10^2, 1), 3, 2) /
  marginal_distribution(y[1:5], X[1:5,], rep(0, 3), diag(10^2, 3), 3, 2)

marginal_distribution(y, X[,-c(2,3)], rep(0, 1), diag(10^2, 1), 3, 2) /
  marginal_distribution(y, X[,], rep(0, 3), diag(10^2, 3), 3, 2)

# -----------------------------
# Predictive distribution
# -----------------------------


# -----------------------------
# Regularization
# 
# we compare the regression model of before (Gaussian prior) with the lasso
# regularization (Laplace prior)
# 
# note that by setting Sigma0 = diag(1 / lambda) we have the same smoothing
# of the slides with parameter lambda
# -----------------------------

data_STAN_reg_gaus <-list(n = 20, 
                          p = 3,  
                          y = y, 
                          X = X, 
                          beta0 = rep(0, 3), 
                          a0 = 2, 
                          b0 = 1)

out_STAN_reg_gaus <- stan(file = "regular_gaussian_prior.stan", 
                          data = data_STAN_reg_gaus,
                          chains = 1, 
                          iter = 1000, 
                          warmup = 200, 
                          seed = 42)

params_STAN_reg_gaus <- as.mcmc(extract(out_STAN_reg_gaus, 
                                        pars = c("beta", "sigma2"), permuted = TRUE))

data_STAN_reg_lap <-list(n = 20, 
                          p = 3,  
                          y = y, 
                          X = X, 
                          beta0 = rep(0, 3), 
                          a0 = 2, 
                          b0 = 1)

out_STAN_reg_lap <- stan(file = "regular_laplace_prior.stan", 
                          data = data_STAN_reg_lap,
                          chains = 1, 
                          iter = 1000, 
                          warmup = 200, 
                          seed = 42)

params_STAN_reg_lap <- as.mcmc(extract(out_STAN_reg_lap, 
                                       pars = c("beta", "sigma2"), permuted = TRUE))

# we generate here a dataframe with the samples of interest, divided by model
# and set model as a factor variable 

df_plot <- as.data.frame(cbind(rbind(params_STAN_reg_gaus[[1]], 
                 params_STAN_reg_lap[[1]]), 
                 c(params_STAN_reg_gaus[[2]], 
                   params_STAN_reg_lap[[2]]), 
                 rep(c(1,2), each = 800)))
colnames(df_plot) <- c("beta1", "beta2", "beta3", "sigma2", "model")
df_plot$model <- factor(df_plot$model, labels = c("gaus", "lap"))

# now some plots to compare the estimates

ggplot(df_plot) + 
  geom_histogram(aes(x = beta1, after_stat(density), group = model, fill = model), 
                 col = 1, alpha = 0.5, position="identity") + 
  theme_bw() + ylab("") + ggtitle("beta 1")
ggplot(df_plot) + 
  geom_histogram(aes(x = beta2, after_stat(density), group = model, fill = model), 
                 col = 1, alpha = 0.5, position="identity") + 
  theme_bw() + ylab("") + ggtitle("beta 2")
ggplot(df_plot) + 
  geom_histogram(aes(x = beta3, after_stat(density), group = model, fill = model), 
                 col = 1, alpha = 0.5, position="identity") + 
  theme_bw() + ylab("") + ggtitle("beta 3")
ggplot(df_plot) + 
  geom_histogram(aes(x = sigma2, after_stat(density), group = model, fill = model), 
                 col = 1, alpha = 0.5, position="identity") + 
  theme_bw() + ylab("") + ggtitle("sigma2")

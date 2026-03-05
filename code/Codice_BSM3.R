# -----------------------------
# We need the following libraries
# -----------------------------

library(mvtnorm)
library(truncnorm)
library(rstan)
library(coda)
library(BayesLogit)
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
  
  an <- (nrow(X) - ncol(X)) / 2
  bn <- 0.5 * t(y - X %*% beta_ML) %*% (y - X %*% beta_ML)
  Sigman <- solve(t(X) %*% X)
  
  for(i in 1:niter){
    out_sigma2[i] <- 1 / rgamma(1, shape = an, scale = 1 / bn)
    out_beta[i,] <- rmvnorm(n = 1, mean = beta_ML, 
                            sigma = out_sigma2[i] * Sigman)
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
                   iter = 1200, 
                   warmup = 200, 
                   seed = 42)

rstan::traceplot(out_STAN_5, pars = c("beta", "sigma2"))
params_STAN_5 <- As.mcmc.list(out_STAN_5, pars = c("beta", "sigma2"))
summary(params_STAN_5)
effectiveSize(params_STAN_5)
geweke.diag(params_STAN_5)

df_plot <- data.frame(x = as.vector(params_STAN_5[[1]]), 
                      group = factor(rep(c("beta1", "beta2", "beta3", "sigma2"), each = 1000)))
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
                    iter = 1200, 
                    warmup = 200, 
                    seed = 42)

rstan::traceplot(out_STAN_20, pars = c("beta", "sigma2"))
params_STAN_20 <- As.mcmc.list(out_STAN_20, pars = c("beta", "sigma2"))
summary(params_STAN_20)
effectiveSize(params_STAN_20)
geweke.diag(params_STAN_20)

df_plot <- data.frame(x = as.vector(params_STAN_20[[1]]), 
                      group = factor(rep(c("beta1", "beta2", "beta3", "sigma2"), each = 1000)))
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

zn1 <- -1.5
xn1 <- c(1, zn1, zn1^2)

sample_pred <- function(y, X, beta0, Sigma0, a0, b0){
  Sigman <- solve(solve(Sigma0) + t(X) %*% X)
  betan <- Sigman %*% (solve(Sigma0) %*% beta0 + t(X) %*% y)
  an <- a0 + length(y) / 2
  bn <- b0 + 0.5 * (t(y) %*% y - t(betan) %*% solve(Sigman) %*% betan + 
                      t(beta0) %*% solve(Sigma0) %*% beta0)
  
  sample <- rt(n = 1000, df = 2 * an) * 
    sqrt(as.numeric((1 + t(xn1) %*% Sigman %*% xn1) * bn / an)) + 
    as.numeric(xn1 %*% betan)
  
  return(sample)
}

sample11 <- sample_pred(y[1:5], X[1:5,], rep(0, 3), diag(25, 3), 2, 5)
sample12 <- sample_pred(y, X, rep(0, 3), diag(25, 3), 2, 5)

sample21 <- sample_pred(y[1:5], X[1:5,], rep(0, 3), diag(5, 3), 2, 5)
sample22 <- sample_pred(y, X, rep(0, 3), diag(5, 3), 2, 5)

data_plt <- data.frame(x = c(sample11, sample12, sample21, sample22), 
                       group = factor(rep(c("large var - small n", 
                                            "large var - large n", 
                                            "small var - small n", 
                                            "small var - large n"), each = 1000)))
ggplot(data_plt) + 
  geom_histogram(aes(x = x), alpha = 0.2, col = 1, fill = 1) + 
  theme_bw() + 
  facet_wrap(.~group, nrow = 2)

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
  theme_bw() + ylab("") + ggtitle("beta 3")2
ggplot(df_plot) + 
  geom_histogram(aes(x = sigma2, after_stat(density), group = model, fill = model), 
                 col = 1, alpha = 0.5, position="identity") + 
  theme_bw() + ylab("") + ggtitle("sigma2")

# -----------------------------
# PROBIT REGRESSION MODEL
# -----------------------------

# -----------------------------
# The data reported in the slides for the probit model
# and the corresponding design matrix
# -----------------------------

set.seed(123); betatrue <- c(-1, -2, 2)
z <- round(rnorm(100, 0, 1), digits = 1)
X <- cbind(rep(1, 100), z, z^2)
tempprobs <- pnorm(X %*% betatrue)
y <- sapply(tempprobs[,1], function(x) rbinom(1,1,x))

# -----------------------------
# In this first part, we build a Gibbs sampler to perform
# inference with the probit model, with the augmentation
# strategy discussed in the slides
# 
# The function is called with four arguments: 
#    - the response variable "y"
#    - the design matrix "X"
#    - total number of iterations "niter"
#    - number of burn-in iterations "nburn"
# -----------------------------

gibbs_augmented_probit <- function(y, X, niter, nburn, b0, Sigma0){
  
  out_beta <- matrix(0, ncol = ncol(X), nrow = niter)
  v <- rep(0, length(y))
  beta <- rep(0, ncol(X))
  
  # Sigman does not depend on the augmented variable
  Sigman <- solve(solve(Sigma0) + t(X) %*% X)
  
  for(i in 1:niter){
    
    # update the augmented variables 
    for(j in 1:length(y)){
      if(y[j] == 1){
        v[j] <- rtruncnorm(1, a = 0, b = Inf, mean = X[j,] %*% beta, sd = 1)  
      } else {
        v[j] <- rtruncnorm(1, a = -Inf, b = 0, mean = X[j,] %*% beta, sd = 1)  
      }
    }
    
    # sample the regression coefficients
    bn <- Sigman %*% (solve(Sigma0) %*% b0 + t(X) %*% v)
    beta <- rmvnorm(n = 1, mean = bn, sigma = Sigman)[1,]
    out_beta[i,] <- beta
  }
  
  return(out_beta[-c(1:nburn),])
}

# -----------------------------
# sample from the posterior distribution
# -----------------------------

out_augmented_probit <- 
  gibbs_augmented_probit(y, X, 2000, 1000, rep(0, 3), diag(10^3, 3))

params_STAN_probit <- as.mcmc(out_augmented_probit)
summary(params_STAN_probit)
geweke.diag(params_STAN_probit)

# plotting the densities
df_plot <- data.frame(label = factor(rep(c("beta1", "beta2", "beta3"), each = 1000)), 
                      x = as.vector(out_augmented_probit))

# -----------------------------
# now some plots to compare the estimates
# -----------------------------

ggplot(df_plot) + 
  geom_histogram(aes(x = x, after_stat(density), fill = label), 
                 col = 1, alpha = 0.5, position="identity") + 
  geom_vline(aes(xintercept = 0), lty = 2) + 
  theme_bw() + 
  ylab("") + 
  facet_wrap(~label) +
  theme(legend.position = "null")

ggplot(df_plot) + 
  geom_violin(aes(x = label, y = x, fill = label), trim=FALSE) +
  geom_boxplot(aes(x = label, y = x), width = 0.1) +
  geom_hline(aes(yintercept = 0), lty = 2) + 
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  theme(legend.position = "null")

# -----------------------------
# test if beta2 is significantly greater than 0
# ... the following is a Bayes factor, do you agree?
# -----------------------------

mean(out_augmented_probit[,2] > 0) / 
  mean(out_augmented_probit[,2] < 0)

# -----------------------------
# predicting the furture!!
# -----------------------------

xn1 <- c(1, 2, 4)
pn1 <- pnorm(out_augmented_probit %*% xn1)
yn1 <- sapply(pn1, function(x) sample(c(1,0), size = 1, prob = c(x, 1 - x)))
table(yn1)

# -----------------------------
# POISSON REGRESSION MODEL
# -----------------------------

# -----------------------------
# The data reported in the slides for the Poisson model
# and the corresponding design matrix
# -----------------------------

set.seed(123); betatrue <- c(-1, -2, 2)
z1 <- round(rnorm(100, 0, 1), digits = 1)
z2 <- round(rnorm(100, 0, 1), digits = 1)
X <- cbind(rep(1, 100), z1, z2)
tempprobs <- exp(X %*% betatrue)
y <- sapply(tempprobs[,1], function(x) rpois(1,x))

# -----------------------------
# sample from the posterior distribution
# -----------------------------

data_STAN_poi <-list(n = 100, 
                    p = 3,  
                    y = y, 
                    X = X, 
                    beta0 = rep(0, 3), 
                    Sigma0 = diag(10^3, 3))

out_STAN_poi <- stan(file = "poisson_glm.stan", 
                    data = data_STAN_poi,
                    chains = 1, 
                    iter = 2000, 
                    warmup = 1000, 
                    seed = 42)

rstan::traceplot(out_STAN_poi, pars = c("beta"))
params_STAN_poi <- As.mcmc.list(out_STAN_poi, pars = c("beta"))
summary(params_STAN_poi)
geweke.diag(params_STAN_poi)

# plotting the densities
df_plot <- data.frame(label = factor(rep(c("beta1", "beta2", "beta3"), each = 1000)), 
                      x = as.vector(params_STAN_poi[[1]]))

# -----------------------------
# now some plots to compare the estimates
# -----------------------------

ggplot(df_plot) + 
  geom_histogram(aes(x = x, after_stat(density), fill = label), 
                 col = 1, alpha = 0.5, position="identity") + 
  geom_vline(aes(xintercept = 0), lty = 2) + 
  theme_bw() + 
  ylab("") + 
  facet_wrap(~label) +
  theme(legend.position = "null")

ggplot(df_plot) + 
  geom_violin(aes(x = label, y = x, fill = label), trim=FALSE) +
  geom_boxplot(aes(x = label, y = x), width = 0.1) +
  geom_hline(aes(yintercept = 0), lty = 2) + 
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  theme(legend.position = "null")


# -----------------------------
# GAMMA REGRESSION MODEL
# -----------------------------

# -----------------------------
# The data reported in the slides for the Gamma model
# and the corresponding design matrix
# -----------------------------

set.seed(123); betatrue <- c(-1, -2, 2)
z1 <- round(rnorm(100, 0, 1), digits = 1)
z2 <- round(rnorm(100, 0, 1), digits = 1)
X <- cbind(rep(1, 100), z1, z2)
tempmeans <- exp(X %*% betatrue)
y <- sapply(tempmeans[,1], function(x) rgamma(1, shape = 1, rate = 1 / x))

# -----------------------------
# sample from the posterior distribution
# -----------------------------

data_STAN_gamma <-list(n = 100, 
                       p = 3,  
                       y = y, 
                       X = X, 
                       alpha = 1,
                       beta0 = rep(0, 3), 
                       Sigma0 = diag(10^3, 3))

out_STAN_gamma <- stan(file = "gamma_glm.stan", 
                       data = data_STAN_gamma,
                       chains = 1, 
                       iter = 2000, 
                       warmup = 1000, 
                       seed = 42)

rstan::traceplot(out_STAN_gamma, pars = c("beta"))
params_STAN_gamma <- As.mcmc.list(out_STAN_gamma, pars = c("beta"))
summary(params_STAN_gamma)
geweke.diag(params_STAN_gamma)

# plotting the densities
df_plot <- data.frame(label = factor(rep(c("beta1", "beta2", "beta3"), each = 1000)), 
                      x = as.vector(params_STAN_gamma[[1]]))

# -----------------------------
# now some plots to compare the estimates
# -----------------------------

ggplot(df_plot) + 
  geom_histogram(aes(x = x, after_stat(density), fill = label), 
                 col = 1, alpha = 0.5, position="identity") + 
  geom_vline(aes(xintercept = 0), lty = 2) + 
  theme_bw() + 
  ylab("") + 
  facet_wrap(~label) +
  theme(legend.position = "null")

ggplot(df_plot) + 
  geom_violin(aes(x = label, y = x, fill = label), trim=FALSE) +
  geom_boxplot(aes(x = label, y = x), width = 0.1) +
  geom_hline(aes(yintercept = 0), lty = 2) + 
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  theme(legend.position = "null")

# -----------------------------
# LINEAR REGRESSION MIXED MODEL
# -----------------------------

# -----------------------------
# The data reported in the slides for the LMM
# and the corresponding design matrices
# -----------------------------

set.seed(123); betatrue <- c(-1, -2, 2); 
gammatrue <- rbind(c(2, 4), c(-2, -4))
z1 <- round(rnorm(100, 0, 1), digits = 1)
z2 <- round(rnorm(100, 0, 1), digits = 1)
z3 <- round(rnorm(100, 0, 1), digits = 1)
z4 <- round(rnorm(100, 0, 1), digits = 1)
c <- rep(c(1, 2), each = 50)
X <- cbind(rep(1, 100), z1, z2)
U <- cbind(z3, z4)
tempmeans <- as.vector(X %*% betatrue) + 
  apply(cbind(U,gammatrue[c,]), 1, function(x) x[1:2] %*% x[3:4]) 
y <- sapply(tempmeans, function(x) rnorm(1, x, 1))

# -----------------------------
# In this first part, we build a Gibbs sampler to perform
# inference with the LMM
# 
# The function is called with four arguments: 
#    - the response variable "y"
#    - the design matrix "X"
#    - the group specific covariates "U"
#    - the group allocation "c"
#    - total number of iterations "niter"
#    - number of burn-in iterations "nburn" 
#    - prior parameters
# -----------------------------

gibbs_LMM <- function(y, X, U, c, niter, nburn, 
                      beta0, Sigma0, tau0, Psi0, a0, b0){
  
  out_beta <- matrix(0, ncol = ncol(X), nrow = niter)
  out_gamma <- array(0, dim = c(niter, max(c), ncol(U)))
  out_sigma2 <- rep(0, niter)
  beta <- rep(0, ncol(X))
  gamma <- matrix(0, nrow = max(c), ncol = ncol(U))
  sigma2 <- 1
  
  q <- ncol(U)
  k <- max(c)
  an <- a0 + length(y) / 2
  
  for(i in 1:niter){
    # update sigma2
    tvec <- (y - X %*% beta - apply(cbind(U,gamma[c,]), 1, 
                                    function(x) x[1:q] %*% x[(q + 1):(2 * q)]))
    bn <- b0 + 0.5 * sum(tvec^2)
    sigma2 <- 1 / rgamma(1, shape = an, scale = 1 / bn)
    
    # update beta
    ygamma <- y - apply(cbind(U,gamma[c,]), 1, 
                        function(x) x[1:q] %*% x[(q + 1):(2 * q)])
    Sigman <- solve(solve(Sigma0) + t(X) %*% X / sigma2)
    betan <- Sigman %*% (solve(Sigma0) %*% beta0 + t(X) %*% ygamma / sigma2)
    beta <- rmvnorm(1, mean = betan, sigma = Sigman)[1,]
    
    # update gamma
    ybeta <- y - X %*% beta
    for(j in 1:k){
      Uj <- U[c == j,]
      Psin <- solve(solve(Psi0) + t(Uj) %*% Uj / sigma2)
      taun <- Psin %*% (Psi0 %*% tau0 + t(Uj) %*% ybeta[c == j] / sigma2)
      gamma[j,] <- rmvnorm(1, mean = taun, sigma = Psin)[1,]
    }
    
    # save the current state
    out_beta[i, ] <- beta
    out_gamma[i,,] <- gamma
    out_sigma2[i] <- sigma2
  }
  return(list(out_beta[-c(1:nburn),], out_gamma[-c(1:nburn),,], out_sigma2[-c(1:nburn)]))
}

# sample from the posterior

out_LMM <- gibbs_LMM(y, X, U, c, 2000, 1000, rep(0, 3), 
                     diag(10^3, 3), rep(0, 2), diag(10^3, 2), 2, 2)

beta_LMM <- as.mcmc(out_LMM[[1]])
gamma1_LMM <- as.mcmc(out_LMM[[2]][,1,])
gamma2_LMM <- as.mcmc(out_LMM[[2]][,2,])
sigma2_LMM <- as.mcmc(out_LMM[[3]])

summary(beta_LMM)
summary(gamma1_LMM)
summary(gamma2_LMM)
summary(sigma2_LMM)

geweke.diag(beta_LMM)
geweke.diag(gamma1_LMM)
geweke.diag(gamma2_LMM)
geweke.diag(sigma2_LMM)

# -----------------------------
# now some plots to compare the estimates
# -----------------------------

df_plot_beta <- data.frame(label = factor(rep(c("beta1", "beta2", "beta3"), each = 1000)), 
                      x = as.vector(beta_LMM))

ggplot(df_plot_beta) + 
  geom_histogram(aes(x = x, after_stat(density), fill = label), 
                 col = 1, alpha = 0.5, position="identity") + 
  geom_vline(aes(xintercept = 0), lty = 2) + 
  theme_bw() + 
  ylab("") + 
  facet_wrap(~label) +
  theme(legend.position = "null")

ggplot(df_plot_beta) + 
  geom_violin(aes(x = label, y = x, fill = label), trim=FALSE) +
  geom_boxplot(aes(x = label, y = x), width = 0.1) +
  geom_hline(aes(yintercept = 0), lty = 2) + 
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  theme(legend.position = "null")

df_plot_gamma <- data.frame(label = factor(rep(c("gamma1", "gamma2","gamma1", "gamma2"), each = 1000)), 
                            group = factor(rep(c("group1", "group2"), each = 2000)), 
                            x = c(as.vector(gamma1_LMM), as.vector(gamma2_LMM)))

ggplot(df_plot_gamma) + 
  geom_histogram(aes(x = x, after_stat(density), group = group, fill = group), 
                 col = 1, alpha = 0.5, position="identity") + 
  geom_vline(aes(xintercept = 0), lty = 2) + 
  theme_bw() + 
  ylab("") + 
  facet_wrap(~label) +
  theme(legend.position = "null")

ggplot(df_plot_gamma) + 
  geom_violin(aes(x = group, y = x, fill = group), trim=FALSE) +
  geom_boxplot(aes(x = group, y = x), width = 0.1) +
  geom_hline(aes(yintercept = 0), lty = 2) + 
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  facet_wrap(~label) +
  theme(legend.position = "null")

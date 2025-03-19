# -----------------------------
# We need the following libraries
# -----------------------------

library(mvtnorm)
library(truncnorm)
library(rstan)
library(coda)
library(BayesLogit)
library(ggplot2)

logsumexp <- function(x){
  M <- max(x)
  lse <- M + log(sum(exp(x - M)))
  return(lse)
}


# -----------------------------
# The data reported in the slides for the LMM
# and the corresponding design matrices
# -----------------------------

set.seed(123); betatrue <- c(-2, 2); 
gammatrue <- rbind(c(-1, 2, 4), c(3, -2, -4))
z1 <- round(rnorm(100, 0, 1), digits = 1)
z2 <- round(rnorm(100, 0, 1), digits = 1)
z3 <- round(rnorm(100, 0, 1), digits = 1)
z4 <- round(rnorm(100, 0, 1), digits = 1)
c <- rep(c(1, 2), each = 50)
X1 <- cbind(z1, z2)
U1 <- cbind(rep(1, 100), z3, z4)
X2 <- cbind(rep(1, 100), z1, z2)
U2 <- cbind(z3, z4)
tempmeans <- as.vector(X1 %*% betatrue) + 
  apply(cbind(U1,gammatrue[c,]), 1, function(x) x[1:2] %*% x[3:4]) 
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
  out_pred_dens <- matrix(0, ncol = length(y), nrow = niter)
  beta <- rep(0, ncol(X))
  gamma <- matrix(0, nrow = max(c), ncol = ncol(U))
  sigma2 <- 1
  
  q <- ncol(U)
  k <- max(c)
  n <- length(y)
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
    
    for(l in 1:n){
      temp_pred <- X[l,] %*% beta + U[l,] %*% gamma[c[l],]
      out_pred_dens[i,l] <- dnorm(y[l], mean = temp_pred, sd = sqrt(sigma2))
    }
  }
  return(list(out_beta[-c(1:nburn),], out_gamma[-c(1:nburn),,], 
              out_sigma2[-c(1:nburn)], out_pred_dens[-c(1:nburn),]))
}

# sample from the posterior

out_LMM1 <- gibbs_LMM(y, X1, U1, c, 2000, 1000, rep(0, 2), 
                     diag(10^3, 2), rep(0, 3), diag(10^3, 3), 2, 2)

out_LMM2 <- gibbs_LMM(y, X2, U2, c, 2000, 1000, rep(0, 3), 
                      diag(10^3, 3), rep(0, 2), diag(10^3, 2), 2, 2)

# check the fit
# model 1

beta_LMM1 <- as.mcmc(out_LMM1[[1]])
gamma1_LMM1 <- as.mcmc(out_LMM1[[2]][,1,])
gamma2_LMM1 <- as.mcmc(out_LMM1[[2]][,2,])
sigma2_LMM1 <- as.mcmc(out_LMM1[[3]])

summary(beta_LMM1)
summary(gamma1_LMM1)
summary(gamma2_LMM1)
summary(sigma2_LMM1)

geweke.diag(beta_LMM1)
geweke.diag(gamma1_LMM1)
geweke.diag(gamma2_LMM1)
geweke.diag(sigma2_LMM1)

# check the fit
# model 2

beta_LMM2 <- as.mcmc(out_LMM2[[1]])
gamma1_LMM2 <- as.mcmc(out_LMM2[[2]][,1,])
gamma2_LMM2 <- as.mcmc(out_LMM2[[2]][,2,])
sigma2_LMM2 <- as.mcmc(out_LMM2[[3]])

summary(beta_LMM2)
summary(gamma1_LMM2)
summary(gamma2_LMM2)
summary(sigma2_LMM2)

geweke.diag(beta_LMM2)
geweke.diag(gamma1_LMM2)
geweke.diag(gamma2_LMM2)
geweke.diag(sigma2_LMM2)

# define the functions to compute LMPL and WAIC
# the functions take as argument
# predictive density values and return
# the corresponding fit measure
# according to the expression of the slides

my_lpml <- function(out_pred_dens){
  out <- sum(log(1 / colMeans(1 / out_pred_dens)))
  return(out)
}

my_waic <- function(out_pred_dens){
  LPPD <-   mean(log(colMeans(out_pred_dens)))
  p_waic <- mean(apply(out_pred_dens, 2, function(x) var(log(x))))
  return(-2 * LPPD + 2 * p_waic)
}

# looking at the models we have
# that the previous measure give the following
# results

out_pred_dens_LMM1 <- out_LMM1[[4]]
out_pred_dens_LMM2 <- out_LMM2[[4]]

my_lpml(out_pred_dens_LMM1)
my_lpml(out_pred_dens_LMM2)

my_waic(out_pred_dens_LMM1)
my_waic(out_pred_dens_LMM2)


# ---------------------------------
# LM comparison
# ---------------------------------

# comparing with a linear model (SLIDE BLOCK 2)
# We have the following gibbs sampler

gibbs_LM <- function(y, X, niter, nburn, beta0, Sigma0, a0, b0){
  
  out_beta <- matrix(0, ncol = ncol(X), nrow = niter)
  out_sigma2 <- rep(0, niter)
  out_pred_dens <- matrix(0, ncol = length(y), nrow = niter)
  beta <- rep(0, ncol(X))
  sigma2 <- 1
  
  n <- length(y)
  an <- a0 + length(y) / 2
  
  
  for(i in 1:niter){
    # update beta
    Sigman <- solve(solve(Sigma0) + t(X) %*% X / sigma2)
    betan <- Sigman %*% (solve(Sigma0) %*% beta0 + t(X) %*% y / sigma2)
    beta <- rmvnorm(1, mean = betan, sigma = Sigman)[1,]
    
    # update sigma2
    bn <- b0 + 0.5 * (t(y) %*% y - t(betan) %*% solve(Sigman) %*% betan + t(beta0) %*% Sigma0 %*% beta0)
    sigma2 <- 1 / rgamma(1, shape = an, scale = 1 / bn)
    
    # save the current state
    out_beta[i, ] <- beta
    out_sigma2[i] <- sigma2
    
    for(l in 1:n){
      temp_pred <- X[l,] %*% beta 
      out_pred_dens[i,l] <- dnorm(y[l], mean = temp_pred, sd = sqrt(sigma2))
    }
  }
  return(list(out_beta[-c(1:nburn),], out_sigma2[-c(1:nburn)], out_pred_dens[-c(1:nburn),]))
}

# sample from the posterior 

out_LM1 <- gibbs_LM(y, X2, 2000, 1000, rep(0, 3), diag(10^3, 3), 2, 2)

# check the fit
# linear model

beta_LM1 <- as.mcmc(out_LM1[[1]])
sigma2_LM1 <- as.mcmc(out_LM1[[2]])

summary(beta_LM1)
summary(sigma2_LM1)

geweke.diag(beta_LM1)
geweke.diag(sigma2_LM1)

# compare LPML and WAIC

out_pred_dens_LM <- out_LM1[[3]]
my_lpml(out_pred_dens_LM)
my_waic(out_pred_dens_LM)

# ---------------------------------
# LM but in STAN 
# ---------------------------------

data_LM_STAN <-list(n = 100, 
                    p = 3,  
                    y = y, 
                    X = X2, 
                    beta0 = rep(0, 3), 
                    Sigma0 = diag(10^3, 3), 
                    a0 = 2, 
                    b0 = 2)

out_STAN <- stan(file = "regression_model.stan", 
                 data = data_LM_STAN,
                 chains = 1, 
                 iter = 2000, 
                 warmup = 1000, 
                 seed = 123)

rstan::traceplot(out_STAN, pars = c("beta", "sigma2"))
param_STAN <- extract(out_STAN, pars = c("beta", "sigma2"), permuted = TRUE)

beta_STAN <- as.mcmc(param_STAN[[1]])
summary(beta_STAN)
geweke.diag(beta_STAN)

sigma2_STAN <- as.mcmc(as.vector(param_STAN[[2]]))
summary(sigma2_STAN)
geweke.diag(sigma2_STAN)

out_pred_dens_STAN <- exp(extract(out_STAN, pars = c("log_lik"))[[1]])
my_lpml(out_pred_dens_STAN)
my_waic(out_pred_dens_STAN)

# ----------------------------
# TEST for logistic regression example
# ----------------------------

set.seed(123); betatrue <- c(-4, 2, -4, 0)
z1 <- round(rnorm(100, 0, 1), digits = 1)
z2 <- round(rnorm(100, 0, 1), digits = 1)
z3 <- round(rnorm(100, 0, 1), digits = 1)
X <- cbind(rep(1, 100), z1, z2, z3)
tempprobs <- exp(X %*% betatrue) / (1 + exp(X %*% betatrue))
y <- sapply(tempprobs[,1], function(x) rbinom(1,1,x))

# -----------------------------
# In this first part, we build a Gibbs sampler to perform
# inference with the logistic model, with the augmentation
# strategy discussed in the slides
# 
# The function is called with four arguments: 
#    - the response variable "y"
#    - the design matrix "X"
#    - total number of iterations "niter"
#    - number of burn-in iterations "nburn"
# -----------------------------

gibbs_augmented_logistic <- function(y, X, niter, nburn, b0, Sigma0){
  
  out_beta <- matrix(0, ncol = ncol(X), nrow = niter)
  out_probs <- matrix(0, ncol = length(y), nrow = niter)
  v <- rep(0, length(y))
  beta <- rep(0, ncol(X))
  n <- length(y)
  
  for(i in 1:niter){
    
    # update the augmented variables 
    for(j in 1:length(y)){
      v[j] <- rpg.devroye(1, 1, X[j,] %*% beta)
    }
    
    # sample the regression coefficients
    Vn <- diag(v)
    Sigman <- solve(solve(Sigma0) + t(X) %*% Vn %*% X)
    bn <- Sigman %*% (solve(Sigma0) %*% b0 + t(X) %*% (y - 0.5))
    beta <- rmvnorm(n = 1, mean = bn, sigma = Sigman)[1,]
    out_beta[i,] <- beta
    
    # evaluate the pmf for each observation
    for(l in 1:n){
      temp_probs <- exp(X[l,] %*% beta) / (1 + exp(X[l,] %*% beta))
      out_probs[i,l] <- ifelse(y[l] == 1, temp_probs, 1 - temp_probs)
    }
    
  }
  
  return(list(out_beta[-c(1:nburn),], out_probs[-c(1:nburn),]))
}

# -----------------------------
# sample from the posterior distribution
# -----------------------------

logistic_M1 <- gibbs_augmented_logistic(y, X, 2500, 500, rep(0, 4), diag(10^3, 4))
logistic_M2 <- gibbs_augmented_logistic(y, X[,1:3], 2500, 500, rep(0, 3), diag(10^3, 3))
logistic_M3 <- gibbs_augmented_logistic(y, X[,1:2], 2500, 500, rep(0, 2), diag(10^3, 2))

mean(apply(logistic_M1[[2]], 1, prod))
mean(apply(logistic_M2[[2]], 1, prod))
mean(apply(logistic_M3[[2]], 1, prod))

log(nrow(logistic_M1[[2]])) - logsumexp(-apply(logistic_M1[[2]], 1, sum))
log(nrow(logistic_M2[[2]])) - logsumexp(-apply(logistic_M2[[2]], 1, sum))
log(nrow(logistic_M3[[2]])) - logsumexp(-apply(logistic_M3[[2]], 1, sum))

# best model is M2, test against the simpler one (M3)

log_marginal_M2 <- log(nrow(logistic_M2[[2]])) - logsumexp(-apply(logistic_M2[[2]], 1, sum))
log_marginal_M3 <- log(nrow(logistic_M3[[2]])) - logsumexp(-apply(logistic_M3[[2]], 1, sum))

BF23 <- exp(log_marginal_M2 - log_marginal_M3)
BF23

# there is a really strong empirical evidence in support of M2

# -----------------------------
# Repeat the same in STAN
# -----------------------------

data_LOGISTIC_STAN1 <-list(n = 100, 
                           p = 4,  
                           y = y, 
                           X = X, 
                           beta0 = rep(0, 4), 
                           Sigma0 = diag(10^3, 4))

data_LOGISTIC_STAN2 <-list(n = 100, 
                           p = 3,  
                           y = y, 
                           X = X[,1:3], 
                           beta0 = rep(0, 3), 
                           Sigma0 = diag(10^3, 3))


data_LOGISTIC_STAN3 <-list(n = 100, 
                           p = 2,  
                           y = y, 
                           X = X[,1:2], 
                           beta0 = rep(0, 2), 
                           Sigma0 = diag(10^3, 2))

out_STAN1 <- stan(file = "logistic_model.stan", 
                  data = data_LOGISTIC_STAN1,
                  chains = 1, 
                  iter = 2000, 
                  warmup = 1000, 
                  seed = 123)

out_STAN2 <- stan(file = "logistic_model.stan", 
                  data = data_LOGISTIC_STAN2,
                  chains = 1, 
                  iter = 2000, 
                  warmup = 1000, 
                  seed = 123)

out_STAN3 <- stan(file = "logistic_model.stan", 
                  data = data_LOGISTIC_STAN3,
                  chains = 1, 
                  iter = 2000, 
                  warmup = 1000, 
                  seed = 123)

rstan::traceplot(out_STAN1, pars = c("beta"))
rstan::traceplot(out_STAN2, pars = c("beta"))
rstan::traceplot(out_STAN3, pars = c("beta"))

param_STAN1 <- as.mcmc(extract(out_STAN1, pars = c("beta"), permuted = TRUE)[[1]])
param_STAN2 <- as.mcmc(extract(out_STAN2, pars = c("beta"), permuted = TRUE)[[1]])
param_STAN3 <- as.mcmc(extract(out_STAN3, pars = c("beta"), permuted = TRUE)[[1]])

summary(param_STAN1)
summary(param_STAN2)
summary(param_STAN3)

geweke.diag(param_STAN1)
geweke.diag(param_STAN2)
geweke.diag(param_STAN3)

out_pred_dens_STAN1 <- exp(extract(out_STAN1, pars = c("log_lik"))[[1]])
out_pred_dens_STAN2 <- exp(extract(out_STAN2, pars = c("log_lik"))[[1]])
out_pred_dens_STAN3 <- exp(extract(out_STAN3, pars = c("log_lik"))[[1]])

mean(apply(out_pred_dens_STAN1, 1, prod))
mean(apply(out_pred_dens_STAN2, 1, prod))
mean(apply(out_pred_dens_STAN3, 1, prod))

log(nrow(out_pred_dens_STAN1)) - logsumexp(-apply(out_pred_dens_STAN1, 1, sum))
log(nrow(out_pred_dens_STAN2)) - logsumexp(-apply(out_pred_dens_STAN2, 1, sum))
log(nrow(out_pred_dens_STAN3)) - logsumexp(-apply(out_pred_dens_STAN3, 1, sum))

# -----------------------------
# VARIABLE SELECTION
# -----------------------------

# -----------------------------
# The data reported in the slides for the variable selection example
# -----------------------------

set.seed(123); betatrue <- c(10, 0, 0, -5, 4, 0, -8, -3, rep(0, 12))
X <- cbind(rep(1, 100), round(matrix(rnorm(1900), ncol = 19), digits = 2))
tempmeans <- X %*% betatrue
y <- sapply(tempmeans, function(x) rnorm(1, x, 1))

# ---------------------------------
# SSVS in STAN 
# ---------------------------------

data_SSVS_STAN <-list(n = 100, 
                    p = 20,  
                    y = y, 
                    X = X, 
                    beta0 = rep(0, 20), 
                    a0 = 2, 
                    b0 = 2, 
                    alpha1 = 1, 
                    alpha2 = 1, 
                    c = 0.001, 
                    tau2 = 10^3)

out_SSVS_STAN <- stan(file = "regression_SSVS_model.stan", 
                      data = data_SSVS_STAN,
                      chains = 1, 
                      iter = 2000, 
                      warmup = 1000, 
                      seed = 123)

rstan::traceplot(out_SSVS_STAN, pars = c("sigma2"))
rstan::traceplot(out_SSVS_STAN, pars = c("beta"))

param_SSVS1 <- as.mcmc(as.vector(extract(out_SSVS_STAN, pars = c("sigma2"), permuted = TRUE)[[1]]))
param_SSVS2 <- as.mcmc(extract(out_SSVS_STAN, pars = c("beta"), permuted = TRUE)[[1]])

summary(param_SSVS1)
summary(param_SSVS2)

geweke.diag(param_SSVS1)
geweke.diag(param_SSVS2)

# compute the bounds

eps <- sqrt(2 * (log(0.001) * 0.001^2) / (0.001^2 - 1))
kappa <- sqrt(10^3) * eps

# reconstruct the binary matrix

gamma_matrix <- ifelse(abs(param_SSVS2) > kappa, 1, 0)

# HPD estimate

unique_model <- unique(gamma_matrix, MARGIN  = 1)
freq <- apply(unique_model, 1, function(b) sum(apply(gamma_matrix, MARGIN = 1, function(a) all(a == b))))
HPD_model <- unique_model[which.max(freq),]

# MPM estimate

MPM_model <- as.numeric(colMeans(gamma_matrix) > 0.5)

# HS estimate

HS_model <- as.numeric(colMeans(gamma_matrix) == 1)

# with the true coefficients

cbind(HPD_model, MPM_model, HS_model, betatrue)

# more plots

ggplot(data.frame(value = colMeans(gamma_matrix), idx = 1:20,
                  var = factor(paste0("covariate ", c(1:20)), levels = paste0("covariate ", c(20:1))), 
                  inc_true = factor(betatrue != 0))) + 
  geom_bar(aes(y = value, x = var, fill = inc_true), stat="identity", alpha = 0.5, col = 1) + 
  geom_hline(mapping = aes(yintercept = .5), col = 2, lty = 2) +
  coord_flip() + 
  theme_minimal() + 
  theme(legend.position="none") + 
  ylab("posterior inclusion probabilities") + 
  xlab("")

# -----------------------------
# One more hierarchy in the model
# -----------------------------

set.seed(123); betatrue <- c(2, -3, 4, -2)
X <- cbind(rep(1, 100), round(matrix(rnorm(300), ncol = 3), digits = 2))
tempmeans <- X %*% betatrue
y <- sapply(tempmeans, function(x) rnorm(1, x, 1))

# estimating the model with no hyperpriors

data_NO_HYPER_STAN <-list(n = 100, 
                          p = 4,  
                          y = y, 
                          X = X, 
                          beta0 = rep(25, 4), 
                          Sigma0 = diag(0.1, 4), 
                          a0 = 100, 
                          b0 = 1)

out_NO_HYPER_STAN <- stan(file = "regression_model.stan", 
                          data = data_NO_HYPER_STAN,
                          chains = 1, 
                          iter = 2000, 
                          warmup = 1000, 
                          seed = 123)

# looking at the estimate

rstan::traceplot(out_NO_HYPER_STAN, pars = c("sigma2"))
rstan::traceplot(out_NO_HYPER_STAN, pars = c("beta"))

param_NO_HYPER1 <- as.mcmc(as.vector(extract(out_NO_HYPER_STAN, pars = c("sigma2"), permuted = TRUE)[[1]]))
param_NO_HYPER2 <- as.mcmc(extract(out_NO_HYPER_STAN, pars = c("beta"), permuted = TRUE)[[1]])

summary(param_NO_HYPER1)
summary(param_NO_HYPER2)

geweke.diag(param_NO_HYPER1)
geweke.diag(param_NO_HYPER2)

# estimating the model with hyperpriors

data_HYPER_STAN <-list(n = 100, 
                       p = 4,  
                       y = y, 
                       X = X, 
                       a1 = 0.1, 
                       b1 = 0.1, 
                       c1 = 0.1, 
                       d1 = 0.1,
                       beta1 = rep(0, 4), 
                       Sigma1 = diag(10^3, 4), 
                       nu1 = 6,
                       Phi1 = diag(1, 4))

out_HYPER_STAN <- stan(file = "regression_model_hyper.stan", 
                       data = data_HYPER_STAN,
                       chains = 1, 
                       iter = 2000, 
                       warmup = 1000, 
                       seed = 123)

# looking at the estimate

rstan::traceplot(out_HYPER_STAN, pars = c("sigma2"))
rstan::traceplot(out_HYPER_STAN, pars = c("beta"))

param_HYPER1 <- as.mcmc(as.vector(extract(out_HYPER_STAN, pars = c("sigma2"), permuted = TRUE)[[1]]))
param_HYPER2 <- as.mcmc(extract(out_HYPER_STAN, pars = c("beta"), permuted = TRUE)[[1]])

summary(param_HYPER1)
summary(param_HYPER2)

geweke.diag(param_HYPER1)
geweke.diag(param_HYPER2)

# comparison plot

data_plt <- data.frame(x = c(as.vector(param_NO_HYPER2), as.vector(param_HYPER2)), 
              param = factor(c(rep(c("beta1", "beta2", "beta3", "beta4"), each = 1000), 
                                  rep(c("beta1", "beta2", "beta3", "beta4"), each = 1000))), 
              group = factor(rep(c("model 1", "model 2"), each = 4000)))

ggplot(data_plt) + 
  geom_violin(aes(x = param, y = x, fill = group), alpha = 1, trim=FALSE) +
  geom_boxplot(aes(x = param, y = x), width = 0.1, position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept = 0), lty = 2) + 
  facet_wrap(.~group) +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  theme(legend.position = "null")

data_plt2 <- data.frame(x = c(param_NO_HYPER1, param_HYPER1), 
                        group = factor(rep(c("model 1", "model 2"), each = 1000)))

ggplot(data_plt2) + 
  geom_histogram(aes(x = x, after_stat(density), fill = group), 
                 col = 1, alpha = 0.5, position="identity") + 
  theme_bw() + 
  ylab("") + 
  theme(legend.position = "null")
ggplot(data_plt2) + 
  geom_histogram(aes(x = log(x), after_stat(density), fill = group), 
                 col = 1, alpha = 0.5, position="identity") + 
  theme_bw() + 
  ylab("") + 
  theme(legend.position = "null")

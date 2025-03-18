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
# test if beta2 is significantly greater than 0
# ... the following is a Bayes factor, do you agree?
# -----------------------------

xn1 <- c(1, 2, 4)
pn1 <- pnorm(out_augmented_probit %*% xn1)
yn1 <- sapply(pn1, function(x) sample(c(1,0), size = 1, prob = c(x, 1 - x)))
table(yn1)

# -----------------------------
# LOGISTIC REGRESSION MODEL
# -----------------------------

# -----------------------------
# The data reported in the slides for the probit model
# and the corresponding design matrix
# -----------------------------

set.seed(123); betatrue <- c(-1, -2, 2)
z <- round(rnorm(100, 0, 1), digits = 1)
X <- cbind(rep(1, 100), z, z^2)
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
  v <- rep(0, length(y))
  beta <- rep(0, ncol(X))
  
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
  }
  
  return(out_beta[-c(1:nburn),])
}

# -----------------------------
# sample from the posterior distribution
# -----------------------------

out_augmented_logistic <- 
  gibbs_augmented_logistic(y, X, 1200, 200, rep(0, 3), diag(10^3, 3))

params_STAN_logistic <- as.mcmc(out_augmented_logistic)
summary(params_STAN_logistic)
geweke.diag(params_STAN_logistic)

# plotting the densities
df_plot <- data.frame(label = factor(rep(c("beta1", "beta2", "beta3"), each = 1000)), 
                      x = as.vector(out_augmented_logistic))

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
# extra plot 
# -----------------------------

z_grid <- seq(-3, 3, by = 0.1)
X_grid <- cbind(rep(1, length(z_grid)), z_grid, z_grid^2)
out_plot <- t(apply(out_augmented_logistic, 1, function(x) exp(X_grid %*% x) / (1 + exp(X_grid %*% x))))
out_median <- apply(out_plot, 2, median)
out_low <- apply(out_plot, 2, quantile, p = 0.1)
out_up <- apply(out_plot, 2, quantile, p = 0.9)

# generating the fitted line for the ML estimate and the true model

m1 <- glm(y ~ 0 + X, family = binomial(link = "logit"))
freq_line <- exp(X_grid %*% m1$coefficients) / (1 + exp(X_grid %*% m1$coefficients))
true_line <- exp(X_grid %*% betatrue) / (1 + exp(X_grid %*% betatrue))

ggplot(data.frame(x = z_grid, y = out_median, 
                  ylow = out_low, yup = out_up)) + 
  geom_ribbon(aes(x = x, ymin = ylow, ymax = yup), fill = "#100f7a", alpha = 0.2) +
  geom_line(aes(x = x, y = y), col = "#100f7a") + 
  geom_line(data = data.frame(x = z_grid, y = freq_line), 
            aes(x = x, y = y), col = "#7a0f10") +
  geom_line(data = data.frame(x = z_grid, y = true_line), 
            aes(x = x, y = y), col = "#7a0f10", lty = 2) +
  theme_bw() + ylab("") + ggtitle("Posterior distribution")

# filled blue line - posterior median
# filled red line - ML estimate
# dashed red line - true model

ggplot(data.frame(x = - out_augmented_logistic[,3] / (2 * out_augmented_logistic[,2]))) + 
  geom_histogram(aes(x = x, after_stat(density)), 
                 alpha = 0.2, position="identity", col = "#100f7a", fill = "#100f7a") + 
  geom_vline(aes(xintercept = mean(- out_augmented_logistic[,3] / (2 * out_augmented_logistic[,2]))), col = "#7a0f10", lty = 2) +
  theme_bw() + ylab("") + ggtitle("Regime change")

# -----------------------------
# test if beta2 is significantly greater than 0
# ... the following is a Bayes factor, do you agree?
# -----------------------------

mean(out_augmented_logistic[,2] > 0) / 
  mean(out_augmented_logistic[,2] < 0)

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
params_STAN_poi <- as.mcmc(extract(out_STAN_poi, pars = c("beta"), permuted = TRUE))
summary(params_STAN_poi$beta)
geweke.diag(params_STAN_poi$beta)

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
params_STAN_gamma <- as.mcmc(extract(out_STAN_gamma, pars = c("beta"), permuted = TRUE))
summary(params_STAN_gamma$beta)
geweke.diag(params_STAN_gamma$beta)

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

df_plot_gamma <- data.frame(label = factor(rep(c("beta1", "beta2","beta1", "beta2"), each = 1000)), 
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

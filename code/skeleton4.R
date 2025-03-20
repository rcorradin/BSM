# -----------------------------
# We need the following libraries
# -----------------------------

library(mvtnorm)
library(truncnorm)
library(rstan)
library(coda)
library(BayesLogit)
library(ggplot2)

# ----------------------------
# EXAMPLE 1
# ----------------------------

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
# LMM
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
    
  }
  return(list(out_beta[-c(1:nburn),], out_gamma[-c(1:nburn),,], 
              out_sigma2[-c(1:nburn)]))
}

# -----------------------------
# sample from the posterior distribution and compute LPML / WAIC
# -----------------------------

out_LMM1 <- gibbs_LMM(y, X1, U1, c, 2000, 1000, rep(0, 2), 
                      diag(10^3, 2), rep(0, 3), diag(10^3, 3), 2, 2)

out_LMM2 <- gibbs_LMM(y, X2, U2, c, 2000, 1000, rep(0, 3), 
                      diag(10^3, 3), rep(0, 2), diag(10^3, 2), 2, 2)


# ---------------------------------
# LM 
# ---------------------------------

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
    
  }
  return(list(out_beta[-c(1:nburn),], out_sigma2[-c(1:nburn)]))
}

# sample from the posterior 

out_LM1 <- gibbs_LM(y, X2, 2000, 1000, rep(0, 3), diag(10^3, 3), 2, 2)

# ----------------------------
# EXAMPLE 2
# ----------------------------
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

gibbs_augmented_logistic <- function(y, X, niter, nburn, b0, Sigma0){
  
  out_beta <- matrix(0, ncol = ncol(X), nrow = niter)
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
  }
  
  return(list(out_beta[-c(1:nburn),]))
}

# -----------------------------
# sample from the posterior distribution and compute the BF
# -----------------------------


# -----------------------------
# Repeat the same in STAN
# -----------------------------

data_LOGISTIC_STAN1 <-list(n = 100, 
                           p = 4,  
                           y = y, 
                           X = X, 
                           beta0 = rep(0, 4), 
                           Sigma0 = diag(10^3, 4))

data_LOGISTIC_STAN2 <-list()


data_LOGISTIC_STAN3 <-list()

out_STAN1 <- stan(file = "logistic_model.stan", 
                  data = data_LOGISTIC_STAN1,
                  chains = 1, 
                  iter = 2000, 
                  warmup = 1000, 
                  seed = 123)

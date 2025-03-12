# -----------------------------
# We need the following libraries
# -----------------------------

library(mvtnorm)
library(ggplot2)
library(bayestestR)

# ----------------------------

data <- read.csv("frogs.csv", head = T)
head(data)
y <- log(data$BoW)
X <- data.frame(intercept = rep(1, length(y)), 
                family = as.numeric(as.factor(data[,2])) - 1, 
                logSVL = log(data[,6]), 
                logBrW = log(data[,7]))

# remove rows with missing data
cond <- apply(X, 1, function(z) sum(is.na(z)) == 0)
X <- as.matrix(X[cond,])
y <- y[cond]

# we want to model the "brain weights" as function of 
# the other quantities
#   - Family (2 categories)
#   - BoW / body weight
#   - SVL / length nose to tail

# implement a function to sample from the posterior distribution of 
# interest, using a Gibbs sampler 

sample_linear_regression <- function(y, X, niter, beta0, 
                                     Sigma0, a0, b0){
  
  # define the parameters of the posterior distribution
  an <- a0 + length(y) / 2
  Sigman <- solve(solve(Sigma0) + t(X) %*% X)
  betan <- as.vector(Sigman %*% (solve(Sigma0) %*% beta0 + t(X) %*% y))
  bn <- b0 + 0.5 * (sum(y^2) - betan %*% solve(Sigman) %*% betan + 
                      beta0 %*% solve(Sigma0) %*% beta0)
  
  # sample from the posterior distribution
  sigma2_out <- 1 / rgamma(niter, shape = an, scale = 1 / bn)
  beta_out <- t(sapply(sigma2_out, function(z) 
    rmvnorm(1, mean = betan, sigma = z * Sigman)))
  
  # return the output
  return(list(beta_out, sigma2_out))
}

# sample from the posterior
sample_post <- sample_linear_regression(y, as.matrix(X), 1000, 
                                        rep(0, ncol(X)), 
                                        diag(100, ncol(X)), 2, 10)

# plot the result
df_plot_beta <- data.frame(x = as.vector((sample_post[[1]])), 
                           coef = rep(paste0("beta", 1:ncol(X)), each = 1000))

ggplot(df_plot_beta) + 
  geom_histogram(aes(x = x, y = after_stat(density), fill = coef), col = 1, alpha = 0.5) + 
  theme_bw() + 
  theme(legend.position = "null") + 
  geom_vline(aes(xintercept = 0), col = 2, lty = 2) + 
  facet_wrap(.~coef, scales = "free")

# point estimates

colMeans(sample_post[[1]])

# credible intervals

ci_hdi_beta1 <- ci(sample_post[[1]][,1], method = "HDI")
ci_hdi_beta2 <- ci(sample_post[[1]][,2], method = "HDI")
ci_hdi_beta3 <- ci(sample_post[[1]][,3], method = "HDI")
ci_hdi_beta4 <- ci(sample_post[[1]][,4], method = "HDI")
ci_hdi_beta1
ci_hdi_beta2
ci_hdi_beta3
ci_hdi_beta4

ci_hdi_beta1 <- ci(sample_post[[1]][,1], method = "ETI")
ci_hdi_beta2 <- ci(sample_post[[1]][,2], method = "ETI")
ci_hdi_beta3 <- ci(sample_post[[1]][,3], method = "ETI")
ci_hdi_beta4 <- ci(sample_post[[1]][,4], method = "ETI")
ci_hdi_beta1
ci_hdi_beta2
ci_hdi_beta3
ci_hdi_beta4

# new model estimate

X2 <- X[,-4]

# sample from the posterior
sample_post2 <- sample_linear_regression(y, as.matrix(X2), 1000, rep(0, ncol(X2)), 
                                        diag(100, ncol(X2)), 2, 10)

# plot the result
df_plot_beta2 <- data.frame(x = as.vector((sample_post2[[1]])), 
                           coef = rep(paste0("beta", 1:ncol(X2)), each = 1000))

ggplot(df_plot_beta2) + 
  geom_histogram(aes(x = x, y = after_stat(density), fill = coef), col = 1, alpha = 0.5) + 
  theme_bw() + 
  theme(legend.position = "null") + 
  geom_vline(aes(xintercept = 0), col = 2, lty = 2) + 
  facet_wrap(.~coef, scales = "free")

# BF function

BF_linear_model <- function(y, X1, X2, Sigma01, Sigma02, beta01, beta02, a0, b0){
  
  # define the parameters of the posterior distribution
  an <- a0 + length(y) / 2
  
  Sigman1 <- solve(solve(Sigma01) + t(X1) %*% X1)
  Sigman2 <- solve(solve(Sigma02) + t(X2) %*% X2)
  
  betan1 <- as.vector(Sigman1 %*% (solve(Sigma01) %*% beta01 + t(X1) %*% y))
  betan2 <- as.vector(Sigman2 %*% (solve(Sigma02) %*% beta02 + t(X2) %*% y))
  
  bn1 <- b0 + 0.5 * (sum(y^2) - betan1 %*% solve(Sigman1) %*% betan1 + 
                      beta01 %*% solve(Sigma01) %*% beta01)
  bn2 <- b0 + 0.5 * (sum(y^2) - betan2 %*% solve(Sigman2) %*% betan2 + 
                       beta02 %*% solve(Sigma02) %*% beta02)
  out <- sqrt(det(Sigman2) / det(Sigman1)) * (bn1 / bn2)^an
  return(as.numeric(out))
}

BF_linear_model(y, X, X2, diag(100, 4), diag(100, 3), 
                rep(0, 4), rep(0, 3), 2, 10)

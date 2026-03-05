# -----------------------------
# We need the following libraries
# -----------------------------

library(mvtnorm)
library(ggplot2)
library(bayestestR)
library(rstan)
library(coda)

logsumexp <- function(x){
  M <- max(x)
  lse <- M + log(sum(exp(x - M)))
  return(lse)
}

# ----------------------------
# CASE STUDY 1 
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

ci_eti_beta1 <- ci(sample_post[[1]][,1], method = "ETI")
ci_eti_beta2 <- ci(sample_post[[1]][,2], method = "ETI")
ci_eti_beta3 <- ci(sample_post[[1]][,3], method = "ETI")
ci_eti_beta4 <- ci(sample_post[[1]][,4], method = "ETI")
ci_eti_beta1
ci_eti_beta2
ci_eti_beta3
ci_eti_beta4

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

# ----------------------------
# CASE STUDY 2 
# ----------------------------

data <- read.csv("grouse.csv", head = T)
y <- data$TICKS
X <- data.frame(x1 = rep(1, length(y)),  
                x2 = data$HEIGHT, 
                g2 = ifelse(data$YEAR == 96, 1, 0), 
                g3 = ifelse(data$YEAR == 97, 1, 0)) 

# ----------------------------

data_POI <-list(N = length(y), 
                p = 4, 
                y = y, 
                X = X, 
                Sigma0 = diag(10^2, 4), 
                b0 = rep(0, 4))

#----- FIT THE MODEL

Poi_GLM <- stan(file = "Poi_GLM.stan", 
                 data = data_POI,
                 chains = 1, 
                 iter = 5000, 
                 warmup = 2000, 
                 seed = 42)

# Show the traceplots

rstan::traceplot(Poi_GLM, pars = "beta", inc_warmup = TRUE)
rstan::traceplot(Poi_GLM, inc_warmup = TRUE)

param_GLM <- As.mcmc.list(Poi_GLM, pars = c("beta"))
summary(param_GLM)
geweke.diag(param_GLM)

# CI coefficients

ci(param_GLM, method = "ETI")
ci(param_GLM, method = "HDI")

# PLOT

params <- rstan::extract(Poi_GLM, c("beta"), perm = T)
x_seq <- seq(380, 570, length.out = 100)
jpar <- cbind(params[[1]][,1:2], rep(0, length(par)), params[[1]][,3:4])
tres <- array(0, dim = c(nrow(jpar), length(x_seq), 3))
for(i in 1:nrow(jpar)){
  for(j in 1:3){
    tres[i,,j] <- exp(jpar[i,1] + jpar[i,2] * x_seq + jpar[i,2 + j])
  }
}

df_plot <- data.frame(x = rep(x_seq, 3),
                      ymean = as.vector(apply(tres, c(2,3), mean)),
                      year = as.factor(c(rep(95, 100), rep(96, 100), rep(97, 100))))
pdf_plot <- data.frame(x = X[,2], y = y, year = as.factor(data$YEAR))

ggplot() + 
  geom_point(data = pdf_plot, mapping = aes(x = x, y = y, col = year)) + 
  geom_line(data = df_plot, mapping = aes(x = x, y = ymean, col = year)) +
  theme_minimal() + 
  theme(legend.position="bottom") + 
  ylab("Numbers of ticks on red grouse") + 
  xlab("") 

# ----------------------------

data_POI <-list(N = length(y), 
                p = 2, 
                y = y, 
                X = X[,1:2], 
                Sigma0 = diag(10^2, 2), 
                b0 = rep(0, 2))

#----- FIT THE MODEL

Poi_GLM2 <- stan(file = "Poi_GLM.stan", 
                data = data_POI,
                chains = 1, 
                iter = 5000, 
                warmup = 2000, 
                seed = 42)

# Show the traceplots

rstan::traceplot(Poi_GLM, pars = "beta", inc_warmup = TRUE)
rstan::traceplot(Poi_GLM, pars = "beta", inc_warmup = FALSE)

param_GLM <- As.mcmc.list(Poi_GLM, pars = c("beta"))
summary(param_GLM)
geweke.diag(param_GLM)

# TEST --------------

llik_GLM <- extract(Poi_GLM, "log_lik")[[1]]
marginal_GLM <- log(nrow(llik_GLM)) - logsumexp(-apply(llik_GLM, 1, sum))

llik_GLM2 <- extract(Poi_GLM2, "log_lik")[[1]]
marginal_GLM2 <- log(nrow(llik_GLM2)) - logsumexp(-apply(llik_GLM2, 1, sum))

# BAYES FACTOR
exp(marginal_GLM - marginal_GLM2)

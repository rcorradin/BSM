# -----------------------------
# We need the following libraries
# -----------------------------

library(ggplot2)
library(rstan)
library(coda)
library(ggpubr)
library(bridgesampling)
library(bayestestR)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)
require(gplots)
require(ggpubr)

logsumexp <- function(x){
  M <- max(x)
  lse <- M + log(sum(exp(x - M)))
  return(lse)
}


# ----------------------------
# frogs data (again)
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

# -------------------
# possible models
# -------------------

models <- as.matrix(expand.grid(1, c(0,1), c(0,1), c(0,1)))

# -------------------
# functions for WAIC / LPML
# -------------------

my_waic <- function(out_pred_dens){
  LPPD <-   mean(log(colMeans(out_pred_dens)))
  p_waic <- mean(apply(out_pred_dens, 2, function(x) var(log(x))))
  return(-2 * LPPD + 2 * p_waic)
}

my_lpml <- function(out_pred_dens){
  out <- sum(log(1 / colMeans(1 / out_pred_dens)))
  return(out)
}

# -------------------
# run the models and compute WAIC / LPML
# -------------------

WAIC_LPML_out <- matrix(0, ncol = 2, nrow = nrow(models))
out_dens_values <- list()
out_coef_values <- list()
model_list <- list()

for(m in 1:nrow(models)){
  if(m == 1){
    X_temp <- matrix(X[,as.logical(models[m,])])
  } else {
    X_temp <- X[,as.logical(models[m,])]
  }
  data_temp <-list(n = nrow(X_temp), 
                   p = ncol(X_temp),  
                   y = y, 
                   X = X_temp, 
                   beta0 = array(0, dim = ncol(X_temp)),
                   Sigma0 = diag(10^3, ncol(X_temp)), 
                   a0 = 2, 
                   b0 = 10)
  temp_regression <- stan(file = "regression.stan", 
                          data = data_temp,
                          chains = 1, 
                          iter = 2500, 
                          warmup = 500, 
                          seed = 42)
  model_list[[m]] <- temp_regression
  out_coef_values[[m]] <- extract(temp_regression, pars = c("beta"))
  out_dens_values[[m]] <- exp(extract(temp_regression, pars = c("log_lik"))[[1]])
  WAIC_LPML_out[m,1] <- my_waic(out_dens_values[[m]])
  WAIC_LPML_out[m,2] <- my_lpml(out_dens_values[[m]])
}

print(cbind(WAIC_LPML_out, models))

# ----------------------------
# test VS simple model
# ----------------------------

marginal_BEST <- bridge_sampler(model_list[[4]])$logml
marginal_SIMPLE <- bridge_sampler(model_list[[1]])$logml
marginal_COMPLEX <- bridge_sampler(model_list[[8]])$logml

wrong_marginal_BEST <- log(nrow(out_dens_values[[4]])) - logsumexp(-apply(out_dens_values[[4]], 1, sum))
wrong_marginal_SIMPLE <- log(nrow(out_dens_values[[1]])) - logsumexp(-apply(out_dens_values[[1]], 1, sum))
wrong_marginal_COMPLEX <- log(nrow(out_dens_values[[8]])) - logsumexp(-apply(out_dens_values[[8]], 1, sum))

# BAYES FACTOR
exp(marginal_BEST - marginal_SIMPLE)
exp(marginal_BEST - marginal_COMPLEX)

# BAYES FACTOR WRONG
exp(wrong_marginal_BEST - wrong_marginal_SIMPLE)
exp(wrong_marginal_BEST - wrong_marginal_COMPLEX)

# ----------------------------
# MODEL AVERAGING
# ----------------------------

# construct the weights
# same prior prob, posterior prop to marginal

wei <- c()
for(i in 1:nrow(models)){
  wei[i] <- bridge_sampler(model_list[[i]])$logml
}
wei <- exp(wei) / sum(exp(wei))

# compute the weights
coef <- matrix(0, nrow = nrow(models), ncol = 4)
for(i in 1:nrow(models)){
  temp_coef <- rep(0, 4)
  temp_coef[which(models[i,] == 1)] <- colMeans(out_coef_values[[i]]$beta)
  coef[i,] <- temp_coef
}

# BMA estiamte
BMA_est <- apply(coef, 2, function(x) x %*% wei)
BMA_est 

# computations
out <- array(0, dim = c(nrow(models), 4, 2))
for(i in 1:nrow(models)){
  temp_var <- rep(0, 4)
  temp_var[which(models[i,] == 1)] <- apply(out_coef_values[[i]]$beta, 2, var)
  out[i,,1] <- temp_var
  
  
  temp_coef <- rep(0, 4)
  temp_coef[which(models[i,] == 1)] <- colMeans(out_coef_values[[i]]$beta)
  out[i,,2] <- (temp_coef - BMA_est)^2
}
# marginal variances
rowSums(apply(out, c(2,3), function(x) x%*% wei))


# ----------------------------
# reach data
# ----------------------------

data <- read.table("reach.txt", head = T)
y <- data$outcome
data$gender <- data$gender - 1
X <- as.matrix(cbind(intercept = rep(1, length(y)), data[,-13]))

# select the best model using spike-and-slab priors

c = 100
kappa = 0.05
tau = kappa / sqrt(2 * log(c) * c^2 / (c^2 - 1))

data_STAN <-list(n = nrow(X), 
                 p = ncol(X),  
                 y = y, 
                 X = X, 
                 beta0 = rep(0, ncol(X)), 
                 alpha1 = 1, 
                 alpha2 = 1, 
                 c2 = c^2, 
                 tau2 = tau^2)

probit_variable_selection <- stan(file = "probit_SpSl.stan", 
                                  data = data_STAN,
                                  chains = 1, 
                                  iter = 12000, 
                                  warmup = 2000, 
                                  seed = 123)

rstan::traceplot(probit_variable_selection, pars = c("beta"))
params <- as.mcmc(extract(probit_variable_selection, pars = c("beta"))[[1]])
summary(params)
geweke.diag(params)

# finding the best model

gamma_matrix <- ifelse(abs(params) > kappa, 1, 0)

# HPD estimate

unique_model <- unique(gamma_matrix, MARGIN  = 1)
freq <- apply(unique_model, 1, function(b) sum(apply(gamma_matrix, MARGIN = 1, function(a) all(a == b))))
HPD_model <- unique_model[which.max(freq),]

# MPM estimate

MPM_model <- as.numeric(colMeans(gamma_matrix) > 0.5)

# HS estimate

HS_model <- as.numeric(colMeans(gamma_matrix) == 1)

# plot everything

p1 <- ggplot(data.frame(value = colMeans(gamma_matrix), idx = 1:13,
                        var = factor(colnames(X)),
                        HPD_model_inc = factor(HPD_model))) + 
  geom_bar(aes(y = value, x = var, fill = HPD_model_inc), stat="identity", alpha = 0.5, col = 1) + 
  geom_hline(mapping = aes(yintercept = .5), col = 2, lty = 2) +
  coord_flip() + 
  theme_minimal() + 
  theme(legend.position="none") + 
  ylab("posterior inclusion probabilities") + 
  xlab("") + 
  ggtitle("Highest posterior probability")

p2 <- ggplot(data.frame(value = colMeans(gamma_matrix), idx = 1:13,
                        var = factor(colnames(X)),
                        MPM_model_inc = factor(MPM_model))) + 
  geom_bar(aes(y = value, x = var, fill = MPM_model_inc), stat="identity", alpha = 0.5, col = 1) + 
  geom_hline(mapping = aes(yintercept = .5), col = 2, lty = 2) +
  coord_flip() + 
  theme_minimal() + 
  theme(legend.position="none") + 
  ylab("posterior inclusion probabilities") + 
  xlab("") + 
  ggtitle("Median probability model")

p3 <- ggplot(data.frame(value = colMeans(gamma_matrix), idx = 1:13,
                        var = factor(colnames(X)),
                        HS_model_inc = factor(HS_model))) + 
  geom_bar(aes(y = value, x = var, fill = HS_model_inc), stat="identity", alpha = 0.5, col = 1) + 
  geom_hline(mapping = aes(yintercept = .5), col = 2, lty = 2) +
  coord_flip() + 
  theme_minimal() + 
  theme(legend.position="none") + 
  ylab("posterior inclusion probabilities") + 
  xlab("") + 
  ggtitle("Hard Shrinkage")

ggarrange(p1, p2, p3, nrow = 1)

# ----------------------------
# Ant data
# ----------------------------

data <- read.csv("ant.csv", head = T)
y <- data$num.fpodz
X <- as.vector(data$trap.days)
Z <- data.frame(z1 = data$summer.precip, 
                z2 = data$summer.temp)
group <- ifelse(data$elevation == "High", 1, 2)
group_nested <- rep(1:6, each = 10)

# ----------------------------

data_POI <-list(N = length(y), 
                p_fix = 1, 
                p_rand = 2,
                y = y, 
                X = X, 
                Z = Z, 
                ngr = 2,
                ngr_nested = 6,
                group = group,
                group_nested = group_nested,
                psi20 = 10^2,
                sigma20 = 10^2, 
                tau20 = 10^2, 
                eta20 = 10^2,
                lambda20 = 10^2)

#----- FIT THE MODEL

Poi_GLMM <- stan(file = "Poi_GLMM.stan", 
                  data = data_POI,
                  chains = 1, 
                  iter = 5000, 
                  warmup = 2000, 
                  seed = 42)

# Show the traceplots

rstan::traceplot(Poi_GLMM, pars = "xi", inc_warmup = TRUE)
rstan::traceplot(Poi_GLMM, pars = "beta", inc_warmup = TRUE)
rstan::traceplot(Poi_GLMM, pars= "theta", inc_warmup = TRUE)
rstan::traceplot(Poi_GLMM, pars= "gamma", inc_warmup = TRUE)
rstan::traceplot(Poi_GLMM, pars = "alpha", inc_warmup = TRUE)

rstan::traceplot(Poi_GLMM, pars = "xi", inc_warmup = FALSE)
rstan::traceplot(Poi_GLMM, pars = "beta", inc_warmup = FALSE)
rstan::traceplot(Poi_GLMM, pars= "theta", inc_warmup = FALSE)
rstan::traceplot(Poi_GLMM, pars= "gamma", inc_warmup = FALSE)
rstan::traceplot(Poi_GLMM, pars = "alpha", inc_warmup = FALSE)

param_GLMM <- As.mcmc.list(Poi_GLMM, pars = c("xi", "beta", "theta", "gamma", "alpha"))
summary(param_GLMM)
geweke.diag(param_GLMM)

ci(param_GLMM, method = "ETI")
ci(param_GLMM, method = "HDI")

#----- FIT THE MODEL 2

Poi_GLMM2 <- stan(file = "Poi_GLMM2.stan", 
                 data = data_POI,
                 chains = 1, 
                 iter = 5000, 
                 warmup = 2000, 
                 seed = 42)

# Show the traceplots

rstan::traceplot(Poi_GLMM2, pars = "beta", inc_warmup = TRUE)
rstan::traceplot(Poi_GLMM2, pars= "gamma", inc_warmup = TRUE)
rstan::traceplot(Poi_GLMM2, pars = "alpha", inc_warmup = TRUE)

rstan::traceplot(Poi_GLMM2, pars = "beta", inc_warmup = FALSE)
rstan::traceplot(Poi_GLMM2, pars= "gamma", inc_warmup = FALSE)
rstan::traceplot(Poi_GLMM2, pars = "alpha", inc_warmup = FALSE)

param_GLMM2 <- As.mcmc.list(Poi_GLMM2, pars = c("beta", "gamma", "alpha"))
summary(param_GLMM2)
geweke.diag(param_GLMM2)

# CI coefficients

ci(param_GLMM2, method = "ETI")
ci(param_GLMM2, method = "HDI")

# show coefs

plot_post <- Poi_GLMM2 %>% 
  rstan::extract("beta") %>% 
  as.data.frame() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density(fill = "gray") + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() + 
  geom_vline(aes(xintercept = 0), lty = 2, col = 2) +
  # xlim(c(-1, 12)) +
  theme(legend.position="none")

plot_post <- Poi_GLMM2 %>% 
  rstan::extract("gamma") %>% 
  as.data.frame() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param) + 
  scale_fill_locuszoom() + 
  theme_minimal() + 
  xlim(c(-0.25, 2)) +
  geom_vline(aes(xintercept = 0), lty = 2, col = 2) +
  theme(legend.position="none")

plot_post <- Poi_GLMM2 %>% 
  rstan::extract("alpha") %>% 
  as.data.frame() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param) + 
  scale_fill_locuszoom() + 
  theme_minimal() + 
  geom_vline(aes(xintercept = 0), lty = 2, col = 2) +
  theme(legend.position="none")

# ----------------------------

data_POI <-list(N = length(y), 
                p_fix = 4, 
                y = y, 
                X = cbind(rep(1, 60), X, Z), 
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
rstan::traceplot(Poi_GLM, pars = "beta", inc_warmup = FALSE)

param_GLM <- As.mcmc.list(Poi_GLM, pars = c("beta"))
summary(param_GLM)
geweke.diag(param_GLM)

# TEST --------------

marginal_GLMM2 <- bridge_sampler(Poi_GLMM2)$logml
marginal_GLM <- bridge_sampler(Poi_GLM)$logml

llik_GLMM2 <- extract(Poi_GLMM2, "log_lik")[[1]]
llik_GLM <- extract(Poi_GLM, "log_lik")[[1]]

wrong_marginal_GLMM2 <- log(nrow(llik_GLMM2)) + logsumexp(rowSums(llik_GLMM2))
wrong_marginal_GLM <- log(nrow(llik_GLM)) + logsumexp(rowSums(llik_GLM))

# BAYES FACTOR
exp(wrong_marginal_GLMM2 - wrong_marginal_GLM)
exp(marginal_GLMM2 - marginal_GLM)


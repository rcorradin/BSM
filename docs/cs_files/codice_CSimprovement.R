# -----------------------------
# We need the following libraries
# -----------------------------

library(ggpubr)
library(ggplot2)
library(rstan)
library(coda)

logsumexp <- function(x){
  M <- max(x)
  lse <- M + log(sum(exp(x - M)))
  return(lse)
}

# -----------------------------
# -----------------------------

data <- read.csv("frogs.csv", head = T)
y <- log(data$BoW)
X <- data.frame(intercept = rep(1, length(y)), 
                family = as.numeric(as.factor(data$Family)) - 1,
                logSVL = log(data$SVL), 
                logBrW = log(data$BrW))

cond <- apply(X, 1, function(x) sum(is.na(x)) == 0)
y <- y[cond]
X <- X[cond,]

# -----------------------------
# -----------------------------

models <- expand.grid(c(1), c(0,1), c(0,1), c(0,1))

# -----------------------------
# -----------------------------

my_waic <- function(out_pred_dens){
  LPPD <- mean(log(colMeans(out_pred_dens)))
  p_waic <- mean(apply(out_pred_dens, 2, function(x) var(log(x))))
  return(- 2 * LPPD + 2 * p_waic)
}

my_lpml <- function(out_pred_dens){
  out <- sum(log(1 / colMeans(1 / out_pred_dens)))
  return(out)
}

# -----------------------------
# -----------------------------

WAIC_LPML_out <- matrix(0, ncol = 2, nrow = 8)
out_dens_values <- list()

for(m in 1:nrow(models)){
  if(m == 1){
    X_temp <- matrix(X[, as.logical(models[m,])])
  } else {
    X_temp <- X[, as.logical(models[m,])]
  }
  
  data_temp <- list(n = nrow(X_temp), 
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
                          iter = 1500, 
                          warmup = 500, 
                          seed = 123)
  out_dens_values[[m]] <- exp(extract(temp_regression, pars = c("log_lik"))[[1]])
  WAIC_LPML_out[m,1] <- my_waic(out_dens_values[[m]])
  WAIC_LPML_out[m,2] <- my_lpml(out_dens_values[[m]])
}

print(cbind(WAIC_LPML_out, models))

marginal_best <- log(nrow(out_dens_values[[4]])) - 
  logsumexp(-apply(out_dens_values[[4]], 1, function(z) sum(log(z))))
marginal_saturated <- log(nrow(out_dens_values[[8]])) - 
  logsumexp(-apply(out_dens_values[[8]], 1, function(z) sum(log(z))))
marginal_intercept <- log(nrow(out_dens_values[[1]])) - 
  logsumexp(-apply(out_dens_values[[1]], 1, function(z) sum(log(z))))

BF_best_saturated <- exp(marginal_best - marginal_saturated)
BF_best_intercept <- exp(marginal_best - marginal_intercept)

# -----------------------------
# -----------------------------

data <- read.table("reach.txt", head = T)
data$gender <- data$gender - 1

y <- data$outcome
X <- data.frame(intercept = rep(1, length(y)), 
                data[,1:12])

# -----------------------------
# -----------------------------

c <- 100
k <- 0.05
tau <- k / sqrt(2 * log(c) * c^2 / (c^2 - 1))

data_stan <- list(n = nrow(X), 
                  p = ncol(X), 
                  y = y, 
                  X = X,
                  beta0 = rep(0, ncol(X)), 
                  alpha1 = 1, 
                  alpha2 = 1,
                  tau2 = tau^2, 
                  c2 = c^2)

probit_variable_selection <- stan(file = "probit_SpSl.stan", 
                                  data = data_stan, 
                                  chain = 1,
                                  iter = 12000, 
                                  warmup = 2000, 
                                  seed = 42)

rstan::traceplot(probit_variable_selection, pars = c("beta"))
params <- as.mcmc(extract(probit_variable_selection, pars = c("beta"))[[1]])
summary(params)
geweke.diag(params)

# -----------------------------
# -----------------------------

gamma_matrix <- ifelse(abs(params) > k, 1, 0)

# HPD 

unique_model <- unique(gamma_matrix, MARGIN = 1)
freq <- apply(unique_model, 1, function(b)
  sum(apply(gamma_matrix, MARGIN = 1, function(a) mean(a == b) == 1)))
HPD_model <- unique_model[which.max(freq),]

# MPM

MPM_model <- as.numeric(colMeans(gamma_matrix) > 0.5)

# HS

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

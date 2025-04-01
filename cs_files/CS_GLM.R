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
# Grouse data
# ----------------------------

data <- read.csv("grouse.csv", head = T)
y <- data$TICKS
X <- data.frame(x1 = rep(1, length(y)),  
                x2 = data$HEIGHT) 
Z <- data.frame(z1 = ifelse(data$YEAR == 95, 1, 0),
                z2 = ifelse(data$YEAR == 96, 1, 0), 
                z3 = ifelse(data$YEAR == 97, 1, 0))

# ----------------------------

data_POI <-list(N = length(y), 
                p_fix = 2, 
                p_rand = 3,
                y = y, 
                X = X, 
                Z = Z, 
                Sigma0 = diag(10^2, 2), 
                b0 = rep(0, 2),
                sigma20 = 10^2)

#----- FIT THE MODEL

Poi_GLMM <- stan(file = "Poi_GLMM.stan", 
                 data = data_POI,
                 chains = 1, 
                 iter = 5000, 
                 warmup = 2000, 
                 seed = 42)

# Show the traceplots

rstan::traceplot(Poi_GLMM, pars = "beta", inc_warmup = TRUE)
rstan::traceplot(Poi_GLMM, pars= "gamma", inc_warmup = TRUE)

rstan::traceplot(Poi_GLMM, pars = "beta", inc_warmup = FALSE)
rstan::traceplot(Poi_GLMM, pars= "gamma", inc_warmup = FALSE)

param_GLMM <- As.mcmc.list(Poi_GLMM, pars = c("beta", "gamma"))
summary(param_GLMM)
geweke.diag(param_GLMM)

# CI coefficients

ci(param_GLMM, method = "ETI")
ci(param_GLMM, method = "HDI")

# --------------------------

data <- read.csv("grouse.csv", head = T)
y <- data$TICKS
X <- data.frame(x1 = rep(1, length(y)),  
                x2 = data$HEIGHT) 
Z <- data.frame(z1 = ifelse(data$YEAR == 96, 1, 0), 
                z2 = ifelse(data$YEAR == 97, 1, 0))

# ----------------------------

data_POI <-list(N = length(y), 
                p_fix = 2, 
                p_rand = 2,  
                y = y, 
                X = X, 
                Z = Z, 
                Sigma0 = diag(10^2, 2), 
                b0 = rep(0, 2),
                sigma20 = 10^2)

#----- FIT THE MODEL

Poi_GLMM <- stan(file = "Poi_GLMM.stan", 
                 data = data_POI,
                 chains = 1, 
                 iter = 5000, 
                 warmup = 2000, 
                 seed = 42)

# Show the traceplots

rstan::traceplot(Poi_GLMM, pars = "beta", inc_warmup = TRUE)
rstan::traceplot(Poi_GLMM, pars= "gamma", inc_warmup = TRUE)

rstan::traceplot(Poi_GLMM, pars = "beta", inc_warmup = FALSE)
rstan::traceplot(Poi_GLMM, pars= "gamma", inc_warmup = FALSE)

param_GLMM <- As.mcmc.list(Poi_GLMM, pars = c("beta", "gamma"))
summary(param_GLMM)
geweke.diag(param_GLMM)

# CI coefficients

ci(param_GLMM, method = "ETI")
ci(param_GLMM, method = "HDI")

# PLOT

params <- rstan::extract(Poi_GLMM, c("beta", "gamma"), perm = T)
x_seq <- seq(380, 570, length.out = 100)
jpar <- cbind(params[[1]], rep(0, length(par)), params[[2]])
tres <- array(0, dim = c(nrow(jpar), length(x_seq), 3))
for(i in 1:nrow(jpar)){
  for(j in 1:3){
    tres[i,,j] <- exp(jpar[i,1] + jpar[i,2] * x_seq + jpar[i,2 + j])
  }
}

df_plot <- data.frame(x = rep(x_seq, 3),
                      ymean = as.vector(apply(tres, c(2,3), mean)),
                      group = as.factor(c(rep(95, 100), rep(96, 100), rep(97, 100))))
pdf_plot <- data.frame(x = X[,2], y = y, group = as.factor(data$YEAR))

ggplot() + 
  geom_point(data = pdf_plot, mapping = aes(x = x, y = y, col = group)) + 
  geom_line(data = df_plot, mapping = aes(x = x, y = ymean, col = group)) +
  theme_minimal() + 
  theme(legend.position="none") + 
  ylab("Numbers of ticks on red grouse") + 
  xlab("") 

# ----------------------------

data_POI <-list(N = length(y), 
                p_fix = 2, 
                y = y, 
                X = X, 
                Sigma0 = diag(10^2, 2), 
                b0 = rep(0, 2))

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

llik_GLMM <- extract(Poi_GLMM, "log_lik")[[1]]
marginal_GLMM <- log(nrow(llik_GLMM)) - logsumexp(-apply(llik_GLMM, 1, sum))

llik_GLM <- extract(Poi_GLM, "log_lik")[[1]]
marginal_GLM <- log(nrow(llik_GLM)) - logsumexp(-apply(llik_GLM, 1, sum))

# BAYES FACTOR
exp(marginal_GLMM - marginal_GLM)

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
                sigma20 = 10^2, 
                tau20 = 10^2, 
                lambda20 = 10^2)

#----- FIT THE MODEL

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

ci(param_GLMM, method = "ETI")
ci(param_GLMM, method = "HDI")

# ----------------------------

data_POI <-list(N = length(y), 
                p_fix = 4, 
                y = y, 
                X = cbind(rep(1, 60), X, Z), 
                Sigma0 = diag(10^2, 4), 
                b0 = rep(0, 4))

#----- FIT THE MODEL

Poi_GLM2 <- stan(file = "Poi_GLM.stan", 
                data = data_POI,
                chains = 1, 
                iter = 5000, 
                warmup = 2000, 
                seed = 42)

# Show the traceplots

rstan::traceplot(Poi_GLM2, pars = "beta", inc_warmup = TRUE)
rstan::traceplot(Poi_GLM2, pars = "beta", inc_warmup = FALSE)

param_GLM <- As.mcmc.list(Poi_GLM2, pars = c("beta"))
summary(param_GLM)
geweke.diag(param_GLM)

# TEST --------------

llik_GLMM2 <- extract(Poi_GLMM2, "log_lik")[[1]]
marginal_GLMM2 <- logsumexp(rowSums(llik_GLMM2))

llik_GLM2 <- extract(Poi_GLM2, "log_lik")[[1]]
marginal_GLM <- logsumexp(rowSums(llik_GLM2))

# BAYES FACTOR
exp(marginal_GLMM - marginal_GLM)

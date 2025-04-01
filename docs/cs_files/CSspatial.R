# -----------------------------
# We need the following libraries
# -----------------------------

library(rstan)
library(coda)
library(ggplot2)
library(mvtnorm)

logsumexp <- function(x){
  M <- max(x)
  lse <- M + log(sum(exp(x - M)))
  return(lse)
}

# -----------------------------
# -----------------------------


load("river.rda")

x0 <- X[1,]
s0 <- unlist(S[1,])
y0 <- y[1]

y <- y[-1]
X <- X[-1,]
S <- S[-1,]

n <- length(y)

# -----------------------------

data_STAN <-list(n = n, 
                 p = 3,  
                 y = y, 
                 X = X, 
                 S = S,
                 x0 = x0, 
                 s0 = s0,
                 ord = 2,
                 
                 beta0 = rep(0, 3), 
                 Lambda0 = diag(10^3, 3), 
                 asigma = 2.1, 
                 bsigma = 10, 
                 atau = 2.1, 
                 btau = 10, 
                 aphi = 2.1, 
                 bphi = 10)

out_point_STAN <- stan(file = "point_referenced.stan", 
                       data = data_STAN,
                       chains = 1, 
                       iter = 7000, 
                       warmup = 2000, 
                       seed = 123)

# -----------------------------
# looking at the estimate

plot(out_point_STAN, pars = c("sigma2", "tau2", "phi"))
plot(out_point_STAN, pars = c("beta"))

param1 <- As.mcmc.list(out_point_STAN, pars = c("sigma2", "tau2", "phi"))
param2 <- As.mcmc.list(out_point_STAN, pars = c("beta"))

summary(param1)
summary(param2)

geweke.diag(param1)
geweke.diag(param2)

plot(param1)
plot(param2)

# -----------------------------

params_pred <- rstan::extract(out_point_STAN, pars = c("new_pred_mean_var"))[[1]]
pred_sample <- apply(params_pred, 1, function(x) rnorm(1, x[1], sqrt(x[2])))
ggplot(data.frame(x = pred_sample)) + 
  geom_histogram(aes(x = x, y = after_stat(density)), alpha = 0.3, col = 1) + 
  theme_bw() +
  xlab("") +
  ylab("") +
  geom_vline(aes(xintercept = y0), col = 2, lwd = 0.8, lty = 2)

# -----------------------------
# no spatial effect

data_STAN <-list(n = n, 
                 p = 3,  
                 y = y, 
                 X = X, 
                 x0 = x0,
                 
                 beta0 = rep(0, 3), 
                 Lambda0 = diag(10^3, 3), 
                 atau = 2.1, 
                 btau = 10)

out_point_null_STAN <- stan(file = "point_referenced_null.stan", 
                            data = data_STAN,
                            chains = 1, 
                            iter = 7000, 
                            warmup = 2000, 
                            seed = 123)

# -----------------------------
# test against no spatial association

llik1 <- rstan::extract(out_point_STAN, pars = "llik")[[1]]
log_marginal_M1 <- log(length(llik1)) - logsumexp(-llik1)


llik2 <- rstan::extract(out_point_null_STAN, pars = "llik")[[1]]
log_marginal_M2 <- log(length(llik2)) - logsumexp(-llik2)

exp(log_marginal_M1 - log_marginal_M2)

# -----------------------------
# -----------------------------

load("sids.rda")

data_STAN <-list(n = length(y), 
                 p = ncol(X),  
                 y = y, 
                 X = X, 
                 W = W,
                 O = O, 
                 
                 beta0 = array(0, dim = ncol(X)), 
                 Lambda0 = diag(10^3, ncol(X)), 
                 atau = 2, 
                 btau = 10, 
                 arho = 1, 
                 brho = 1)

out_areal_STAN <- stan(file = "areal.stan", 
                       data = data_STAN,
                       chains = 1, 
                       iter = 7000, 
                       warmup = 2000, 
                       seed = 123)

plot(out_areal_STAN, pars = c("tau2", "rho"))
plot(out_areal_STAN, pars = c("beta"))

param1 <- As.mcmc.list(out_areal_STAN, pars = c("tau2", "rho"))
param2 <- As.mcmc.list(out_areal_STAN, pars = c("beta"))

summary(param1)
summary(param2)

geweke.diag(param1)
geweke.diag(param2)

plot(param1)
plot(param2)

# random eff

library(sf)
        
nc.sids <- sf::st_read(system.file("shapes/sids.gpkg", package="spData")[1])
row.names(nc.sids) <- as.character(nc.sids$FIPS)
rn <- row.names(nc.sids)

ggplot(nc.sids) + 
  geom_sf(aes(fill = SID74)) + 
  scale_fill_gradient2(low="blue", mid="white",
                       high="red", space ="Lab") + 
  theme_bw() +
  theme(legend.position = "null")

r_eff <- As.mcmc.list(out_areal_STAN, pars = c("phi"))
ggplot(nc.sids) + 
  geom_sf(aes(fill = colMeans(r_eff[[1]]))) + 
  scale_fill_gradient2(low="blue", mid="white",
                       high="red", space ="Lab") + 
  theme_bw() +
  theme(legend.position = "null")

# -----------------------------
# null model

data_STAN <-list(n = length(y), 
                 p = ncol(X) + 1,  
                 y = y, 
                 X = cbind(rep(1, length(y)), X),
                 O = O, 
                 
                 beta0 = array(0, dim = ncol(X) + 1), 
                 Lambda0 = diag(10^3, ncol(X) + 1))

out_areal_null_STAN <- stan(file = "areal_null.stan", 
                            data = data_STAN,
                            chains = 1, 
                            iter = 7000, 
                            warmup = 2000, 
                            seed = 123)

# -----------------------------
# test against no spatial association

llik1 <- rstan::extract(out_areal_STAN, pars = "llik")[[1]]
log_marginal_M1 <- log(nrow(llik1)) - logsumexp(-rowSums(llik1))


llik2 <- rstan::extract(out_areal_null_STAN, pars = "llik")[[1]]
log_marginal_M2 <- log(nrow(llik2)) - logsumexp(-rowSums(llik2))

exp(log_marginal_M1 - log_marginal_M2)

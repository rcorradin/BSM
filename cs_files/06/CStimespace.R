# -----------------------------
# We need the following libraries
# -----------------------------

library(rstan)
library(coda)
library(ggplot2)
library(mvtnorm)
library(bridgesampling)
library(sf)
library(spdep)
library(matrixsampling)
library(reshape2)

my_lpml <- function(out_pred_dens){
  out <- sum(log(1 / colMeans(1 / out_pred_dens)))
  return(out)
}

my_waic <- function(out_pred_dens){
  LPPD <-   mean(log(colMeans(out_pred_dens)))
  p_waic <- mean(apply(out_pred_dens, 2, function(x) var(log(x))))
  return(-2 * LPPD + 2 * p_waic)
}

# -----------------------------

time_series <- read.csv("time_series.csv")
amazon <- time_series[,1]

# -----------------------------
# TIME SERIES - univariate
# -----------------------------

sample_univ <- function(y, y0, niter, nburn, m0, S0, a0, b0){
  out_coef <- matrix(0, nrow = niter, ncol = length(m0))
  out_s2 <- rep(0, niter)
  out_lik <- matrix(0, nrow = niter, ncol = length(y))
  
  p <- length(y0)
  n <- length(y)
  
  y_temp <- c(y0, y)
  Z <- matrix(1, nrow = length(y), ncol = length(y0) + 1)
  for(i in 1:length(y)){
    Z[i,-1] <- y_temp[(i + p - 1):i]
  }
  
  pb <- txtProgressBar(min = 1, max = niter, style = 3)
  for(iter in 1:niter){
    Sn = solve(solve(S0) + t(Z) %*% Z)
    mn = (Sn %*% (solve(S0) %*% m0 + t(Z) %*% y))
    an = a0 + n / 2
    bn = b0 + 0.5 *(t(y) %*% y - t(mn) %*% solve(Sn) %*% mn + t(m0) %*% solve(S0) %*% m0)
    
    out_s2[iter] <- 1 / rgamma(1, shape = an, rate = bn)
    out_coef[iter,] <- rmvnorm(1, mean = mn, sigma = Sn)[1,]
    for(i in 1:length(y)){
      out_lik[iter, i] <- dnorm(y[i], Z[i,] %*% out_coef[iter,], sqrt(out_s2[iter]))
    }
    setTxtProgressBar(pb, iter)
  }
  close(pb)
  return(list(out_coef[-c(1:nburn),], out_s2[-c(1:nburn)], out_lik[-c(1:nburn),]))
}

# -----------------------------

amazons <- scale(amazon[-1])
model_univ_3 <- sample_univ(amazons[-c(1:7)], amazons[5:7], 7000, 2000, 
                            rep(0, 4), diag(1, 4), 2, 1)
out_mcmc <- as.mcmc(model_univ_3[[1]])
plot(out_mcmc)
summary(out_mcmc)
geweke.diag(out_mcmc)

# -----------------------------

model_univ_7 <- sample_univ(amazons[-c(1:7)], amazons[1:7], 7000, 2000, 
                            rep(0, 8), diag(1, 8), 2, 1)
out_mcmc <- as.mcmc(model_univ_7[[1]])
plot(out_mcmc)
summary(out_mcmc)
geweke.diag(out_mcmc)

# -----------------------------

my_lpml(model_univ_3[[3]])
my_lpml(model_univ_7[[3]])

my_waic(model_univ_3[[3]])
my_waic(model_univ_7[[3]])

# -----------------------------
# TIME SERIES - multivariate
# -----------------------------

sample_multiv <- function(Y, niter, nburn, m0, k0, nu0, Lambda0, M0, Psi0){
  out_drift  <- matrix(0, nrow = niter, ncol = ncol(Y))
  out_effect <- array(0, dim = c(ncol(Y), ncol(Y), niter))
  out_Sigma  <- array(0, dim = c(ncol(Y), ncol(Y), niter))
  out_lik    <- matrix(0, nrow = niter, ncol = nrow(Y))
  
  d <- ncol(Y)
  n <- nrow(Y)
  
  out_effect[,,1] <- matrix(rnorm(d * d, 0, 0.1), ncol = d)
  
  Z <- as.matrix(Y[-n,])
  Ytemp <- as.matrix(Y[-1,])
  
  pb <- txtProgressBar(min = 1, max = niter, style = 3)
  for(iter in 1:niter){
    X <- Ytemp - Z %*% out_effect[,,1]
    nun <- nu0 + n + d - 1
    kn <- k0 + n - 1
    Lambdan <- Lambda0 + t(out_effect[,,1] - M0) %*% solve(Psi0) %*% (out_effect[,,1] - M0) +
      t(X) %*% X + k0 * (n - 1) / kn * (m0 - colMeans(X)) %*% t(m0 - colMeans(X))
    mn <- (k0 * m0 + colSums(Z %*% out_effect[,,1])) / kn
    
    out_Sigma[,,iter] <- rinvwishart(1, nun, Lambdan)[,,1]
    out_drift[iter,] <- rmvnorm(1, mn, out_Sigma[,,iter])[1,]
    
    U <- Ytemp - rep(1, n-1) %*% t(out_drift[iter,])
    Psin <- solve(Psi0 + t(Z) %*% Z)
    Mn <- Psin %*% (solve(Psi0) %*% M0 + t(Z) %*% U)
    out_effect[,,iter] <- rmatrixnormal(n = 1, M = Mn, V = round(Psin, digits = 10), U = out_Sigma[,,iter])[,,1]
    
    for(i in 1:length(y)){
      out_lik[iter, i] <- dmvnorm(x = Y[i,], mean = out_drift[iter,] + Z[i,] %*% out_effect[,,iter], 
                                  sigma = out_Sigma[,,iter])
    }
    setTxtProgressBar(pb, iter)
  }
  close(pb)
  return(list(out_drift[-c(1:nburn),], out_effect[,,-c(1:nburn)],
              out_Sigma[,,-c(1:nburn)], out_lik[-c(1:nburn),]))
}

# -----------------------------

time_seriess <- scale(time_series)
model_multiv_1 <- sample_multiv(time_seriess, 7000, 2000, rep(0, 3), 
                                1, 5, diag(1, 3), matrix(0, 3, 3), diag(1, 3))

colMeans(model_multiv_1[[1]] > 0)
apply(model_multiv_1[[2]] > 0, c(1,2), mean)
coefs <- apply(model_multiv_1[[2]], c(1,2), mean)
rownames(coefs) <- colnames(coefs) <- colnames(time_series)
diag(coefs) <- rep(0, 3)

matrix_plot <- melt(abs(coefs))
plot_heat <- ggplot(matrix_plot) +
  geom_tile(aes(x = Var2, y = Var1, fill = value), na.rm = TRUE, size = 0.0,show.legend = FALSE) +
  scale_fill_gradient(
    low = "white",
    high = "black"
  ) +
  labs(x = "", y = "", fill = "") +
  theme_linedraw() +
  coord_cartesian(xlim = c(0, ncol(coefs) + 1), ylim = c(0, ncol(coefs) + 1), expand = FALSE) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(size=0.1))

# -----------------------------
# POINT REFERENCED - river data
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

log_marginal_M1 <- bridge_sampler(out_point_STAN)$logml
log_marginal_M2 <- bridge_sampler(out_point_null_STAN)$logml
exp(log_marginal_M1 - log_marginal_M2)

# -----------------------------
# AREAL DATA - sids
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

log_marginal_M1 <- bridge_sampler(out_areal_STAN)$logml
log_marginal_M2 <- bridge_sampler(out_areal_null_STAN)$logml
exp(log_marginal_M1 - log_marginal_M2)

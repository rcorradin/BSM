library(mvtnorm)
library(ggplot2)
library(GGally)
library(extraDistr)
library(microbenchmark)

# ------------------
# SIMULATE THE DATA
# ------------------

n <- 150
true_clust <- sample(1:3, size = n, replace = T, prob = c(0.25, 0.25, 0.5))
means <- c(-3, 0, 3)
vars <- c(0.5, 0.5, 0.5)
Y <- rnorm(n, mean = means[true_clust], sd = sqrt(vars[true_clust]))

# ------------------
# Functions
# ------------------

gibbs_univariate_full <- function(Y, k, niter, nburn, m0, k0, a0, b0, alpha){
  
  n <- length(Y)
  
  # clusters output matrix 
  out_cluster <- matrix(0, nrow = niter, ncol = n)
  out_weights <- matrix(0, nrow = niter, ncol = k)
  out_params  <- array(0, dim = c(k, 2, niter))
  
  # group-specific parameters (scales and locations), weights and partition
  params <- matrix(nrow = k, ncol = 2)
  for(j in 1:k){
    params[j,2] <- 1 / rgamma(1, shape = a0, rate = b0)
    params[j,1] <- rnorm(1, mean = m0, sd = sqrt(params[j,2] / k0))
  }
  
  w <- rdirichlet(1, rep(alpha, k))[1,]
  allocations <- sample(1:k, size = n, replace = T)
  temp_prob <- rep(0, k)
  pb <- txtProgressBar(min = 1, max = niter, style = 3)
  
  # loop over the number of iterations
  for(iter in 1:niter){
    
    # ------------------
    # update allocations 
    # ------------------
    
    for(i in 1:n){
      for(j in 1:k){
        temp_prob[j] <- w[j] * dnorm(Y[i], params[j,1], sqrt(params[j,2]))
      }
      allocations[i] <- sample(1:k, size = 1, prob = temp_prob)
    }  
    
    # ------------------    
    # update parameters
    # ------------------
    
    for(j in 1:k){
      Y_temp <- Y[allocations == j]
      n_temp <- sum(allocations == j)
      if(n_temp > 0){
        y_bar_temp <- mean(Y_temp)  
      } else {
        y_bar_temp <- 0
      }
      
      an = a0 + n_temp/2
      kn = k0 + n_temp
      
      mn = (k0 * m0 + n_temp * y_bar_temp) / kn
      bn = b0 + 0.5 * (sum((Y_temp - y_bar_temp)^2) + 
                         k0 * n_temp / kn * (y_bar_temp- m0)^2)
      
      params[j,2] <- 1 / rgamma(1, shape = an, rate = bn)
      params[j,1] <- rnorm(1, mean = mn, sd = sqrt(params[j,2] / kn))
    }
    
    # ------------------
    # update the weights
    # ------------------
      
    alphan <- rep(alpha, k) + sapply(1:k, function(x) sum(allocations == x))
    w <- rdirichlet(1, alphan)[1,]
    
    # save the results 
    out_cluster[iter, ] <- allocations
    out_weights[iter, ] <- w
    out_params[,,iter] <- params
    setTxtProgressBar(pb, iter)
  }
  close(pb)
  return(list(out_cluster[-c(1:nburn),],
              out_weights[-c(1:nburn),],
              out_params[,,-c(1:nburn)]))
}

# ------------------
# TEST THE FUNCTION
# ------------------

test1 <- gibbs_univariate_full(Y = Y, k = 6, niter = 2500, nburn = 500, m0 = 0, 
                               k0 = 0.1, a0 = 2, b0 = 1, alpha = 1)
sampsize <- ncol(test1[[1]])
entropy <- apply(test1[[1]], 1, function(x) -sum(table(x)/sampsize * log(table(x)/sampsize)))
active_comp <- apply(test1[[1]], 1, function(x) length(table(x)))

ggplot(data.frame(x = rep(1:length(entropy), 2), 
                  y = c(entropy, active_comp), 
                  z = factor(rep(c("entropy", "active"), each = length(entropy))))) +
  geom_line(aes(x = x, y = y)) + 
  theme_bw() + 
  facet_wrap(.~z, scales = "free", ncol = 1)

# ------------------
# POINT ESTIMATE
# ------------------

point_Binder <- function(partition_matrix){
  
  n <- ncol(partition_matrix)
  PSM <- matrix(0, ncol = n, nrow = n)
  dist_Binder <- c()
  
  # compute the PSM
  for(i in 1:n){
    for(j in 1:i){
      PSM[i,j] <- PSM[j,i] <- mean(partition_matrix[,i] == partition_matrix[,j])
    }
  }
  
  # find the best matrix
  for(i in 1:nrow(partition_matrix)){
    temp_mat <- matrix(as.numeric(sapply(partition_matrix[i,], 
                                         function(z) z == partition_matrix[i,])), ncol = n)
    dist_Binder[i] <- sum((temp_mat - PSM)^2)
  }
  point_estimate <- partition_matrix[which.min(dist_Binder),]
  
  return(as.numeric(as.factor(point_estimate)))
}

est_part <- point_Binder(test1[[1]])
table(true_clust, est_part)


# ------------------
# PLOT THE DATA
# ------------------

pl <- ggplot(data.frame(y = Y)) +
  geom_histogram(aes(x = y, y = ..density..), bins = 50, alpha = 0.25, col = 1) +
  theme_bw() 

xgrid <- seq(-7, +7, by = 0.1)
wei <- test1[[2]]
params <- test1[[3]]
temp_dens_mean <- rep(0, length(xgrid))

for(i in 1:nrow(wei)){
    temp_dens <- sapply(xgrid, function(x) sum(wei[i,] * 
                         dnorm(x, mean = params[,1,i], sd = sqrt(params[,2,i]))))
    temp_dens_mean <- temp_dens_mean + temp_dens
    if(i < 250){
      pl <- pl + geom_line(data = data.frame(x = xgrid, y = temp_dens), 
                         aes(x = x, y = y), lwd = 0.1) 
    }
}

temp_dens_mean <- temp_dens_mean / nrow(wei)
pl <- pl + geom_line(data = data.frame(x = xgrid, y = temp_dens_mean), 
                   aes(x = x, y = y), col = 2, lwd = 1)
pl <- pl + geom_point(data = data.frame(y = rep(0, length(Y)), 
                          x = Y, col = factor(est_part)), 
               aes(x = x, y = y, col = col)) + 
  theme(legend.position = "null")
pl


# ------------------
# TEST THE FUNCTION 2
# ------------------

test2 <- gibbs_univariate_full(Y = Y, k = 6, niter = 2500, nburn = 500, m0 = 0, 
                               k0 = 10, a0 = 2, b0 = 1, alpha = 1)
sampsize <- ncol(test2[[1]])
entropy <- apply(test2[[1]], 1, function(x) -sum(table(x)/sampsize * log(table(x)/sampsize)))
active_comp <- apply(test2[[1]], 1, function(x) length(table(x)))

ggplot(data.frame(x = rep(1:length(entropy), 2), 
                  y = c(entropy, active_comp), 
                  z = factor(rep(c("entropy", "active"), each = length(entropy))))) +
  geom_line(aes(x = x, y = y)) + 
  theme_bw() + 
  facet_wrap(.~z, scales = "free", ncol = 1)

est_part <- point_Binder(test2[[1]])
table(true_clust, est_part)

pl <- ggplot(data.frame(y = Y)) +
  geom_histogram(aes(x = y, y = ..density..), bins = 50, alpha = 0.25, col = 1) +
  theme_bw() 

xgrid <- seq(-7, +7, by = 0.1)
wei <- test2[[2]]
params <- test2[[3]]
temp_dens_mean <- rep(0, length(xgrid))

for(i in 1:nrow(wei)){
  temp_dens <- sapply(xgrid, function(x) sum(wei[i,] * 
                                               dnorm(x, mean = params[,1,i], sd = sqrt(params[,2,i]))))
  temp_dens_mean <- temp_dens_mean + temp_dens
  if(i < 250){
    pl <- pl + geom_line(data = data.frame(x = xgrid, y = temp_dens), 
                         aes(x = x, y = y), lwd = 0.1) 
  }
}

temp_dens_mean <- temp_dens_mean / nrow(wei)
pl <- pl + geom_line(data = data.frame(x = xgrid, y = temp_dens_mean), 
                     aes(x = x, y = y), col = 2, lwd = 1)
pl <- pl + geom_point(data = data.frame(y = rep(0, length(Y)), 
                                        x = Y, col = factor(est_part)), 
                      aes(x = x, y = y, col = col)) + 
  theme(legend.position = "null")
pl

# ------------------
# CHECK THE PERFORMANCES
# ------------------

# DIFFERENT SAMPLE SIZES

ns <- c(150, 250, 1000)  
means <- c(-3, 0, 3)
vars <- c(0.5, 0.5, 0.5)
bench <- list()
out <- list()
binder <- c()

for(i in 1:3){
  true_clust <- sample(1:3, size = ns[i], replace = T, prob = c(0.25, 0.25, 0.5))
  Y <- rnorm(ns[i], mean = means[true_clust], sd = sqrt(vars[true_clust]))
  bench[[i]] <- microbenchmark(out[[i]] <- gibbs_univariate_full(Y = Y, k = 6, niter = 2500, nburn = 500, m0 = 0, 
                                                k0 = 10, a0 = 2, b0 = 1, alpha = 1), times = 1)
  est_clust <- point_Binder(out[[i]][[1]])
  binder[i] <- sum((table(true_clust) / ns[i])^2) + sum((table(est_clust) / ns[i])^2) - 2 * 
    sum((table(true_clust, est_clust) / ns[i])^2)
}

# DIFFERENT DISPERSION

gamma <- c(0.5, 1, 2)  
means <- c(-3, 0, 3)
vars <- c(0.5, 0.5, 0.5)
binder <- c()

for(i in 1:3){
  true_clust <- sample(1:3, size = 150, replace = T, prob = c(0.25, 0.25, 0.5))
  Y <- rnorm(150, mean = gamma[i] * means[true_clust], sd = sqrt(vars[true_clust]))
  est_model <- gibbs_univariate_full(Y = Y, k = 6, niter = 2500, nburn = 500, m0 = 0, 
                        k0 = 10, a0 = 2, b0 = 1, alpha = 1)
  est_clust <- point_Binder(est_model[[1]])
  binder[i] <- sum((table(true_clust) / ns[i])^2) + sum((table(est_clust) / ns[i])^2) - 2 * 
    sum((table(true_clust, est_clust) / ns[i])^2)
}

# ------------------
# STARS DATA
# ------------------

# First we write a function to perform multivariate clustering using a 
# Gaussian kernel function, mixing with respect to both locations and scales

gibbs_multivariate <- function(Y, k, niter, nburn, m0, k0, nu0, Lambda0, alpha){
  
  n <- nrow(Y)
  p <- ncol(Y)
  
  # clusters output matrix 
  out_cluster <- matrix(0, nrow = niter, ncol = n)
  
  # group-specific parameters (scales and locations), weights and partition
  covariances <- array(0, dim = c(p, p, k))
  locations <- matrix(0, ncol = p, nrow = k)
  for(j in 1:k){
    covariances[,,j] <- solve(rWishart(n = 1,df = nu0, Sigma = solve(Lambda0))[,,1])
    locations[j,] <- rmvnorm(1, mean = m0, sigma = Lambda0 / ((nu0 - p - 1) * k0))
  }
  
  w <- rdirichlet(1, rep(alpha, k))[1,]
  allocations <- sample(1:k, size = n, replace = T)
  temp_prob <- rep(0, k)
  
  pb <- txtProgressBar(min = 1, max = niter, style = 3)
  
  # loop over the number of iterations
  for(iter in 1:niter){
    
    # update allocations 
    
    for(i in 1:n){
      for(j in 1:k){
        temp_prob[j] <- w[j] * dmvnorm(Y[i,], locations[j,], covariances[,,j])
      }
      allocations[i] <- sample(1:k, size = 1, prob = temp_prob)
    }  
    
    # update parameters
    
    for(j in 1:k){
      Y_temp <- Y[allocations == j,]
      n_temp <- sum(allocations == j)
      
      if(n_temp > 0){
        nun <- nu0 + n_temp
        kn <- k0 + n_temp
        
        if(n_temp > 1){
          y_bar_temp <- colMeans(Y_temp)
        } else {
          y_bar_temp <- Y_temp
        }
        
        y_bar_mat <- matrix(rep(y_bar_temp, n_temp), nrow = n_temp, byrow = T)
        Lambdan <- Lambda0 + t(Y_temp - y_bar_mat) %*% (Y_temp - y_bar_mat) + 
          (k0 * n_temp) / kn * (y_bar_temp - m0) %*% t(y_bar_temp - m0)
        mn <- (k0 * m0 + n_temp * y_bar_temp) / (k0 + n_temp)
        covariances[,,j] <- solve(rWishart(n = 1,df = nun, Sigma = solve(Lambdan))[,,1])
        locations[j,] <- rmvnorm(1, mn, covariances[,,j] / kn)
      } else {
        covariances[,,j] <- solve(rWishart(n = 1,df = nu0, Sigma = solve(Lambda0))[,,1])
        locations[j,] <- rmvnorm(1, mean = m0, sigma = covariances[,,j] / k0)
      }
    }
    
    # update the weights
    
    alphan <- rep(alpha, k) + sapply(1:k, function(x) sum(allocations == x))
    w <- rdirichlet(1, alphan)[1,]
    
    # save the results 
    out_cluster[iter, ] <- allocations
    setTxtProgressBar(pb, iter)
  }
  close(pb)
  return(out_cluster[-c(1:nburn),])
}

# ------------------
# TEST THE FUCNTION
# ------------------

Y <- rbind(rmvnorm(50, c(-3, -3), diag(1, 2)), rmvnorm(25, c(3, 3), diag(1, 2)))
set.seed(123)
out <- gibbs_multivariate(Y = Y, niter = 2500, nburn = 500, k = 10, m0 = c(0,0), 
                          k0 = 0.1, nu0 = 5, Lambda0 = diag(3, 2), alpha = 0.5)

# ---------------

part_estimate <- point_Binder(out)
data_plt <- data.frame(x = Y[,1], y = Y[,2], group = as.factor(part_estimate))
ggplot(data_plt) + 
  geom_point(aes(x = x, y = y, color = group)) + 
  theme_bw() + 
  theme(legend.position = "null")

# ------------------
# RUN EVERYTHING ON THE DATA
# ------------------

data0 <- read.table("globular_cluster.txt", header = T)
data <- as.matrix(scale(data0))

out <- gibbs_multivariate(Y = data, niter = 2500, nburn = 500, k = 10, m0 = rep(0,4), 
                          k0 = 0.1, nu0 = 7, Lambda0 = diag(3, 4), alpha = 0.5)
part_estimate <- point_Binder(out)

ggpairs(data0, aes(color = factor(part_estimate)), alpha = 0.5, 
        upper = list(continuous = "points"), diag = list(continuous = "blankDiag")) + 
  theme_bw()

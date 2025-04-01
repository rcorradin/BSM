library(mvtnorm)
library(LaplacesDemon)
library(ggplot2)
library(GGally)

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
  
  return(point_estimate)
}

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

plot(data0, col = as.numeric(as.factor(part_estimate)), pch = 20)

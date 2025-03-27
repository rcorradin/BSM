# -----------------------------
# We need the following libraries
# -----------------------------

library(coda)
library(LaplacesDemon)
library(ggplot2)

# -----------------------------
# data from the example
# -----------------------------

set.seed(42)
y <- c(rnorm(50, -3, 1), rnorm(25, 3, 1))

# plot the data generating process and distribution

y_seq <- seq(-6, 6, by = 0.1)
dens <- sapply(y_seq, function(x) 0.67 * dnorm(x, -3, 1) + 0.33 * dnorm(x, 3, 1))

# with standard plots
plot(x = y_seq, y = dens, type = "l", xlab = "y", ylab = "dens")
points(x = y[1:50], y = rep(0, 50), col = 2, pch = 20)
points(x = y[51:75], y = rep(0, 25), col = 3, pch = 20)

#  with ggplot
ggplot(data.frame(x = y_seq, y = dens)) + 
  geom_line(aes(x = x, y = y)) + 
  theme_bw() + 
  xlab("y") + 
  ylab("dens") + 
  geom_point(data = data.frame(x = y, group = factor(c(rep(1, 50), rep(2, 25)))), 
             aes(x = x, y = 0, colour = group)) +
  theme(legend.position = "null")

# -----------------------------
# implementing the gibbs sampler for univariate data
# with a simmetric prior for the weights
# -----------------------------

gibbs_univariate <- function(y, k, niter, nburn, m0, k0, a0, b0, alpha){
  
  n <- length(y)
  
  # clusters output matrix 
  out_cluster <- matrix(0, nrow = niter, ncol = n)
  
  # group-specific parameters (scales and locations), weights and partition
  variances <- 1 / rgamma(k, shape = a0, scale = 1 / b0)
  locations <- rnorm(k, m0, sqrt(b0 / ((a0 - 1) * k0)))
  w <- rdirichlet(1, rep(alpha, k))[1,]
  allocations <- sample(1:k, size = n, replace = T)
  temp_prob <- rep(0, k)
  
  pb <- txtProgressBar(min = 1, max = niter, style = 3)
  
  # loop over the number of iterations
  for(iter in 1:niter){
    
    # update allocations 
    
    for(i in 1:n){
      for(j in 1:k){
        temp_prob[j] <- w[j] * dnorm(y[i], locations[j], sqrt(variances[j]))
      }
      allocations[i] <- sample(1:k, size = 1, prob = temp_prob)
    }  
    
    # update parameters
    
    for(j in 1:k){
      y_temp <- y[allocations == j]
      n_temp <- length(y_temp)
      
      if(n_temp > 0){
        an <- a0 + n_temp / 2
        kn <- k0 + n_temp
        bn <- b0 + 0.5 * (sum((y_temp - mean(y_temp))^2) + 
                            k0 * n_temp / kn * (mean(y_temp) - m0)^2)
        mn <- (k0 * m0 + sum(y_temp)) / (k0 + n_temp)
        variances[j] <- 1 / rgamma(1, shape = an, scale = 1 / bn)
        locations[j] <- rnorm(1, mn, sqrt(variances[j] / kn))
      } else {
        variances[j] <- 1 / rgamma(1, shape = a0, scale = 1 / b0)
        locations[j] <- rnorm(1, m0, sqrt(variances[j] / k0))
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

out_gibbs_clustering <- gibbs_univariate(y = y, k = 10, niter = 1500, nburn = 500, 
                                         m0 = 0, k0 = 0.1, a0 = 2.01, b0 = 1.01, alpha = 0.1)

# -----------------------------
# check the chain
# -----------------------------

n <- ncol(out_gibbs_clustering)
entropies <- apply(out_gibbs_clustering, 1, function(x) 
  sum(table(x) / n * log(table(x) / n)))
mcmc_entropies <- as.mcmc(entropies)
summary(mcmc_entropies)
geweke.diag(mcmc_entropies)

# -----------------------------
# use the Binder loss function and provide a point estimate
# of the latent partition in the previous example
# -----------------------------

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

point_Binder(out_gibbs_clustering)

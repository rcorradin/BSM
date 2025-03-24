data{
  int<lower = 0> n;
  int<lower = 0> p;
  
  vector[n] y;
  matrix[n,p] X;
  
  vector[p] beta0;
  matrix[p,p] Sigma0;
  real a0; 
  real b0;
}

parameters{
  real<lower = 0> sigma2; 
  vector[p] beta;
}

transformed parameters{
  vector[n] mu;
  for(i in 1:n){
    mu[i] = row(X, i) * beta;
  }
}

model{
  beta ~ multi_normal(beta0, sigma2 * Sigma0);
  sigma2 ~ inv_gamma(a0, b0);
  y ~ normal(mu, pow(sigma2, 0.5));
}

generated quantities{
  vector[n] log_lik;
  for(i in 1:n){
    log_lik[i] = normal_lpdf(y[i] | mu[i], pow(sigma2, 0.5));
  }
}
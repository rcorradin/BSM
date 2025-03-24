data{
  int<lower = 0> n;
  int<lower = 0> p;
  int y[n];
  matrix[n,p] X;
  
  vector[p] beta0;
  real<lower = 0> c2;
  real<lower = 0> tau2;
  real<lower = 0> alpha1;
  real<lower = 0> alpha2;
}

parameters{
  vector[p] beta;
  real<lower = 0, upper = 1> probs[p];
}

transformed parameters{
  vector[n] mu; 
  for(i in 1:n){
    mu[i] = Phi(row(X, i) * beta);
  }
}

model{
  // likelihood 
  y ~ bernoulli(mu);
  // prior
  for(j in 1:p){
    probs[j] ~ beta(alpha1, alpha2);
    target += log_mix(probs[j], 
              normal_lpdf(beta[j]| beta0[j], pow(c2 * tau2, 0.5)),
              normal_lpdf(beta[j]| beta0[j], pow(tau2, 0.5)));
  }
}

















data
{
	int<lower = 0> n;       
	int<lower = 0> p;     
	int<lower = 0> y[n];
	matrix[n,p] X;
	matrix[n,n] W;
	vector[n] O;
	
	vector[p] beta0;
	matrix[p,p] Lambda0;
	real atau;
	real btau;
	real arho; 
	real brho;
}

transformed data
{
  vector[n] w_plus;
  for(i in 1:n){
    w_plus[i] = sum(row(W, i));
  }
}

parameters
{
	vector[n] phi;
	vector[p] beta;
	real<lower = 0> tau2; 
	real<lower = 0, upper = 1> rho;
}

transformed parameters 
{
	vector[n] mu;
	for(i in 1:n){
	  mu[i] = exp(row(X,i) * beta + log(O[i]) + phi[i]);
	}
}

model
{
	beta ~ multi_normal(beta0, Lambda0);
	tau2 ~ inv_gamma(atau, btau);
	rho ~ beta(arho, brho);
	
	for(i in 1:n){
	  phi[i] ~ normal(rho * (row(W,i) * phi) / w_plus[i], tau2 / w_plus[i]);
	}
	
  y ~ poisson(mu);
}

generated quantities 
{
  vector[n] llik;
  for(i in 1:n){
    llik[i] = poisson_lpmf(y[i] | mu[i]);
  }
}

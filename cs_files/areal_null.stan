data
{
	int<lower = 0> n;       
	int<lower = 0> p;     
	int<lower = 0> y[n];
	matrix[n,p] X;
	vector[n] O;
	
	vector[p] beta0;
	matrix[p,p] Lambda0;
}

parameters
{
	vector[p] beta;
}

transformed parameters 
{
	vector[n] mu;
	for(i in 1:n){
	  mu[i] = exp(row(X,i) * beta + log(O[i]));
	}
}

model
{
	beta ~ multi_normal(beta0, Lambda0);
  y ~ poisson(mu);
}

generated quantities 
{
  vector[n] llik;
  for(i in 1:n){
    llik[i] = poisson_lpmf(y[i] | mu[i]);
  }
}
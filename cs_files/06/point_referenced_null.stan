data
{
	int<lower = 0> n;       
	int<lower = 0> p;     
	vector[n] y;
	matrix[n,p] X;
	
	vector[p] beta0;
	matrix[p,p] Lambda0;
	real atau;
	real btau;
}

parameters
{
	vector[p] beta;
	real tau2; 
}

transformed parameters 
{
	vector[n] mu;
	for(i in 1:n){
	  mu[i] = row(X,i) * beta;
	}
}

model
{
	beta ~ multi_normal(beta0, Lambda0);
	tau2 ~ inv_gamma(atau, btau);
	y ~ normal(mu, pow(tau2, 0.5));
}

generated quantities 
{
  real llik;
  llik = multi_normal_lpdf(y | mu, diag_matrix(rep_vector(tau2, n)));
}
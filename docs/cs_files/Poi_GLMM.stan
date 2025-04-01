///////////////////////// DATA /////////////////////////////////////
data {
	int<lower = 0> N;       // number of data
	int<lower = 0> p_fix;   // number of covariates, fixed effects
	int<lower = 0> p_rand;  // number of covariates, fixed effects
	
	int<lower = 0> y[N];  	// response vector
	matrix[N, p_fix] X;   	// design matrix (fixed effects)
	matrix[N, p_rand] Z;    // random effects
	
	matrix[p_fix, p_fix] Sigma0;
	vector[p_fix] b0;
	real<lower = 0> sigma20;
}

//////////////////// PARAMETERS /////////////////////////////////
parameters {
	vector[p_fix] beta;        	// regression coefficients 
	vector[p_rand] gamma;	// variances for the prior on beta
}

//////////////////// TRANSFORMED PARAMETERS /////////////////////////////////
transformed parameters 
{
	vector[N] mu;
	for(i in 1:N){
	  mu[i] = exp(row(X, i) * beta + row(Z, i) * gamma);
	}
}

////////////////// MODEL ////////////////////////
model {

	// Likelihood     
	for (s in 1:N)
	{
		y[s] ~ poisson(mu[s]);  
	} 
  
  beta ~ multi_normal(b0, Sigma0);
	for (j in 1:p_rand) 
	{
	 	gamma[j] ~ normal(0.0, pow(sigma20, 0.5));
	}
}

////////////////// GENERATED QUANTITIES ////////////////////////
generated quantities 
{
  	vector[N] log_lik;
  	for (j in 1:N) 
	{
    		log_lik[j] = poisson_lpmf(y[j] | mu[j]);
  	}
}
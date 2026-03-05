///////////////////////// DATA /////////////////////////////////////
data {
	int<lower = 0> N;       // number of data
	int<lower = 0> p;   // number of covariates, fixed effects
	
	int<lower = 0> y[N];  	// response vector
	matrix[N, p] X;   	// design matrix (fixed effects)
	
	matrix[p, p] Sigma0;
	vector[p] b0;
}

//////////////////// PARAMETERS /////////////////////////////////
parameters {
	vector[p] beta;        	// regression coefficients 
}

//////////////////// TRANSFORMED PARAMETERS /////////////////////////////////
transformed parameters 
{
	vector[N] mu;
	for(i in 1:N){
	  mu[i] = exp(row(X, i) * beta);
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
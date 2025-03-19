// in the data block we can also parse the 
// hyperparameters of the model

data
{
	int<lower = 0> n;       
	int<lower = 0> p;       
	vector[n] y;     	
	matrix[n,p] X;		
	
	vector[p] beta0;
	matrix[p,p] Sigma0;
	real a0;
	real b0;
}

// here we have the parameters on which we set a prior distributional assumption

parameters
{
	real<lower = 0> sigma2;
	vector[p] beta;
}

// a possible strategy is to produce as transformation 
// the linear predictor for each observation

transformed parameters 
{
	vector[n] mu;
	for(i in 1:n){
	  mu[i] = row(X, i) * beta;
	}
}

// in the model block we specify the distributional assumption for the prior
// and the likelihood

model
{
	// Prior:
  beta ~ multi_normal(beta0, Sigma0 * sigma2);
	sigma2 ~ inv_gamma(2., 10.);
	
	// Likelihood:
	y ~ normal(mu, pow(sigma2, 0.5));
}

// return the likelihood values

generated quantities 
{
  vector[n] log_lik;
  for (l in 1:n) 
	{
    log_lik[l] = normal_lpdf(y[l] | mu[l], pow(sigma2, 0.5));
  }
}

// data block

data
{
	int<lower = 0> n;       
	int<lower = 0> p;       
	int y[n];     	
	matrix[n,p] X;		
	
	vector[p] beta0;
	real<lower = 0> alpha1; 
	real<lower = 0> alpha2;
	real<lower = 0> c2;
	real<lower = 0> tau2;
}

// parameters

parameters
{
	vector[p] beta;
	real<lower = 0, upper = 1> probs[p];
}

// transformed linear predictor

transformed parameters 
{
	vector[n] mu;
	for(i in 1:n){
	  mu[i] = Phi(row(X, i) * beta);
	}
}

// in the model block we specify the distributional assumption for the prior
// and the likelihood

model
{
	// Prior:
	for(j in 1:p){
	  probs[j] ~ beta(alpha1, alpha2);
	  target += log_mix(probs[j], normal_lpdf(beta[j] | beta0[j], pow(c2 * tau2, 0.5)),
                    normal_lpdf(beta[j] | beta0[j], pow(tau2, 0.5)));
	}
	
	// Likelihood:
	y ~ bernoulli(mu);
}

// return the log-likelihood values

generated quantities 
{
  vector[n] log_lik;
  for (l in 1:n) 
	{
    log_lik[l] = bernoulli_lpmf(y[l] | mu[l]);
  }
}
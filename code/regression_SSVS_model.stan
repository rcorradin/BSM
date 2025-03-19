// in the data block we can also parse the 
// hyperparameters of the model

data
{
	int<lower = 0> n;       
	int<lower = 0> p;       
	vector[n] y;     	
	matrix[n,p] X;		
	
	vector[p] beta0;
	real<lower = 0> a0;
	real<lower = 0> b0;
	real<lower = 0> alpha1; 
	real<lower = 0> alpha2;
	real<lower = 0> c;
	real<lower = 0> tau2;
}

// here we have the parameters on which we set a prior distributional assumption

parameters
{
	real<lower = 0> sigma2;
	vector[p] beta;
	real<lower = 0, upper = 1> probs[p];
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
	sigma2 ~ inv_gamma(a0, b0);
	for(j in 1:p){
	  probs[j] ~ beta(alpha1, alpha2);
	  target += log_mix(probs[j], normal_lpdf(beta[j] | beta0[j], pow(tau2, 0.5)),
                    normal_lpdf(beta[j] | beta0[j], pow(c * tau2, 0.5)));
	}
	
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

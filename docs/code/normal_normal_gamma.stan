// in the data block we can also parse the 
// hyperparameters of the model

data
{
	int<lower = 0> n;       
	vector[n] y;     	
	real m0; 
	real s20;
	real a0;
	real b0;
}

// here we have the parameters on which we set a prior distributional assumption

parameters
{
	real<lower = 0> sigma2;
	real mu;
}


// in the model block we specify the distributional assumption for the prior
// and the likelihood

model
{
	// Prior:
	sigma2 ~ inv_gamma(a0, b0);
	mu ~ normal(m0, pow(s20, 0.5));
	
	// Likelihood:
	y ~ normal(mu, pow(sigma2, 0.5));
}

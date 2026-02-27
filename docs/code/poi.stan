// in the data block we can also parse the 
// hyperparameters of the model

data
{
	int<lower = 0> n;       
	int<lower = 0> y[n];     	
	real a0;
	real b0;
}

// here we have the parameters on which we set a prior distributional assumption

parameters
{
	real lambda;
}


// in the model block we specify the distributional assumption for the prior
// and the likelihood

model
{
	// Prior:
	lambda ~ gamma(a0, b0);

	// Likelihood:
	y ~ poisson(lambda);
}

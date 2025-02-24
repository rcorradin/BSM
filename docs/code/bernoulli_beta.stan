// in the data block we can also parse the 
// hyperparameters of the model

data
{
	int<lower = 0> n;       
	int y[n];     	
	real a0;
	real b0;
}

// here we have the parameters on which we set a prior distributional assumption

parameters
{
	real<lower = 0, upper = 1> theta;
}

// in the model block we specify the distributional assumption for the prior
// and the likelihood

model
{
	// Prior:
	theta ~ beta(a0, b0);
	
	// Likelihood:
	y ~ bernoulli(theta);
	
	// for(i in 1:n){
	//   y[i] ~ bernoulli(theta);
	// }
}

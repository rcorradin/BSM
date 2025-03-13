// in the data block we can also parse the 
// hyperparameters of the model

data
{
	int<lower = 0> n;       
	int<lower = 0> p;       
	int<lower = 0> y[n];     	
	matrix[n,p] X;		
	
	vector[p] beta0;
	matrix[p,p] Sigma0;
}

// here we have the parameters on which we set a prior distributional assumption

parameters
{
	vector[p] beta;
}

// a possible strategy is to produce as transformation 
// the exp of linear predictor for each observation

transformed parameters 
{
	vector[n] lambda;
	for(i in 1:n){
	  lambda[i] = exp(row(X, i) * beta);
	}
}

// in the model block we specify the distributional assumption for the prior
// and the likelihood

model
{
	// Prior:
  beta ~ multi_normal(beta0, Sigma0);
  
	// Likelihood:
	for(i in 1:n){
	 y[i] ~ poisson(lambda[i]);
	}
}

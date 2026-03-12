///////////////////////// DATA /////////////////////////////////////
data {
	int<lower = 0> N;       
	int<lower = 0> p_fix;   
	int<lower = 0> p_rand;  
	
	int<lower = 0> y[N];  	
	real X[N];  	          
	matrix[N, p_rand] Z;    
	
	int ngr;
	int ngr_nested;
	
	int group[N];
	int group_nested[N];
	
	real psi20;
	real sigma20;
	real tau20;
	real eta20;
	real lambda20;
}

//////////////////// PARAMETERS /////////////////////////////////
parameters {
  real xi;
	real beta;        	
	vector[ngr] theta;
	matrix[ngr,p_rand] gamma;
	vector[ngr_nested] alpha;
}

//////////////////// TRANSFORMED PARAMETERS /////////////////////////////////
transformed parameters 
{
	vector[N] mu;
	for(i in 1:N){
	  mu[i] = exp(xi + 
	              X[i] * beta + 
	              theta[group[i]] + 
	              row(Z, i) * (row(gamma, group[i]))' + 
	              alpha[group_nested[i]]);
	}
}

////////////////// MODEL ////////////////////////
model {
	for (s in 1:N)
	{
		y[s] ~ poisson(mu[s]);  
	} 
  
  beta ~ normal(0, pow(sigma20, 0.5));
	theta ~ normal(0, pow(eta20, 0.5));
	alpha ~ normal(0, pow(lambda20, 0.5));
	for(i in 1:ngr){
	  for(j in 1:p_rand){
	    gamma[i,j] ~ normal(0, pow(tau20, 0.5));
	  }
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
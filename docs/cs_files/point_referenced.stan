data
{
	int<lower = 0> n;       
	int<lower = 0> p;     
	int<lower = 0> ord;
	vector[n] y;
	matrix[n,p] X;
	matrix[n,2] S;
	vector[p] x0;
	vector[2] s0;
	
	vector[p] beta0;
	matrix[p,p] Lambda0;
	real asigma;
	real bsigma;
	real atau;
	real btau;
	real aphi; 
	real bphi;
}

transformed data
{
  matrix[n,n] dist_matrix;
  vector[n] dist_new;
  for(i in 1:n){
    for(j in 1:n){
      dist_matrix[i,j] = sqrt(sum(pow(row(S,i) - row(S,j), 2)));
    }
    dist_new[i] = sqrt(sum(pow(row(S,i) - s0', 2)));
  }
}

parameters
{
	vector[p] beta;
	real<lower = 0> sigma2; 
	real<lower = 0> tau2; 
	real<lower = 0> phi;
}

transformed parameters 
{
	vector[n] mu;
	vector[n] gamma;
	for(i in 1:n){
	  mu[i] = row(X,i) * beta;
	}
	matrix[n,n] H;
	for(i in 1:n){
    for(j in 1:n){
      H[i,j] = exp(- pow(phi * dist_matrix[i,j], ord));
    }
    gamma[i] = sigma2 * exp(- pow(phi * dist_new[i], ord));
  }
}

model
{
	beta ~ multi_normal(beta0, Lambda0);
	sigma2 ~ inv_gamma(asigma, bsigma);
	tau2 ~ inv_gamma(atau, btau);
	phi ~ inv_gamma(aphi, bphi);
  y ~ multi_normal(mu, add_diag(sigma2 * H, tau2));
}

generated quantities 
{
  real llik;
  llik = multi_normal_lpdf(y | mu, add_diag(sigma2 * H, tau2));
  
  vector[2] new_pred_mean_var;
  new_pred_mean_var[1] = x0' * beta + gamma' * inverse(add_diag(sigma2 * H, tau2)) * (y - mu);
  new_pred_mean_var[2] = tau2 + sigma2 - gamma' * inverse(add_diag(sigma2 * H, tau2)) * gamma;
}

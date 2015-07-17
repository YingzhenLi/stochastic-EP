// Multivariate GMM model with diagonal matrix prior for the means
// now the sampler knows the true label of each datapoint


data {
	int<lower=1> D; 		// number of dimensions
	int<lower=1> N; 		// number of samples
	int<lower=1> J; 		// number of mixture components
	vector[D] X[N]; 		// data to train
	int<lower=1> y[N];		// label of datapoints
	simplex[J] weights; 	// mixture weights
	vector[D] mu[J];		// mean of prior on mixture means
	matrix[D, D] sigma[J];		// covariance (diagonal)
	real<lower=0.0> std_noise;	// std of noise
}
parameters {
	vector[D] phi[J]; // means of mixture components
}
transformed parameters {
    vector[D] beta[J];
    beta <- phi;
}
model {
	for (j in 1:J){
		phi[j] ~ multi_normal_prec(mu[j], sigma[j]);	// sample the cluster means
	}
	// sample datapoints
	for (n in 1:N){
		for (d in 1:D){
			X[n, d] ~ normal(beta[y[n], d], std_noise);
		}
	}
} 

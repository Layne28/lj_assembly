/* 
   Functions to compute the sample mean and variance of a collection of samples.
   Can use an array or a vector. 
   Also includes a minimum variance estimator in the case where each sample is independent 
   and comes with its own mean and variance. 
   
*/


void sampleStats(double* X, int N, double& M, double& V) {
	//return the sample mean and variance of the array X with N elements

	//init the mean and variance at 0, get sample size
	M = 0; V = 0;

	//compute the mean
	for (int i = 0; i < N; i++) {
		M += X[i];
	}
	M /= float(N);

	//compute the variance 
	for (int i = 0; i < N; i++) {
		V += (X[i]-M) * (X[i]-M);
	}
	V /= float(N-1.0);

}

void sampleStats(std::vector<double> X, double& M, double& V) {
	//return the sample mean of the vector X

	//init the mean and variance at 0, get sample size
	M = 0; V = 0;
	int N = X.size();

	//compute the mean
	for (int i = 0; i < N; i++) M += X[i];
	M /= float(N);

	//compute the variance 
	for (int i = 0; i < N; i++) V += (X[i]-M) * (X[i]-M);
	V /= (N-1);
}

void minVarEstimate(int sampleSize, double* means, double* variances, double& M, double& V) {
  	/* Compute the minimum variance weighted estimator, given a collection of sample means and variances. 
     	Assume each sample mean and variance comes from an IID source (independent and identically distributed).
    	 
     	The weighting is given by M = 1/S * sum(m_j/V_j), where S = 1/sum(1/V_j). Derivable via Lagrange multipliers.
     	The resulting variance is 1/S.
 	*/

	//init the output mean and variance, and normalizer
	M = 0; V = 0;
	double S = 0;

	//compute the normalizing term, S, and sum the means over the variances
	for (int i = 0; i < sampleSize; i++) {
		S += 1.0/variances[i];
    		M += means[i]/variances[i];
	}

	V = 1.0/S;
	M /= S;
}

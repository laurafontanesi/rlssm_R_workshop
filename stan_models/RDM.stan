functions {
     real race_pdf(real t, real b, real v){
          real pdf;
          pdf = b/sqrt(2 * pi() * pow(t, 3)) * exp(-pow(v*t-b, 2) / (2*t));
          return pdf;
     }

     real race_cdf(real t, real b, real v){
          real cdf;
          cdf = Phi((v*t-b)/sqrt(t)) + exp(2*v*b) * Phi(-(v*t+b)/sqrt(t));
          return cdf;
     }

     real race_lpdf(matrix RT, vector  ndt, vector b, vector drift_cor, vector drift_inc){

          real t;
          vector[rows(RT)] prob;
          real cdf;
          real pdf;
          real out;

          for (i in 1:rows(RT)){
               t = RT[i,1] - ndt[i];
               if(t > 0){
                  if(RT[i,2] == 1){
                    pdf = race_pdf(t, b[i], drift_cor[i]);
                    cdf = 1 - race_cdf(t, b[i], drift_inc[i]);
                  }
                  else{
                    pdf = race_pdf(t, b[i], drift_inc[i]);
                    cdf = 1 - race_cdf(t, b[i], drift_cor[i]);
                  }
                  prob[i] = pdf*cdf;

                if(prob[i] < 1e-10){
                    prob[i] = 1e-10;
                }
               }
               else{
                    prob[i] = 1e-10;
               }
          }
          out = sum(log(prob));
          return out;
     }
}
data {
	int<lower=1> N;									// number of data items
	int<lower=1,upper=2> choices[N];				// coded as 1 and 2
	real<lower=0> rt[N];							// rt
}
transformed data {
	matrix [N, 2] RT;
	for (n in 1:N){
		RT[n, 1] = rt[n];
		RT[n, 2] = choices[n];
	}
}
parameters {
	real ndt;
	real threshold;
	real drift1;
	real drift2;
}
transformed parameters {
	vector<lower=0> [N] drift1_t;				// trial-by-trial drift rate for predictions
	vector<lower=0> [N] drift2_t;				// trial-by-trial drift rate for predictions
	vector<lower=0> [N] threshold_t;				// trial-by-trial threshold
	vector<lower=0> [N] ndt_t;						// trial-by-trial ndt

	real<lower=0> transf_drift1;
	real<lower=0> transf_drift2;
	real<lower=0> transf_threshold;
	real<lower=0> transf_ndt;
	transf_drift1 = log(1 + exp(drift1));
	transf_drift2 = log(1 + exp(drift2));
	transf_threshold = log(1 + exp(threshold));
	transf_ndt = log(1 + exp(ndt));

	for (n in 1:N) {
		drift1_t[n] = transf_drift1;
		drift2_t[n] = transf_drift2;
		threshold_t[n] = transf_threshold;
		ndt_t[n] = transf_ndt;
	}
}
model {
	ndt ~  normal(-1, 3);
	threshold ~ normal(0, 3);
	drift1 ~ normal(1, 5);
	drift2 ~ normal(1, 5);
	RT ~ race(ndt_t, threshold_t, drift1_t, drift2_t);
}
generated quantities {
	vector[N] log_lik;
	{
	for (n in 1:N){
		log_lik[n] = race_lpdf(block(RT, n, 1, 1, 2)| segment(ndt_t, n, 1), segment(threshold_t, n, 1), segment(drift1_t, n, 1), segment(drift2_t, n, 1));
	}
	}
}
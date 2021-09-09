data {
	int<lower=1> N;									// number of data items
	int<lower=-1,upper=1> accuracy[N];				// accuracy (-1, 1)
	real<lower=0> rt[N];							// rt
	real<lower=0, upper=1> starting_point;			// starting point diffusion model not to estimate
}
parameters {
	real drift_trialmu;
	real threshold;
	real ndt;
	real<lower=0> drift_trialsd;
	real z_drift_trial[N];
}
transformed parameters {
	real drift_ll[N];								// trial-by-trial drift rate for likelihood (incorporates accuracy)
	real drift_t[N];								// trial-by-trial drift rate for predictions
	real<lower=0> threshold_t[N];					// trial-by-trial threshold
	real<lower=0> ndt_t[N];							// trial-by-trial ndt

	real transf_threshold;
	real transf_ndt;
	real<lower=0> transf_drift_trialsd;

	transf_threshold = log(1 + exp(threshold));
	transf_ndt = log(1 + exp(ndt));
	transf_drift_trialsd = drift_trialsd;

	for (n in 1:N) {
		drift_t[n] = drift_trialmu + z_drift_trial[n]*drift_trialsd;
		drift_ll[n] = drift_t[n]*accuracy[n];
		threshold_t[n] = transf_threshold;
		ndt_t[n] = transf_ndt;
	}
}
model {
	drift_trialmu ~ normal(0, 3);
	threshold ~ normal(0, 5);
	ndt ~ normal(0, 1);

	drift_trialsd ~ normal(0, 3);
	z_drift_trial ~ normal(0, 1);

	rt ~ wiener(threshold_t, ndt_t, starting_point, drift_ll);
}
generated quantities {
	vector[N] log_lik;

	{for (n in 1:N) {
		log_lik[n] = wiener_lpdf(rt[n] | threshold_t[n], ndt_t[n], starting_point, drift_ll[n]);
	}
}
}
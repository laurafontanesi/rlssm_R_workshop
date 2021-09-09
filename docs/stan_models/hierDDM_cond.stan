data {
	int<lower=1> N;									// number of data items
	int<lower=1> L;									// number of levels
	int<lower=1> C;									// number of conditions
	int<lower=1, upper=L> participant[N];			// level (participant)
	int<lower=1, upper=C> condition[N];				// condition
	int<lower=-1,upper=1> accuracy[N];				// accuracy (-1, 1)
	real<lower=0> rt[N];							// rt
	real<lower=0, upper=1> starting_point;			// starting point diffusion model not to estimate
}
parameters {
	real mu_drift[C];
	real mu_threshold[C];
	real mu_ndt;

	real<lower=0> sd_drift[C];
	real<lower=0> sd_threshold[C];
	real<lower=0> sd_ndt;

	vector[C] z_drift[L];
	vector[C] z_threshold[L];
	real z_ndt[L];
}
transformed parameters {
	real drift_ll[N];								// trial-by-trial drift rate for likelihood (incorporates accuracy)
	real drift_t[N];								// trial-by-trial drift rate for predictions
	real<lower=0> threshold_t[N];					// trial-by-trial threshold
	real<lower=0> ndt_t[N];							// trial-by-trial ndt

	vector[C] drift_sbj[L];
	vector[C] threshold_sbj[L];
	real<lower=0> ndt_sbj[L];
	
	for (l in 1:L) {
		for (c in 1:C) {
			drift_sbj[l][c] = mu_drift[c] + z_drift[l][c]*sd_drift[c];
			threshold_sbj[l][c] = log(1 + exp(mu_threshold[c] + z_threshold[l][c]*sd_threshold[c]));
		}
		ndt_sbj[l] = log(1 + exp(mu_ndt + z_ndt[l]*sd_ndt));
	}

	for (n in 1:N) {
		drift_t[n] = drift_sbj[participant[n]][condition[n]];
		drift_ll[n] = drift_t[n]*accuracy[n];
		threshold_t[n] = threshold_sbj[participant[n]][condition[n]];
		ndt_t[n] = ndt_sbj[participant[n]];
	}
}
model {
	mu_drift ~ normal(0, 3);
	mu_threshold ~ normal(0, 5);
	mu_ndt ~ normal(-1, 3);

	sd_drift ~ normal(0, 3);
	sd_threshold ~ normal(0, 5);
	sd_ndt ~ normal(0, 1);

	for (l in 1:L) {
		z_drift[l] ~ normal(0, 1);
		z_threshold[l] ~ normal(0, 1);
		z_ndt[l] ~ normal(0, 1);
	}

	rt ~ wiener(threshold_t, ndt_t, starting_point, drift_ll);
}
generated quantities {
	vector[N] log_lik;

	{for (n in 1:N) {
		log_lik[n] = wiener_lpdf(rt[n] | threshold_t[n], ndt_t[n], starting_point, drift_ll[n]);
	}
}
}
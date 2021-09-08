data {
	int<lower=1> N;									// number of data items
	int<lower=1> K;									// number of options
	vector[N] f_cor;								// feedback correct option
	vector[N] f_inc;								// feedback incorrect option
	vector[N] trial;								// trial number
	int<lower=1, upper=K> cor_option[N];			// correct option
	int<lower=1, upper=K> inc_option[N];			// incorrect option
	int<lower=-1,upper=1> accuracy[N];				// accuracy (-1, 1)
	real<lower=0> rt[N];							// rt

	real initial_value;								// intial value for learning in the first block

	real<lower=0, upper=1> starting_point;			// starting point diffusion model not to estimate
}
transformed data {
	vector[K] Q0;
	Q0 = rep_vector(initial_value, K);
}
parameters {
	real alpha;
	real drift_scaling;
	real threshold;
	real ndt;
}
transformed parameters {
	real drift_ll[N];								// trial-by-trial drift rate for likelihood (incorporates accuracy)
	real drift_t[N];								// trial-by-trial drift rate for predictions
	real<lower=0> threshold_t[N];					// trial-by-trial threshold
	real<lower=0> ndt_t[N];							// trial-by-trial ndt

	vector[K] Q;									// Q state values

	real transf_alpha;
	real transf_drift_scaling;
	real transf_threshold;
	real transf_ndt;

	transf_alpha = Phi(alpha);						// for the output
	transf_drift_scaling = log(1 + exp(drift_scaling));
	transf_threshold = log(1 + exp(threshold));
	transf_ndt = log(1 + exp(ndt));

	Q = Q0;
	for (n in 1:N) {
		drift_t[n] = transf_drift_scaling*(Q[cor_option[n]] - Q[inc_option[n]]);
		drift_ll[n] = drift_t[n]*accuracy[n];
		threshold_t[n] = transf_threshold;
		ndt_t[n] = transf_ndt;

		Q[cor_option[n]] = Q[cor_option[n]] + transf_alpha*(f_cor[n] - Q[cor_option[n]]);
		Q[inc_option[n]] = Q[inc_option[n]] + transf_alpha*(f_inc[n] - Q[inc_option[n]]);
	}
}
model {
	alpha ~ normal(0, .8);
	drift_scaling ~ normal(0, 5);
	threshold ~ normal(0, 5);
	ndt ~ normal(-1, 3);

	rt ~ wiener(threshold_t, ndt_t, starting_point, drift_ll);
}
generated quantities {
	vector[N] log_lik;

	{for (n in 1:N) {
		log_lik[n] = wiener_lpdf(rt[n] | threshold_t[n], ndt_t[n], starting_point, drift_ll[n]);
	}
	}
}
data {
    // declare variables
    int n_sp;
    int n_obs;
    int map[n_obs];
    real y[n_obs];
    real x[n_obs];
}
parameters {
    // declare variables
    real b0_mean;
    real b1_mean;
    real z_b0[n_sp];
    real z_b1[n_sp];
    real<lower=0> sig_b0;
    real<lower=0> sig_b1;
    real<lower=0> sig_res;
}
transformed parameters {
    // delcare variables
    real b0[n_sp];
    real b1[n_sp];
    // intercepts and slopes by group
    for(j in 1:n_sp){
        b0[j] = b0_mean + sig_b0*z_b0[j];
        b1[j] = b1_mean + sig_b1*z_b1[j];
    }
}
model {
    // priors
    b0_mean ~ normal(0, 1);
    b1_mean ~ normal(0, 1);
    z_b0 ~ normal(0, 1);
    z_b1 ~ normal(0, 1);
    sig_b0 ~ normal(0, 1) T[0, ];
    sig_b1 ~ normal(0, 1) T[0, ];
    sig_res ~ normal(0, 1) T[0, ];
    // observations
    for(i in 1:n_obs){
        y[i] ~ normal(b0[map[i]] + b1[map[i]]*x[i], sig_res);
    }
}

data {
    // indices
    int n_obs;                          // # observations
    int n_ts;                           // # time series
    int obs_per[n_ts];                  // # observations per time series
    int n_coef;                         // # coefficients (fixed effects + intercepts)
    int g_per_ff[n_coef];               // # groups per fixed effect
    int lev_per_g[sum(g_per_ff)];       // # levels per group (repeated by fixed effect)
    int b_groups[n_ts, sum(g_per_ff)];  // grouping structure for betas
    int p_groups[n_ts];                 // grouping structure for phis
    // data
    real y[n_obs];                      // response variables
    real x[n_obs, n_coef];              // predictor variables
    real<lower=0> time[n_obs];          // response variable times
    real<lower=0> p_bound;              // upper bound for phis
}
parameters {
    real alpha[n_coef];                         // fixed effects and intercepts
    real z[sum(lev_per_g)];                     // standardized variates for group levels
    real<lower=0> sig_beta[sum(g_per_ff)];      // group standard deviations
    real<lower=0> sig_res;                      // residual standard deviation
}
transformed parameters {
    real beta[n_ts,n_coef]; // coefficients
    real y_pred[n_obs];     // predicted values
    {
        int xy_pos = 1;         // position in x and y vectors
        // loop over time series:
        for (ts in 1:n_ts){
            int beta_pos = 1;   // position in sig_beta vector
            // loop over coefficients:
            for (c in 1:n_coef){
                if(g_per_ff[c]==0) {
                    beta[ts,c] = alpha[c];
                } else {
                    real sigs[g_per_ff[c]];
                    real zs[g_per_ff[c]];
                    // loop over groups:
                    for(i in beta_pos:(beta_pos + g_per_ff[c] - 1)){
                        sigs[i - beta_pos + 1] = sig_beta[i];
                        zs[i - beta_pos + 1] = z[b_groups[ts, i]];
                    } // i
                    // calculate coefficients (fixed + random effect):
                    beta[ts,c] = alpha[c] + dot_product(sigs, zs);
                    beta_pos += g_per_ff[c];
                }
            } // c
            // predicted values:
            for (t in xy_pos:(xy_pos + obs_per[ts] - 1)) {
                    y_pred[t] = dot_product(beta[ts,], x[t,]);
           } // t
           xy_pos += obs_per[ts];
        } // ts
    }
}
model {
    // priors:
    alpha ~ normal(0, 1);
    z ~ normal(0, 1);
    for (i in 1:sum(g_per_ff)){
        sig_beta[i] ~ gamma(1.5, 3);
    }
    sig_res ~ gamma(1.5, 3);
    // model:
    y ~ normal(y_pred, sig_res);
}
generated quantities {
  real log_lik[n_obs];
  real log_lik_sum;
  for(i in 1:n_obs){
    log_lik[i] = normal_lpdf(y[i] | y_pred[i], sig_res);
  }
  log_lik_sum = sum(log_lik);
}

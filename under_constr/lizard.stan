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
}
parameters {
    real alpha[n_coef];                         // fixed effects and intercepts
    real z[sum(lev_per_g)];                     // standardized variates for group levels
    real<lower=0, upper=1> phi[max(p_groups)];  // autoregressive parameter for each
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
            y_pred[xy_pos] = dot_product(beta[ts,], x[xy_pos,]);
            for (t in (xy_pos + 1):(xy_pos + obs_per[ts] - 1)) {
                    y_pred[t] = dot_product(beta[ts,], x[t,]) +
                        phi[p_groups[ts]] * (y[t-1] - dot_product(beta[ts,], x[t-1,]));
           } // t
        } // ts
    }
}
model {
    // priors:
    alpha ~ normal(0, 1);
    z ~ normal(0, 1);
    phi ~ uniform(0, 1);
    for (i in 1:sum(g_per_ff)){
        sig_beta[i] ~ normal(0, 1) T[0, ];
    }
    sig_res ~ normal(0, 1) T[0, ];
    // model:
    y ~ normal(y_pred, sig_res);
}

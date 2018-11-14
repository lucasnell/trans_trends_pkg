data {
    // indeces
    int n_obs; // number of observations
    int ns; // number of time series
    int obs_per[ns]; // number of observations per time series
    int nq; // number of coefficients (same as number of fixed effects + intercepts)
    int k[nq]; // number of groups per fixed effect;
    int l[sum(k)]; // number of levels per group (repeated for each fixed effect)
    int g[ns, sum(k)]; // grouping structure for betas
    int h[ns]; // grouping structure for phis
    // data
    real y[n_obs]; // response variables
    real x[n_obs,nq]; // predictor variables
}
parameters {
    real alpha[nq]; // fixed effects and intercepts
    real z[sum(l)]; // standardized variates for group levels
    real<lower=0, upper=1> phi[max(h)]; // autoregressive paramter for each
    real<lower=0> sig_beta[sum(k)]; // group standard deviations
    real<lower=0> sig_res; // residual standard deviation
}
transformed parameters {
    real beta[ns,nq]; // coefficients
    real y_pred[n_obs]; // predicted values
    {
        int pos1 = 1;
        // loop over time series
        for (s in 1:ns){
            int pos2 = 1;
            // loop over coefficients
            for (q in 1:nq){
                if(k[q]==0) {
                    beta[s,q] = alpha[q];
                } else {
                    real sigs[k[q]];
                    real zs[k[q]];
                    // loop over groups
                    for(i in pos2:(pos2 + k[q] - 1)){
                        sigs[i - pos2 + 1] = sig_beta[i];
                        zs[i - pos2 + 1] = z[g[s, i]];
                    } // i
                    // calculate coefficients (fixed + random effect)
                    beta[s,q] = alpha[q] + dot_product(sigs, zs);
                    pos2 += k[q];
                }
            } // q
            // predicted values
            y_pred[pos1] = dot_product(beta[s,], x[pos1,]);
            for (t in (pos1 + 1):(pos1 + obs_per[s] - 1)) {
                    y_pred[t] = dot_product(beta[s,], x[t,]) + phi[h[s]]*(y[t-1] - dot_product(beta[s,], x[t-1,]));
           } // t
        } // s
    }
}
model {
    // priors
    alpha ~ normal(0, 1);
    z ~ normal(0, 1);
    phi ~ uniform(0, 1);
    for (i in 1:sum(k)){
        sig_beta[i] ~ normal(0, 1) T[0, ];
    }
    sig_res ~ normal(0, 1) T[0, ];
    // model
    y ~ normal(y_pred, sig_res);
}

data {
    // indeces
    int Y // number of observations
    int S; // number of time series
    int Y_S; // number of observations per time series
    int Q; // number of fixed effects
    int K[Q]; // number of groups per fixed effect;
    int L[sum[K]]; // number of levels per group (repeated for each fixed effect)
    // data
    real y[Y]; // number of observations
}
parameters {
    real alpha[Q]; // fixed effect slope
    real z[sum(L)]; // standardized variates for group levels
    real<lower=0> sig_beta[sum(K)]; // group standard deviations
    real<lower=0> sig_res; // residual standard deviation
}
transformed parameters {
    real beta[S,Q];
    {   // loop over time series
        for (s in 1:S){
            int pos = 1;
            // loop over fixed effects
            for (q in 1:Q){
                real sigs[K[q]];
                real zs[K[q]];
                // loop over groups
                for(i in pos:(pos + K[q] - 1)){
                    sigs[i - pos + 1] = sig_beta[i];
                    zs[i - pos + 1] = z[G[s, i]];
                } // i
                // calculate fixed + group effect
                beta[s,q] = alpha[q] + sum(sigs.*zs);
                pos = pos + K[q];
            } // q
        } // s
    }
}
model {
    // priors
    z ~ normal(0, 1);
    for (i in 1:sum(K)){
        sig_beta ~ normal(0, 1) T[0, ];
    }
    sig_res ~ normal(0, 1) T[0, ];
    // model
    y ~ normal(y_pred, sig_res)
}

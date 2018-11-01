data {
    // indeces
    int Y // number of observations
    int S; // number of time series
    int Y_S; // number of observations per time series
    int Q; // number of coefficients (same as number of fixed effects + intercepts)
    int K[Q]; // number of groups per fixed effect;
    int L[sum[K]]; // number of levels per group (repeated for each fixed effect)
    int G[S, sum(Q)]; // grouping structure
    // data
    real y[Y]; // response variables
    real X[Y,QQ]; // predictor variables
}
parameters {
    real alpha[Q]; // fixed effects and intercepts
    real z[sum(L)]; // standardized variates for group levels
    real<lower=0> sig_beta[sum(K)]; // group standard deviations
    real<lower=0> sig_res; // residual standard deviation
}
transformed parameters {
    real beta[S,Q]; // coefficients
    real y_pred[Y]; // predicted values
    {
        int pos1 = 1;
        // loop over time series
        for (s in 1:S){
            int pos2 = 1;
            // loop over coefficients
            for (q in 1:Q){
                real sigs[K[q]];
                real zs[K[q]];
                // loop over groups
                for(i in pos2:(pos2 + K[q] - 1)){
                    sigs[i - pos2 + 1] = sig_beta[i];
                    zs[i - pos2 + 1] = z[G[s, i]];
                } // i
                // calculate coefficients (fixed + random effect)
                beta[s,q] = alpha[q] + dot_product(sigs, zs);
                pos2 = pos2 + K[q];
            } // q
            // predicted values
            y_pred[pos1] = dot_product(beta[s,], X[pos1,]);
            for (t in (pos1 + 1):(pos1 + Y_S[s] -1)) {
                    y_pred[t] = dot_product(beta[s,], X[t,]) + phi[s]*(y[t-1] - dot_product(beta[s,], X[t-1,]));
           } // t
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

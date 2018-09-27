data{
  int n_ts; // number of time series
  int n_per[n_ts]; // number of observations per time series
  int n_obs; // total number of observations
  real pop_obs[n_obs]; // scaled population abundances
}
parameters{
  // declare variables
  real mu[n_ts];
  real logit_rho[n_ts];
  real pop_init[n_ts];
  real z[n_obs - n_ts];
}
transformed parameters{
  // declare variables
  real rho[n_ts];
  real pop[n_obs];
  int pos;
  // interate over each time series
  pos = 1;
  for(i in 1:n_ts){
    rho[i] = exp(logit_rho[i])/(exp(logit_rho[i]) + 1);
    // intial values
    pop[pos] = pop_init[i];
    for(t in (pos + 1):(pos + n_per[n_ts] - 1)){
      // state process
      n[t] = mu[i]*x[t] + rho[i]*(pop[t-1] - mu[i]*x[t-1]) + sig_proc*z_proc[t - (i - 1)];
    }
    pos = pos + n_per[i];
  }
}
model{
  // priors
  mu ~ normal(mu_prior[1], mu_prior[2]);
  logit_rho ~ normal(logit_rho_prior[1], logit_rho_prior[2]);
  pop_init ~ normal(pop_prior[1], pop_prior[2]);
  // z for process error
  z ~ normal(0, 1);
  // observation process
  n_obs ~ normal(n, sig_obs);
}
data{
  int n_ts; // number of time series
  int n_per[n_ts]; // number of observations per time series
  int n_obs; // total number of observations
  real pop_obs[n_obs]; // scaled population abundances
  real x[n_obs]; // scaled environmental variables
}
parameters{
  // declare variables
  real beta[n_ts]; // coefficient for environmental response
  real logit_rho[n_ts]; // autoregressive parameter (on logit scale)
  real pop_init[n_ts]; // initial population size
  real z[n_obs - n_ts]; // z-scores for process error
}
transformed parameters{
  // declare variables
  real rho[n_ts]; // autoregressive parameter
  real pop[n_obs]; // "true" population size 
  int pos; // index for iterating 
  // interate over each time series
  pos = 1;
  for(i in 1:n_ts){
      // inverse logit transform
    rho[i] = exp(logit_rho[i])/(exp(logit_rho[i]) + 1);
    // intial values
    pop[pos] = pop_init[i];
    for(t in (pos + 1):(pos + n_per[n_ts] - 1)){
      // state process
      n[t] = beta[i]*x[t] + rho[i]*(pop[t-1] - beta[i]*x[t-1]) + sig_proc*z_proc[t - (i - 1)];
    }
    pos = pos + n_per[i];
  }
}
model{
  // priors
  beta ~ normal(beta_prior[1], beta_prior[2]);
  logit_rho ~ normal(logit_rho_prior[1], logit_rho_prior[2]);
  pop_init ~ normal(pop_prior[1], pop_prior[2]);
  // z for process error
  z ~ normal(0, 1);
  // observation process
  n_obs ~ normal(n, sig_obs);
}
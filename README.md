[![R build status](https://github.com/lucasnell/trans_trends_pkg/workflows/R-CMD-check/badge.svg)](https://github.com/lucasnell/trans_trends_pkg/actions)

This repository contains an R package that fits the models in
<https://github.com/lucasnell/trans_trends>.


Most of the files here are dedicated to building and compiling the stan files
and making them available to the R package.
If you want to dig into the meat of this package, you should check out 
the following files:

- `data/myv_arth.rda` is the cleaned version of data used in the analysis
- `data-raw/` contains raw data and code used to clean it
    - `clean_data.R` cleans the raw data
    - `myvatn_infalls.csv` is raw data from infall traps that estimate
      midge deposition
    - `myvatn_pitfalls.csv` is raw data from predatory arthropod pitfall traps
- `exec/` contains the stan files defining models
    - `armm.stan` defines an autoregressive mixed model, 
      with no observation error
      and Gaussian response distributions
    - `armm_ss.stan` defines an autoregressive mixed model, 
      with observation error 
      and Gaussian response distributions
      (this version was used in the paper)
    - `armm_ss_lnp.stan` defines an autoregressive mixed model, 
      with observation error 
      and lognormal Poisson response distributions
    - `mm.stan` defines a non-autoregressive mixed model, 
      with no observation error
      and Gaussian response distributions
- `R/armmr.R` contains the `armmr` function [at the end of the file],
  which fits autoregressive mixed models.
  This is the function you'll use from R if you want to use this package.
- `R/mode.R` contains the `posterior_mode` function, which estimates 
  the marginal parameter mode using kernel density estimation
  (not used in the paper)
- `src/hpdi.cpp` contains the `hpdi` function, which calculates 
  Highest Posterior Density intervals (not used in the paper)

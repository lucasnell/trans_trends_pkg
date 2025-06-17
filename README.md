[![R build status](https://github.com/lucasnell/trans_trends_pkg/workflows/R-CMD-check/badge.svg)](https://github.com/lucasnell/trans_trends_pkg/actions)

[![DOI](https://zenodo.org/badge/314584876.svg)](https://zenodo.org/badge/latestdoi/314584876)


This repository contains an R package that fits the models in
<https://github.com/lucasnell/trans_trends>.


Most of the files here are dedicated to building and compiling the stan files
and making them available to the R package.
If you want to dig into the meat of this package, you should check out 
the following files:


- `inst/stan/`: the stan files defining models
    - `armm.stan` defines an autoregressive mixed model, 
      with no observation error
      and Gaussian response distributions
    - `armm_ss.stan` defines an autoregressive mixed model, 
      with observation error 
      and Gaussian response distributions
    - `armm_ss_lnp.stan` defines an autoregressive mixed model, 
      with observation error 
      and lognormal Poisson response distributions
      (this version was used in the paper)
    - `mm.stan` defines a non-autoregressive mixed model, 
      with no observation error
      and Gaussian response distributions
- `R/armmr.R`: defines the `armmr` function [at the end of the file],
  which fits autoregressive mixed models.
  This is the function you'll use from R if you want to use this package.
- `R/armmMod_class.R`: class functions (e.g., `print`, `summary`, `ranef`,
  `fixef`) for the `armmMod` class, which is output from `armmr`.
  Many of these are used throughout `trans_trends`.
- `R/mode.R`: the `posterior_mode` function, which estimates 
  the marginal parameter mode using kernel density estimation
  (not used in the paper)
- `src/hpdi.cpp`: the `hpdi` function, which calculates 
  Highest Posterior Density intervals (not used in the paper)

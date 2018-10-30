#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)
library(GGally)

# load data
load("data/myv_arth.rda")





#==========
#========== Prepare data
#==========

# log1p transform, then z-score midges
myv_arth2 = myv_arth %>%
    tbl_df() %>%
    mutate(log_count = log1p(count),
           log_midges = log1p(midges),
           z_midges = (log_midges - mean(log_midges, na.rm=T))/sd(log_midges, na.rm=T),
           group = paste0(trans, dist, taxon) %>% as.factor() %>% as.integer()) %>%
    # group by taxon and them z-score count within taxon
    group_by(taxon) %>%
    mutate(z_count = (log_count - mean(log_count, na.rm=T))/sd(log_count, na.rm=T)) %>%
    ungroup()

# plot: vs. midges
myv_arth2 %>%
    ggplot(aes(z_midges, z_count))+
    facet_wrap(~taxon)+
    geom_point(alpha = 0.4, size = 1)+
    geom_smooth(method="lm", formula = y ~ x, se=F, color="firebrick")+
    theme_bw()

# plot: vs year
myv_arth2 %>%
    ggplot(aes(year, z_count, group = group))+
    facet_wrap(~taxon)+
    geom_line(alpha = 0.5, size = 0.5)+
    theme_bw()

# package data
data_list = myv_arth2 %>%
    arrange(trans, dist, taxon, year) %>%
{list(n_obs = nrow(.),
      n_sp = .$taxon %>% unique %>% length(),
      map = .$taxon %>% as.factor() %>% as.integer(),
      taxon = .$taxon %>% as.factor(),
      n_groups = length(.$group %>% unique()),
      obs_per = {myv_arth2 %>% group_by(group) %>% summarize(n = length(z_count))}$n,
      x = .$z_midges,
      y = .$z_count)}





#==========
#========== Fit model: LMM
#==========

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

# specificaionts
chains = 4
iter = 2000

# fit_lmm model
fit_lmm = stan(file = "under_constr/lmm.stan", data = data_list, seed=1, chains = chains, iter = iter)

# summary of fit_lmm
fit_lmm_summary = summary(fit_lmm, probs=c(0.16, 0.5, 0.84))$summary %>%
{as_data_frame(.) %>%
        mutate(var = rownames(summary(fit_lmm)$summary))}

# check Rhat & n_eff
# note that initial values beta0, alpha, and rho sample for an extra year
# these extra values do not contribute to the likelihood
fit_lmm_summary %>% filter(Rhat > 1.05) %>% select(Rhat, n_eff, var)
fit_lmm_summary %>% filter(n_eff < 0.5*(chains*iter/2)) %>% select(Rhat, n_eff, var) %>% arrange(n_eff) %>%
    mutate(eff_frac = n_eff/(chains*iter/2))




#==========
#========== Plot output: LMM
#==========

# extract coefficients
b0 = fit_lmm_summary %>%
    filter(var %in% paste0("b0[",1:7,"]")) %>%
    mutate(spec = 1:7) %>%
    select(var, spec, `16%`,`50%`,`84%`) %>%
    left_join(data_frame(spec = data_list$map %>% unique, taxon = data_list$taxon %>% unique))

b1 = fit_lmm_summary %>%
    filter(var %in% paste0("b1[",1:7,"]")) %>%
    mutate(spec = 1:7) %>%
    select(var, spec, `16%`,`50%`,`84%`) %>%
    left_join(data_frame(spec = data_list$map %>% unique, taxon = data_list$taxon %>% unique))

# plot coefs
b0 %>%
    ggplot(aes(taxon, `50%`))+
    geom_hline(yintercept = 0, size = 0.3, alpha = 0.5)+
    geom_point(size = 2)+
    geom_errorbar(aes(ymin = `16%`, ymax = `84%`), width = 0)+
    theme_bw()

b1 %>%
    ggplot(aes(taxon, `50%`))+
    geom_hline(yintercept = 0, size = 0.3, alpha = 0.5)+
    geom_point(size = 2)+
    geom_errorbar(aes(ymin = `16%`, ymax = `84%`), width = 0)+
    theme_bw()

# plot predictions
preds = b0 %>%
    expand(taxon, z_midges = -3:1.5) %>%
    left_join(b0 %>% rename(b0 = `50%`) %>% select(taxon, b0)) %>%
    left_join(b1 %>% rename(b1 = `50%`) %>% select(taxon, b1)) %>%
    mutate(z_count = b0 + b1*z_midges)
myv_arth2 %>%
    ggplot(aes(z_midges, z_count))+
    facet_wrap(~taxon)+
    geom_point(alpha = 0.4, size = 1)+
    geom_line(data = preds, color = "red", size = 0.8)+
    geom_smooth(method="lm", formula = y ~ x, se=F, color="blue", linetype = 2)+
    theme_bw()

# fixed parameters by step
fixed_par_v = c("sig_b0","sig_b1","sig_res")
fixed_pars = rstan::extract(fit_lmm, pars=fixed_par_v) %>%
    lapply(as_data_frame) %>%
    bind_cols() %>%
    mutate(chain = rep(1:chains, each = iter/2), step = rep(c(1:(iter/2)), chains))
names(fixed_pars) = c(fixed_par_v,"chain","step")

# examine chains for parameters
fixed_pars %>%
    gather(par, value, -chain, -step) %>%
    filter(par != "lp__") %>%
    ggplot(aes(step, value, color=factor(chain)))+
    facet_wrap(~par, scales="free_y")+
    geom_line(alpha=0.5)+
    theme_bw()

# pairs plot for parameters
ggpairs(fixed_pars %>% select(-chain, -step))

# posterior densities
fixed_pars %>%
    gather(par, value, -chain, -step) %>%
    filter(par != "lp__") %>%
    ggplot(aes(value))+
    facet_wrap(~par, ncol=1)+
    stat_density(alpha=0.5, geom = "line")+
    theme_bw()




#==========
#========== Fit model: AR
#==========

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

# specificaionts
chains = 4
iter = 2000

# fit_AR model
fit_AR = stan(file = "under_constr/lmm_AR.stan", data = data_list, seed=1, chains = chains, iter = iter)

# summary of fit_AR
fit_AR_summary = summary(fit_AR, probs=c(0.16, 0.5, 0.84))$summary %>%
{as_data_frame(.) %>%
        mutate(var = rownames(summary(fit_AR)$summary))}

# check Rhat & n_eff
# note that initial values beta0, alpha, and rho sample for an extra year
# these extra values do not contribute to the likelihood
fit_AR_summary %>% filter(Rhat > 1.05) %>% select(Rhat, n_eff, var)
fit_AR_summary %>% filter(n_eff < 0.5*(chains*iter/2)) %>% select(Rhat, n_eff, var) %>% arrange(n_eff) %>%
    mutate(eff_frac = n_eff/(chains*iter/2))




#==========
#========== Plot output: LMM
#==========

# extract coefficients
b0 = fit_AR_summary %>%
    filter(var %in% paste0("b0[",1:7,"]")) %>%
    mutate(spec = 1:7) %>%
    select(var, spec, `16%`,`50%`,`84%`) %>%
    left_join(data_frame(spec = data_list$map %>% unique, taxon = data_list$taxon %>% unique))

b1 = fit_AR_summary %>%
    filter(var %in% paste0("b1[",1:7,"]")) %>%
    mutate(spec = 1:7) %>%
    select(var, spec, `16%`,`50%`,`84%`) %>%
    left_join(data_frame(spec = data_list$map %>% unique, taxon = data_list$taxon %>% unique))

phi = fit_AR_summary %>%
    filter(var %in% paste0("phi[",1:7,"]")) %>%
    mutate(spec = 1:7) %>%
    select(var, spec, `16%`,`50%`,`84%`) %>%
    left_join(data_frame(spec = data_list$map %>% unique, taxon = data_list$taxon %>% unique))

# plot coefs
b0 %>%
    ggplot(aes(taxon, `50%`))+
    geom_hline(yintercept = 0, size = 0.3, alpha = 0.5)+
    geom_point(size = 2)+
    geom_errorbar(aes(ymin = `16%`, ymax = `84%`), width = 0)+
    theme_bw()

b1 %>%
    ggplot(aes(taxon, `50%`))+
    geom_hline(yintercept = 0, size = 0.3, alpha = 0.5)+
    geom_point(size = 2)+
    geom_errorbar(aes(ymin = `16%`, ymax = `84%`), width = 0)+
    theme_bw()

phi %>%
    ggplot(aes(taxon, `50%`))+
    geom_hline(yintercept = 0, size = 0.3, alpha = 0.5)+
    geom_point(size = 2)+
    geom_errorbar(aes(ymin = `16%`, ymax = `84%`), width = 0)+
    theme_bw()


# plot predictions
preds = b0 %>%
    expand(taxon, z_midges = -3:1.5) %>%
    left_join(b0 %>% rename(b0 = `50%`) %>% select(taxon, b0)) %>%
    left_join(b1 %>% rename(b1 = `50%`) %>% select(taxon, b1)) %>%
    mutate(z_count = b0 + b1*z_midges)
myv_arth2 %>%
    ggplot(aes(z_midges, z_count))+
    facet_wrap(~taxon)+
    geom_point(alpha = 0.4, size = 1)+
    geom_line(data = preds, color = "red", size = 0.8)+
    geom_smooth(method="lm", formula = y ~ x, se=F, color="blue", linetype = 2)+
    theme_bw()

# fixed parameters by step
fixed_par_v = c("sig_b0","sig_b1","sig_res")
fixed_pars = rstan::extract(fit_AR, pars=fixed_par_v) %>%
    lapply(as_data_frame) %>%
    bind_cols() %>%
    mutate(chain = rep(1:chains, each = iter/2), step = rep(c(1:(iter/2)), chains))
names(fixed_pars) = c(fixed_par_v,"chain","step")

# examine chains for parameters
fixed_pars %>%
    gather(par, value, -chain, -step) %>%
    filter(par != "lp__") %>%
    ggplot(aes(step, value, color=factor(chain)))+
    facet_wrap(~par, scales="free_y")+
    geom_line(alpha=0.5)+
    theme_bw()

# pairs plot for parameters
ggpairs(fixed_pars %>% select(-chain, -step))

# posterior densities
fixed_pars %>%
    gather(par, value, -chain, -step) %>%
    filter(par != "lp__") %>%
    ggplot(aes(value))+
    facet_wrap(~par, ncol=1)+
    stat_density(alpha=0.5, geom = "line")+
    theme_bw()




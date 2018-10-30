#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)

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

# plot
myv_arth2 %>%
    ggplot(aes(z_midges, z_count))+
    facet_wrap(~taxon)+
    geom_point(alpha = 0.4, size = 1)+
    geom_smooth(method="lm", formula = y ~ x, se=F, color="firebrick")+
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
#========== Fit model
#==========

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

# specificaionts
chains = 4
iter = 2000

# fit model
fit = stan(file = "under_constr/lmm_AR.stan", data = data_list, seed=1, chains = chains, iter = iter)

# summary of fit
fit_summary = summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>%
{as_data_frame(.) %>%
        mutate(var = rownames(summary(fit)$summary))}

# check Rhat & n_eff
# note that initial values beta0, alpha, and rho sample for an extra year
# these extra values do not contribute to the likelihood
fit_summary %>% filter(Rhat > 1.05) %>% select(Rhat, n_eff, var)
fit_summary %>% filter(n_eff < 0.5*(chains*iter/2)) %>% select(Rhat, n_eff, var) %>% arrange(n_eff) %>%
    mutate(eff_frac = n_eff/(chains*iter/2))

# extract coefficients
coefs = fit_summary %>%
    filter(var %in% paste0("b1[",1:7,"]")) %>%
    mutate(spec = 1:7) %>%
    select(var, spec, `16%`,`50%`,`84%`) %>%
    left_join(data_frame(spec = data_list$map %>% unique, taxon = data_list$taxon %>% unique))



# plot coefs
coefs %>%
    ggplot(aes(taxon, `50%`))+
    geom_hline(yintercept = 0, size = 0.3, alpha = 0.5)+
    geom_point(size = 2)+
    geom_errorbar(aes(ymin = `16%`, ymax = `84%`), width = 0)+
    scale_y_continuous(limits=c(-0.25, 0.25))+
    theme_bw()

fit_summary %>%
    filter(var %in% paste0("phi[",1:7,"]")) %>%
    mutate(spec = 1:7) %>%
    select(var, spec, `16%`,`50%`,`84%`) %>%
    left_join(data_frame(spec = data_list$map %>% unique, taxon = data_list$taxon %>% unique)) %>%
    ggplot(aes(taxon, `50%`))+
    geom_hline(yintercept = 0, size = 0.3, alpha = 0.5)+
    geom_point(size = 2)+
    geom_errorbar(aes(ymin = `16%`, ymax = `84%`), width = 0)+
    theme_bw()

test = coefs %>%
    expand(taxon, z_midges = -3:1.5) %>%
    left_join(coefs %>% select(taxon, `50%`)) %>%
    mutate(z_count = `50%`*z_midges)

myv_arth2 %>%
    ggplot(aes(z_midges, z_count))+
    facet_wrap(~taxon)+
    geom_point(alpha = 0.4, size = 1)+
    geom_line(data = test, color = "red", size = 0.8)+
    geom_smooth(method="lm", formula = y ~ x, se=F, color="blue", linetype = 2)+
    theme_bw()





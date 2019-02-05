#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(truncnorm)
library(rstan)

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)




#==========
#========== Simulate data
#==========

# number of taxa
n_taxa = 5

# number of time series
n_ts = 10*n_taxa

# obserations per time series
obs_per = rep(10, n_ts)

# ovreall intercept
alpha0 = 0

# overall slope
alpha1 = 0.25

# standard deviations
sig_res = 0.5
sig_beta0 = 0.1
sig_beta1 = 0.5
sig_phi = 0.75

# random slopes
beta0 = alpha0 + rnorm(n = n_taxa, mean = 0, sd = sig_beta0)
beta1 = alpha1 + rnorm(n = n_taxa, mean = 0, sd = sig_beta1)

# phi's
phis = rtruncnorm(n = n_taxa, a = 0, b = 1, mean = 0.5, sd = sig_phi)
beta1 = alpha1 + phis - 0.5

# taxon
taxon = rep(rep(c(1:n_taxa), each = 10), 10)

# times
time = rep(rep(c(0:9), n_taxa), 10)

# environmental variable
Xx = rep(c(1:10), n_ts) + rnorm(n_ts*10, mean = 0, sd = 1)
xx = (Xx - mean(Xx))/(2*sd(Xx))

y = rep(0, sum(obs_per))
j = 1
for(i in 1:n_ts){
    y[j] = beta0[taxon[j]] + rnorm(n = 1, mean = beta1[taxon[j]]*xx[j], sd = sig_res)

    for(t in (j + 1):(j + obs_per[i] - 1)) {
        y[t] = rnorm(n = 1, sd = sig_res,
                     mean = beta0[taxon[j]] + beta1[taxon[j]]*xx[t] +
                         phis[taxon[j]]^(time[t] - time[t-1])*(y[t-1] -
                                                                   (beta0[taxon[j]] + beta1[taxon[j]]*xx[t-1])))
    }

    j = j + obs_per[i]
}


data = data_frame(taxon = taxon, time = time, series = rep(1:n_ts, each = 10), y = y, xx = xx)

data %>%
    ggplot(aes(time, y, group = series))+
    facet_wrap(~taxon)+
    geom_line()+
    theme_bw()

data %>%
    ggplot(aes(xx, y))+
    facet_wrap(~taxon)+
    geom_point(alpha = 0.5, size = 2.5)+
    theme_bw()



#==========
#========== Fit model
#==========

n_obs = length(y)

x = cbind(rep(1, n_obs), xx)

n_coef = ncol(x)

g_per_ff = c(1,1)

b_groups = cbind(rep(1:n_taxa, 10),rep(1:n_taxa, 10))
b_groups2 = apply(b_groups, 1, function(x){x + cumsum(c(0,lev_per_g[1:(length(lev_per_g )-1)]))}) %>% t()

lev_per_g = c(5, 5)

p_groups = b_groups[,1]
# p_groups = rep(1, 50)



# package data
data_list = list(n_obs = length(y),
                 n_ts = n_ts,
                 obs_per = obs_per,
                 n_coef = n_coef,
                 g_per_ff = g_per_ff,
                 lev_per_g = lev_per_g,
                 b_groups = b_groups2,
                 p_groups = p_groups,
                 y = y,
                 x = x,
                 time = time)

# specificaionts
chains = 1
iter = 1000

# fit_lizard model
fit_lizard = stan(file = "under_constr/lizard.stan", data = data_list, seed=1, chains = chains, iter = iter)

# summary of fit_lizard
fit_lizard_summary = summary(fit_lizard, probs=c(0.025, 0.5, 0.975))$summary %>%
{as_data_frame(.) %>%
        mutate(var = rownames(summary(fit_lizard)$summary))} %>%
    rename(lower = `2.5%`,
           middle = `50%`,
           upper = `97.5%`)

# check Rhat & n_eff
# note that initial values beta0, alpha, and rho sample for an extra year
# these extra values do not contribute to the likelihood
fit_lizard_summary %>% filter(Rhat > 1.05) %>% select(Rhat, n_eff, var)
fit_lizard_summary %>% filter(n_eff < 0.5*(chains*iter/2)) %>% select(Rhat, n_eff, var) %>% arrange(n_eff) %>%
    mutate(eff_frac = n_eff/(chains*iter/2))





#==========
#========== Plot output
#==========


# extract coefficients
coefs = lapply(1:n_coef, function(b){
    df = fit_lizard_summary %>%
        filter(var %in% paste0("beta[",1:n_ts,",",b,"]")) %>%
        mutate(series = 1:n_ts)
}) %>%
    bind_rows() %>%
    mutate(taxon = rep(rep(c(1:n_taxa), 10),2),
           beta = rep(1:2, each = 50)) %>%
    group_by(taxon, beta) %>%
    summarize(lower = unique(lower),
              middle = unique(middle),
              upper = unique(upper))

phi = fit_lizard_summary %>%
    filter(var %in% paste0("phi[",1:n_taxa,"]")) %>%
    mutate(group = 1:n_taxa) %>%
    select(group, lower, middle, upper)

alphas = fit_lizard_summary %>%
    filter(var %in% paste0("alpha[",1:n_coef,"]")) %>%
    mutate(beta = 1:n_coef) %>%
    select(beta, lower,middle,upper)

# plot coefs
coefs %>%
    ggplot(aes(taxon, middle, color))+
    facet_wrap(~beta)+
    geom_point(data = data_frame(beta = rep(1:2, each = 5), taxon = rep(1:5, 2), real = c(beta0, beta1)),
               aes(y= real),
               color = "red", size = 2)+
    geom_point(size = 2)+
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)+
    scale_y_continuous(limits=c(-1.1,1.1))+
    theme_bw()


phi %>%
    ggplot(aes(group, middle))+
    geom_point(data = data_frame(group = 1:5, real = phis),
               aes(y= real),
               color = "red", size = 2)+
    geom_point(size = 2)+
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)+
    theme_bw()

alphas %>%
    ggplot(aes(beta, middle))+
    geom_hline(yintercept = 0, size = 0.3, alpha = 0.5)+
    geom_point(data = data_frame(beta = 1:2, real = c(alpha0,alpha1)),
               aes(y= real),
               color = "red", size = 2)+
    geom_point(size = 2)+
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)+
    scale_y_continuous(limits=c(-1,1))+
    theme_bw()









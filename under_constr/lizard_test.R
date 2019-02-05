#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)

# load data
load("data/myv_arth.rda")
myv_arth = myv_arth %>% tbl_df





#==========
#========== Prepare data
#==========

# transform variables
# define time series
myv_arth2 = myv_arth %>%
    arrange(trans, dist, taxon, year) %>%
    group_by(taxon) %>%
    mutate(y = log1p(count),
           y = (y - mean(y))/(2*sd(y))) %>%
    ungroup() %>%
    mutate(series = paste0(trans, dist, taxon) %>% as.factor() %>% as.integer(),
           midges_z = log1p(midges),
           midges_z = (midges_z - mean(midges_z))/(2*sd(midges_z)),
           time = factor(year, levels = c(1:max(year))) %>% as.numeric(),
           time = time - min(time),
           time_z = (time - mean(time))/(2*sd(time)),
           dist_z = (dist - mean(dist))/(2*sd(dist))) %>%
    arrange(series, year)

# reponse variable
y = myv_arth2$y

# predictor variables
x = model.matrix( ~ 1 + midges_z + time_z + dist_z, myv_arth2)

# observation times
time = myv_arth2$time

# number of observations
n_obs = length(y)

# number of time series
n_ts = myv_arth2$series %>% unique() %>% length()

# observations per time series
obs_per = {myv_arth2 %>% group_by(series) %>% summarize(n = length(y))}$n

# number of coefficients
n_coef = ncol(x)

# groups per fixed effect
g_per_ff = c(2, 1, 1, 1)

# grouping structure for beta's
b_groups_verbose = myv_arth2 %>%
    group_by(series) %>%
    summarize(taxon = unique(taxon),
              trans = unique(trans)) %>%
    ungroup() %>%
    mutate(group1 = as.factor(taxon) %>% as.numeric(),
           group2 = as.factor(trans) %>% as.numeric(),
           group3 = as.factor(taxon) %>% as.numeric(),
           group4 = as.factor(taxon) %>% as.numeric(),
           group5 = as.factor(taxon) %>% as.numeric())
b_groups = b_groups_verbose %>%
    select(group1, group2, group3, group4, group5) %>%
    as.matrix()

# levels per fixed effect
lev_per_g = apply(b_groups, 2, function(x){x %>% unique() %>% length()})

# b_groups2
b_groups2 = apply(b_groups, 1, function(x){x + cumsum(c(0,lev_per_g[1:(length(lev_per_g )-1)]))}) %>% t()

# grouping for phi's
p_groups_verbose = myv_arth2 %>%
    group_by(series) %>%
    summarize(taxon = unique(taxon)) %>%
    ungroup() %>%
    mutate(group = as.factor(taxon) %>% as.numeric())
p_groups = {b_groups_verbose %>%
    select(group1)}$group1

# package data
data_list = list(n_obs = n_obs,
            n_ts = n_ts,
            obs_per = obs_per,
            n_coef = n_coef,
            g_per_ff = g_per_ff,
            lev_per_g = lev_per_g,
            b_groups = b_groups2,
            p_groups = p_groups,
            p_bound = 1,
            y = y,
            x = x,
            time = time)





#==========
#========== Fit model: LMM
#==========

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

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
         mutate(series = 1:n_ts) %>%
         select(var, series, lower, middle, upper) %>%
         left_join(data_frame(series = b_groups_verbose$series,
                             taxon = b_groups_verbose$taxon,
                             trans = b_groups_verbose$trans)) %>%
         select(-var, -series) %>%
         group_by(taxon, trans) %>%
         summarize_all(unique) %>%
         mutate(beta = colnames(x)[b])
}) %>%
    bind_rows()

phi = fit_lizard_summary %>%
    filter(var %in% paste0("phi[",1:7,"]")) %>%
    mutate(group = 1:7) %>%
    select(var, group, lower,middle,upper) %>%
    left_join(data_frame(group = b_groups_verbose$group1,
                         taxon = b_groups_verbose$taxon)) %>%
    select(-var, -group) %>%
    group_by(taxon) %>%
    summarize_all(unique)

alphas = fit_lizard_summary %>%
    filter(var %in% paste0("alpha[",1:4,"]")) %>%
    mutate(group = 1:4) %>%
    select(var, group, lower,middle,upper) %>%
    mutate(beta = colnames(x)[group])

# plot coefs
coefs %>%
    ggplot(aes(taxon, middle, color))+
    facet_wrap(~beta)+
    geom_hline(yintercept = 0, size = 0.3, alpha = 0.5)+
    geom_point(size = 2)+
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)+
    theme_bw()

phi %>%
    ggplot(aes(taxon, middle))+
    geom_hline(yintercept = 0, size = 0.3, alpha = 0.5)+
    geom_point(size = 2)+
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)+
    theme_bw()

alphas %>%
    ggplot(aes(beta, middle))+
    geom_hline(yintercept = 0, size = 0.3, alpha = 0.5)+
    geom_point(size = 2)+
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)+
    theme_bw()

# plot predictions
preds = coefs %>%
    tidyr::expand(taxon = unique(taxon), midges_z = -1.5:0.5) %>%
    left_join(coefs %>% filter(beta == "(Intercept)") %>% rename(b1 = middle) %>% select(taxon, b1)) %>%
    left_join(coefs %>% filter(beta == "midges_z") %>% rename(b2 = middle) %>% select(taxon, b2)) %>%
    mutate(y = b1 + b2*midges_z)

myv_arth2 %>%
    ggplot(aes(midges_z, y))+
    facet_wrap(~taxon)+
    geom_point(alpha = 0.4, size = 1)+
    geom_line(data = preds, color = "red", size = 0.8)+
    geom_smooth(method="lm", formula = y ~ x, se=F, color="blue", linetype = 2)+
    theme_bw()



# plot predictions
preds = coefs %>%
    tidyr::expand(taxon = unique(taxon), time_z = seq(-0.75,0.75,0.1)) %>%
    left_join(coefs %>% filter(beta == "(Intercept)") %>% rename(b1 = middle) %>% select(taxon, b1)) %>%
    left_join(coefs %>% filter(beta == "time_z") %>% rename(b3 = middle) %>% select(taxon, b3)) %>%
    mutate(y = b1 + b3*time_z)
myv_arth2 %>%
    ggplot(aes(time_z, y))+
    facet_wrap(~taxon)+
    geom_point(alpha = 0.4, size = 1)+
    geom_line(data = preds, color = "red", size = 0.8)+
    geom_smooth(method="lm", formula = y ~ x, se=F, color="blue", linetype = 2)+
    theme_bw()

# plot predictions
preds = coefs %>%
    tidyr::expand(taxon = unique(taxon), dist_z = seq(-0.75,1.75,0.1)) %>%
    left_join(coefs %>% filter(beta == "(Intercept)") %>% rename(b1 = middle) %>% select(taxon, b1)) %>%
    left_join(coefs %>% filter(beta == "dist_z") %>% rename(b4 = middle) %>% select(taxon, b4)) %>%
    mutate(y = b1 + b4*dist_z)
myv_arth2 %>%
    ggplot(aes(dist_z, y))+
    facet_wrap(~taxon)+
    geom_point(alpha = 0.4, size = 1)+
    geom_line(data = preds, color = "red", size = 0.8)+
    geom_smooth(method="lm", formula = y ~ x, se=F, color="blue", linetype = 2)+
    theme_bw()

# sd's
fixed_par_v = paste0("sig_beta[",1:5,"]")
fixed_pars = rstan::extract(fit_lizard, pars=fixed_par_v) %>%
    lapply(as_data_frame) %>%
    bind_cols() %>%
    mutate(chain = rep(1:chains, each = iter/2), step = rep(c(1:(iter/2)), chains))
names(fixed_pars) = c(fixed_par_v,"chain","step")
fixed_pars %>%
    gather(par, value, -chain, -step) %>%
    filter(par != "lp__") %>%
    ggplot(aes(value))+
    facet_wrap(~par, ncol=1, scales = "free_y")+
    stat_density(alpha=0.5, geom = "line")+
    scale_x_continuous(limits=c(0,0.7))+
    theme_bw()

# alpha's
fixed_par_v = paste0("alpha[",2:4,"]")
fixed_pars = rstan::extract(fit_lizard, pars=fixed_par_v) %>%
    lapply(as_data_frame) %>%
    bind_cols() %>%
    mutate(chain = rep(1:chains, each = iter/2), step = rep(c(1:(iter/2)), chains))
names(fixed_pars) = c(fixed_par_v,"chain","step")
fixed_pars %>%
    gather(par, value, -chain, -step) %>%
    filter(par != "lp__") %>%
    ggplot(aes(value))+
    facet_wrap(~par, ncol=1)+
    stat_density(alpha=0.5, geom = "line")+
    theme_bw()

library(lme4)
m = lmer(y ~ midges_z + dist_z + time_z + (midges_z + dist_z + time_z|taxon), myv_arth2, REML = F,
         lmerControl(optimizer ="bobyqa"))
coef(m)$taxon[,"midges_z"]

vars = c("midges_z","time_z","dist_z")
var = vars[1]
comp = coefs %>% filter(beta == var) %>% group_by(taxon) %>% summarize(middle = unique(middle)) %>%
mutate(lmer = coef(m)$taxon[,var])
comp
plot(lmer ~ middle, comp)
abline(a = 0, b = 1)

library(car)
Anova(m)
summary(m)

m = lm(y ~ midges_z*taxon + dist_z*taxon + time_z*taxon, myv_arth2)
Anova(m, type = 2)

fixed_pars$`sig_beta[3]` %>% median()
fixed_pars$`sig_beta[4]` %>% median()
fixed_pars$`sig_beta[5]` %>% median()

summary(m)

apply(fixed_pars %>% select(-chain, -step), 2, sd)

rstan::extract(fit_lizard, pars="sig_phi")$sig_res %>% median

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

# prepare data
myv_arth = myv_arth %>% tbl_df
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
           dist_z = (dist - mean(dist))/(2*sd(dist)),
           trans = factor(trans),
           distf = factor(dist),
           taxon = factor(taxon)) %>%
    arrange(series, year)





#==========
#========== Full model: Fit
#==========

# fit model
m = lizfit(formula = y ~ midges_z + time_z + dist_z + (1 + midges_z + time_z + dist_z | taxon),
           time_form = ~  time | trans + distf + taxon,
           ar_form = ~ taxon,
           data = myv_arth2)

# clean fit
m_clean = summary(m$stan, probs=c(0.16, 0.5, 0.84))$summary %>%
    {as_data_frame(.) %>%
            mutate(var = rownames(summary(m$stan)$summary))} %>%
    rename(lower = `16%`,
           middle = `50%`,
           upper = `84%`) %>%
    select(var, lower, middle, upper) %>%
    mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]))

# extract alphas
alphas = m_clean %>%
    filter(name == "alpha") %>%
    mutate(coef = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
           coef = ifelse(coef==1, "intercept",
                         ifelse(coef==2,"midges_z",
                                ifelse(coef==3,"time_z",
                                       ifelse(coef==4,"dist_z",NA)))),
           coef = factor(coef, levels=c("intercept","midges_z","time_z","dist_z")))

# extract coefficients
coefs = m_clean %>%
    filter(name == "beta") %>%
    mutate(series = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
           coef = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
           coef = ifelse(coef==1, "intercept",
                         ifelse(coef==2,"midges_z",
                                ifelse(coef==3,"time_z",
                                       ifelse(coef==4,"dist_z",NA))))) %>%
    left_join(myv_arth2 %>%
                  group_by(series) %>%
                  summarize(taxon = unique(taxon))) %>%
    mutate(coef = factor(coef, levels=c("intercept","midges_z","time_z","dist_z")))

# extract phis
phis = m_clean %>%
    filter(name == "phi") %>%
    mutate(series = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
    left_join(myv_arth2 %>%
                  group_by(series) %>%
                  summarize(taxon = unique(taxon)))





#==========
#========== Full model: Plot
#==========

# coefs
alphas %>%
    ggplot(aes(coef, middle, color))+
    geom_hline(yintercept = 0, size = 0.3, alpha = 0.5)+
    geom_point(size = 2)+
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)+
    theme_bw()

# coefs
coefs %>%
    ggplot(aes(taxon, middle, color))+
    facet_wrap(~coef)+
    geom_hline(yintercept = 0, size = 0.3, alpha = 0.5)+
    geom_point(size = 2)+
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)+
    theme_bw()

# phis
phis %>%
    ggplot(aes(taxon, middle, color))+
    geom_hline(yintercept = 0, size = 0.3, alpha = 0.5)+
    geom_point(size = 2)+
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)+
    theme_bw()





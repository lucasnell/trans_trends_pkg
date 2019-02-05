# load packages and data
library(lizard)
load("data/myv_arth.rda")

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

# fit model
m = lizfit(formula = y ~ midges_z + time_z + dist_z + (1 + midges_z | taxon),
           time_form = ~  time | trans + distf + taxon,
           ar_form = ~ taxon,
           data = myv_arth2)


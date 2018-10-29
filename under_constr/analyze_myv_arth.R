#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
load("data/myv_arth.rda")





#==========
#========== Prepare data
#==========

# log1p transform, then z-score midges
myv_arth2 = myv_arth %>%
    tbl_df() %>%
    mutate(log_count = log1p(count),
           log_midges = log1p(midges),
           z_midges = (log_midges - mean(log_midges, na.rm=T))/sd(log_midges, na.rm=T)) %>%
    # group by taxon and them z-score count within taxon
    group_by(taxon) %>%
    mutate(z_count = (log_count - mean(log_count, na.rm=T))/sd(log_count, na.rm=T))

# plot
myv_arth2 %>%
    ggplot(aes(z_midges, z_count))+
    facet_wrap(~taxon)+
    geom_point(alpha = 0.4, size = 1)+
    geom_smooth(method="lm", formula = y ~ x, se=F, color="firebrick")+
    theme_bw()

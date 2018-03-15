#==================================================
#Preliminaries
#==================================================

# Build: Install and Restart

# library(repts)

# source(".Rprofile")

library(tidyverse)




# ==================================================
# Simulations
# ==================================================

# Simulate time seires
cube <- do.call(sim_pops_ar,
                generate_data(n_time = 30, n_loc = 16, n_spp = 5, corr_method = "phylo"))


expr_ <- setNames(rlang::syms(paste0("V", 1:(ncol(cube)+1))),
                  c("loc", paste0("sp", 1:ncol(cube))))

cube_df <- melt_cube(cube) %>% 
    as_data_frame() %>% 
    select(!!!expr_) %>% 
    mutate(loc = factor(as.integer(loc)), 
           time = rep(1:sum(loc == 1), length(levels(loc))))


cube_df %>%
    gather('species', 'N', -loc, -time, factor_key = TRUE) %>% 
    ggplot(aes(time, N, color = species)) +
    theme_classic() +
    geom_line() +
    facet_wrap(~ loc, nrow = 2)






#==================================================
#Preliminaries
#==================================================

# Build: Install and Restart

# library(lizard)

# source(".Rprofile")

library(tidyverse)




# ==================================================
# Simulations
# ==================================================

# Generate input data
gen_data = generate_data(n_time = 30, n_loc = 10, n_spp = 5,
              corr_method = "none", sigma_obs = 1.5,
              mean_b0 = 0, sigma_b0 = 1, mean_b1 = 0,
              sigma_b1 = 4, mean_rho = 0, sigma_rho = 3,
              sigma_eps = 1)

# Simulate communities
cube <- do.call(sim_pops_ar, gen_data)


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






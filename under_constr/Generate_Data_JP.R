#==================================================
#Preliminaries
#==================================================

# Build: Install and Restart

# library(ReplicateTimeseries)
# devtools::load_all()

# source(".Rprofile")

library(tidyverse)





#==================================================
#Define parameters
#==================================================

# Dimensions
n_time = 30
n_loc = 16
n_spp = 5

# Means
mean_b0 = 5
mean_b1 = 0.5
mean_rho = 0 #inverse logit scale

# Standard Deviations
sigma_b0 = 0.1
sigma_b1 = 0.1
sigma_rho = 0.5 #inverse logit scale
sigma_eps = 0.1
sigma_obs = 0





#==================================================
#Simulations
#==================================================

# Set up matries
# Environmental variables
X_ <- matrix(rnorm(n = n_time, mean = 0 , sd = 1), 
             nrow = n_time, ncol = n_loc)
# Mean abundance in average environment
b0_mat_ <- matrix(rnorm(n = n_spp, mean = mean_b0, sd  = sigma_b0), 
                  nrow = n_spp, ncol = n_loc)
# Response to environmental variables
b1_mat_ <- matrix(rnorm(n = n_spp, mean = mean_b1, sd = sigma_b1), 
                  nrow = n_spp, ncol = n_loc)
# Autoregressive paramer
rho_mat_inv_ <- matrix(rnorm(n = n_spp, mean = mean_rho, sd = sigma_rho), 
                       nrow = n_spp, ncol = n_loc)
rho_mat_ = exp(rho_mat_inv_)/(1 + exp(rho_mat_inv_)) 
# Observation error
obs_sigma_ <- rep(sigma_obs, n_spp)
# Initial population size (set to b0)
N0_mat_ <- b0_mat_


# Specify epsilons: phyogenetic correlation
# vcv1_ <- ape::vcv(ape::rcoal(n_spp))
# dimnames(vcv1_) <- NULL
# vcv_ = (sigma_eps^2)*vcv1_/diag(vcv1_)[1]
# vcv_cube_ <- replicate(n_loc, vcv_)

# Specify epsilons: idenpendent responses
vcv_ = matrix(0, nrow=n_spp, ncol=n_spp)
diag(vcv_) = sigma_eps^2 
vcv_cube_ <- replicate(n_loc, vcv_)

# Simulate time seires
cube <- sim_pops_ar(X = X_,
                    N0_mat = N0_mat_,
                    b0_mat = b0_mat_,
                    b1_mat = b1_mat_,
                    rho_mat = rho_mat_,
                    vcv_cube = vcv_cube_,
                    obs_sigma = obs_sigma_)


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






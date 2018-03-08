# library(ReplicateTimeseries)
# devtools::load_all()

# source(".Rprofile")
z <- sim_pops(n_gen = 100, N0 = rep(100, 10), r = rep(0.5, 10), alpha = rep(1e-3, 10),
              sigma = 0.01, seed = 496)
# Should max out at r / alpha
matplot(z, type = 'l')
max(z)

# log(N) for one population using an AR process
ar_N <- sim_pop_ar(X = 0.1 * 1:100, N0 = 5, b0 = 5, b1 = 0.5, rho = 0.2, 
                   sigma = 1, seed = 1)
plot(ar_N, type = 'l', ylab = 'N')




# 100 time points, 10 locations, 5 species
n_time = 100; n_locs = 10; n_spp = 5
set.seed(1)
X_ <- matrix(rnorm(n_time * n_locs), n_time, n_locs)
N0_mat_ <- matrix(log(100), n_spp, n_locs)
b0_mat_ <- matrix(log(100), n_spp, n_locs)
b1_mat_ <- matrix(0.5, n_spp, n_locs)
rho_mat_ <- matrix(0.2, n_spp, n_locs)
vcv_ <- ape::vcv(ape::rcoal(n_spp))
dimnames(vcv_) <- NULL
vcv_cube_ <- replicate(n_locs, vcv_)
obs_sigma_ <- runif(n_spp)

cube <- sim_pops_ar(X = X_,
                    N0_mat = N0_mat_,
                    b0_mat = b0_mat_,
                    b1_mat = b1_mat_,
                    rho_mat = rho_mat_,
                    vcv_cube = vcv_cube_,
                    obs_sigma = obs_sigma_)

library(tidyverse)

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


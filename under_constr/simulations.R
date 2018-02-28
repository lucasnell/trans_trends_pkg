# library(ReplicateTimeseries)
devtools::load_all()

# source(".Rprofile")
z <- sim_pops(n_gen = 100, N0 = rep(100, 10), r = rep(0.5, 10), alpha = rep(1e-3, 10),
              sigma = 0.01, seed = 496)
# Should max out at r / alpha
matplot(z, type = 'l')
max(z)

#' Melt a cube into a single data frame.
#'
#' @param C Three-dimensional array that you want to melt into a two-dimensional
#'     data frame.
#'
#' @return
#'
#' @noRd
#'
melt_cube <- function(C) {
    C_rows = nrow(C)
    M = matrix(0, C_rows * dim(C)[3], ncol(C) + 1)
    for (i in 1:(dim(C)[3])) {
        rows_ = ((i-1) * C_rows + 1):(i * C_rows)
        M[rows_, 1] = i
        M[rows_, 2:ncol(C)] = C[,,i]
    }
    return(as.data.frame(M))
}


#' Generate data for simulations.
#'
#' Generate multi-location, multi-species time series data.
#'
#'
#' @param n_time Number of time steps.
#' @param n_loc Number of locations.
#' @param n_spp Number of species.
#' @param mean_b0 Mean for the b0 parameter relating X to N.
#' @param mean_b1 Mean for the b1 parameter relating X to N.
#' @param mean_rho Mean for the rho parameter relating X to N.
#'     This parameter is on the inverse logit scale.
#' @param sigma_b0 Standard deviation for the b0 parameter relating X to N.
#' @param sigma_b1 Standard deviation for the b1 parameter relating X to N.
#' @param sigma_rho Standard deviation for the rho parameter relating X to N.
#'     This parameter is on the inverse logit scale.
#' @param sigma_eps Standard deviation for the epsilon parameter.
#' @param sigma_obs Standard deviation for observation error. Defaults to 0.
#' @param corr_method Method for determining correlations between species.
#'     Options are "none", "phylo", or "random". Defaults to "none".
#'
#'
#' @export
#'
#'
generate_data <- function(n_time,
                          n_loc,
                          n_spp,
                          mean_b0,
                          mean_b1,
                          mean_rho,
                          sigma_b0,
                          sigma_b1,
                          sigma_rho,
                          sigma_eps,
                          sigma_obs = 0,
                          corr_method = "none") {

    corr_method = match.arg(corr_method, c("none", "phylo", "random"))

    # Set up matrices
    # Environmental variables
    X_ = matrix(0, n_time, n_loc)
    rnd = cbind(rnorm(n_time))
    for (i in 1:n_loc) X_[,i] = rnd
    # Mean abundance in average environment
    b0_mat_ = matrix(0, n_spp, n_loc)
    rnd = cbind(rnorm(n_spp, mean_b0, sigma_b0))
    for (i in 1:n_loc) b0_mat_[,i] = rnd
    # Response to environmental variables
    b1_mat_ = matrix(0, n_spp, n_loc)
    rnd = cbind(rnorm(n_spp, mean_b1, sigma_b1))
    for (i in 1:n_loc) b1_mat_[,i] = rnd
    # Autoregressive parameter
    rho_mat_ = matrix(0, n_spp, n_loc)
    rnd = cbind(rnorm(n_spp, mean_rho, sigma_rho))
    rnd = exp(rnd) / (1 + exp(rnd))
    for (i in 1:n_loc) rho_mat_[,i] = rnd
    # Observation error
    obs_sigma_ = cbind(rep(sigma_obs, n_spp))

    # Initial population size (set to b0)
    N0_mat_ = b0_mat_

    vcv_ = diag(n_spp)
    vcv_cube_ = array(0, dim = c(n_spp, n_spp, n_loc))

    # LEFT OFF -----


    # if (corr_method == "phylo") {
    #     Environment ape = Environment::namespace_env("ape");
    #     Function rcoal = ape["rcoal"];
    #     Function vcv_phylo = ape["vcv.phylo"];
    #     SEXP phylo = rcoal(n_spp);
    #     SEXP vcv_SEXP_ = vcv_phylo(phylo);
    #     vcv_ = as<arma::mat>(vcv_SEXP_);
    #     vcv_ /= vcv_.diag()(0);
    # } else if (corr_method == "random") {
    #     uint n_corrs = ((n_spp - 1) / 2) * (1 + (n_spp - 1));
    #     arma::vec rnd_phy = as<arma::vec>(Rcpp::runif(n_corrs, -1, 1));
    #     for (uint i = 0, rnd_i = 0; i < n_spp; i++) {
    #         for (uint j = i+1; j < n_spp; j++, rnd_i++) {
    #             vcv_(i,j) = rnd_phy(rnd_i);
    #             vcv_(j,i) = rnd_phy(rnd_i);
    #         }
    #     }
    # }
    # # Going from correlations to covariances:
    # vcv_ *= (sigma_eps * sigma_eps);
    # vcv_cube_.each_slice() += vcv_;
    #
    # List out = List::create(_["X"] = X_,
    # _["N0_mat"] = N0_mat_,
    # _["b0_mat"] = b0_mat_,
    # _["b1_mat"] = b1_mat_,
    # _["rho_mat"] = rho_mat_,
    # _["vcv_cube"] = vcv_cube_,
    # _["obs_sigma"] = obs_sigma_);
    #
    # return out;

    return(0)
}

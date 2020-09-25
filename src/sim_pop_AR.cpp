#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <random>  // normal distribution
#include <vector>  // vector class
#ifdef _OPENMP
#include <omp.h>  // OpenMP
#endif

#include "armmr_types.h"
#include "pcg.h"  // pcg seeding

using namespace Rcpp;




// /*
//  One population simulated using AR1 process.
//  N in this case is already logged
//  */
// //[[Rcpp::export]]
// arma::vec sim_pop_ar(const arma::vec& X, const double& N0,
//                      const double& b0, const double& b1, const double& rho,
//                      const double& sigma) {
//     uint32 n_gen = X.n_elem - 1;
//     arma::vec N(n_gen + 1);
//     N(0) = N0;
//
//     std::normal_distribution<double> rnorm_distr(0.0, sigma);
//     pcg32 engine = seeded_pcg();
//
//     for (uint32 t = 0; t < n_gen; t++) {
//         N(t+1) = b0 + b1 * X(t+1) + rho * ( N(t) - b0 - b1 * X(t) );
//         N(t+1) += rnorm_distr(engine);
//     }
//
//     return N;
// }






/*
 Check that input parameters have the correct dimensionality.
 For more info, see description of `sim_pops_ar` below.
 */
void check_dims(const arma::mat& X, const arma::mat& N0_mat,
                const arma::mat& b0_mat, const arma::mat& b1_mat,
                const arma::mat& rho_mat, const arma::cube& vcv_cube) {

    uint32 n_locs = X.n_cols;
    if (N0_mat.n_cols != n_locs) {
        stop("N0_mat should have the same number of columns as does X");
    }
    if (b0_mat.n_cols != n_locs) {
        stop("b0_mat should have the same number of columns as does X");
    }
    if (b1_mat.n_cols != n_locs) {
        stop("b1_mat should have the same number of columns as does X");
    }
    if (rho_mat.n_cols != n_locs) {
        stop("rho_mat should have the same number of columns as does X");
    }
    if (vcv_cube.n_slices != n_locs) {
        stop("vcv_cube should have the same number of slices as X has columns");
    }

    uint32 n_spp = N0_mat.n_rows;
    if (N0_mat.n_rows != n_spp) {
        stop("N0_mat should have the same number of rows as does N0");
    }
    if (b0_mat.n_rows != n_spp) {
        stop("b0_mat should have the same number of rows as does N0");
    }
    if (b1_mat.n_rows != n_spp) {
        stop("b1_mat should have the same number of rows as does N0");
    }
    if (rho_mat.n_rows != n_spp) {
        stop("rho_mat should have the same number of rows as does N0");
    }
    if (vcv_cube.n_rows != n_spp || vcv_cube.n_cols != n_spp) {
        stop("vcv_cube should have the same number of rows and columns as N0 has rows");
    }

    return;
}


/*
 Turn variance-covariance matrix into matrix that can be used to generate random
 normal deviates with same covariance structure

 "a vector of independent normal random variables,
  when multiplied by the transpose of the Cholesky deposition of [vcv] will
  have covariance matrix equal to [vcv]."
 */
arma::cube make_chol_decomp(const arma::cube& vcv_cube) {

    arma::cube chol_decomp(vcv_cube);
    for (uint32 i = 0; i < vcv_cube.n_slices; i++) {
        chol_decomp.slice(i) = arma::chol(vcv_cube.slice(i)).t();
    }

    return chol_decomp;
}




//' Multiple populations simulated using AR1 process.
//'
//' Input and output N values are logged.
//' All input matrices other than `X` and `vcv_cube` should have rows associated
//' with a given species and columns associated with a given location.
//' See descriptions for `X` and `vcv_cube`.
//'
//' @param X Matrix of environmental variable.
//'     It should have rows associated with a given time point and
//'     columns associated with a given location.
//' @param N0_mat Matrix of starting population abundances (`log(# individuals)`)
//'     by species and location.
//' @param b0_mat Matrix of \eqn{\beta_0} (the population-abundance intercept) values
//'     by species and location.
//' @param b1_mat Matrix of \eqn{\beta_1} (the effect of \eqn{X} on \eqn{N}) values
//'     by species and location.
//' @param rho_mat Matrix of growth rates by species and location.
//' @param vcv_cube Cube representing variance-covariance matrices for process error
//'     among species, one matrix for each location.
//'     It should have rows and columns associated with a given species,
//'     and slices associated with a given location.
//' @param obs_sigma Vector of standard deviations of observation error for each species.
//' @param n_cores Number of cores to use. Defaults to 1.
//'
//'
//'
//' @return A 3-dimensional array.
//' The output will have rows associated with a given time point,
//' columns associated with a given species, and
//' slices associated with a given location.
//'
//'
//'
//'
//' @export
//'
//' @examples
//' X <- matrix(rlnorm(20), 10)
//' N0 <- matrix(rep(log(10), 6), 3, 2)
//' b0 <- matrix(rep(log(100), 6), 3, 2)
//' b1 <- matrix(rep(0.1, 6), 3, 2)
//' rho <- matrix(rep(0.2, 6), 3, 2)
//' vcv <- diag(3)
//' vcv[lower.tri(vcv)] <- vcv[upper.tri(vcv)] <- 0.1
//' vcv <- array(vcv, dim = c(3, 3, 2))
//' obs <- rep(0.1, 3)
//' sim_pops_ar(X, N0, b0, b1, rho, vcv, obs)
//'
//[[Rcpp::export]]
arma::cube sim_pops_ar(const arma::mat& X, const arma::mat& N0_mat,
                       const arma::mat& b0_mat, const arma::mat& b1_mat,
                       const arma::mat& rho_mat, const arma::cube& vcv_cube,
                       const arma::vec& obs_sigma, const uint32& n_cores = 1) {

    check_dims(X, N0_mat, b0_mat, b1_mat, rho_mat, vcv_cube);

    // For random number generator
    const std::vector<std::vector<uint64>> seeds = mc_seeds(n_cores);
    // For turning ~N(0,1) to multivariate with given covariance matrix
    const arma::cube chol_decomp(make_chol_decomp(vcv_cube));

    const uint32 n_time(X.n_rows);
    const uint32 n_locs(X.n_cols);
    const uint32 n_spp(N0_mat.n_rows);

    arma::cube N(n_time, n_spp, n_locs);
    N.tube(arma::span(0), arma::span()) = N0_mat;

    #ifdef _OPENMP
    #pragma omp parallel default(shared) num_threads(n_cores) if(n_cores > 1)
    {
    #endif

    std::vector<uint64> active_seeds;

    // Write the active seed per core or just write one of the seeds.
    #ifdef _OPENMP
    uint32 active_thread = omp_get_thread_num();
    active_seeds = seeds[active_thread];
    #else
    active_seeds = seeds[0];
    #endif

    pcg32 engine = seeded_pcg(active_seeds);
    // Random normal distribution:
    std::normal_distribution<double> rnorm_distr(0.0, 1.0);

    // Parallelize the Loop
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (uint32 loc = 0; loc < n_locs; loc++) {
        const arma::mat& cd_loc(chol_decomp.slice(loc));
        const arma::vec& b0s(b0_mat.col(loc));
        const arma::vec& b1s(b1_mat.col(loc));
        const arma::vec& rhos(rho_mat.col(loc));
        const arma::vec& Xs(X.col(loc));
        arma::mat& Ns(N.slice(loc));
        for (uint32 t = 0; t < n_time - 1; t++) {
            arma::vec proc_rnd(n_spp);
            arma::vec obs_rnd(n_spp);
            for (uint32 sp = 0; sp < n_spp; sp++) {
                Ns(t+1,sp) = b0s(sp) + b1s(sp) * Xs(t+1) + rhos(sp) * (
                    Ns(t,sp) - b0s(sp) - b1s(sp) * Xs(t));
                proc_rnd(sp) = rnorm_distr(engine);
                obs_rnd(sp) = rnorm_distr(engine) * obs_sigma(sp);
            }
            // Generate variance and covariance:
            proc_rnd = cd_loc * proc_rnd;
            Ns.row(t+1) += proc_rnd.t();
            Ns.row(t+1) += obs_rnd.t();
        }
    }

    #ifdef _OPENMP
    }
    #endif


    return N;
}



//' Melt a cube into a single data frame.
//'
//' @param C Three-dimensional array that you want to melt into a two-dimensional
//'     data frame.
//'
//' @return A melted data frame.
//'
//' @noRd
//'
//[[Rcpp::export]]
DataFrame melt_cube(const arma::cube& C) {

    uint32 C_rows = C.n_rows;
    arma::mat M(C.n_rows * C.n_slices, C.n_cols + 1);
    for (uint32 i = 0; i < C.n_slices; i++) {
        M(arma::span(i * C_rows, (i+1) * C_rows - 1), arma::span(0)).fill(i+1);
        M(arma::span(i * C_rows, (i+1) * C_rows - 1),
          arma::span(1, C.n_cols)) = C.slice(i);
    }

    return M;
}


//' Generate parameter values for simulations.
//'
//' Generate parameter values to simulate multi-location, multi-species time series data.
//'
//'
//' @param n_time Number of time steps.
//' @param n_loc Number of locations.
//' @param n_spp Number of species.
//' @param mean_b0 Mean for the b0 parameter relating X to N.
//' @param mean_b1 Mean for the b1 parameter relating X to N.
//' @param mean_rho Mean for the rho parameter relating X to N.
//'     This parameter is on the inverse logit scale.
//' @param sigma_b0 Standard deviation for the b0 parameter relating X to N.
//' @param sigma_b1 Standard deviation for the b1 parameter relating X to N.
//' @param sigma_rho Standard deviation for the rho parameter relating X to N.
//'     This parameter is on the inverse logit scale.
//' @param sigma_eps Standard deviation for the epsilon parameter.
//' @param sigma_obs Standard deviation for observation error. Defaults to 0.
//' @param corr_method Method for determining correlations between species.
//'     Options are "none", "phylo", or "random". Defaults to "none".
//'
//'
//' @export
//'
//' @examples
//' generate_pars(10, 2, 3,
//'               mean_b0 = log(100),
//'               mean_b1 = 0.1,
//'               mean_rho = 0.25,
//'               sigma_b0 = 0.1,
//'               sigma_b1 = 0.1,
//'               sigma_rho = 0.1,
//'               sigma_eps = 0.1)
//'
//'
//[[Rcpp::export]]
List generate_pars(const uint32& n_time,
                   const uint32& n_loc,
                   const uint32& n_spp,
                   const double& mean_b0,
                   const double& mean_b1,
                   const double& mean_rho,
                   const double& sigma_b0,
                   const double& sigma_b1,
                   const double& sigma_rho,
                   const double& sigma_eps,
                   const double& sigma_obs = 0,
                   const std::string& corr_method = "none") {

    std::normal_distribution<double> rnorm_distr(0.0, 1.0);
    pcg32 eng = seeded_pcg();

    // Set up matrices
    // Environmental variables
    arma::mat X_(n_time, n_loc, arma::fill::zeros);
    arma::vec rnd(n_time);
    for (double& r : rnd) r = rnorm_distr(eng);
    X_.each_col() += rnd;
    // Mean abundance in average environment
    arma::mat b0_mat_(n_spp, n_loc, arma::fill::zeros);
    rnd.set_size(n_spp);
    for (double& r : rnd) r = rnorm_distr(eng) * sigma_b0 + mean_b0;
    b0_mat_.each_col() += rnd;
    // Response to environmental variables
    arma::mat b1_mat_(n_spp, n_loc, arma::fill::zeros);
    for (double& r : rnd) r = rnorm_distr(eng) * sigma_b1 + mean_b1;
    b1_mat_.each_col() += rnd;
    // Autoregressive parameter
    arma::mat rho_mat_(n_spp, n_loc, arma::fill::zeros);
    for (double& r : rnd) r = rnorm_distr(eng) * sigma_rho + mean_rho;
    rnd = arma::exp(rnd) / (1 + arma::exp(rnd));
    rho_mat_.each_col() += rnd;
    // Observation error
    arma::vec obs_sigma_(n_spp);
    obs_sigma_.fill(sigma_obs);

    // Initial population size (set to b0)
    arma::mat N0_mat_(b0_mat_);

    arma::mat vcv_(n_spp, n_spp, arma::fill::eye);
    arma::cube vcv_cube_(n_spp, n_spp, n_loc, arma::fill::zeros);
    if (corr_method == "phylo") {
        Environment ape = Environment::namespace_env("ape");
        Function rcoal = ape["rcoal"];
        Function vcv_phylo = ape["vcv.phylo"];
        SEXP phylo = rcoal(n_spp);
        SEXP vcv_SEXP_ = vcv_phylo(phylo);
        vcv_ = as<arma::mat>(vcv_SEXP_);
        vcv_ /= vcv_.diag()(0);
    } else if (corr_method == "random") {
        uint32 n_corrs = ((n_spp - 1) / 2) * (1 + (n_spp - 1));
        arma::vec rnd_phy(n_corrs);
        for (double& r : rnd_phy) r = runif_ab(eng, -1, 1);
        for (uint32 i = 0, rnd_i = 0; i < n_spp; i++) {
            for (uint32 j = i+1; j < n_spp; j++, rnd_i++) {
                vcv_(i,j) = rnd_phy(rnd_i);
                vcv_(j,i) = rnd_phy(rnd_i);
            }
        }
    } else if (corr_method == "none") {
        ;
    } else {
        stop("corr_method must be none, phylo, or random.");
    }
    // Going from correlations to covariances:
    vcv_ *= (sigma_eps * sigma_eps);
    vcv_cube_.each_slice() += vcv_;

    List out = List::create(_["X"] = X_,
                            _["N0_mat"] = N0_mat_,
                            _["b0_mat"] = b0_mat_,
                            _["b1_mat"] = b1_mat_,
                            _["rho_mat"] = rho_mat_,
                            _["vcv_cube"] = vcv_cube_,
                            _["obs_sigma"] = obs_sigma_);

    return out;
}

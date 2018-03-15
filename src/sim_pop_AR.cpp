#include <RcppArmadillo.h>
#include <sitmo.h>  // sitmo::prng_engine
#include <random>  // normal distribution
#include <vector>  // vector class
#ifdef _OPENMP
#include <omp.h>  // OpenMP
#endif

#include "repts_types.h"

using namespace Rcpp;


/*
 One population simulated using AR1 process.
 N in this case is already logged
 */
//[[Rcpp::export]]
arma::vec sim_pop_ar(const arma::vec& X, const double& N0,
                     const double& b0, const double& b1, const double& rho,
                     const double& sigma, const uint& seed) {
    uint n_gen = X.n_elem - 1;
    arma::vec N(n_gen + 1);
    N(0) = N0;
    
    std::normal_distribution<double> rnorm_distr(0.0, sigma);
    sitmo::prng_engine engine(seed);
    
    for (uint t = 0; t < n_gen; t++) {
        N(t+1) = b0 + b1 * X(t+1) + rho * ( N(t) - b0 - b1 * X(t) );
        N(t+1) += rnorm_distr(engine);
    }
    
    return N;
}





/*
 Get one N value based on input parameters.
 Does NOT include process or observation error.
 */
inline double N_t1(const double& b0, const double& b1,
                   const double& rho,
                   const double& Xt, const double& Xt1,
                   const double& Nt) {
    double Nt1 = b0 + b1 * Xt1 + rho * ( Nt - b0 - b1 * Xt );
    return Nt1;
}


/*
 Check that input parameters have the correct dimensionality.
 For more info, see description of `sim_pops_ar` below.
 */
void check_dims(const arma::mat& X, const arma::mat& N0_mat,
                const arma::mat& b0_mat, const arma::mat& b1_mat,
                const arma::mat& rho_mat, const arma::cube& vcv_cube) {
    
    uint n_locs = X.n_cols;
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
    
    uint n_spp = N0_mat.n_rows;
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
    for (uint i = 0; i < vcv_cube.n_slices; i++) {
        chol_decomp.slice(i) = arma::chol(vcv_cube.slice(i)).t();
    }
    
    return chol_decomp;
}



/*
 Multiple populations simulated using AR1 process.
 Input and output N values are logged.
 The `X` matrix should have rows associated with a given time point and
 columns associated with a given location.
 The `vcv_cube` cube should have rows and columns associated with a given species,
 and slices associated with a given location.
 All other input matrices should have rows associated with a given species and
 columns associated with a given location.
 `obs_sigma` contains the standard deviations of observation error for each species.
 The output array will have rows associated with a given time point,
 columns associated with a given species, and 
 slices associated with a given location.
 */
//[[Rcpp::export]]
arma::cube sim_pops_ar(const arma::mat& X, const arma::mat& N0_mat,
                       const arma::mat& b0_mat, const arma::mat& b1_mat,
                       const arma::mat& rho_mat, const arma::cube& vcv_cube,
                       const arma::vec& obs_sigma, const uint& n_cores = 1) {
    
    check_dims(X, N0_mat, b0_mat, b1_mat, rho_mat, vcv_cube);
    
    // For random number generator
    const std::vector<uint> seeds(as<std::vector<uint>>(runif(n_cores, 0, 2147483647)));
    // For turning ~N(0,1) to multivariate with given covariance matrix
    const arma::cube chol_decomp(make_chol_decomp(vcv_cube));

    const uint n_time(X.n_rows);
    const uint n_locs(X.n_cols);
    const uint n_spp(N0_mat.n_rows);
    
    arma::cube N(n_time, n_spp, n_locs);
    N.tube(arma::span(0), arma::span()) = N0_mat;
    
    #ifdef _OPENMP
    #pragma omp parallel default(shared) num_threads(n_cores) if(n_cores > 1)
    {
    #endif
    
    uint active_seed;
    
    // Write the active seed per core or just write one of the seeds.
    #ifdef _OPENMP
    uint active_thread = omp_get_thread_num();
    active_seed = seeds[active_thread];
    #else
    active_seed = seeds[0];
    #endif
    
    sitmo::prng_engine engine(active_seed);
    // Random normal distribution:
    std::normal_distribution<double> rnorm_distr(0.0, 1.0);
    
    // Parallelize the Loop
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (uint loc = 0; loc < n_locs; loc++) {
        const arma::mat& cd_loc(chol_decomp.slice(loc));
        const arma::vec& b0s(b0_mat.col(loc));
        const arma::vec& b1s(b1_mat.col(loc));
        const arma::vec& rhos(rho_mat.col(loc));
        const arma::vec& Xs(X.col(loc));
        arma::mat& Ns(N.slice(loc));
        for (uint t = 0; t < n_time - 1; t++) {
            arma::vec rnd(n_spp);
            arma::vec obs_rnd(n_spp);
            for (uint sp = 0; sp < n_spp; sp++) {
                Ns(t+1,sp) = N_t1(b0s(sp), b1s(sp), rhos(sp), Xs(t), Xs(t+1), Ns(t,sp));
                rnd(sp) = rnorm_distr(engine);
                obs_rnd(sp) = rnorm_distr(engine) * obs_sigma(sp);
            }
            // Generate variance and covariance:
            rnd = cd_loc * rnd;
            Ns.row(t+1) += rnd.t();
            Ns.row(t+1) += obs_rnd.t();
        }
    }
    
    #ifdef _OPENMP
    }
    #endif

    
    return N;
}


/*
 Melt a cube into a single data frame.
 */
//[[Rcpp::export]]
DataFrame melt_cube(arma::cube C) {
    
    uint C_rows = C.n_rows;
    arma::mat M(C.n_rows * C.n_slices, C.n_cols + 1);
    for (uint i = 0; i < C.n_slices; i++) {
        M(arma::span(i * C_rows, (i+1) * C_rows - 1), arma::span(0)).fill(i+1);
        M(arma::span(i * C_rows, (i+1) * C_rows - 1),
          arma::span(1, C.n_cols)) = C.slice(i);
    }
    
    return M;
}

#include <RcppArmadillo.h>
#include <sitmo.h>  // sitmo::prng_engine
#include <random>  // normal distribution
#include <vector>  // vector class
#ifdef _OPENMP
#include <omp.h>  // OpenMP
#endif

#include "ReplicateTimeseries_types.h"

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
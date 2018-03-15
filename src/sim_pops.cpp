#include <RcppArmadillo.h>
#include <sitmo.h>  // sitmo::prng_engine
#include <random>  // normal distribution
#include <vector>  // vector class
#ifdef _OPENMP
#include <omp.h>  // OpenMP
#endif

#include "repts_types.h"

using namespace Rcpp;



// Draw from normal distribution and fill row
inline void rnd_row(arma::rowvec& rnd_vec,
                    std::normal_distribution<double>& distr,
                    sitmo::prng_engine& eng) {
    for (uint i = 0; i < rnd_vec.n_elem; i++) {
        rnd_vec(i) = distr(eng);
    }
}


//[[Rcpp::export]]
arma::mat sim_pops(const uint& n_gen, const arma::rowvec& N0, const arma::rowvec& r,
                   const arma::mat& alpha, const double& sigma, const uint& seed) {
    uint n_spp = N0.n_elem;
    if (r.n_elem != n_spp || alpha.n_rows != n_spp || alpha.n_cols != n_spp) {
        stop("N0, r, and alpha lengths must be the same.");
    }
    arma::mat N(n_gen + 1, n_spp);
    N.row(0) = N0;
    
    std::normal_distribution<double> rnorm_distr(0.0, sigma);
    sitmo::prng_engine engine(seed);
    arma::rowvec rnd(n_spp, arma::fill::zeros);
    
    for (uint t = 0; t < n_gen; t++) {
        rnd_row(rnd, rnorm_distr, engine);
        N.row(t+1) = N.row(t) % arma::exp(r - N.row(t) * alpha + rnd);
    }
    
    return N;
    
}


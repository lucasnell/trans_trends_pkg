#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <random>  // normal distribution
#include <vector>  // vector class
#ifdef _OPENMP
#include <omp.h>  // OpenMP
#endif

#include "repts_types.h"
#include "pcg.h"  // pcg seeding

using namespace Rcpp;




//' Simulate populations with competition.
//'
//' This does most of the work for the R-exported function below.
//'
//' @param N Output matrix of population abundances.
//' @param distr Normal distribution generator from C++ STL.
//' @param eng pcg32 object that generates random numbers.
//' @inheritParams sim_pops
//'
//'
//' @noRd
//'
void sim_pops_(arma::mat& N,
               const arma::rowvec& r,
               const arma::mat& alpha,
               std::normal_distribution<double>& distr,
               pcg32& eng) {

    uint32 n_spp = N.n_cols;
    uint32 max_t = N.n_rows - 1;

    arma::rowvec rnd(n_spp, arma::fill::zeros);

    for (uint32 t = 0; t < max_t; t++) {

        for (uint32 i = 0; i < rnd.n_elem; i++) rnd(i) = distr(eng);

        N.row(t+1) = N.row(t) % arma::exp(r - N.row(t) * alpha + rnd);
    }

    return;

}



//' Simulate populations with competition.
//'
//' @param max_t Number of time points to simulate.
//' @param N0 Vector of starting population abundances, one for each species.
//' @param r Vector of growth rates, one for each species.
//' @param alpha Matrix of intra- and inter-specific density dependences.
//' @param sigma Standard deviation of process error.
//'
//'
//' @export
//'
//[[Rcpp::export]]
arma::mat sim_pops(const uint32& max_t, const arma::rowvec& N0, const arma::rowvec& r,
                   const arma::mat& alpha, const double& sigma) {

    uint32 n_spp = N0.n_elem;
    if (r.n_elem != n_spp || alpha.n_rows != n_spp || alpha.n_cols != n_spp) {
        stop("N0, r, and alpha lengths must be the same.");
    }
    pcg32 eng = seeded_pcg();

    arma::mat N(max_t + 1, n_spp);
    N.row(0) = N0;

    std::normal_distribution<double> distr(0.0, sigma);

    sim_pops_(N, r, alpha, distr, eng);

    return N;

}


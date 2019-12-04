#include <RcppArmadillo.h>

using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

typedef arma::uword uint_t;


/*

 Compute Highest Posterior Density Intervals for the parameters in an MCMC sample.

 This is a C++ version of \code{\link[coda]{HPDinterval}}.

 @param input Matrix containing the MCMC sample.
 @param prob A number in the interval (0,1) giving the target probability
     content of the intervals.
     The nominal probability content of the intervals is the multiple of
     `1/nrow(input)` nearest to `prob`.

*/
//[[Rcpp::export]]
NumericMatrix hpdi(const arma::mat& input, const double& prob = 0.95) {

    if (input.n_rows <= 1) stop("Number of rows in matrix must be > 1.");

    arma::mat vals = arma::sort(input);

    uint_t n = vals.n_rows;
    uint_t p = vals.n_cols;

    uint_t gap = std::round(n * prob);
    if (gap > (n - 1)) gap = n - 1;
    if (gap < 1) gap = 1;

    arma::uvec init = arma::regspace<arma::uvec>(0, n - gap - 1);

    arma::mat init_mat = vals.rows(init + gap) - vals.rows(init);

    arma::urowvec inds = arma::index_min(init_mat);

    NumericMatrix out(p, 2);
    for (uint_t i = 0; i < p; i++) {
        out(i, 0) = vals(inds(i), i);
        out(i, 1) = vals(inds(i) + gap, i);
    }
    colnames(out) = CharacterVector::create("lower", "upper");

    out.attr("Probability") = static_cast<double>(gap) / static_cast<double>(n);

    return out;

}

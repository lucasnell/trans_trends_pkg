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
NumericVector hpdi(NumericVector input, const double& prob = 0.95) {

    if (input.size() <= 1) stop("Number of items in vector must be > 1.");

    NumericVector vals = input.sort();

    uint_t n = input.size();

    uint_t gap = std::round(static_cast<double>(n) * prob);
    if (gap > (n - 1)) gap = n - 1;
    if (gap < 1) gap = 1;

    arma::uvec init = arma::regspace<arma::uvec>(0, n - gap - 1);

    double min_d = vals[init(0) + gap] - vals[init(0)];
    uint_t ind = 0;

    for (uint_t i = 1; i < init.n_elem; i++) {
        double d = vals[init(i) + gap] - vals[init(i)];
        if (d < min_d) {
            min_d = d;
            ind = i;
        }
    }

    NumericVector out = {vals[ind], vals[ind + gap]};

    out.names() = CharacterVector({"lower", "upper"});

    out.attr("Probability") = static_cast<double>(gap) / static_cast<double>(n);

    return out;

}

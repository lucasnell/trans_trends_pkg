#include <RcppEigen.h>


using namespace Rcpp;



using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;

typedef uint_fast32_t uint32;
typedef uint_fast64_t uint64;



//' Compute Highest Posterior Density Intervals for the parameters in an MCMC sample.
//'
//' This is a C++ version of \code{\link[coda]{HPDinterval}}.
//'
//' @param input Vector containing the MCMC sample.
//' @param prob A number in the interval (0,1) giving the target probability
//'     content of the intervals.
//'     The nominal probability content of the intervals is the multiple of
//'     `1/nrow(input)` nearest to `prob`.
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector hpdi(NumericVector input, const double& prob = 0.95) {

    if (input.size() <= 1) stop("Number of items in vector must be > 1.");

    NumericVector vals = input.sort();

    uint32 n = input.size();

    uint32 gap = std::round(static_cast<double>(n) * prob);
    if (gap > (n - 1)) gap = n - 1;
    if (gap < 1) gap = 1;

    NumericVector init(n - gap);
    for (uint32 i = 0; i < init.size(); i++) init(i) = i;

    double min_d = vals[init(0) + gap] - vals[init(0)];
    uint32 ind = 0;

    for (uint32 i = 1; i < init.size(); i++) {
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

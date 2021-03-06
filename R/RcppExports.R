# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Compute Highest Posterior Density Intervals for the parameters in an MCMC sample.
#'
#' This is a C++ version of \code{\link[coda]{HPDinterval}}.
#'
#' @param input Vector containing the MCMC sample.
#' @param prob A number in the interval (0,1) giving the target probability
#'     content of the intervals.
#'     The nominal probability content of the intervals is the multiple of
#'     `1/nrow(input)` nearest to `prob`.
#'
#' @export
#'
hpdi <- function(input, prob = 0.95) {
    .Call(`_TransTrendsPkg_hpdi`, input, prob)
}


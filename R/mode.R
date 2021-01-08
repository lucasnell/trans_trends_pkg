#' Estimates the marginal parameter mode using kernel density estimation
#'
#' @param x numeric vector of posteriors
#' @param adjust numeric, passed to density to adjust the bandwidth of the
#'     kernal density
#' @param ... other parameters to pass to `density`
#'
#' @export
#'
posterior_mode <- function (x, adjust = 0.1, ...) {
    if (is.numeric(x) == FALSE & !is.null(dim(x))) {
        warning("posterior_mode expects a numeric vector")
    }
    dx <- density(x, adjust = adjust, ...)
    .mode <- dx$x[which.max(dx$y)]
    return(.mode)
}


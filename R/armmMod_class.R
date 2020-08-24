

new_armmMod <- function(.stan, .call, .hmc, .x_means_sds, .y_means_sds, .stan_data) {

    stopifnot(inherits(.call, "call"))
    stopifnot(inherits(.hmc, "logical"))
    stopifnot(is.null(.x_means_sds) || inherits(.x_means_sds, "data.frame"))
    stopifnot(is.null(.y_means_sds) || inherits(.y_means_sds, "data.frame"))

    # So it doesn't show the whole function if using do.call:
    if (.call[1] != as.call(quote(armm()))) {
        .call[1] <- as.call(quote(armm()))
    }

    liz_obj <- structure(list(stan = .stan, call = .call,
                              hmc = .hmc,
                              x_means_sds = .x_means_sds,
                              y_means_sds = .y_means_sds,
                              stan_data = .stan_data),
                         class = "armmMod")

    return(liz_obj)

}


#' Print a `armmMod` object.
#'
#' @export
#'
#' @noRd
#'
print.armmMod <- function(x, digits = max(3, getOption("digits") - 3), ...) {

    cat("\nCall to armm:\n")
    cat(paste(trimws(deparse(x$call)), collapse = " "), "\n\n")
    cat(sprintf("  * Standardized X: %s\n", !is.null(x$x_means_sds)))
    cat(sprintf("  * Standardized Y: %s\n", !is.null(x$y_means_sds)))

}

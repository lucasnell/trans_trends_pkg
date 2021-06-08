

new_armmMod <- function(.stan, .call, .options, .x_means_sds, .y_means_sds,
                        .stan_data, .orig_data) {

    stopifnot(inherits(.call, "call"))
    stopifnot(is.null(.x_means_sds) || inherits(.x_means_sds, "data.frame"))
    stopifnot(is.null(.y_means_sds) || inherits(.y_means_sds, "data.frame"))
    stopifnot(inherits(.orig_data, "environment"))

    # Convert from environment to data.frame for storage:
    .orig_data <- as.data.frame(as.list(.orig_data))

    # So it doesn't show the whole function if using do.call:
    if (.call[1] != as.call(quote(armm()))) {
        .call[1] <- as.call(quote(armm()))
    }

    armm_obj <- structure(list(stan = .stan, call = .call,
                               options = .options,
                               x_means_sds = .x_means_sds,
                               y_means_sds = .y_means_sds,
                               stan_data = .stan_data,
                               orig_data = .orig_data),
                          class = "armmMod")

    return(armm_obj)

}


#' Print an `armmMod` object.
#'
#' @method print armmMod
#'
#' @export
#'
#' @noRd
#'
print.armmMod <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("Autoregressive mixed model\n")
    cat(sprintf("  method:         %s\n",
                ifelse(x$options$hmc, "Hamiltonian Monte Carlo",
                       "Direct optimization")))
    cat(sprintf("  family:         %s\n", family(x)))
    form <- x$options$formula
    cat("  formula:       ", paste(trimws(deparse(form)), collapse = " "), "\n")
    cat("  data:          ", paste(trimws(deparse(x$call$data)),
                                   collapse = " "))
    scale_strs <- c()
    if (!is.null(x$x_means_sds)) scale_strs <- c(scale_strs, "scaled X")
    if (!is.null(x$y_means_sds)) scale_strs <- c(scale_strs, "scaled Y")
    if (length(scale_strs) > 0) {
        cat(sprintf(" << %s >>\n", paste(scale_strs, collapse = " & ")))
    } else cat("\n")
    cat("  observations:  ", nobs(x), "\n")
    if (length(x$options$rstan_control) > 0) {
        cat("  rstan options: ", gsub("list\\(|\\)", "",
                                      deparse(x$options$rstan_control)), "\n")
    }
    cat(sprintf("  obs. error:     %s\n", "sig_obs" %in% names(x$stan)))
    cat("------\n")

    LL <- median(rstan::extract(x$stan, "log_lik_sum")[[1]])
    cat("Median posterior logLik:", LL, "\n")
    loo_obj <- tryCatch(loo(x),
                        warning = function(w) {
                            return(list(loo = suppressWarnings(loo(x)),
                                        warn = w))
                        },
                        error = function(e) {
                            return(e)
                        })
    if (inherits(loo_obj, "list")) {
        loo_est <- loo_obj[["loo"]][["estimates"]]["looic","Estimate"]
        loo_warn <- paste0("** loo warning: ",
                           trimws(loo_obj[["warn"]][["message"]]), "\n")
    } else if (inherits(loo_obj, "error")) {
        loo_est <- NA_real_
        loo_warn <- paste0("*** loo error: ", paste(loo_obj))
    } else {
        loo_est <- loo_obj[["estimates"]]["looic","Estimate"]
        loo_warn <- ""
    }

    waic_obj <- tryCatch(waic(x),
                         warning = function(w) {
                             return(list(waic = suppressWarnings(waic(x)),
                                         warn = w))
                         },
                         error = function(e) {
                             return(e)
                         })
    if (inherits(waic_obj, "list")) {
        waic_est <- waic_obj[["waic"]][["estimates"]]["waic","Estimate"]
        waic_warn <- paste0("** waic warning: ",
                            trimws(waic_obj[["warn"]][["message"]]), "\n")
    } else if (inherits(waic_obj, "error")) {
        waic_est <- NA_real_
        waic_warn <- paste0("*** waic error: ", paste(waic_obj))
    } else {
        waic_est <- waic_obj[["estimates"]]["waic","Estimate"]
        waic_warn <- ""
    }

    print(c(LOO = loo_est, WAIC = waic_est), digits = digits)
    cat(loo_warn, waic_warn, sep = "")

    invisible(NULL)
}

#' Summary for an `armmMod` object.
#'
#' Standard errors are based on either quantiles or HPDIs,
#' and are half the width of the 68% uncertainty interval.
#'
#' @param object an object of class `armmMod`, a result of a call to `armm`.
#' @param se_method A single string, for either using quantiles (`"quantile"`)
#'     or HPDI (`"hpdi"`) to compute the standard errors.
#'     Defaults to `"quantile"`.
#' @param digits the number of significant digits to use when printing.
#' @param ... Not used.
#'
#' @export
#'
#' @method summary armmMod
#'
#'
summary.armmMod <- function(object,
                            se_method = c("quantile", "hpdi"),
                            digits = max(3, getOption("digits") - 3),
                            ...) {


    print(object, digits = digits)

    cat("------\n")

    AR <- autoreg(object)
    if (!is.null(AR)) {
        cat("Autoregressive parameters:\n")
        rownames(AR) <- paste0(" ", rownames(AR))
        print(AR, digits = digits)
        rownames(AR) <- substring(rownames(AR), 2)
    } else cat("No autoregressive parameters\n")

    RE <- NULL
    if (any(grepl("^sig_beta", names(object$stan)))) {
        cat("------\n")
        cat("Random effects (Standard deviation):\n")
        RE <- print_sigma_betas(object, se_method, digits)
    } else cat("------\nNo random effects\n")

    cat("------\n")
    cat("Fixed effects:\n")
    FF <- print_fixef_w_se(object, se_method, digits)

    invisible(list(autoreg = AR, ranef = RE, fixef = FF))

}






#' @name loo
#' @title Extract LOO estimates from an `armmMod` object
#' @aliases loo loo.armmMod
#' @docType methods
#' @param x A fitted model with class `armmMod`.
#' @inheritParams rstan::loo
#' @importFrom loo loo
#' @export loo
#' @method loo armmMod
#' @seealso \code{\link[loo]{loo}} \code{\link[rstan]{loo.stanfit}}
#' @export
loo.armmMod <- function(x,
                        save_psis = FALSE,
                        cores = getOption("mc.cores", 1),
                        moment_match = FALSE,
                        k_threshold = 0.7,
                        ...) {
    if (!is.null(x$options$rstan_control$iter) &&
        x$options$rstan_control$iter < 3) {
        stop("loo cannot be calculated when the number of iterations is < 3")
    }

    .loo <- rstan::loo(x$stan, pars = "log_lik",
                       save_psis = save_psis, cores = cores,
                       moment_match = moment_match,
                       k_threshold = k_threshold, ...)
    return(.loo)
}


#' @name waic
#' @title Extract WAIC from an `armmMod` object
#' @aliases waic waic.armmMod
#' @docType methods
#' @param x A fitted model with class `armmMod`.
#' @param \dots Ignored.
#' @importFrom loo waic
#' @export waic
#' @method waic armmMod
#' @seealso \code{\link[loo]{waic}}
#' @export
waic.armmMod <- function(x, ...) {
    if (!is.null(x$options$rstan_control$iter) &&
        x$options$rstan_control$iter < 3) {
        stop("waic cannot be calculated when the number of iterations is < 3")
    }
    ll_m <- loo::extract_log_lik(x$stan, parameter_name = "log_lik")
    .waic <- loo::waic.matrix(ll_m)
    return(.waic)
}



#
# Print fixed effects coefficients and standard errors.
# This is used in the summary method
#
print_fixef_w_se <- function(object, se_method, digits) {

    A <- rstan::extract(object$stan, "alpha")[[1]]

    SEs <- bayesian_se(A, se_method)

    fixef_df <- data.frame(Median = apply(A, 2, median),
                           `Std.Error` = SEs)

    rownames(fixef_df) <- paste0(" ", colnames(object$stan_data$x))
    print(fixef_df, digits = digits)

    rownames(fixef_df) <- substring(rownames(fixef_df), 2)

    return(fixef_df)
}


#
# This is for use in the summary method.
#
print_sigma_betas <- function(object, se_method, digits) {

    if (!any(grepl("^sig_beta", names(object$stan)))) invisible(NULL)

    B <- rstan::extract(object$stan, "sig_beta")[[1]]

    sigmaB_df <- cbind(object$rnd_names[,c("Groups", "Name")],
                       data.frame(Median = apply(B, 2, median)))
    sigmaB_df$`Std.Error` <- bayesian_se(B, se_method)
    sigmaB_df <- sigmaB_df[order(sigmaB_df$Groups),]
    rownames(sigmaB_df) <- NULL

    sigmaB_df2 <- sigmaB_df
    sigmaB_df2$Groups[sigmaB_df2$Groups ==
                         c("", sigmaB_df2$Groups[-nrow(sigmaB_df2)])] <- ""

    print(sigmaB_df2, row.names = FALSE, right = FALSE, digits = digits)

    return(sigmaB_df)
}



# Estimate Bayesian standard errors by computing half the width of the
# 68% CI.
# This CI can be calculated using quantiles or HPDI
# Used multiple times below.
bayesian_se <- function(x, se_method) {
    se_method <- match.arg(tolower(se_method), c("quantile", "hpdi"))
    if (!isTRUE(length(dim(x)) == 2)) {
        stop("\n`x` arg to `bayesian_se` should be 2-dimensional.")
    }
    if (se_method == "quantile") {
        upper <- unname(sapply(1:ncol(x), function(i) quantile(x[,i], 0.84)))
        lower <- unname(sapply(1:ncol(x), function(i) quantile(x[,i], 0.16)))
    } else if (se_method == "hpdi") {
        ints <- lapply(1:ncol(x), function(i) hpdi(x[,i], 0.68))
        upper <- sapply(ints, function(x) x[["upper"]])
        lower <- sapply(ints, function(x) x[["lower"]])
    } else stop("\n`se_method` arg to `bayesian_se` should be `\"quantile\"` or",
                "`\"hpdi\"`")
    SEs <- 0.5 * (upper - lower)
    return(SEs)
}







#' @name ranef
#' @title Extract random-effects estimates from an `armmMod` object
#' @aliases ranef random.effects ranef.armmMod
#' @docType methods
#' @inheritParams fixef
#' @return A data frame of random-effects estimates and standard errors.
#' @importFrom lme4 ranef
#' @export ranef
#' @method ranef armmMod
#' @export
ranef.armmMod <- function(object,
                          se_method = c("quantile", "hpdi"),
                          ...) {

    if (sum(object$stan_data$g_per_ff) == 0) return(NULL)
    sigma_names <- names(object$stan)[grepl("^sig_beta\\[", names(object$stan))]
    z_names <- names(object$stan)[grepl("^z\\[", names(object$stan))]

    S <- rstan::extract(object$stan, sigma_names)
    S <- lapply(1:length(S),
                function(i) {
                    matrix(as.numeric(S[[i]]), length(S[[i]]),
                           object$stan_data$lev_per_g[i])
                })
    S <- do.call(cbind, S)
    Z <- do.call(cbind, rstan::extract(object$stan, z_names))
    Z <- S * Z

    ests <- unname(apply(Z, 2, median))
    ses <- bayesian_se(Z, se_method)

    ranef_df <- cbind(object$rnd_lvl_names,
                      data.frame(Median = ests,
                                 `Std.Error` = ses),
                      stringsAsFactors = FALSE)

    return(ranef_df)

}




#' Extract autoregressive parameters
#'
#' @param object A fitted model with class `armmMod`
#' @param \dots Additional arguments, ignored for method compatibility.
#'
#' @return AR parameters.
#'
#' @export
autoreg <- function(object, ...) {
    UseMethod("autoreg")
}
#' @export
#' @param se_method A single string, for either using quantiles (`"quantile"`)
#'     or HPDI (`"hpdi"`) to compute the standard errors.
#'     Defaults to `"quantile"`.
#' @method autoreg armmMod
#' @describeIn autoreg Autoregressive parameters for an `armmMod` object
autoreg.armmMod <- function(object,
                            se_method = c("quantile", "hpdi"),
                            ...) {

    if (is.null(eval(object$call$ar_form))) return(NULL)

    phis <- rstan::extract(object$stan, "phi")[[1]]

    SEs <- bayesian_se(phis, se_method)

    autoreg_df <- data.frame(Median = apply(phis, 2, median),
                             `Std.Error` = SEs)
    rownames(autoreg_df) <- object$ar_names
    return(autoreg_df)
}





#' @name fixef
#' @title Extract fixed-effects estimates from an `armmMod` object
#' @aliases fixef fixed.effects fixef.armmMod
#' @docType methods
#' @inheritParams summary.armmMod
#' @return A data frame of fixed-effects estimates and standard errors.
#' @importFrom lme4 fixef
#' @export fixef
#' @method fixef armmMod
#' @export
fixef.armmMod <- function(object,
                          se_method = c("quantile", "hpdi"),
                          ...) {

    fef <- rstan::extract(object$stan, "alpha")[[1]]

    A <- data.frame(Median = apply(fef, 2, median),
                    `Std.Error` = bayesian_se(fef, se_method))

    rownames(A) <- colnames(object$stan_data$x)

    return(A)
}






#' Coefficients from an `armmMod` object
#'
#' @inheritParams summary.armmMod
#' @method coef armmMod
#'
#' @return A data frame of fixed and (if present) random effects estimates
#'     with standard errors.
#'
#' @export
coef.armmMod <- function(object,
                         se_method = c("quantile", "hpdi"),
                         ...) {

    if (sum(object$stan_data$g_per_ff) == 0) return(fixef(object, se_method))

    fef <- rstan::extract(object$stan, "alpha")[[1]]
    colnames(fef) <- colnames(object$stan_data$x)


    sigma_names <- names(object$stan)[grepl("^sig_beta\\[", names(object$stan))]
    z_names <- names(object$stan)[grepl("^z\\[", names(object$stan))]

    S <- rstan::extract(object$stan, sigma_names)
    S <- lapply(1:length(S),
                function(i) {
                    matrix(as.numeric(S[[i]]), length(S[[i]]),
                           object$stan_data$lev_per_g[i])
                })
    S <- do.call(cbind, S)
    Z <- do.call(cbind, rstan::extract(object$stan, z_names))
    R <- S * Z


    # Combining random and fixed effects:

    coefs <- cbind(object$rnd_lvl_names,
                   data.frame(Median = NA_real_,
                              `Std.Error` = NA_real_),
                   stringsAsFactors = FALSE)
    if (!all(coefs$Name %in% colnames(fef))) {
        stop(paste("\nINTERNAL ERROR: Not all coefficient names from",
                   "`object$rnd_lvl_names$Name` are present in",
                   "`colnames(object$stan_data$x)`.",
                   "Make sure you don't edit `stan_data` or `rnd_lvl_names`",
                   "fields manually."))
    }
    if (nrow(coefs) != ncol(R)) {
        stop(paste("\nINTERNAL ERROR: `nrow(object$rnd_lvl_names) !=",
                   "names(object$stan)[grepl(\"^z\\[\", names(object$stan))]`.",
                   "Make sure you don't edit `stan` or `rnd_lvl_names`",
                   "fields manually."))
    }
    E <- R
    for (i in 1:nrow(coefs)) {
        .f <- fef[,coefs$Name[i]]
        E[,i] <- E[,i] + .f
    }

    coefs$Median <- unname(apply(E, 2, median))
    coefs$`Std.Error` <- bayesian_se(E, se_method)


    return(coefs)

}


#' Residuals of `armmMod` objects
#'
#' Getting different types of residuals for `armmMod` objects.
#'
#' @param object A fitted model with class `armmMod`.
#' @param type Type of residuals, currently only `"response"` is programmed.
#' @param \dots Additional arguments, ignored for method compatibility.
#' @return A vector of residuals.
#'
#' @method residuals armmMod
#'
#' @export
residuals.armmMod <- function(object,
                              type = "response",
                              ...) {
    if (family(object) == "normal"){
        y <- object$stan_data$y
        mu <- fitted(object)
        res <- y - mu
        return(res)
    }

    stop("\nresiduals only for normal so far.")
}


#' Fitted values for an `armmMod` object
#'
#' @param object A fitted model with class `armmMod`
#' @param \dots Additional arguments, ignored for method compatibility.
#' @return Fitted values.
#'
#' @method fitted armmMod
#' @export
fitted.armmMod <- function(object, ...) {
    if (!object$options$hmc) {
        stop("\n`fitted.armmMod` method not yet implemented for direct ",
             "optimization", call. = FALSE)
    }
    apply(rstan::extract(object$stan, "y_pred")[[1]], 2, median)
}


#' Extracting the model frame from an `armmMod` object
#'
#' @inheritParams stats::model.frame
#' @method model.frame armmMod
#'
#'
#' @export
model.frame.armmMod <- function(formula, ...) {
    model.frame(formula$formula, formula$data)
}


#' Number of observations in an `armmMod` object
#'
#' @inheritParams stats::nobs
#' @method nobs armmMod
#' @export
nobs.armmMod <- function(object, ...) {
    return(object$stan_data$n_obs)
}


#' Family for an `armmMod` object
#'
#' @inheritParams stats::family
#' @method family armmMod
#'
#' @export
family.armmMod <- function(object, ...) {
    fam <- if (is.null(object$call$family)) {
        formals(armm)[["family"]]
    } else object$call$family
    return(fam)
}




#' @name pp_check
#' @title Posterior predictive checks for an `armmMod` object
#' @aliases pp_check pp_check.armmMod
#' @docType methods
#' @param object A fitted model with class `armmMod`.
#' @param ... Arguments for function \code{\link[bayesplot]{ppc_dens_overlay}}.
#' @return A `ggplot2` object containing the plot.
#' @importFrom bayesplot pp_check
#' @export pp_check
#' @method pp_check armmMod
#' @seealso \code{\link[bayesplot]{ppc_dens_overlay}}
#' @export
pp_check.armmMod <- function(object, ...) {

    yrep <- rstan::extract(object$stan, "y_pred")[[1]]
    y <- object$stan_data$y

    p <- bayesplot::ppc_dens_overlay(y, yrep, ...)

    return(p)

}




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

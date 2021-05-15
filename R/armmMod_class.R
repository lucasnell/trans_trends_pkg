


new_armmMod <- function(.stan, .call, .hmc, .x_means_sds, .y_means_sds,
                        .stan_data, .orig_data) {

    stopifnot(inherits(.call, "call"))
    stopifnot(inherits(.hmc, "logical"))
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
                               hmc = .hmc,
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
                ifelse(x$hmc, "Hamiltonian Monte Carlo",
                       "Direct optimization")))
    cat(sprintf("  family:         %s\n", family(x)))
    form <- if (inherits(x$call$formula, "formula")) {
        x$call$formula
    } else {
        eval(x$call$formula, parent.frame(1L))
    }
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
    if (!is.null(x$call$rstan_control)) {
        cat("  rstan options: ", gsub("list\\(|\\)", "",
                                      deparse(x$call$rstan_control)), "\n")
    }
    cat(sprintf("  obs. error:     %s\n", "sig_obs" %in% names(x$stan)))
    cat("------\n")

    LL <- median(rstan::extract(x$stan, "log_lik_sum")[[1]])
    cat("Median posterior logLik:", LL, "\n")
    loo_obj <- tryCatch(loo(x),
                        warning = function(w) {
                            return(list(loo = suppressWarnings(loo(x)),
                                        warn = w))
                        })
    if (inherits(loo_obj, "list")) {
        loo_est <- loo_obj[["loo"]][["estimates"]]["looic","Estimate"]
        loo_warn <- paste0("** loo warning: ",
                           trimws(loo_obj[["warn"]][["message"]]), "\n")
    } else {
        loo_est <- loo_obj[["estimates"]]["looic","Estimate"]
        loo_warn <- ""
    }

    waic_obj <- tryCatch(waic(x),
                         warning = function(w) {
                             return(list(waic = suppressWarnings(waic(x)),
                                         warn = w))
                         })
    if (inherits(waic_obj, "list")) {
        waic_est <- waic_obj[["waic"]][["estimates"]]["waic","Estimate"]
        waic_warn <- paste0("** waic warning: ",
                            trimws(waic_obj[["warn"]][["message"]]), "\n")

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

    se_method <- match.arg(tolower(se_method), c("quantile", "hpdi"))

    cat("------\n")

    AR <- autoreg(object)
    if (!is.null(AR)) {
        cat("Autoregressive parameters:\n")
        rownames(AR) <- paste0(" ", rownames(AR))
        print(AR, digits = digits)
    } else cat("No autoregressive parameters\n")

    if (any(grepl("^sig_beta", names(object$stan)))) {
        cat("------\n")
        cat("Random effects:\n")
        print_sigma_betas(object, digits = digits)
    } else cat("------\nNo random effects\n")

    cat("------\n")
    cat("Fixed effects:\n")
    print_fixef_w_se(object, se_method, digits)

    invisible(NULL)

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
    .loo <- rstan::loo(x$stan, pars = "log_lik_sum",
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
    ll_m <- loo::extract_log_lik(x$stan, parameter_name = "log_lik_sum")
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

    invisible(NULL)
}


#
# This is for use in the summary method.
#
print_sigma_betas <- function(object, digits) {

    if (!any(grepl("^sig_beta", names(object$stan)))) invisible(NULL)

    B <- rstan::extract(object$stan, "sig_beta")[[1]]

    sigmaB_df <- cbind(object$rnd_names[,c("Groups", "Name")],
                       data.frame(`Std.Dev.` = apply(B, 2, median)))
    sigmaB_df <- sigmaB_df[order(sigmaB_df$Groups),]
    rownames(sigmaB_df) <- NULL

    sigmaB_df$Groups[sigmaB_df$Groups ==
                         c("", sigmaB_df$Groups[-nrow(sigmaB_df)])] <- ""

    print(sigmaB_df, row.names = FALSE, right = FALSE, digits = digits)

    invisible(NULL)
}



# Estimate Bayesian standard errors by computing half the width of the
# 68% CI.
# This CI can be calculated using quantiles or HPDI
# Used multiple times below.
bayesian_se <- function(x, se_method) {
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
#' @return A list of random-effects estimates.
#' @importFrom lme4 ranef
#' @export ranef
#' @method ranef armmMod
#' @export
ranef.armmMod <- function(object, ...) {

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

    ests <- unname(apply(S * Z, 2, median))

    ranef_df <- cbind(object$rnd_lvl_names, data.frame(Estimate = ests),
                      stringsAsFactors = FALSE)

    ranef_list <- lapply(split(ranef_df, ranef_df$Groups),
                         function(x) {
                             x <- x[,c("Name", "Level", "Estimate")]
                             xdf <- lapply(split(x, x$Name), function(y) {
                                 ydf <- data.frame(y$Estimate)
                                 colnames(ydf) <- y$Name[1]
                                 rownames(ydf) <- y$Level
                                 return(ydf)
                             })
                             # Make sure all rownames are in same order
                             if (length(xdf) > 1) {
                                 rnames <- rownames(xdf[[1]])
                                 for (i in 2:length(xdf)) {
                                     if (!all(rownames(xdf[[i]]) == rnames)) {
                                         xdf[[i]] <- xdf[[i]][rnames,,drop=F]
                                     }
                                 }
                             }
                             xdf <- do.call(cbind, xdf)
                             return(xdf)
                         })

    return(ranef_list)

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

    se_method <- match.arg(tolower(se_method), c("quantile", "hpdi"))
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
#' @param object A fitted model with class `armmMod`.
#' @param ... Ignored.
#' @return A vector of fixed-effects estimates.
#' @importFrom lme4 fixef
#' @export fixef
#' @method fixef armmMod
#' @export
fixef.armmMod <- function(object, ...) {

    A <- apply(rstan::extract(object$stan, "alpha")[[1]], 2, median)

    names(A) <- colnames(object$stan_data$x)

    return(A)
}





#' Coefficients from an `armmMod` object
#'
#' @param object An object of class `armmMod`.
#' @param ... Ignored.
#' @method coef armmMod
#'
#' @return If random variables present, a list of data frames (1 for each
#'     random-effect grouping variable).
#'     Each data frame contains coefficients for each level of the variable.
#'     If the model doesn't contain random variables, then it returns a
#'     vector of coefficients.
#'
#' @export
coef.armmMod <- function(object, ...) {

    fixed <- fixef(object)
    random <- ranef(object)
    if (is.null(random)) return(fixed)

    coef_obj <- lapply(random,
                       function(x) {
                           m <- matrix(fixed, nrow(x), length(fixed),
                                       byrow = TRUE)
                           colnames(m) <- names(fixed)
                           rownames(m) <- rownames(x)
                           for (n in colnames(x)) m[,n] <- m[,n] + x[,n]
                           return(m)
                       })

    return(coef_obj)

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
    if (!object$hmc) {
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

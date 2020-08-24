


#' Version of lapply that returns the list flattened to a vector or array.
#'
#' @noRd
#'
f_apply <- function(X, FUN, bind_fxn = base::c, ...) {
    x <- lapply(X, FUN, ...)
    x <- do.call(bind_fxn, x)
    return(x)
}

#' Deparse a formula to a single string.
#'
#' @noRd
#'
f_deparse <- function(form) {
    ss <- deparse(form, 500L)
    if (length(ss) > 1) ss <- paste(ss, collapse = "")
    return(ss)
}


#' Expand interactions for a formula, part of a formula, or a formula string.
#'
#' @param x The input formula, part of a formula, or a formula string.
#'
#' @noRd
#'
expand_inters <- function(x) {
    if (inherits(x, "call") || inherits(x, "name") || inherits(x, "formula")) {
        fs <- f_deparse(x)
    } else if (isTRUE(x == 1)) {
        return("1")
    } else if (!inherits(x, "character")) {
        stop("Must be character, formula, name, call, or 1.")
    }
    if (!grepl("~", fs)) fs <- paste("DEPENDENT_VARIABLE ~", fs)
    f <- as.formula(fs)
    fo <- attr(terms.formula(f), "term.labels")
    # Add back intercept if necessary:
    if ("1" %in% trimws(strsplit(fs, "[\\+|\\*|\\~]")[[1]])) {
        fo <- c("1", fo)
    }
    return(fo)
}


#' Make info for random effects and check inputs.
#'
#'
#' @noRd
#'
make_check_rand <- function(fixed, rand_chunks, data) {

    # Check strings for weird characters:
    # ----------------------*
    if (any(f_apply(rand_chunks, function(x)
        sum(gregexpr("|", f_deparse(x), fixed=TRUE)[[1]] > 0) > 1))) {
        stop("\nWhen using bars to represent random effects, use parentheses ",
             "around each random effect and don't use multiple bars.",
             call. = FALSE)
    }
    # The variables being grouped into random effects:
    rand <- f_apply(rand_chunks, function(x) f_deparse(x[[2]]))
    # The variables grouping the random effects:
    rand_g <- f_apply(rand_chunks, function(x) f_deparse(x[[3]]))
    # Checking for weird chars in LHS (`rand`)
    if (any(grepl("[\\||\\(|\\)|~]", rand))) {
        stop("\nNone of the following characters should be on the left side of ",
             "a random effect-specifying chunk: |, (, ), ~",
             call. = FALSE)
    }
    # Check for weird characters in the grouping sides:
    if (any(grepl("[\\||\\(|\\)|~]", rand_g))) {
        stop("\nNone of the following characters should be on the right side of ",
             "a random effect-specifying chunk: |, (, ), ~",
             call. = FALSE)
    }
    # Check for any interactions in grouping side:
    if (any(grepl("\\*", rand_g))) {
        stop("\nDo not include interactions in the grouping variables for random ",
             "effects. Instead create a new factor to group by.",
             call. = FALSE)
    }

    # Split up, then group random effects:
    # ----------------------*
    rand <- lapply(rand_chunks, function(x) expand_inters(x[[2]]))
    rand_g <- f_apply(1:length(rand_chunks),
                      function(i) replicate(length(rand[[i]]),
                                            all.vars(rand_chunks[[i]][[3]]),
                                            simplify = FALSE))
    rand <- c(rand, recursive = TRUE)
    names(rand_g) <- rand
    if ("0" %in% rand) {
        stop("\n`0` cannot be a covariate grouped for random effects.", call. = FALSE)
    }
    if (sum(duplicated(rand)) > 0) {
        stop("\nYou have duplicates in the terms used on the left-hand side of the ",
             "bar in your random effects.",
             call. = FALSE)
    }
    # Make sure all grouping variables are factors
    g_factors <- sapply(unique(c(rand_g, recursive = TRUE)),
                        function(x) inherits(get(x, envir = data), "factor"))
    if (!all(g_factors)) {
        stop("\nAll mixed-effects grouping variables must be factors. ",
             "Please change the following variable(s) to factor(s): `",
             paste(names(g_factors)[!g_factors], collapse = ", "), "`.",
             call. = FALSE)
    }
    # Take care of weird scenario where a 0 is in `fixed` and a 1 in `rand`
    if ("0" %in% fixed && "1" %in% rand) {
        stop("\nIf including the intercept as a random effect (using `+ (1 | ...)`) ",
             "you cannot also add `+ 0` to the formula.",
             call. = TRUE)
    }
    # Make sure all random-effects variables have an intercept:
    if (!all(rand %in% fixed)) {
        stop("\nAll random effects are required to have an intercept. ",
             "So `y ~ x1 + (x1 | g1)` works, but `y ~ (x1 | g1)` does not.",
             call. = TRUE)
    }

    return(list(rand = rand, rand_g = rand_g))
}



#' Create vector of the number of columns in model matrix per covariate.
#'
#' This will always be 1 for non-factors.
#'
#'
#' @noRd
#'
n_cols_per_cov <- function(fixed, formula, data) {

    ncpc <- f_apply(fixed,
                             function(v) {
                                 if (v == "0") return(0L)
                                 if (v == "1") return(1L)
                                 f <- reformulate(v, f_deparse(formula[[2]]),
                                                  intercept = FALSE)
                                 z <- model.matrix(f, data = data)
                                 return(ncol(z))
                             })
    names(ncpc) <- fixed
    if (sum(is_fctr <- ncpc > 1) > 0) {
        has_inter <- "1" %in% fixed
        first_fctr <- TRUE
        for (i in which(is_fctr)) {
            if (first_fctr) {
                first_fctr <- FALSE
                if (has_inter) ncpc[i] <- ncpc[i] - 1
            } else ncpc[i] <- ncpc[i] - 1
        }
    }

    return(ncpc)
}


#' Get information from some object by time series.
#'
#' This returns an error if the time series returns >1 unique value.
#'
#'
#' @noRd
#'
get_ts_info <- function(x, start_end_mat, err_msg_arg) {
    Z <- NULL
    err_msg <- sprintf("\n%s must not differ within any time series.", err_msg_arg)
    if (!is.null(dim(x))) {
        if (length(dim(x)) != 2) stop("get_ts_info not for non-2D arrays");
        get_info <- function(i) {
            z <- x[start_end_mat[i,1]:start_end_mat[i,2],,drop=FALSE]
            if (nrow(z) > 1 && any(apply(z, 2, sd) != 0)) stop(err_msg, call. = FALSE)
            return(z[1,,drop=FALSE])
        }
        Z <- f_apply(1:nrow(start_end_mat), get_info, base::rbind)
    } else {
        get_info <- function(i) {
            z <- unique(x[start_end_mat[i,1]:start_end_mat[i,2]])
            if (length(z) != 1) stop(err_msg, call. = FALSE)
            return(z)
        }
        Z <- f_apply(1:nrow(start_end_mat), get_info)
    }
    return(Z)
}



#' Make information list when formula has random effects.
#'
#'
#' @noRd
#'
form_info_w_rand <- function(formula, fixed, rand_chunks, data, start_end_mat) {

    rand_info <- make_check_rand(fixed, rand_chunks, data)
    rand <- rand_info$rand
    rand_g <- rand_info$rand_g

    # Number of columns in model matrix for each covariate.
    # This will always be 1 unless the covariate is a factor.
    ncpc <- n_cols_per_cov(fixed, formula, data)

    # number of coefficients (fixed effects + intercepts)
    # number of groups per fixed effect:
    g_per_ff <- f_apply(fixed,
                        function(x) {
                            z <- ifelse(!x %in% rand, 0L, length(rand_g[[x]]))
                            return(rep(z, ncpc[[x]]))
                        })
    # levels per group (repeated by fixed effect):
    lev_per_g <- f_apply(fixed,
                         function(x) {
                             if (!x %in% rand) return(integer(0))
                             z <- f_apply(rand_g[[x]],
                                          function(xx) {
                                              length(levels(get(xx, envir = data)))
                                          })
                             return(rep(z, ncpc[[x]]))
                         })
    # grouping structure for betas:
    # Getting dataset of grouping variables
    g_mat <- f_apply(fixed,
                     function(x) {
                         if (!x %in% rand) return(NULL)
                         z <- lapply(rand_g[[x]],
                                     function(xx) cbind(get(xx, envir = data)))
                         return(do.call(cbind, rep(z, ncpc[[x]])))
                     }, bind_fxn = base::cbind)
    # Because we're mapping to a vector of all groups combined, we need to
    # add the max level from previous columns for columns >1
    if (ncol(g_mat) > 1) {
        for (i in 2:ncol(g_mat)) {
            g_mat[,i] <- g_mat[,i] + max(g_mat[,(i-1)])
        }
    }
    b_groups <- get_ts_info(g_mat, start_end_mat,
                            "Random effects-grouping variables")

    info_list <- list(
        g_per_ff = g_per_ff,
        lev_per_g = lev_per_g,
        b_groups = b_groups
    )

    return(info_list)
}


#' Check that formulas are of proper structure.
#'
#' This is only used to create error messages, so this function invisibly returns `NULL`.
#'
#' @param formula Formula
#' @param arg Which argument in `liz_fit` the formula is for.
#'
#'
#' @noRd
#'
proper_formula <- function(formula, arg) {

    arg <- match.arg(arg, c("formula", "time_form", "ar_form", "y_scale"))

    if (arg == "ar_form" && is.null(formula)) return(NULL)

    err_msg <- NULL

    if (arg == "formula") {
        allow_bars <- TRUE
        one_sided <- FALSE
        allow_inter <- TRUE
    } else if (arg == "time_form") {
        allow_bars <- TRUE
        one_sided <- TRUE
        allow_inter <- FALSE
    } else {  # ar_form and y_scale
        allow_bars <- FALSE
        one_sided <- TRUE
        allow_inter <- FALSE
    }

    if (!inherits(formula, "formula") || !identical(quote(`~`), formula[[1]])) {
        err_msg <- paste("\nArgument", arg, "is not a formula.")
        stop(err_msg, call. = FALSE)
    }
    if (sum(all.names(formula) == "~") > 1) {
        stop("\nYou should never include > 1 tilde (`~`) in any `liz_fit` argument.",
             call. = TRUE)
    }

    if (one_sided && length(formula) != 2) err_msg <- "one-sided"
    if (!one_sided && length(formula) != 3) err_msg <- "two-sided"
    if (!is.null(err_msg)) {
        stop(sprintf("\nThe provided `%s` argument is not %s.", arg, err_msg),
             call. = FALSE)
    }

    # First check for colons to provide more useful error message for this case:
    if (grepl("\\:", f_deparse(formula))) {
        stop("\nIn the `", arg, "` argument, you've included a colon. ",
             "This is not allowed in `liz_fit`, so please just use an asterisk to ",
             "specify interactive effects.", call. = FALSE)
    }
    # Do the same for double bars:
    if (sum(all.names(formula) == "||") > 0) {
        stop("\nDouble bars (`||`) are not allowed in formulas.", call. = FALSE)
    }
    allowed_chars <- c("\\.", "\\_", "\\~", "\\+")
    if (allow_inter) allowed_chars <- c(allowed_chars, "\\*")
    if (allow_bars) allowed_chars <- c(allowed_chars, "\\(", "\\)", "\\|")
    grep_str <- paste0(c(paste0("(?!", allowed_chars, ")"), "[[:punct:]]"), collapse = "")
    if (grepl(grep_str, f_deparse(formula), perl = TRUE)) {
        stop("\nIn the `", arg, "` argument, you're only allowed the following ",
             "characters: ", paste(gsub("\\\\", "", allowed_chars), collapse = ", "), ".",
             call. = FALSE)
    }

    if (arg == "time_form") {
        if (length(formula[[2]]) != 3 || !identical(quote(`|`), formula[[2]][[1]]) ||
            length(formula[[2]][[2]]) != 1 ||
            any(grepl("(?!\\.)(?!\\_)(?!\\+)[[:punct:]]",
                      f_deparse(formula[[2]][[3]]), perl = TRUE))) {
        stop("\nThe `time_form` argument must be a one-sided formula ",
             "with a bar separating the time variable from the variable(s) ",
             "separating time series (e.g., `~ time | species + site`). ",
             "Only one variable can be to the left of the bar, and variables to ",
             "the right of the bar must only be separated by `+`.",
             call. = FALSE)
        }
    }

    if (arg == "formula" && length(formula[[2]]) != 1) {
        stop("\nThe `formula` argument should only have one variable on the ",
             "left-hand side of the formula.",
             call. = FALSE)
    }

    if (arg == "y_scale") {
        if (length(all.vars(formula)) > 1) {
            stop("\nThe `y_scale` argument, if a formula, should only include ",
                 "one variable.", call. = FALSE)
        }
    }

    invisible(NULL)

}



#' Do initial checks for inputs to `liz_fit`.
#'
#' @inheritParams liz_fit
#'
#' @noRd
#'
initial_input_checks <- function(formula,
                                 time_form,
                                 ar_form,
                                 y_scale,
                                 data,
                                 ar_bound,
                                 x_scale) {

    if (!inherits(data, "environment")) {
        stop("\nIn `liz_fit`, the `data` argument must be a list, data frame, or ",
             "environment.", call. = FALSE)
    }

    # Checking structure of the formulas:
    proper_formula(formula, "formula")
    proper_formula(time_form, "time_form")
    proper_formula(ar_form, "ar_form")
    if (!is.null(y_scale)) {
        if (!inherits(y_scale, c("formula", "character")) ||
            (inherits(y_scale, "character") && length(y_scale) != 1)) {
            stop("\nIn `liz_fit`, argument `y_scale` should be NULL, a formula, or a ",
                 "single string.", call. = FALSE)
        }
        if (inherits(y_scale, "formula")) proper_formula(y_scale, "y_scale")
    }
    if (!inherits(ar_bound, "logical") || length(ar_bound) != 1) {
        stop("\nIn `liz_fit`, the `ar_bound` argument must be a logical ",
             "of length 1.", call. = FALSE)
    }
    if (!inherits(x_scale, "logical") || length(x_scale) != 1) {
        stop("\nIn `liz_fit`, the `x_scale` argument must be a logical ",
             "of length 1.", call. = FALSE)
    }

    # Check that all variables specified in formulas are found in `data`:
    if (!all(all.vars(formula) %in% ls(envir = data))) {
        stop("\nThe following variables in `formula` are not present in the `data` ",
             "object:\n  ",
             paste(all.vars(formula)[!all.vars(formula) %in% ls(envir = data)],
                   collapse = ", "),
             call. = FALSE)
    }
    if (!all(all.vars(time_form) %in% ls(envir = data))) {
        stop("\nThe following variables in `time_form` are not present in the `data` ",
             "object:\n  ",
             paste(all.vars(time_form)[!all.vars(time_form) %in% ls(envir = data)],
                   collapse = ", "),
             call. = FALSE)
    }
    if (!is.null(ar_form)) {
        if (!all(all.vars(ar_form) %in% ls(envir = data))) {
            stop("\nThe following variables in `ar_form` are not present in the `data` ",
                 "object:\n  ",
                 paste(all.vars(ar_form)[!all.vars(ar_form) %in% ls(envir = data)],
                       collapse = ", "),
                 call. = FALSE)
        }
        # Check that all items in `ar_form` are in `time_form`
        if (length(all.vars(ar_form)) != 0 &&
            !all(all.vars(ar_form) %in% all.vars(time_form[[2]][[3]]))) {
            stop("\nAll variables in `ar_form` must be present in the part of `time_form` ",
                 "to the right of the bar.", call. = FALSE)
        }
    }

    invisible(NULL)

}



#' Check data lengths, sort data, return vector of observations per time series.
#'
#' This function also checks for duplicate times and for "bad" factors (those with
#' missing levels or just one level).
#'
#' It also checks for NAs.
#'
#' Changes made to `data` are in place.
#'
#' @noRd
#'
check_len_sort_data <- function(formula,
                                time_form,
                                ar_form,
                                data) {

    time_var_names <- all.vars(time_form)
    if (length(time_var_names) > 1) {
        time_var_names <- c(time_var_names[-1], time_var_names[1])
    }
    time_vars <- lapply(time_var_names, get, envir = data)

    len <- sapply(time_vars, length)
    if (length(unique(len)) != 1) {
        stop("\nThe variables in `time_form` aren't the same length.",
             call. = FALSE)
    }
    len <- len[1]

    order_ <- do.call(order, time_vars)
    do_reorder <- !isTRUE(all.equal(order_, 1L:len))

    # Checking for duplicate times and creating `obs_per`:
    tv_df <- do.call(cbind, time_vars)
    tv_df <- as.data.frame(tv_df[order_,])
    # split by grouping variables to make looking for repeats easier:
    tv_list <- split(tv_df, interaction(tv_df[, (ncol(tv_df)-1):1], drop = TRUE))
    obs_per <- integer(length(tv_list))
    for (i in 1:length(tv_list)) {
        tv <- tv_list[[i]]
        obs_per[i] <- nrow(tv)
        if (nrow(tv) > length(unique(tv[,ncol(tv)]))) {
            stop("\nDuplicate times detected in your data. ",
                 "Please verify that you're grouping time series properly.",
                 call. = FALSE)
        }
    }

    # Check for proper lengths and reorder if necessary:
    all_vars <- f_apply(list(formula, time_form, ar_form), all.vars)
    all_vars <- unique(all_vars)
    for (v in all_vars) {
        x <- get(v, envir = data)
        if (length(x) != len) {
            stop("\nThe following variable in one of the input formulas isn't the ",
                 "same length as others: `", v, "`.",
                 call. = FALSE)
        }
        if (any(is.na(x))) {
            stop("\nNAs are not allowed in `liz_fit`. The following variable in one of ",
                 "the input formulas has them: `", v, "`.", call. = FALSE)
        }
        if (inherits(x, "factor")) {
            if (length(levels(x)) < length(unique(x))) {
                stop("\nInput factors should not contain missing levels. ",
                     "Use the `droplevels` function on the following column to ",
                     "fix this problem: `", v, "`.",
                     call. = FALSE)
            }
            if (length(levels(x)) < 2) {
                stop("\nInput factors should contain at least two levels. ",
                     "Remove or fix the following variable: `", v, "`.",
                     call. = FALSE)
            }
        }
        if (do_reorder) {
            assign(v, x[order_], envir = data)
        }
    }

    return(obs_per)
}








#' Make objects related to coefficients.
#'
#' Output list has the following items:
#'
#' - `n_coef`:      number of coefficients (fixed effects + intercepts)
#' - `g_per_ff`:    number of groups per fixed effect
#' - `lev_per_g`:   number of levels per group (repeated by fixed effect)
#' - `b_groups`:    grouping structure for betas
#' - `x`:           predictor variables
#' - `x_means_sds`: (optional) means and SDs of x variables if they were scaled
#'
#' @noRd
#'
make_coef_objects <- function(formula, time_form, ar_form, data, obs_per,
                              ar_bound, x_scale) {

    # Starting and ending positions for each time series
    # (this will be useful for later functions):
    start_end_mat <- cbind(c(1, cumsum(head(obs_per, -1)) + 1),
                           cumsum(obs_per))

    rand_chunks <- findbars(formula)
    fixed <- expand_inters(nobars(formula)[[3]])

    # Make sure intercept is explicitly included in fixed if not already:
    if (!any(c("0","1") %in% fixed)) fixed <- c("1", fixed)
    # If the intercept isn't the first item, switch it in the names so that it is:
    if (sum(fixed == "1") > 0 && (i <- which(fixed == "1")) != 1) {
        fixed <- c("1", fixed[-i])
    }

    # Scale x variables if desired:
    x_means_sds <- NULL
    if (x_scale) {
        x_vars <- unique(fixed[!(grepl(":", fixed) | fixed %in% 0:1)])
        x_means_sds <- data.frame(name = x_vars, stringsAsFactors = FALSE)
        x_means_sds$mean  <- NA_real_
        x_means_sds$sd  <- NA_real_
        for (i in 1:nrow(x_means_sds)) {
            x_name  <- x_means_sds$name[i]
            x <- get(x_name, data)
            # Don't scale factors:
            if (inherits(x, "factor")) next
            .mean <- mean(x, na.rm = TRUE)
            .sd <- sd(x, na.rm = TRUE)
            x <- (x - .mean) / .sd
            assign(x_name, x, envir = data)
            x_means_sds$mean[i]  <- .mean
            x_means_sds$sd[i]  <- .sd
        }

    }


    out <- list()

    if (is.null(rand_chunks)) {

        out$x <- model.matrix(formula, data = data)
        out$n_coef <- ncol(out$x)
        out$g_per_ff <- rep(0, out$n_coef)
        out$lev_per_g <- structure(integer(0), .Dim = 0)
        out$b_groups <- matrix(0L, nrow = nrow(start_end_mat), ncol = 0)

    } else {

        ff_form <- reformulate(fixed, f_deparse(formula[[2]]))

        out <- form_info_w_rand(formula, fixed, rand_chunks, data, start_end_mat)
        out$x <- model.matrix(ff_form, data = data)
        out$n_coef <- ncol(out$x)


    }

    # Set up groups of autoregressive parameters
    if (length(all.vars(ar_form)) == 0) {
        out$p_groups <- rep(1, nrow(start_end_mat))
    } else {
        out$p_groups <- do.call(interaction,
                                      c(lapply(all.vars(ar_form), get, envir = data),
                                        drop = TRUE))
        out$p_groups <- as.integer(out$p_groups)
        out$p_groups <- get_ts_info(out$p_groups, start_end_mat,
                                          "Autoregressive term-grouping variables")
    }

    # Adding other info:
    out$y <- eval(formula[[2]], envir = data)
    out$n_obs <- length(out$y)
    out$n_ts <- length(obs_per)
    out$obs_per <- obs_per
    out$time <- eval(time_form[[2]][[2]], envir = data)
    out$x_means_sds <- x_means_sds

    if (ar_bound) {
        out$p_bound <- 1
    } else out$p_bound <- Inf

    # This fixes the issue of when a parameter is defined in stan as an array,
    # but only a single number is input.
    # Even though this is a vector in R, stan returns an error in this case.
    for (x in c("g_per_ff", "lev_per_g")) {
        if (length(out[[x]]) == 1) out[[x]] <- structure(out[[x]], .Dim = 1)
    }

    return(out)

}










set_priors <- function(stan_data, priors, x_scale, y_scale) {

    # Defaults:
    prior_list <- list(alpha = cbind(0, 1),
                       sig_beta = do.call(rbind, rep(list(cbind(0, 1)),
                                                     sum(stan_data$g_per_ff))),
                       phi = do.call(rbind, rep(list(cbind(0, 0.5)),
                                                max(stan_data$p_groups))),
                       sig_res = cbind(0, 1))

    # Use defaults if no priors provided, but return error if no scaling done:
    if (is.null(priors)) {
        if (!x_scale || is.null(y_scale)) {
            stop("\nYou must specify all priors when you decide not to scale the ",
                 "x or y variables. The `priors` argument must therefore contain ",
                 "all of the following names: ",
                 paste(sprintf("\"%s\"", names(prior_list)), collapse = ", "), ".",
                 call. = FALSE)
        }
        return(prior_list)
    }

    # Basic checks on priors:
    if (!inherits(priors, "list") || is.null(names(priors)) ||
        !all(sapply(priors, inherits, what = "matrix")) ||
        !all(sapply(priors, is.numeric))) {
        stop("\nThe `priors` argument should be a named list of numeric matrices.",
             call. = FALSE)
    }
    if (any(names(priors) == "")) {
        stop("\nThe `priors` argument should have names for all items within.",
             call. = FALSE)
    }
    if (any(duplicated(names(priors)))) {
        stop("\nThe `priors` argument should not have duplicate names.",
             call. = FALSE)
    }
    if (!all(names(priors) %in% names(prior_list))) {
        stop("\nThe `priors` argument should only contain the following names: ",
             paste(sprintf("\"%s\"", names(prior_list)), collapse = ", "), ".",
             call. = FALSE)
    }

    # We are forcing users to provide all priors if they are not scaling
    # x or y variables:
    if (!x_scale || is.null(y_scale)) {
        if (!all(names(prior_list) %in% names(priors))) {
            stop("\nYou must specify all priors when you decide not to scale the ",
                 "x or y variables. The `priors` argument must therefore contain ",
                 "all of the following names: ",
                 paste(sprintf("\"%s\"", names(prior_list)), collapse = ", "), ".",
                 call. = FALSE)
        }
    }


    # Make sure dimensions are correct, then replace defaults
    for (n in names(priors)) {
        if (any(dim(prior_list[[n]]) != dim(priors[[n]]))) {
            stop("\nDimensions for the priors matrix for \"", n, "\" should be ",
                 nrow(prior_list[[n]]), "x", ncol(prior_list[[n]]), ", not ",
                 nrow(priors[[n]]), "x", ncol(priors[[n]]), ".",
                 call. = FALSE)
        }
        prior_list[[n]] <- priors[[n]]
    }

    # So that they don't conflict with current objects in stan file:
    names(prior_list) <- paste0("prior_", names(prior_list))

    return(prior_list)

}




# start doc ------
#' Mixed-effects, autoregressive model.
#'
#'
#' @section Setting priors:
#' Priors are for the hyperparameters for the normal distributions of
#' fixed-effect coefficients and intercepts, random-effect standard deviations,
#' autoregressive parameters, and residual standard deviations.
#' Autoregressive parametes and standard deviation sampling distributions are
#' truncated above zero.
#' By default, most have priors of \eqn{\mu = 0} and \eqn{\sigma = 1}, except for
#' the autoregressive parameters that have priors of \eqn{\mu = 0} and \eqn{\sigma = 0.5}.
#'
#' To pass priors, you must provide a named list.
#' Each item in the list must be a 2-column matrix, with the first column
#' containing priors for \eqn{\mu} and the second containing those for \eqn{\sigma}.
#' The possible names in the list are the following:
#' \describe{
#'     \item{alpha}{ Fixed effects and intercepts }
#'     \item{phi}{ Autoregressive parameters }
#'     \item{sig_beta}{ Random effect group standard deviations }
#'     \item{sig_res}{ Residual standard deviation }
#' }
#'
#'
#' If scaling is not done on either the x or y variables, then the user must pass
#' *all* priors to replace these defaults.
#'
#'
#'
#' @param formula A required, two-sided, linear formula specifying both fixed
#'     and random effects of the model.
#'     The structure is similar to that in the `lme4` package, but more restrictive.
#'     Random effects should be contained inside parentheses, and a bar (`|`)
#'     should separate the name of the fixed effect from the grouping variables.
#'     All grouping variables must be factors without any missing levels.
#'     Example: `y ~ x1 + (x2 | g1 + g2) + (x3 | g1 + g3)`.
#'     _Note:_ If you do not explicitly specify a random effect for the intercept,
#'     `liz_fit` will not include one. This differs from `lmer`.
#'     In the above example, the intercept will not have a random effect.
#' @param time_form A required, one-sided formula specifying the structure of
#'     the time series (e.g., `~ time | species + site + rep`).
#'     On the left side of the bar there should be the object indicating the actual
#'     time point (e.g., day, hour), and on the right side there should be the
#'     variables that, together, separate all time series.
#'     No time series should ever span multiples of these variables.
#' @param ar_form A required, one-sided formula specifying the grouping to use for
#'     the autoregressive parameter(s).
#'     All groups present here should also be present in the grouping part of the
#'     `time_form` argument; an error is thrown otherwise.
#'     If you want to use a grouping variable here that a variable in `time_form`
#'     is nested within (e.g., using genus here and species in `time_form`),
#'     just insert the higher-level variable (genus in the example) in `time_form`,
#'     as it won't effect the results.
#'     Providing `NULL` for this argument causes `liz_fit` to use a model that does
#'     not account for temporal autocorrelation. This can be useful for testing.
#'
#' @param y_scale
#'
#' @param data An optional list, data frame, or environment that contains
#'     the dependent, independent, and grouping variables.
#'     By default, it uses the environment the function was executed in.
#'
#' @param x_scale
#'
#' @param ar_bound An optional logical for whether to bound the autoregressive
#'     parameter(s) <= 1. Defaults to `FALSE`.
#' @param hmc An optional logical for whether to use Hamiltonian Monte Carlo
#'     sampling for the model fit.
#'     This is `stan`'s default, and it gives samples from a posterior distribution
#'     as output.
#'     When `FALSE`, `liz_fit` will obtain point estimates by maximizing the
#'     joint posterior from the model;
#'     it also returns standard errors based on the Hessian.
#'     Defaults to `FALSE`.
#' @param change An optional logical for whether predictors model the change
#'     in the response variable between time points.
#'     The alternative is for predictors to model the mean in the stationary
#'     distribution of the response variable.
#'     Defaults to `TRUE`.
#' @param priors Named list specifying priors.
#'     `NULL` results in the default priors being used.
#'     An error is returned if this argument is not provided when not scaling x or y
#'     variables.
#'     See "Setting Priors" section for more information.
#'     Defaults to `NULL`.
#' @param rstan_control A list of arguments passed to `rstan::sampling`
#'     or `rstan::optimizing`
#'     (e.g., `iter`, `chains`, `cores`, `algorithm`).
#'     See \code{\link[rstan]{sampling}} or \code{\link[rstan]{optimizing}}.
#'
#' @return A list containing, among other things, a `stanfit` object with the model fit.
#'
#' @export
#'
#' @importFrom rstan sampling
#' @importFrom rstan optimizing
#'
#' @examples
#' formula <- y ~ x1 * x2 + x3 + (x1 * x2 | g1 + g2) + (x3 | g1) + (1 | g2)
#' time_form <- ~ t | tg + g1
#' ar_form <- ~ g1
#' y_scale <- ~ g1
#'
#' set.seed(1)
#' x1_coef <- 1.5
#' x2_coefs <- runif(10, 1, 5)
#' data <- data.frame(g1 = factor(rep(1:5, each = 20)),
#'                    g2 = factor(rep(2:1, each = 50)),
#'                    y = rnorm(100),
#'                    x1 = runif(100),
#'                    x2 = rnorm(100),
#'                    x3 = factor(rep(1:4, 25)),
#'                    t = rep(1:10, 10),
#'                    tg = factor(rep(1:10, each = 10)))
#' data$y <- data$y + data$x1 * x1_coef
#' data$y <- data$y + data$x2 *
#'     sapply(as.integer(interaction(data$g1, data$g2)),
#'            function(i) x2_coefs[i])
#'
#' liz <- liz_fit(formula, time_form, ar_form, y_scale, data,
#'               rstan_control = list(chains = 1, iter = 100))
#'
liz_fit <- function(formula,
                    time_form,
                    ar_form,
                    y_scale,
                    data = parent.frame(1L),
                    x_scale = TRUE,
                    ar_bound = FALSE,
                    hmc = FALSE,
                    priors = NULL,
                    change = TRUE,
                    rstan_control = list()) {

    call_ = match.call()

    if (missing(formula)) {
        stop("\nThe `liz_fit` function requires the `formula` argument.",
             call. = FALSE)
    }
    if (missing(time_form)) {
        stop("\nThe `liz_fit` function requires the `time_form` argument.",
             call. = FALSE)
    }
    if (missing(ar_form)) {
        stop("\nThe `liz_fit` function requires the `ar_form` argument.",
             call. = FALSE)
    }
    if (missing(y_scale)) {
        stop("\nThe `liz_fit` function requires the `y_scale` argument.",
             call. = FALSE)
    }
    if (inherits(data, c("data.frame", "list"))) {
        data <- list2env(data)
    }

    if (!inherits(hmc, "logical") || length(hmc) != 1) {
        stop("\nThe `liz_fit` argument `hmc` must be a single logical.", call. = FALSE)
    }
    if (!inherits(change, "logical") || length(change) != 1) {
        stop("\nThe `liz_fit` argument `change` must be a single logical.",
             call. = FALSE)
    }
    if (!inherits(rstan_control, "list")) {
        stop("\n`rstan_control` must be a list.", call. = FALSE)
    }

    # Check for proper inputs:
    initial_input_checks(formula, time_form, ar_form, y_scale, data, ar_bound, x_scale)
    # Checks for variables being same length and reorders by time if necessary
    obs_per <- check_len_sort_data(formula, time_form, ar_form, data)

    # Create the data to input to the stan model:
    stan_data <- make_coef_objects(formula, time_form, ar_form, data, obs_per,
                                   ar_bound, x_scale)
    stan_data$change <- as.integer(change)

    # Deal with potential scaling of x variables that may or may not have been done
    # inside `make_coef_objects`:
    x_means_sds <- stan_data$x_means_sds
    stan_data$x_means_sds <- NULL

    # Now scale y variables if desired:
    y_means_sds <- NULL
    if (!is.null(y_scale)) {
        if (inherits(y_scale, "formula")) y_scale <- all.vars(y_scale)
        y_scale_var_vec <- tryCatch(
            eval(parse(text = y_scale), data),
            error = function(e) stop("\nIn `liz_fit`, the `y_scale` argument ",
                                     "provided is not found in the `data` argument."))
        if (!inherits(y_scale_var_vec, "factor")) {
            stop("\nIn `liz_fit`, the `y_scale` argument is being parsed to ",
                 "an object of class \"", class(y_scale_var_vec), "\", but it needs ",
                 "to be a factor.")
        }
        y_means_sds <- data.frame(level = sort(unique(y_scale_var_vec)))
        y_means_sds$mean  <- NA_real_
        y_means_sds$sd  <- NA_real_
        for (i in 1:nrow(y_means_sds)) {
            l <- y_means_sds$level[i]
            yy <- stan_data$y[y_scale_var_vec == l]
            if (length(yy) == 0) {
                stop("\nIn `liz_fit`, there is at least one level in the factor ",
                     "that's been input to the `y_scale` argument that is missing.")
            }
            .mean <- mean(yy, na.rm = TRUE)
            .sd <- sd(yy, na.rm = TRUE)
            stan_data$y[y_scale_var_vec == l] <- (yy - .mean) / .sd
            y_means_sds$mean[y_means_sds$level == l] <- .mean
            y_means_sds$sd[y_means_sds$level == l] <- .sd
        }

    }


    # UNCOMMENT BELOW ONCE STAN FILES CAN ACCEPT THE PRIORS
    # # Adding priors to stan_data
    # prior_list <- set_priors(priors, stan_data, x_scale, y_scale)
    # for (n in names(prior_list)) stan_data[[n]] <- prior_list[[n]]

    if (!is.null(ar_form)) {
        rstan_control <- c(rstan_control, list(object = stanmodels[["lizard"]],
                                               data = stan_data))
    } else {
        rstan_control <- c(rstan_control, list(object = stanmodels[["snake"]],
                                               data = stan_data))
    }

    if (hmc) {
        stan_fit <- do.call(sampling, rstan_control)
    } else {
        if (!isTRUE(rstan_control$hessian)) rstan_control$hessian <- TRUE
        stan_fit <- do.call(optimizing, rstan_control)
    }

    lizard_obj <- new_lizard(stan_fit, call_, hmc, x_means_sds, y_means_sds, stan_data)

    return(lizard_obj)
}

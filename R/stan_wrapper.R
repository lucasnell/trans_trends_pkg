
#' Version of lapply that returns the list flattened to a vector or array.
#'
#' @noRd
#'
f_apply <- function(X, FUN, bind_fxn = base::c, ...) {
    x <- lapply(X, FUN, ...)
    x <- do.call(bind_fxn, x)
    return(x)
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
            if (any(apply(z, 2, sd) != 0)) stop(err_msg, call. = FALSE)
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


#' Check that formulas are of proper structure.
#'
#' This is only used to create error messages, so this function invisibly returns `NULL`.
#'
#' @param form Formula
#' @param arg Which argument in `lizard` the formula is for.
#'
#'
#' @noRd
#'
proper_formula <- function(form, arg) {

    arg <- match.arg(arg, c("formula", "time_form", "ar_form"))

    err_msg <- NULL

    if (arg == "formula") {
        allow_bars <- TRUE
        one_sided <- FALSE
    } else if (arg == "time_form") {
        allow_bars <- TRUE
        one_sided <- TRUE
    } else {
        allow_bars <- FALSE
        one_sided <- TRUE
    }

    if (!inherits(form, "formula") || !identical(quote(`~`), form[[1]])) {
        err_msg <- paste("\nArgument", arg, "is not a formula.")
        stop(err_msg, call. = FALSE)
    }

    if (one_sided && length(form) != 2) err_msg <- "one-sided"
    if (!one_sided && length(form) != 3) err_msg <- "two-sided"
    if (!is.null(err_msg)) {
        stop(sprintf("\nThe provided `%s` argument is not %s.", arg, err_msg),
             call. = FALSE)
    }

    allowed_chars <- c("\\.", "\\_", "\\~", "\\+", "\\(", "\\)", "\\|")
    if (!allow_bars) allowed_chars <- allowed_chars[1:4]
    grep_str <- paste0(c(paste0("(?!", allowed_chars, ")"), "[[:punct:]]"), collapse = "")
    if (any(grepl(grep_str, deparse(form), perl = TRUE))) {
        stop("\nIn the `", arg, "` argument, you're only allowed the following ",
             "characters: ", paste(gsub("\\\\", "", allowed_chars), collapse = ", "), ".",
             call. = FALSE)
    }

    if (arg == "time_form") {
        if (length(form[[2]]) != 3 || !identical(quote(`|`), form[[2]][[1]]) ||
            length(form[[2]][[2]]) != 1 ||
            any(grepl("(?!\\+)[[:punct:]]", deparse(form[[2]][[3]]), perl = TRUE))) {
        stop("\nThe `time_form` argument must be a one-sided formula ",
             "with a bar separating the time variable from the variable(s) ",
             "separating time series (e.g., `~ time | species + site`). ",
             "Only one variable can be to the left of the bar, and variables to ",
             "the right of the bar must only be separated by `+`.",
             call. = FALSE)
        }
    }

    if (arg == "formula" && length(form[[2]]) != 1) {
        stop("\nThe `formula` argument should only have one variable on the ",
             "left-hand side of the formula.",
             call. = FALSE)
    }

    if (grepl("\\|\\|", deparse(form))) {
        stop("\nDouble bars (`||`) are not allowed in formulas.",
             call. = FALSE)
    }

    invisible(NULL)

}





#' Do initial checks for inputs to `lizard`.
#'
#' @inheritParams lizard
#'
#' @noRd
#'
initial_input_checks <- function(formula,
                         time_form,
                         data,
                         ar_form) {

    if (!inherits(data, "environment")) {
        stop("\nThe `data` argument to the `lizard` function, if provided, must ",
             "be a list, data frame, or environment.",
             call. = FALSE)
    }

    # Checking structure of the formulas:
    proper_formula(formula, "formula")
    proper_formula(time_form, "time_form")
    proper_formula(ar_form, "ar_form")

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

    invisible(NULL)

}



#' Check data lengths, sort data, return vector of observations per time series.
#'
#' This function also checks for duplicate times and for "bad" factors (those with
#' missing levels or just one level).
#' Changes made to `data` are in place.
#'
#' @noRd
#'
check_len_sort_data <- function(formula,
                                time_form,
                                data,
                                ar_form) {

    time_vars <- all.vars(time_form)
    if (length(time_vars) > 1) time_vars <- c(time_vars[-1], time_vars[1])
    time_vars <- lapply(time_vars, get, envir = data)

    len <- sapply(time_vars, length)
    if (length(unique(len)) != 1) {
        stop("\nThe variables in `time_form` aren't the same length.",
             call. = FALSE)
    }
    len <- len[1]

    order_ <- do.call(order, time_vars)
    do_reorder <- !all.equal(order_, 1L:len)

    # Checking for duplicate times and creating `obs_per`:
    tv_df <- do.call(cbind, time_vars)
    tv_df <- as.data.frame(tv_df[order_,])
    # split by grouping variables to make looking for repeats easier:
    tv_list <- split(tv_df, interaction(tv_df[, 1:(ncol(tv_df)-1)], drop = TRUE))
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
#'
#' @noRd
#'
make_coef_objects <- function(formula, data, start_end_mat) {

    rand_chunks <- lme4::findbars(formula)
    fixed <- strsplit(deparse(lme4::nobars(formula)[[3]], 500L), " \\+ ")[[1]]

    if (is.null(rand_chunks)) {
        # predictor variables
        x <- model.matrix(formula, data = data)
        # number of coefficients:
        n_coef <- ncol(x)
        # groups per independent variable
        g_per_ff <- rep(0, n_coef)
        # levels per group
        lev_per_g <- integer(0)
        # grouping structure for betas
        b_groups <- matrix(0L, nrow = nrow(start_end_mat), ncol = 0)
    } else {
        if (any(f_apply(rand_chunks, function(x)
            length(gregexpr("|", deparse(x), fixed=TRUE)[[1]])) > 1)) {
            stop("\nWhen using bars to represent random effects, use parentheses ",
                 "around each random effect and don't use multiple bars.",
                 call. = FALSE)
        }
        # The variables being grouped into random effects:
        rand <- f_apply(rand_chunks, function(x) deparse(x[[2]], 500L))
        if ("0" %in% rand) {
            stop("\n`0` cannot be a covariate grouped for random effects.", call. = FALSE)
        }
        if (any(grepl("\\+", rand))) {
            stop("\nDon't include multiple variables on the left side of a random ",
                 "effect-specifying chunk. ",
                 "Instead, make multiple chunks, one for each variable being grouped.",
                 call. = FALSE)
        }
        if (any(grepl("[\\||\\(|\\)|~]", rand))) {
            stop("\nNone of the following characters should be on the left side of ",
                 "a random effect-specifying chunk: |, (, ), ~",
                 call. = FALSE)
        }
        if (sum(duplicated(rand)) > 0) {
            stop("\nYou have duplicates in the terms used on the left-hand side of the ",
                 "bar in your random effects.",
                 call. = FALSE)
        }
        # The variables grouping random-effects:
        rand_groups <- lapply(rand_chunks,
                              function(x) strsplit(deparse(x[[3]], 500L)," \\+ ")[[1]])
        names(rand_groups) <- rand
        # Check for weird characters in the grouping sides of the random-effect chunks
        if (any(grepl("[\\||\\(|\\)|~]", c(rand_groups, recursive = TRUE)))) {
            stop("\nNone of the following characters should be on the right side of ",
                 "a random effect-specifying chunk: |, (, ), ~",
                 call. = FALSE)
        }
        # Make sure all grouping variables are factors
        g_factors <- vapply(c(rand_groups, recursive = TRUE),
                            function(x) inherits(get(x, envir = data), "factor"),
                            TRUE)
        if (!all(g_factors)) {
            stop("\nAll mixed-effects grouping variables must be factors. ",
                 "Please change the following variable(s) to factor(s): `",
                 paste(names(g_factors)[!g_factors], collapse = ", "), "`.",
                 call. = FALSE)
        }
        # Take care of weirdness with lme4::nobars changing fixed to y~1 when there
        # are no fixed effects:
        if (all(fixed == "1") && "1" %in% rand) fixed <- character(0)
        # Take care of weird scenario where a 0 is in `fixed` and a 1 in `rand`
        if ("0" %in% fixed && "1" %in% rand) {
            stop("\nDo not include intercept as both random (using `+ (1 | ...)`) ",
                 "and fixed (using `+ 0`) effects.",
                 call. = TRUE)
        }

        # Necessary objects for creating output:
        if (!any(c("0","1") %in% fixed) && !"1" %in% rand) fixed <- c("1", fixed)
        ff_names <- c(fixed, rand)
        if (sum(duplicated(ff_names)) > 0) {
            stop("\nThe following covariate(s) is/are duplicated: ",
                 paste(ff_names[duplicated(ff_names)], collapse = ", "), call. = FALSE)
        }
        # If the intercept isn't the first item, switch it in the names so that it is:
        if (sum(ff_names == "1") > 0 && (i <- which(ff_names == "1")) != 1) {
            ff_names <- c("1", ff_names[-i])
        }
        ff_form <- reformulate(ff_names, deparse(formula[[2]]))
        # Number of columns in model matrix for each covariate.
        # This will always be 1 unless the covariate is a factor.
        n_cols_per_ff <- f_apply(ff_names,
                                 function(v) {
                                     if (v == "0") return(0L)
                                     if (v == "1") return(1L)
                                     z <- get(v, envir = data)
                                     if (inherits(z, "factor")) return(length(levels(z)))
                                     return(1L)
                                 })
        names(n_cols_per_ff) <- ff_names
        if (sum(is_fctr <- n_cols_per_ff > 1) > 0) {
            has_inter <- "1" %in% ff_names
            first_fctr <- TRUE
            for (i in which(is_fctr)) {
                if (first_fctr) {
                    first_fctr <- FALSE
                    if (has_inter) n_cols_per_ff[i] <- n_cols_per_ff[i] - 1
                } else n_cols_per_ff[i] <- n_cols_per_ff[i] - 1
            }
        }
        x <- model.matrix(ff_form, data = data)

        # number of coefficients (fixed effects + intercepts)
        n_coef <- ncol(x)
        # number of groups per fixed effect:
        g_per_ff <- f_apply(ff_names,
                            function(x) {
                                if (x %in% fixed) z <- 0L
                                else z <- length(rand_groups[[x]])
                                return(rep(z, n_cols_per_ff[[x]]))
                            })
        # levels per group (repeated by fixed effect):
        lev_per_g <- f_apply(ff_names,
                             function(x) {
                                 if (x %in% fixed) return(integer(0))
                                 z <- f_apply(rand_groups[[x]],
                                              function(xx) {
                                                  length(levels(get(xx, envir = data)))
                                                  })
                                 return(rep(z, n_cols_per_ff[[x]]))
                             })
        # grouping structure for betas:
        # Getting dataset of grouping variables
        g_mat <- f_apply(ff_names,
                         function(x) {
                             if (x %in% fixed || x == "0") return(NULL)
                             z <- lapply(rand_groups[[x]],
                                         function(xx) {
                                             cbind(get(xx, envir = data))
                                         })
                             return(do.call(cbind, rep(z, n_cols_per_ff[[x]])))
                         }, bind_fxn = base::cbind)
        b_groups <- get_ts_info(g_mat, start_end_mat,
                                "Random effects-grouping variables")

    }

    out <- list(
        n_coef = n_coef,        # number of coefficients (fixed effects + intercepts)
        g_per_ff = g_per_ff,    # number of groups per fixed effect
        lev_per_g = lev_per_g,  # number of levels per group (repeated by fixed effect)
        b_groups = b_groups,    # grouping structure for betas
        x = x                   # predictor variables
    )

    return(out)

}





#' Mixed-effects, autoregressive model.
#'
#' @param formula A required, two-sided, linear formula specifying both fixed
#'     and random effects of the model.
#'     The structure is similar to that in the `lme4` package, but more restrictive.
#'     Random effects should be contained inside parentheses, and a bar (`|`)
#'     should separate the name of the fixed effect from the grouping variables.
#'     All grouping variables must be factors without any missing levels.
#'     Example: `y ~ x1 + (x2 | g1 + g2) + (x3 | g1 + g3)`.
#' @param time_form A required, one-sided formula specifying the structure of
#'     the time series (e.g., `~ time | species + site + rep`).
#'     On the left side of the bar there should be the object indicating the actual
#'     time point (e.g., day, hour), and on the right side there should be the
#'     variables that, together, separate all time series.
#'     No time series should ever span multiples of these variables.
#' @param data An optional list, data frame, or environment that contains
#'     the dependent, independent, and grouping variables.
#'     By default, it uses the environment the function was execute in.
#' @param ar_form An optional formula specifying the grouping to use for
#'     the autoregressive parameter(s).
#'     All groups present here should also be present in the grouping part of the
#'     `time_form` argument; an error is thrown otherwise.
#'     If you want to use a grouping variable here that a variable in `time_form`
#'     is nested within (e.g., using genus here and species in `time_form`),
#'     just insert the higher-level variable (genus in the example) in `time_form`,
#'     as it won't effect the results.
#'     The left-hand side of this formula is ignored.
#'     Defaults to no grouping, which results in a single parameter estimate.
#' @param rstan_control A list of arguments passed to `rstan::sampling`
#'     (e.g., iter, chains). See \code{\link[rstan]{sampling}}.
#'
#' @return A `stanfit` object containing the model fit.
#' @export
#'
#'
#'
lizard <- function(formula,
                   time_form,
                   data,
                   ar_form,
                   rstan_control = list()) {

    # formula <- y ~ x1 + (x2 | g1 + g2)
    # time_form <- ~ t | tg + g1
    #
    # set.seed(1)
    # x1_coef <- 1.5
    # x2_coefs <- runif(10, 1, 5)
    # data <- data.frame(g1 = factor(rep(1:5, each = 20)),
    #                 g2 = factor(rep(2:1, each = 50)),
    #                 y = rnorm(100),
    #                 x1 = runif(100),
    #                 x2 = rnorm(100),
    #                 x3 = factor(rep(1:4, 25)),
    #                 t = rep(1:10, 10),
    #                 tg = factor(rep(1:10, each = 10)))
    # data$y <- data$y + data$x1 * x1_coef
    # data$y <- data$y + data$x2 * sapply(as.integer(interaction(data$g1, data$g2)),
    #                                     function(i) x2_coefs[i])
    #
    # ar_form <- ~ g1


    if (missing(formula)) {
        stop("\nThe `lizard` function requires the `formula` argument.",
             call. = FALSE)
    }
    if (missing(time_form)) {
        stop("\nThe `lizard` function requires the `time_form` argument.",
             call. = FALSE)
    }
    if (missing(data)) {
        data <- parent.frame(1L)
    } else if (inherits(data, "list") || inherits(data, "data.frame")) {
        data <- list2env(data)
    }
    if (missing(ar_form)) {
        ar_form <- ~ 1
    }
    if (!inherits(rstan_control, "list")) {
        stop("\n`rstan_control` must be a list.", call. = FALSE)
    }

    # Check for proper inputs:
    initial_input_checks(formula, time_form, data, ar_form)
    # Checks for variables being same length and reorders by time if necessary
    obs_per <- check_len_sort_data(formula, time_form, data, ar_form)
    # Starting and ending positions for each time series
    # (this will be useful for later functions):
    start_end_mat <- cbind(c(1, cumsum(head(obs_per, -1)) + 1),
                           c(cumsum(obs_per[-1]), sum(obs_per)))

    # Most of the info for the stan model:
    stan_data <- make_coef_objects(formula, data, start_end_mat)
    # Adding other info:
    stan_data$n_obs <- length(stan_data$y)
    stan_data$n_ts <- length(obs_per)
    stan_data$obs_per <- obs_per
    stan_data$y <- eval(formula[[2]], envir = data)
    stan_data$time <- eval(time_form[[2]][[2]], envir = data)

    # Set up groups of autoregressive parameters
    if (length(all.vars(ar_form)) == 0) {
        stan_data$p_groups <- rep(1, nrow(start_end_mat))
    } else {
        stan_data$p_groups <- do.call(interaction,
                                      c(lapply(all.vars(ar_form), get, envir = data),
                                        drop = TRUE))
        stan_data$p_groups <- as.integer(stan_data$p_groups)
        stan_data$p_groups <- get_ts_info(stan_data$p_groups, start_end_mat,
                                          "Autoregressive term-grouping variables")
    }

    rstan_control <- c(rstan_control, list(object = stanmodels[["lizard"]],
                                           data = stan_data))

    stan_fit <- do.call(rstan::sampling, rstan_control)

    return(stan_fit)
}


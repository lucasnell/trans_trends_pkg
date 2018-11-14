

#' Mixed-effects, autoregressive model.
#'
#' @param formula A required, two-sided, linear formula specifying both fixed
#'     and random effects of the model.
#'     The structure is similar to \code{\link[lme4]{lmer}}.
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
#' @param ar_groups An optional formula specifying the grouping to use for
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
                   ar_groups,
                   rstan_control = list()) {

    formula <- y ~ x1 + (x2 | g1 + g2)
    time_form <- ~ t | tg

    x1_coef <- 1.5
    x2_coefs <- runif(10, 1, 5)
    set.seed(1)
    data <- data.frame(g1 = factor(rep(1:5, each = 20)),
                    g2 = factor(rep(2:1, 50)),
                    y = rnorm(100),
                    x1 = runif(100),
                    x2 = rnorm(100),
                    t = rep(1:20, 5),
                    tg = rep(1:5, each = 20))
    data$y <- data$y + data$x1 * x1_coef
    data$y <- data$y + data$x2 * sapply(as.integer(interaction(data$g1, data$g2)),
                                        function(i) x2_coefs[i])

    ar_groups <- ~ 1
    model.matrix(ar_groups, data = data)



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
    } else {
        if (!any(sapply(c("list", "data.frame", "environment"), inherits, x = data))) {
            stop("\nThe `data` argument to the `lizard` function, if provided, must ",
                 "be a list, data frame, or environment.",
                 call. = FALSE)
        }
        data <- list2env(data)
    }


    if (!all(all.vars(formula) %in% ls(envir = data))) {
        stop("\nThe following variables in `formula` are not present in the `data` ",
             "object:\n  ",
             paste(all.vars(formula)[!all.vars(formula) %in% ls(envir = data)],
                   collapse = ", "),
             call. = FALSE)
    }



    if (missing(ar_groups)) {
        ar_groups <- ~ 1
    }


    stopifnot(inherits(data_df, "data.frame"))
    if (missing(line)) line <- quote(line)
    if (missing(rep)) rep <- quote(rep)
    if (missing(date)) date <- quote(date)
    if (missing(X)) X <- quote(X)
    # theta is already defined:
    if (missing(theta_)) theta_ <- theta

    line <- substitute(line)
    rep <- substitute(rep)
    date <- substitute(date)
    X <- substitute(X)

    data_df <- data_df %>%
        dplyr::select(!!line, !!rep, !!date, !!X) %>%
        dplyr::arrange(!!line, !!rep, !!date) %>%
        mutate_if(is.factor, as.integer) %>%
        identity()

    # number of observations for each time series:
    n_ts_ <- data_df %>%
        dplyr::group_by(!!line, !!rep) %>%
        summarize() %>%
        nrow()
    n_obs_ <- nrow(data_df)
    n_per_ <- data_df %>%
        dplyr::group_by(!!line, !!rep) %>%
        summarize(n_ = n()) %>%
        .[["n_"]] %>%
        set_names(NULL)
    stopifnot(sum(n_per_) == n_obs_)

    X_ <- data_df[[X]]

    n_lines_ <- length(unique(data_df[[line]]))
    # Line number for each time series:
    L_ <- data_df %>%
        dplyr::group_by(!!line, !!rep) %>%
        summarize() %>%
        .[["line"]] %>%
        set_names(NULL)

    model_data_ <- list(
        n_ts = n_ts_,
        n_obs = n_obs_,
        n_per = n_per_,
        X = X_,
        n_lines = n_lines_,
        L = L_,
        theta = theta_
    )

    # // indices
    # int n_obs;                          // # observations
    # int n_ts;                           // # time series
    # int obs_per[n_ts];                  // # observations per time series
    # int n_coef;                         // # coefficients (fixed effects + intercepts)
    # int g_per_ff[n_coef];               // # groups per fixed effect
    # int lev_per_g[sum(g_per_ff)];       // # levels per group (repeated by fixed effect)
    # int b_groups[n_ts, sum(g_per_ff)];  // grouping structure for betas
    # int p_groups[n_ts];                 // grouping structure for phis
    # // data
    # real y[n_obs];                      // response variables
    # real x[n_obs, n_coef];              // predictor variables

    rstan_control <- c(rstan_control, list(object = stanmodels[["lizard"]],
                                           data = model_data_))

    stan_fit <- do.call(rstan::sampling, rstan_control)

    return(stan_fit)
}


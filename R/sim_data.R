
#' Simulate data for a lizard call.
#'
#'
#' @examples
#' library(dplyr)
#' library(tidyr)
#' set.seed(1279015933)
#' data <- tibble(g1 = factor(rep(1:5, each = 8), levels = 1:5),
#'                g2 = factor(rep(1:8, 5), levels = 1:8)) %>%
#'     mutate(t = rep(list(1:20), n()),
#'            x1 = rep(list(runif(20)), n()),
#'            x2 = rep(list(runif(20)), n())) %>%
#'     unnest()
#'
#' formula <- y ~ x1 + x2 + (x1 | g1 + g2) + (1 | g2)
#' time_form <- ~ t | g1+g2
#' ar_form <- ~ g1
#' ar_bound <- FALSE
#'
#' coef = list(x1 = list(g1 = list(0, 0.5), g2 = runif(8)),
#'             x2 = 2,
#'             `1` = list(g2 = list(3, 1.2)))
#' auto_regr = list(0.1, 0.2)
#' resid_sd = 0.5
#'
#' lizard:::sim_data(formula, time_form, ar_form, coef, auto_regr, resid_sd,
#'                   data, ar_bound)
#'
#'
#'
#' @noRd
#'
sim_data <- function(formula,
                     time_form,
                     ar_form,
                     coef,
                     auto_regr,
                     resid_sd,
                     data,
                     ar_bound = FALSE) {

    call_ = match.call()

    stopifnot(inherits(coef, "list"))
    stopifnot((inherits(auto_regr, "list") && length(auto_regr) == 2) ||
                  inherits(auto_regr, "numeric"))
    stopifnot(is.numeric(resid_sd) && resid_sd >= 0 && length(resid_sd) == 1)

    if (!inherits(data, "data.frame")) stop("data must be a data frame.")
    data[[paste(formula[[2]])]] <- 0
    if (!"1" %in% names(coef)) coef[["1"]] <- 0
    data_ <- list2env(data)

    # Check for proper inputs:
    initial_input_checks(formula, time_form, ar_form, data_, ar_bound)
    # Checks for variables being same length and reorders by time if necessary
    obs_per <- check_len_sort_data(formula, time_form, ar_form, data_)

    stan_data <- make_coef_objects(formula, time_form, ar_form, data_, obs_per, ar_bound)

    # Autoregressive parameters for each time series
    if (!is.null(ar_form)) {
        n <- length(unique(stan_data$p_groups))
        if (inherits(auto_regr, "list")) {
            if (ar_bound && auto_regr[[2]] > 1) auto_regr[[2]] <- 1
            stopifnot(auto_regr[[1]] >= 0)
            ar_phis <- runif(n, auto_regr[[1]], auto_regr[[2]])
        } else {
            if (length(auto_regr) != n) {
                stop("The length of auto_regr, if a numeric vector, must be ", n)
            }
            ar_phis <- auto_regr
            if (ar_bound) ar_phis[ar_phis > 1] <- 1
        }
        phis <- sapply(stan_data$p_groups, function(i) ar_phis[i])
    } else phis <- rep(0, length(stan_data$p_groups))



    # Starting and ending positions for each time series
    # (this will be useful for later functions):
    start_end_mat <- cbind(c(1, cumsum(head(obs_per, -1)) + 1),
                           cumsum(obs_per))

    rand_chunks <- lme4::findbars(formula)
    fixed <- expand_inters(lme4::nobars(formula)[[3]])
    # Make sure intercept is explicitly included in fixed if not already:
    if (!any(c("0","1") %in% fixed)) fixed <- c("1", fixed)
    # If the intercept isn't the first item, switch it in the names so that it is:
    if (sum(fixed == "1") > 0 && (i <- which(fixed == "1")) != 1) {
        fixed <- c("1", fixed[-i])
    }

    rand <- lapply(rand_chunks, function(x) expand_inters(x[[2]]))
    rand_g <- f_apply(1:length(rand_chunks),
                      function(i) replicate(length(rand[[i]]),
                                            all.vars(rand_chunks[[i]][[3]]),
                                            simplify = FALSE))
    rand <- c(rand, recursive = TRUE)
    names(rand_g) <- rand
    if (sum(rand == "1") > 0 && (i <- which(rand == "1")) != 1) {
        rand <- c(rand[i], rand[-i])
        rand_g <- rand_g[rand]
    }


    if (!all(c(fixed, rand) %in% names(coef))) {
        stop("\nNot all fixed and random effects are found in `coef` argument.")
    }


    # Extract fixed-effect betas:
    fixed_betas <- numeric(length(fixed))
    names(fixed_betas) <- fixed
    for (f in fixed) {
        v <- coef[[f]]
        if (is.numeric(v)) {
            fixed_betas[f] <- v
        } else {
            if (is.null(names(v))) {
                stop(paste("\nFor", f, "in the `coef` object, no names are provided."))
            }
            if (!all(rand_g[[f]] %in% names(v))) {
                stop(paste("\nFor", f, "in the `coef` object, it doesn't specify an",
                           "object for each grouping variable."))
            }
            if (!all(sapply(v, inherits, what = c("numeric", "list")))) {
                stop(paste("\nFor", f, "in the `coef` object, not all inner objects",
                           "are of class `numeric` or `list`."))
            }
            for (i in 1:length(v)) {
                if (inherits(coef[[f]][[i]], "list")) {
                    fixed_betas[f] <- fixed_betas[f] + coef[[f]][[i]][[1]]
                    coef[[f]][[i]][[1]] <- 0
                } else if (inherits(coef[[f]][[i]], "numeric")) {
                    fixed_betas[f] <- fixed_betas[f] + mean(coef[[f]][[i]])
                    coef[[f]][[i]] <- coef[[f]][[i]] - mean(coef[[f]][[i]])
                } else {
                    stop("All terminal items inside coef must be numerics or lists")
                }
            }
        }
    }


    # Generate random effects:
    random_betas <- rep(list(NA), length(rand))
    names(random_betas) <- rand
    k <- 1
    for (i in 1:length(rand)) {
        x <- rand[[i]]
        random_betas[[i]] <- rep(list(NA), length(rand_g[[x]]))
        names(random_betas[[i]]) <- rand_g[[x]]
        for (j in 1:length(random_betas[[i]])) {
            g <- rand_g[[x]][j]
            n <- diff(range(stan_data$b_groups[,k])) + 1
            if (inherits(coef[[x]][[g]], "list")) {
                random_betas[[i]][[g]] <- rnorm(n,
                                                coef[[x]][[g]][[1]],
                                                coef[[x]][[g]][[2]])
            } else if (inherits(coef[[x]][[g]], "numeric")) {
                if (length(coef[[x]][[g]]) != n) {
                    stop("\nIn `coef` for coefficient `", x, "`, group `", g, "`, ",
                         "its length must be ", n, ", not ", length(coef[[x]][[g]]), ".")
                }
                random_betas[[i]][[g]] <- coef[[x]][[g]]
            } else {
                stop("\nIn `coef` for coefficient ", x, ", group ", g, ", the object ",
                     "must be of class list or numeric.")
            }
            k <- k + 1
        }
    }

    # Make matrix of final coefficients
    betas <- matrix(0, stan_data$n_obs, stan_data$n_coef)
    colnames(betas) <- fixed

    bg_ind <- 1
    for (i in 1:length(fixed)) {
        p <- fixed[i]
        if (p %in% rand) {
            n_bg <- length(random_betas[[p]])
            re_groups <- stan_data$b_groups[,bg_ind:(bg_ind+n_bg-1),drop=FALSE]
            if ((bg_ind+n_bg-1) > 1) {
                inds <- (bg_ind+n_bg-1):bg_ind
                for (ii in inds[inds > 1]) {
                    rei <- ii - bg_ind + 1
                    re_groups[,rei] <- re_groups[,rei] - max(stan_data$b_groups[,(ii-1)])
                }
            }
            random_effects <- lapply(
                1:n_bg,
                function(k) {
                    sapply(re_groups[,k], function(l) random_betas[[p]][[k]][[l]])
                })
            random_effects <- do.call(cbind, random_effects)
            random_effects <- rowSums(random_effects)
            overall <- fixed_betas[[p]] + random_effects
            betas[,p] <- f_apply(1:length(overall), function(k) {
                rep(overall[k], diff(range(start_end_mat[k,]))+1)
            })
            bg_ind <- bg_ind + n_bg
        } else {
            betas[,p] <- fixed_betas[[p]];
        }
    }


    # Deterministic portion of model:
    determ <- stan_data$x * betas
    determ <- rowSums(determ)
    names(determ) <- NULL

    Y <- f_apply(1:nrow(start_end_mat),
                 function(k) {
                     start <- start_end_mat[k,1]
                     end <- start_end_mat[k,2]
                     N_ <- end - start + 1
                     phi <- phis[k]
                     time <- stan_data$time[start:end]
                     determ_ <- determ[start:end]
                     resids <- rnorm(N_, sd = resid_sd)
                     y <- numeric(N_)
                     y[1] <- determ_[1] + resids[1]
                     if (N_ == 1) return(y)
                     for (t in 2:N_) {
                         y[t] <- determ_[t] +
                             phi^(time[t] - time[t-1]) * (y[t-1] - determ_[t-1]) +
                             resids[t]
                     }
                     return(y)
                 })

    new_data <- as.data.frame(as.list(data_))[,colnames(data)]
    new_data[,paste(formula[[2]])] <- Y

    out_obj <- list(data = new_data, fixed = fixed_betas, random = random_betas,
                    ar_terms = phis, resid_sd = resid_sd,
                    call = call_)

    class(out_obj) <- "lizard_data"

    return(out_obj)
}

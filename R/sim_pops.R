
#' Check that input parameters have the correct dimensionality and types.
#'
#' For more info, see description of `sim_pops_ar` below.
#'
#'
#' @noRd
#'
check_inputs <- function(X, N0_mat, b0_mat, b1_mat, rho_mat, vcv_cube) {

    stopifnot(inherits(X, "matrix"))
    stopifnot(inherits(N0_mat, "matrix"))
    stopifnot(inherits(b0_mat, "matrix"))
    stopifnot(inherits(b1_mat, "matrix"))
    stopifnot(inherits(rho_mat, "matrix"))
    stopifnot(inherits(vcv_cube, "array"))

    n_locs = ncol(X)

    if (ncol(N0_mat) != n_locs) {
        stop("N0_mat should have the same number of columns as does X");
    }
    if (ncol(b0_mat) != n_locs) {
        stop("b0_mat should have the same number of columns as does X");
    }
    if (ncol(b1_mat) != n_locs) {
        stop("b1_mat should have the same number of columns as does X");
    }
    if (ncol(rho_mat) != n_locs) {
        stop("rho_mat should have the same number of columns as does X");
    }
    if (dim(vcv_cube)[3] != n_locs) {
        stop("vcv_cube should have the same number of slices as X has columns");
    }

    n_spp = nrow(N0_mat)
    if (nrow(N0_mat) != n_spp) {
        stop("N0_mat should have the same number of rows as does N0");
    }
    if (nrow(b0_mat) != n_spp) {
        stop("b0_mat should have the same number of rows as does N0");
    }
    if (nrow(b1_mat) != n_spp) {
        stop("b1_mat should have the same number of rows as does N0");
    }
    if (nrow(rho_mat) != n_spp) {
        stop("rho_mat should have the same number of rows as does N0");
    }
    if (nrow(vcv_cube) != n_spp || ncol(vcv_cube) != n_spp) {
        stop("vcv_cube should have the same number of rows and columns as N0 has rows");
    }

    invisible(NULL);
}



#' Convert cube of variance-covariance matrices into matrices for generating deviates.
#'
#' Output is array of matrices that can be used to generate random normal deviates with
#' same covariance structure as input matrices.
#'
#' "a vector of independent normal random variables,
#' when multiplied by the transpose of the Cholesky deposition of [vcv] will
#' have covariance matrix equal to [vcv]."
#'
#'
#' @noRd
#'
make_chol_decomp <- function(vcv_cube) {
    chol_decomp <- vcv_cube
    for (i in 1:(dim(vcv_cube)[3])) {
        chol_decomp[,,i] <- t(chol(vcv_cube[,,i,drop=FALSE]))
    }
    return(chol_decomp)
}





sim_pops_ar <- function(X, N0_mat, b0_mat, b1_mat, rho_mat, vcv_cube, obs_sigma,
                        n_cores = 1) {

    if (n_cores > 1 & .Platform$OS.type == "windows") {
        warning("Not programmed to use multiple cores on Windows. Changing to just one.")
        n_cores <- 1
    }

    check_inputs(X, N0_mat, b0_mat, b1_mat, rho_mat, vcv_cube)

    # For turning ~N(0,1) to multivariate with given covariance matrix
    chol_decomp <- make_chol_decomp(vcv_cube)

    n_time <- nrow(X)
    n_locs <- ncol(X)
    n_spp <- nrow(N0_mat)

    one_location <- function(loc) {
        cd_loc = chol_decomp[,,loc,drop=FALSE]
        b0s = b0_mat[,loc,drop=FALSE]
        b1s = b1_mat[,loc,drop=FALSE]
        rhos = rho_mat[,loc,drop=FALSE]
        Xs = X[,loc,drop=FALSE]
        Ns = N[,,loc,drop=FALSE]
        for (t in 1:n_time) {
            rnd = rbind(rnorm(n_spp))
            obs_rnd = rbind(rnorm(n_spp) * obs_sigma)
            for (sp in 1:n_spp) {
                Ns[t+1,sp] = b0s[sp] + b1s[sp] * Xs[t+1] + rhos[sp] *
                    (Ns[t,sp] - b0s[sp] - b1s[sp] * Xs[t] );
            }
            # Generate variance and covariance:
            rnd = cd_loc * rnd
            Ns[t+1,] = Ns[t+1,] + rnd
            Ns[t+1,] = Ns[t+1,] + obs_rnd
        }
        return(Ns)
    }

    N <- array(NA_real_, dim = c(n_time, n_spp, n_locs))
    N[1,,] <- N0_mat

    if (n_cores <= 1) {
        for (loc in 1:n_locs) N[,,loc] = one_location(loc)
    } else {
        if (!requireNamespace("parallel", quietly = TRUE)) {
            stop("\nPackage \"parallel\" needed for simulating to work with multiple ",
                 "cores. Please install it.",
                 call. = FALSE)
        }
        if (n_cores > parallel::detectCores()) stop("\ntoo many cores requested")
        N_sims <- parallel::mclapply(1:n_locs, one_location, mc.cores = n_cores)
        for (loc in 1:n_locs) N[,,loc] = N_sims[[loc]]
    }

    return(N)
}





# library(testthat)
# library(armmr)


test_that("hpdi produces accurate output",
{
    skip_if_not_installed("coda")
    library(coda)

    r <- cbind(lower = numeric(100), upper = numeric(100),
                prob = numeric(100))
    cpp <- cbind(lower = numeric(100), upper = numeric(100),
                  prob = numeric(100))
    for (i in 1:100) {
        n <- as.integer(runif(1, 100, 500))
        z <- runif(n)
        x <- matrix(z, n, 1)
        class(x) <- "mcmc"
        pr <- runif(1)
        r_i <- HPDinterval(x, prob = pr)
        cpp_i <- armmr:::hpdi(z, pr)
        for (n in c("lower", "upper")) {
            r[i, n] <- r_i[, n]
            cpp[i, n] <- cpp_i[[n]]
        }
        r[i,"prob"] <- attr(r_i, "Probability")
        cpp[i,"prob"] <- attr(cpp_i, "Probability")
    }
    expect_equal(r, cpp)

})





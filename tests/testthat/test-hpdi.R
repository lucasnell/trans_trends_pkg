

# library(testthat)
# library(lizard)


test_that("hpdi produces accurate output",
{
    skip_if_not_installed("coda")
    library(coda)

    for (i in 1:100) {
        n <- as.integer(runif(1, 100, 500))
        p <- as.integer(runif(1, 2, 11))
        x <- matrix(runif(n * p), n, p)
        class(x) <- "mcmc"
        pr <- runif(1)
        r <- HPDinterval(x, prob = pr)
        cpp <- lizard:::hpdi(x, pr)
        expect_equivalent(r, cpp)
    }

})





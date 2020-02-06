

# library(testthat)
# library(lizard)


test_that("hpdi produces accurate output",
{
    skip_if_not_installed("coda")
    library(coda)

    for (i in 1:100) {
        n <- as.integer(runif(1, 100, 500))
        z <- runif(n)
        x <- matrix(z, n, 1)
        class(x) <- "mcmc"
        pr <- runif(1)
        r <- HPDinterval(x, prob = pr)
        cpp <- lizard:::hpdi(z, pr)
        expect_equivalent(r, cpp)
    }

})





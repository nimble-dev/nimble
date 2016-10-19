context("Testing of miscellaneous functionality.")

oldmsg <- geterrmessage()

code <- nimbleCode({
    y ~ dnorm(mu, 1)
    mu ~ dnorm(0, 1)
})

try(m <- nimbleModel(code, data = list(y = 0), check = TRUE))

errmsg <- geterrmessage()    

test_that("Test of full model check", expect_equal(oldmsg, errmsg, info = "found error in running full model check"))

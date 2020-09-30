source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context("Testing of expandNodeNames")

test_that("expandNodeNames works for various cases, including going beyond extent of variable", {
    
   code <- nimbleCode({
    for(i in 1:4)
        mu[i] ~ dnorm(0,1)
    for(i in 1:3)
        for(j in 1:3)
            theta[i,j] ~ dnorm(0,1)
    p[1:4] ~ ddirch(alpha[1:4])
   })
   m <- nimbleModel(code, inits = list(alpha = rep(1, 4)))

   ## vector variable
   expect_equal(m$expandNodeNames("mu"), c("mu[1]","mu[2]","mu[3]","mu[4]"))
   expect_equal(m$expandNodeNames("mu[3:5]"), c("mu[3]","mu[4]"))
   expect_equal(m$expandNodeNames("mu[5:7]"), character(0))

   ## matrix variable
   expect_equal(m$expandNodeNames("theta"), c("theta[1, 1]","theta[2, 1]","theta[3, 1]","theta[1, 2]","theta[2, 2]","theta[3, 2]","theta[1, 3]","theta[2, 3]","theta[3, 3]"))
   expect_equal(m$expandNodeNames("theta[3:5,1:2]"), c("theta[3, 1]","theta[3, 2]"))
   expect_equal(m$expandNodeNames("theta[1:2,3:5]"), c("theta[1, 3]","theta[2, 3]"))
   expect_equal(m$expandNodeNames("theta[4:6,5]"), character(0))

   ## multiple inputs, mixed
   expect_equal(m$expandNodeNames(c("theta[1, 7]", "mu")), c("mu[1]","mu[2]","mu[3]","mu[4]"))
   expect_equal(m$expandNodeNames(c("theta[1, 3:5]", "mu[3:5]")), c("theta[1, 3]", "mu[3]", "mu[4]"))
   expect_equal(m$expandNodeNames(c("theta[1, 7]", "mu[5]")), character(0))
   expect_equal(m$expandNodeNames(c("mu[1, 1]", "theta[1, 1]")), "theta[1, 1]")
                
   ## multiple inputs, mixed, not unique
   expect_equal(m$expandNodeNames(c("mu[3:5]", "mu[3:9]"), unique = FALSE),
                c("mu[3]","mu[4]","mu[3]","mu[4]"))
   expect_equal(m$expandNodeNames(c("theta[1:3, 3:5]", "theta[3:5, 1:3]"), unique = FALSE),
                c("theta[1, 3]","theta[2, 3]","theta[3, 3]","theta[3, 1]","theta[3, 2]","theta[3, 3]"))
                                  
   ## indexing with a variable
   expect_error(m$expandNodeNames("theta[[a]]"), "variable was found in the indexing")

   ## multivariate node
   expect_equal(m$expandNodeNames("p[1:4]"), "p[1:4]")
   expect_equal(m$expandNodeNames("p[1:2]"), "p[1:4]")
   expect_equal(m$expandNodeNames("p[1:5]"), "p[1:4]")
   expect_equal(m$expandNodeNames("p[5:7]"), character(0))
   expect_equal(m$expandNodeNames(c("p[1:5]", "mu[3:5]")), c("p[1:4]", "mu[3]", "mu[4]"))
   expect_equal(m$expandNodeNames(c("p[1:5]", "p"), unique = FALSE), c("p[1:4]", "p[1:4]"))              
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)

   


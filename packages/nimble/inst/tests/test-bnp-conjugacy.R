

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of bnp conjugacy")

#RwarnLevel <- options('warn')$warn
#options(warn = -1)

#source(system.file(file.path('tests', 'dynamicIndexingTestLists.R'), package = 'nimble'))

#nimbleOptions(allowDynamicIndexing = TRUE)

## check variations on use of dynamic indexing in BUGS code including building and compiling,
## dependencies, and valid and invalid dynamic index values

#ans1 <- sapply(testsDynIndex, test_dynamic_indexing_model)
#ans2 <- sapply(testsInvalidDynIndex, test_dynamic_indexing_model)
#ans3 <- sapply(testsInvalidDynIndexValue, test_dynamic_indexing_model)

## check conjugacy detection

test_that("Testing conjugacy detection with bnp models", { 
  
  code = nimbleCode({
    for(i in 1:4) 
      y[i] ~ dnorm(thetatilde[xi[i]], sd = 1)
    xi[1:4] ~ dCRP(1)
    for(j in 1:5) 
      thetatilde[j] ~ dnorm(0, 1)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), thetatilde=rnorm(5)))
  conf <- configureMCMC(m)
  expect_match(conf$getSamplers()[[1]]$name, "dCRP_conjugate_dnorm_dnorm",
               info = "failed to detect normal-normal conjugacy")
  
  
})

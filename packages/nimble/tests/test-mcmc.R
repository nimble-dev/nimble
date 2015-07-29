library(testthat)
library(nimble)

## do all tests in inst/tests
##test_package('nimble')

##test_package('nimble', 'copy')
##test_package('nimble', 'math')
test_package('nimble', 'mcmc')
##test_package('nimble', 'meta')
##test_package('nimble', 'models')
##test_package('nimble', 'trunc')
##test_package('nimble', 'user')

warnings()

library(testthat)
library(nimble)

# this will do all the tests in inst/tests
test_package("nimble")

# to only test some of the files in insts/tests, comment out the line above and uncomment as follows:
# test_package("nimble", "mcmc") # e.g., for test-mcmc.R


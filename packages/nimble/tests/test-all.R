library(testthat)
library(nimble)

# this will do all the tests in inst/tests (if any; as of version 0.1 there were none)
test_package("nimble")

# to only test some of the files in insts/tests, comment out the line above and uncomment as follows:
# test_file(system.file("testFile.R", package = "nimble"))


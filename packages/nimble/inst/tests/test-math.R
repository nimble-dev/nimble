### INSTRUCTIONS:
## got to mathTestLists.R:
## enter each test as a list, with an informative name, NIMBLE expression to evaluate, vector of input dimensions, value of output dimension, and (if NIMBLE expression cannot be directly evaluated in R) the equivalent pure R expression whose result should match the NIMBLE result

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of math functions in NIMBLE code")

source(system.file(file.path('tests', 'mathTestLists.R'), package = 'nimble'))

set.seed(0)
ans1 <- sapply(testsVaried, test_math)    ## 12
ans2 <- sapply(testsBasicMath, test_math) ## 70
if(.Platform$OS.type == 'windows') {
    message("Since you are running on Windows, tests stopped prior to reaching max DLL limit.  Please use test-math2 to continue")
    stop()
}
ans3 <- sapply(testsMoreMath, test_math)  ## 41
ans4 <- sapply(testsReduction, test_math) ## 13
ans5 <- sapply(testsComparison, test_math)## 6
ans6 <- sapply(testsMatrix, test_math)    ## 19




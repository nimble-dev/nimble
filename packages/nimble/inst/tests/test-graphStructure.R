## These tests address construction and querying of the graph structure
## test-getDependencies has some more basic tests.
source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))


correctOutputFilename <- 'graphStructureTestResults_Correct.Rout'
## Use this to regenerate results in the test directory:
## library(testthat)
## library(nimble)
## run code below, then:
## writeOutput(cases, correctOutputFilename)
testOutputFilename <- 'graphStructureTestResults.Rout'

clearOldOutput <- function(filename) {
    if(file.exists(filename)) file.remove(filename)
}

appendOutput <- function(filename, case, caseName, casePrefix = "") {
    outputConnection <- file(filename, open = 'at')
    writeLines(caseName, con = outputConnection)
    outputAns <- lapply(case, function(x) writeLines(paste0(casePrefix, paste(x, collapse = " ")), con = outputConnection))
    close(outputConnection)
}

writeOutput <- function(cases, filename) {
    clearOldOutput(filename)
    for(i in seq_along(cases)) appendOutput(filename, cases[[i]], names(cases)[i], casePrefix = paste0(i,": "))
}



cases <- list()
caseName <- 'graph structure tests case 1'
m <- nimbleModel(
    code = nimbleCode({
        a ~ dnorm(0,1)
        for(i in 1:10) x[i] ~ dnorm(a, 1)
        for(j in 1:10) y[j] ~ dnorm(x[j], 1)
    })
)
cases[[caseName]] <- list(
    m$getDependencies('a'),
    m$getDependencies('a', omit = 'a'),
    m$getDependencies('a', self = FALSE),
    m$getDependencies('a', omit = 'a', self = FALSE),
    m$getDependencies('x[2]'),
    m$getDependencies('x[2]', self = FALSE),
    m$getDependencies('x[2]', omit = 'x[2]'),
    m$getDependencies('x[2]', omit = 'y[2:3]'),
    m$getDependencies('x[2:4]'),
    m$getDependencies('x[2:4]', self = FALSE),
    m$getDependencies('x[2:4]', omit = 'x[4:6]'),
    m$getDependencies('a', downstream = TRUE),
    m$getDependencies('a', downstream = TRUE, self = FALSE),
    m$getDependencies('a', downstream = TRUE, returnScalarComponents = TRUE),
    m$getDependencies('a', downstream = TRUE, self = FALSE, returnScalarComponents = TRUE)
    )

caseName <- 'graph structure tests case 2 (dmnorm fully split)'
m <- nimbleModel(
    code = nimbleCode({
        x[1:5] ~ dmnorm(mu[1:5], cov[1:5, 1:5])                        
        for(i in 1:5) {
            ypred[i] <- x[i] + 1
            y[i] ~ dnorm(ypred[i], 1)
        }
    }),
    constants = list(mu = rep(0,5), cov = diag(5))
)

cases[[caseName]] <- list(
    m$getDependencies('x[1:5]'),
    m$getDependencies('x[1:5]', returnScalarComponents = TRUE),
    m$getDependencies('x[1:5]', self = FALSE),
    m$getDependencies('x[1:5]', self = FALSE, returnScalarComponents = TRUE),
    m$getDependencies('x[1:5]',  omit = c('x[2:3]','y[c(1, 5)]')),
    m$getDependencies('x[2]'),
    m$getDependencies('x[2]', returnScalarComponents = TRUE),
    m$getDependencies('x[2]', self = FALSE),
    m$getDependencies('x[2:4]'),
    m$getDependencies('x[2:4]', self = FALSE),
    m$getDependencies('mu[2:4]'),
    m$getDependencies('mu[2:4]', includeRHSonly = TRUE),
    m$getDependencies('mu[2:4]', includeRHSonly = TRUE, returnScalarComponents = TRUE),
    m$getDependencies('mu', includeRHSonly = TRUE)
)

writeOutput(cases, testOutputFilename)
trialResults <- readLines(testOutputFilename)
##correctResults <- trialResults
correctResults <- readLines(system.file(file.path('tests', correctOutputFilename), package = 'nimble'))

test_that('same number of output lines',
          expect_equal(length(trialResults), length(correctResults)))

linesToTest <- min(length(trialResults), length(correctResults))
mapply(function(lineno, trialLine, correctLine) {
    test_that(paste0("output line #", lineno),
              expect_identical(trialLine, correctLine))
}, 1:linesToTest, trialResults, correctResults)

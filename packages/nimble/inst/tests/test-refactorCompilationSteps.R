require(testthat)

context("Testing of old vs. new generated C++ during refactoring steps")

oldWarnLevel <- options('warn')
options(warn = -1)

compareOldAndNewCompilationRC <- function(input) {
    run <- input$run
    name <- input$name
    require(testthat)
    foo <- nimbleFunction(run = run, name = 'foo')
    nimbleOptions(useRefactoredSizeProcessing = FALSE)
    testProject <- nimble:::nimbleProjectClass(name = 'test')
    ## management of compilation through the control list is crude and leads to error
    ## but the project does contain the result, so we can run this, catch the error
    ## and extract the result from the project
    old <- try(compileNimble(foo, project = testProject, control = list(writeFiles = TRUE, compileCpp = FALSE, loadSO = FALSE)))
    filename <- testProject$RCfunInfos[['foo']][['cppClass']]$filename
    ## we could regenerate it, but might as well read it from the file
    oldFooCppOutput <- readLines(file.path(tempdir(), 'nimble_generatedCode', paste0(filename, '.cpp')))
    oldFooHOutput <- readLines(file.path(tempdir(), 'nimble_generatedCode', paste0(filename, '.h')))
    
    nimbleOptions(useRefactoredSizeProcessing = TRUE)
    testProject <- nimble:::nimbleProjectClass(name = 'test')
    new <- try(compileNimble(foo, project = testProject, control = list(writeFiles = TRUE, compileCpp = FALSE, loadSO = FALSE)))
    filename <- testProject$RCfunInfos[['foo']][['cppClass']]$filename
    ## we could regenerate it, but might as well read it from the file
    newFooCppOutput <- readLines(file.path(tempdir(), 'nimble_generatedCode', paste0(filename, '.cpp')))
    newFooHOutput <- readLines(file.path(tempdir(), 'nimble_generatedCode', paste0(filename, '.h')))
    test_that(paste0(name,': .cpp matches'), expect_identical(oldFooCppOutput, newFooCppOutput))
    test_that(paste0(name,': .h matches'), expect_identical(oldFooHOutput, newFooHOutput))
}

testCases <- list(
    list(name = 'Y <- x',
         run = function(x = double(1)) {
             Y <- x
         }))

ans <- lapply(testCases, compareOldAndNewCompilationRC)
options(warn = as.numeric(oldWarnLevel))

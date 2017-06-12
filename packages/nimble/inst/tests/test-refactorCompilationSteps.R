source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of old vs. new generated C++ during refactoring steps")

oldWarnLevel <- options('warn')
options(warn = -1)

## Known concern: ordering of asRow()/asCol() and intermediates

compareOldAndNewCompilationRC <- function(input) {
    run <- input$run
    name <- input$name
    require(testthat)

    foo <- nimbleFunction(run = run, name = 'foo')
    nimbleOptions(useRefactoredSizeProcessing = FALSE)
    nimble:::resetLabelFunctionCreators() ## sets any generated IDs back to 1
    testProject <- nimble:::nimbleProjectClass(name = 'for_comparison')
    ## management of compilation through the control list is crude and leads to error we don't care about
    ## but the project does contain the result, so we can run this, catch the error to keep running
    ## and extract the result from the project
    old <- try(compileNimble(foo, project = testProject, control = list(writeFiles = TRUE, compileCpp = FALSE, loadSO = FALSE)))
    filename <- testProject$RCfunInfos[['foo']][['cppClass']]$filename
    newfilename <- paste0(filename,'_original')
    pathedfilename <- file.path(tempdir(), 'nimble_generatedCode', filename)
    original_pathedfilename <- file.path(tempdir(), 'nimble_generatedCode', newfilename)
    for(ext in c('.h', '.cpp')) file.copy(paste0(pathedfilename, ext), paste0(original_pathedfilename, ext), overwrite = TRUE)
    
    nimbleOptions(useRefactoredSizeProcessing = TRUE)
    nimble:::resetLabelFunctionCreators() ## sets any generated IDs back to 1
    testProject <- nimble:::nimbleProjectClass(name = 'for_comparison')
    new <- try(compileNimble(foo, project = testProject, control = list(writeFiles = TRUE, compileCpp = FALSE, loadSO = FALSE)))
    filename <- testProject$RCfunInfos[['foo']][['cppClass']]$filename
    ## we could regenerate it, but might as well read it from the file
    refactored_pathedfilename <- file.path(tempdir(), 'nimble_generatedCode', filename)

    for(ext in c('.h','.cpp'))
        compareFilesUsingDiff(paste0(refactored_pathedfilename, ext), paste0(original_pathedfilename, ext),
                              main = paste0(ext, ' files do not match for: ', name))    
}

testCases <- list(
    list(name = 'Y <- x',
         run = function(x = double(1)) {
             Y <- x
         }))

ans <- lapply(testCases, compareOldAndNewCompilationRC)

compareOldAndNewMathTest <- function(input) {
    runFun <- gen_runFun(input)
    input$run <- runFun
    compareOldAndNewCompilationRC(input)
}

source(system.file(file.path('tests', 'mathTestLists.R'), package = 'nimble'))
ans2 <- lapply(testsBasicMath, compareOldAndNewMathTest)

options(warn = as.numeric(oldWarnLevel))

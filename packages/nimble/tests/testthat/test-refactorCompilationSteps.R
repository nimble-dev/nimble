source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

context("Testing of old vs. new generated C++ during refactoring steps")

RwarnLevel <- options('warn')$warn
## Many warnings of "creating a .Call() expression with no DLL information",
## so suppressing them.
options(warn = -1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

## Known concern: ordering of asRow()/asCol() and intermediates.

compareOldAndNewCompilationRC <- function(input) {
    name <- paste0('math: ', input$name, ': compiles')
    test_that(name, {
        wrap_if_matches(input$knownFailure, name, expect_error, {
            run <- input$run
            name <- input$name

            foo <- nimbleFunction(run = run, name = 'foo')
            nimbleOptions(useRefactoredSizeProcessing = FALSE)
            nimble:::resetLabelFunctionCreators() ## This sets any generated IDs back to 1.
            testProject <- nimble:::nimbleProjectClass(name = 'for_comparison')
            ## Management of compilation through the control list is crude and leads to error we don't care about
            ## but the project does contain the result, so we can run this, catch the error to keep running
            ## and extract the result from the project.
            compileNimble(foo, project = testProject, control = list(writeFiles = TRUE, compileCpp = FALSE, loadSO = FALSE))
            filename <- testProject$RCfunInfos[['foo']][['cppClass']]$filename
            newfilename <- paste0(filename,'_original')
            pathedfilename <- file.path(tempdir(), 'nimble_generatedCode', filename)
            original_pathedfilename <- file.path(tempdir(), 'nimble_generatedCode', newfilename)
            for(ext in c('.h', '.cpp')) {
                file.copy(paste0(pathedfilename, ext), paste0(original_pathedfilename, ext), overwrite = TRUE)
            }
            
            nimbleOptions(useRefactoredSizeProcessing = TRUE)
            nimble:::resetLabelFunctionCreators() ## This sets any generated IDs back to 1.
            testProject <- nimble:::nimbleProjectClass(name = 'for_comparison')
            compileNimble(foo, project = testProject, control = list(writeFiles = TRUE, compileCpp = FALSE, loadSO = FALSE))
            filename <- testProject$RCfunInfos[['foo']][['cppClass']]$filename
            ## We could regenerate it [edit: what is it?], but might as well read it from the file.
            refactored_pathedfilename <- file.path(tempdir(), 'nimble_generatedCode', filename)

            for(ext in c('.h','.cpp')) {
                compareFilesUsingDiff(paste0(refactored_pathedfilename, ext), paste0(original_pathedfilename, ext),
                                      main = paste0(ext, ' files do not match for: ', name))    
            }
        })
    })
}

testCases <- list(
    list(name = 'Y <- x',
         run = function(x = double(1)) {
             Y <- x
         }))

ans <- lapply(testCases, compareOldAndNewCompilationRC)

compareOldAndNewMathTest <- function(input) {
    runFun <- gen_runFun(input, logicalArgs = input$logicalArgs,
                         returnType = ifelse(is.null(input$returnType), "double", input$returnType))
    input$run <- runFun
    if('knownFailureReport' %in% names(input) && input$knownFailureReport)
        cat("\nBegin expected error message:\n")
    compareOldAndNewCompilationRC(input)
    if('knownFailureReport' %in% names(input) && input$knownFailureReport)
        cat("End expected error message.\n")
}

source(system.file(file.path('tests', 'testthat', 'mathTestLists.R'), package = 'nimble'))
## compilation of two modulo tests only fails if actual compilation is done,
## not just C++ generation, so this clunkily unsets knownFailure
testsBasicMathModified <- lapply(testsBasicMath, function(x) {
    if(x$name %in% c('modulo of vectors', 'modulo of vector and scalar'))
        x$knownFailure <- NULL
    x
})

ans2 <- lapply(testsBasicMathModified, compareOldAndNewMathTest)

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)

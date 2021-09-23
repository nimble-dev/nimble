source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)


context("Testing of testing tools")

## new way to force anything compiled to be clear-compiled
## (note, no check for windows is done [yet?])
test_that('testing withTempProject',
          withTempProject( {
              ## using two nimbleFunctions to make sure the second one is found during compilation
              foo <- nimbleFunction(
                  run = function(x = double(1)) {return(x + 1); returnType(double(1))}
              )
              temporarilyAssignInGlobalEnv(foo) ## unfortunately this is necessary for bar to find foo
              bar <- nimbleFunction(
                  run = function(x = double(1)) {return(foo(x) + 1); returnType(double(1))}
              )
              cbar <- compileNimble(bar) ## will automatically clear and dyn.unload()
              cbar(1:2)
          } )
          )

test_that('works', {
    A <- 4
    expect_failure( expect_equal(A, 5, info = 'A test') )
    })

## For pull request: demonstrate that withTempProject will register a failed test and continue when run through test_package
test_that('testing withTempProject: expected failure',
          withTempProject( {
              ## using two nimbleFunctions to make sure the second one is found during compilation
              foo <- nimbleFunction(
                  run = function(x = double(1)) {return(x + 1); returnType(double(1))}
              )
              A <- 4
              expect_failure( expect_equal(A, 5, info = 'A test') )
              ## uncomment to see an actual failure
              ##expect_equal(A, 5, info = 'A test') 
              temporarilyAssignInGlobalEnv(foo) ## unfortunately this is necessary for bar to find foo
              bar <- nimbleFunction(
                  run = function(x = double(1)) {return(foo(x) + 1); returnType(double(1))}
              )
              cbar <- compileNimble(bar) ## will automatically clear and dyn.unload()
              cbar(1:2)
          } )
          )

## As of 2020-05-05, ability to stop compilation before linking and to use -O1 is disabled,
## as we need to switch to system2 and handling these options would have been a pain.
## new way to check only if compilation up to (not including) linking works
## this also by default replaces any -O# flags to -O1
test_that('testing expect_compiles', {    
    foo <- nimbleFunction(
        run = function(x = double(1)) {return(x + 1); returnType(double(1))}
    )
    temporarilyAssignInGlobalEnv(foo) ## unfortunately this is necessary for bar to find foo
    bar <- nimbleFunction(
        run = function(x = double(1)) {return(foo(x) + 1); returnType(double(1))}
    )
    ## temporarilyAssignInGlobalEnv(foo) ## why is this being done again?
    ## cat("\nBegin expected error message:\n")
    expect_compiles(bar, dirName = '.', link = TRUE, forceO1 = FALSE) ## arguments as for compileNimble
    ## cat("End expected error message.\n")
    ## will bail out with "Error : safely stopping before linking"
    ## alternative format in test_util would use
    ##    expect_compiles(compileNimble(bar))
})

test_that('testing expect_compiles expected failure', {    
    foo <- nimbleFunction(
        run = function(x = double(1)) {return(x + 1); returnType(double(1))}
    )
    temporarilyAssignInGlobalEnv(foo) ## unfortunately this is necessary for bar to find foo
    expect_message(bar <- nimbleFunction(
        run = function(x = double(1)) {return(foo_oops(x) + 1); returnType(double(1))}
    ), "For this nimbleFunction to compile")
    ## temporarilyAssignInGlobalEnv(foo) ## why is this being done again?
    cat("\nBegin expected error message:\n")
    expect_failure(expect_compiles(bar, dirName = '.', info = 'trying to compile foobar')) ## arguments as for compileNimble
    cat("End expected error message.\n")
    ## uncomment to see an actual failure:
    ##expect_compiles(bar, dirName = '.', info = 'trying to compile foobar')
    ## will bail out with "Error : safely stopping before linking"
    ## alternative format in test_util would use
    ##    expect_compiles(compileNimble(bar))
})

options(warn = RwarnLevel)

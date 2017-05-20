source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of testing tools")

## new way to force anything compiled to be clear-compiled
## (note, no check for windows is done [yet?])
test_that('withTempProjec test',
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

## new way to check only if compilation up to (not including) linking works
## this also by default replaces any -O# flags to -O1
test_that('testing', {    
    foo <- nimbleFunction(
        run = function(x = double(1)) {return(x + 1); returnType(double(1))}
    )
    temporarilyAssignInGlobalEnv(foo) ## unfortunately this is necessary for bar to find foo
    bar <- nimbleFunction(
        run = function(x = double(1)) {return(foo(x) + 1); returnType(double(1))}
    )
    temporarilyAssignInGlobalEnv(foo)
    expect_compiles(bar, dirName = '.' ) ## arguments as for compileNimble
    ## will bail out with "Error : safely stopping before linking"
    ## alternative format in test_util would use
    ##    expect_compiles(compileNimble(bar))
})

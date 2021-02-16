source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
context('Testing return() type error trapping')

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)
nimbleVerboseErrorsSetting <- nimbleOptions('verboseErrors')
nimbleOptions(verboseErrors = FALSE)

## We can use upcoming expect_compiles in the future,
## but for now I'm doing it more crudely

test_that('return type error caught', {
    foo <- nimbleFunction(
        run = function(x = double(1)) {
            return(x + 1)
            returnType(integer(1))
        })
    expect_error(compileNimble(foo))
})

test_that('return nDim error caught', {
    foo <- nimbleFunction(
        run = function(x = double(1)) {
            return(x + 1)
            returnType(double(2))
        })
    expect_error(compileNimble(foo))
})

test_that('return type and nDim error caught', {
    foo <- nimbleFunction(
        run = function(x = double(1)) {
            return(x + 1)
            returnType(integer())
        })
    expect_error(compileNimble(foo))
})

test_that('nimbleList return type error caught', {
    nlConf1 <- nimbleList(nlVector = double(1))
    nlConf2 <- nimbleList(nlVector = double(1)) ## same content, different conf
    temporarilyAssignInGlobalEnv(nlConf1)
    temporarilyAssignInGlobalEnv(nlConf2)
    nf1 <- nimbleFunction(
        run = function(){
            l1 <- nlConf1$new()
            returnType(nlConf2())
            return(l1)
        }
    )    
    expect_error(compileNimble(nf1))
})

test_that('void() return passes return type error trapping', {
    foo <- nimbleFunction(
        run = function(x = double(1)) {
            return()
        })
    cfoo <- compileNimble(foo)
    expect_identical(class(cfoo), "function")
})

test_that('return type error caught when non-void object is returned', {
    foo <- nimbleFunction(
        run = function(x = double(1)) {
            return(x)
            returnType(void())
        })
    expect_error(compileNimble(foo))
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(verboseErrors = nimbleVerboseErrorsSetting)

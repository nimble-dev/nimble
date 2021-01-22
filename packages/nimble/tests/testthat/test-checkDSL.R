source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

# Might want to make a test_dslCheck function.

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context('Testing of nimbleFunction (DSL) code checking')

test_that("Test of DSL check of valid RCfunction", expect_silent(
test1 <- nimbleFunction(
    run = function(x = double(0)) {
        returnType(double())
        y = cos(x)
        tmp <- c(seq(0, 1, by = 0.1), rep(4, 5))
        tmp2 <- logical(5)
        tmp3 <- numeric(3)
        tmp4 <- identityMatrix(3)
        return(y)
    }
    )
))

test_that("Test of DSL check of valid nimbleFunction", expect_silent(
test2 <- nimbleFunction(
    setup = function(model, node) {
        calcNodes <- model$getDependencies(node)
    },
    run = function(x = double(0)) {
        returnType(double())
        y = model$calculate(calcNodes)
        return(y)
    }
    )
))

test_that("Test of DSL check of invalid RCfunction with R code present", expect_warning(
test3 <- nimbleFunction(
    run = function(x = double(1), y = double(1)) {
        returnType(double(1))
        out <- lm(y ~ x)
        return(out$coefficients)
    }
    )
))


test_that("Test of DSL check of invalid RCfunction with check turned off", expect_silent(
test3a <- nimbleFunction(
    run = function(x = double(1), y = double(1)) {
        returnType(double(1))
        out <- lm(y ~ x)
        return(out$coefficients)
    }, check = FALSE
    )
))

# we'd like to catch this but 'aa' could be a nf and don't want to catch that; at the moment the only thing we readily have access to at time of checking is the names of objects in setup, not their types
test_that("Test of DSL check of invalid nimbleFunction with R code present", expect_failure(expect_warning(
test4 <- nimbleFunction(
    setup = function() {
        aa <- function() { return(0) }
    },
    run = function(x = double(0)) {
        returnType(double(0))
        tmp <- aa()
        return(tmp)
    }
)
)))

test_that("Test of DSL check of valid nimbleFunction with other nimbleFunctions present", expect_silent(
test5 <- nimbleFunction(
    setup = function(model, target) {
        calcNodes <- model$getDependencies(target)
        my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
    },
    run = function(par = double(1)) {
        returnType(double(0))
        values(model, target) <<- par
        ans <- model$calculate(calcNodes)
        tmp = my_decideAndJump(0,0,0,0)
        return(tmp)
    }
)
))

test_that("Test of DSL check of valid nimbleFunction with other methods present", expect_silent(
test6 <- nimbleFunction(
    setup = function(model, target) {
        calcNodes <- model$getDependencies(target)
    },
    run = function(par = double(1)) {
        returnType(double(0))
        values(model, target) <<- par
        ans <- model$calculate(calcNodes)
        tmp = auxil()
        return(tmp)
    }, methods = list(
           auxil = function() {
               returnType(double(0))
               return(3)
           })
)
))

test_that("Test of DSL check of valid nimbleFunction using undefined RC function", expect_warning(
test7 <- nimbleFunction(
    setup = function(model, target) {
        calcNodes <- model$getDependencies(target)
    },
    run = function(par = double(1)) {
        returnType(double(0))
        values(model, target) <<- par
        ans <- model$calculate(calcNodes)
        tmp = myRCfun()
        return(tmp)
    })
))

test_that("Test of DSL check of valid nimbleFunction using RC function from setup", expect_silent(
test8 <- nimbleFunction(
    setup = function(model, target) {
        calcNodes <- model$getDependencies(target)
        myRCfun <- nimbleFunction(run = function() {})
    },
    run = function(par = double(1)) {
        returnType(double(0))
        values(model, target) <<- par
        ans <- model$calculate(calcNodes)
        tmp = myRCfun()
        return(tmp)
    })
))

myRCfun <- nimbleFunction(run = function() {})
assign('myRCfun', myRCfun, envir = .GlobalEnv)

test_that("Test of DSL check of valid nimbleFunction using RC function from R", expect_silent( 
test9 <- nimbleFunction(
    setup = function(model, target) {
        calcNodes <- model$getDependencies(target)
    },
    run = function(par = double(1)) {
        returnType(double(0))
        values(model, target) <<- par
        ans <- model$calculate(calcNodes)
        tmp = myRCfun()
        return(tmp)
    })
))

test_that("Test of DSL check of valid nimbleFunction with undefined nimbleList", expect_warning(
test10 <- nimbleFunction(
  setup = function(){    
  },
  run = function(){
    outList <- method1()
    returnType(testListDef())
     return(outList)
  },
  methods = list(
    method1 = function(){
      outList <- testListDef$new(nlCharacter = 'hello world')
       returnType(testListDef())
      return(outList)
    }
  )
)
))

test_that("Test of DSL check of valid nimbleFunction with nimbleList in setup", expect_silent(
test11 <- nimbleFunction(
  setup = function(){
      testListDef <- nimbleList(nlScalar = double(0), nlVector = double(1), nlMatrix = double(2))
  },
  run = function(){
    outList <- method1()
    returnType(testListDef())
     return(outList)
  },
  methods = list(
    method1 = function(){
      outList <- testListDef$new(nlCharacter = 'hello world')
       returnType(testListDef())
      return(outList)
    }
  )
)
))

testListDef <- nimbleList(nlScalar = double(0), nlVector = double(1), nlMatrix = double(2))
assign('testListDef', testListDef, envir = .GlobalEnv)
test_that("Test of DSL check of valid nimbleFunction with nimbleList in R", expect_silent(
test12 <- nimbleFunction(
  setup = function(){
  },
  run = function(){
    outList <- method1()
    returnType(testListDef())
     return(outList)
  },
  methods = list(
    method1 = function(){
      outList <- testListDef$new(nlCharacter = 'hello world')
       returnType(testListDef())
      return(outList)
    }
  )
  )))

test_that("Test of DSL check of valid nimbleFunction with call to undefined nimbleFunction", expect_warning(
test13 <- nimbleFunction(
    setup = function(model, target) {
        calcNodes <- model$getDependencies(target)
    },
    run = function(par = double(1)) {
        returnType(double(0))
        values(model, target) <<- par
        ans <- model$calculate(calcNodes)
        tmp = mynf()
        return(tmp)
    })
))


test_that("Test of DSL check of valid nimbleFunction using nimbleFunction from setup", expect_silent({
    nfGen <- nimbleFunction(
    setup = function() {},
    run = function() {}
    );
test14 <- nimbleFunction(
    setup = function(model, target) {
        calcNodes <- model$getDependencies(target)
        mynf <- nfGen()
    },
    run = function(par = double(1)) {
        returnType(double(0))
        values(model, target) <<- par
        ans <- model$calculate(calcNodes)
        tmp = mynf()
        return(tmp)
    })
}
))


options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)

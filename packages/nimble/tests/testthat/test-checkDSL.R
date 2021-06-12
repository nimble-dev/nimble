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

## Actually not any different than test7.
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

test_that("Test of DSL check of valid nimbleFunction with calls to undefined nimbleFunctions and list elements", {
    nf_base= nimbleFunctionVirtual(
        name = 'nfbase',
        run = function() {returnType(double(0))},
        methods = list(
            foo = function() { }
        )
    )

    nfextend = nimbleFunction(
        contains = nf_base,
        setup = function() {},
        run = function() {
            returnType(double(0))
            return(3)
        },
        methods = list(
            foo2 = function() {},
            foo3 = function() {}
        )
    )

    expect_silent(
        nf_more <- nimbleFunction(
            setup = function() {
                fxns = nimbleFunctionList(nf_base)
                fxns[[1]] <- nf_extend()
            },
            run = function() {
                returnType(double(0))
                fxns[[1]]$foo()
                fxns[[1]]$foo2()
                i <- 1
                fxns[[i]]$foo3()
                out = fxns[[1]]$run()
                fxns[[i]]$foo4()
                return(out)
            }))
})

test_that("Test of detection of nf methods in findMethods", {
    code <- quote({
        m$cc1
        m$cc2()
        m$cc3(a)
        m[[1]]$cc4
        m[[i]]$cc5
        m[[1]]$cc6()
        m[[i]]$cc7()
        m[[x+3]]$cc8
        m[[x+3]]$cc9()
        m[[1]]$cc10(mm$cc11)
        m[[1]]$cc12(mm$cc13())
        m[[mm$cc14]]$cc15()
        m[[mm$cc16()]]$cc17()
    })
    expr <- RparseTree2ExprClasses(code)
    expect_identical(findMethodsInExprClass(expr), c('cc2', 'cc3', 'cc6', 'cc7', 'cc9', 'cc10', 'cc12', 'cc13', 'cc15', 'cc17', 'cc16'))
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)

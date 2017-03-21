source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

# need a test_dslCheck function???

nimbleOptions(showCompilerOutput = TRUE)

test_that("Test of DSL check of valid RCfunction", expect_silent(
test1 <- nimbleFunction(
    run = function(x = double(0)) {
        returnType(double())
        y = cos(x)
        tmp <- c(seq(0, 1, by = 0.1), rep(4, 5))
        tmp2 <- logical(5)
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

test_that("Test of DSL check of invalid nimbleFunction with R code present", expect_warning(
test4 <- nimbleFunction(
    setup = function(model, target, mvSaved) {
        aa <- function() { return(0) }
    },
    run = function(x = double(0)) {
        returnType(double(0))
        tmp <- aa()
        return(tmp)
    }
)
))

test_that("Test of DSL check of valid nimbleFunction with other nimbleFunctions present", expect_warning(
test5 <- nimbleFunction(
    setup = function(model, target, mvSaved) {
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
test5 <- nimbleFunction(
    setup = function(model, target, mvSaved) {
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

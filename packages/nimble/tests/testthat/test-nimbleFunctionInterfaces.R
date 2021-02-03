source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context("Testing of nimbleFunction interfaces")

## goals here:
##   passing every kind of argument
##   populating every kind of member data upon instantiation
##   getting and setting every kind of member data via full interface and multi-interface
##
##  By "every kind of argument", we mean [scalar, non-scalar] x [bool, int, double, string]
##  (So we are not testing nimbleLists or nested nimbleFunctions here)

## passing every kind of argument
test_that("copying of argument passing in nimbleFunctionInterface", {
    f1 <- nimbleFunction(
        run = function(bs = logical(),
                       b1 = logical(1),
                       b2 = logical(2),
                       is = integer(),
                       i1 = integer(1),
                       i2 = integer(2),
                       ds = double(),
                       d1 = double(1),
                       d2 = double(2),
                       si = character(),
                   s2 = character(1)) {
            returnType(double())
            return(bs + sum(b1) + sum(b2) + is + sum(i1) + sum(i2) + ds + sum(d1) + sum(d2))
        })
    
    cf1 <- compileNimble(f1)
    
    ans <- f1(TRUE,
              c(FALSE, TRUE, TRUE),
              matrix(c(FALSE, TRUE, TRUE, TRUE), nrow = 2),
              100,
              10107:10111,
              matrix(1234:1237, nrow = 2),
              4315.3512,
              seq(15043.234, 15096, length = 3),
              matrix(seq(492.1, 593.9, length = 4), nrow = 2),
              'hello',
              c('hw','goodbye'))
    
    cans <- cf1(TRUE,
                c(FALSE, TRUE, TRUE),
                matrix(c(FALSE, TRUE, TRUE, TRUE), nrow = 2),
                100,
                10107:10111,
                matrix(1234:1237, nrow = 2),
                4315.3512,
                seq(15043.234, 15096, length = 3),
                matrix(seq(492.1, 593.9, length = 4), nrow = 2),
                'hello',
                c('hw','goodbye'))
    
    
    expect_identical(ans, cans)
})

## populating and using every kind of argument with a full interface
    
test_that("population and getting of member data in full nimbleFunctionInterface", {
    f2G <- nimbleFunction(
        setup = function() {
            bs <- TRUE
            b1 <- c(FALSE, TRUE, TRUE)
            b2 <- matrix(c(FALSE, TRUE, TRUE, TRUE), nrow = 2)
            is <- as.integer(100)
            i1 <- as.integer(10107:10111)
            i2 <- matrix(as.integer(1234:1237), nrow = 2)
            ds <- 4315.3512
            d1 <- seq(15043.234, 15096, length = 3)
            d2 <- matrix(seq(492.1, 593.9, length = 4), nrow = 2)
            si <- 'hello'
            s1 <- c('hw','goodbye')
            setupOutputs(si, s1)
        },
        run = function() {
            returnType(double())
            return(bs + sum(b1) + sum(b2) + is + sum(i1) + sum(i2) + ds + sum(d1) + sum(d2))
        })
    
    f2 <- f2G()
    ans2 <- f2$run()    
    cf2 <- compileNimble(f2)
    cans2 <- cf2$run()
    expect_identical(ans2, cans2)

    f2$bs <- cf2$bs <- FALSE
    f2$b1 <- cf2$b1 <- c(FALSE, TRUE, FALSE, FALSE)
    f2$b2 <- cf2$b2 <- matrix(rep(c(FALSE, TRUE, FALSE), 3), nrow = 3)
    f2$is <- cf2$is <- as.integer(200)
    f2$i1 <- cf2$i1 <- as.integer(2003:2010)
    f2$i2 <- cf2$i2 <-  matrix(as.integer(4321 + 1:9), nrow = 3)
    f2$ds <- cf2$ds <- 843.1
    f2$d1 <- cf2$d1 <- seq(5834.1, 60654.2, length = 4)
    f2$d2 <- cf2$d2 <-  matrix(seq(4912.1, 5593.9, length = 9), nrow = 3)
    
    f2$si <- cf2$si <- 'hello again'
    f2$s1 <- cf2$s1 <- c('hw again', 'goodbye again', 'and another thing')
    ans2b <- f2$run()
    cans2b <- cf2$run()
    expect_equal(ans2b, cans2b)
    
    ## these test the getters
    f2$bs <- cf2$bs
    f2$b1 <- cf2$b1
    f2$b2 <- cf2$b2
    f2$is <- cf2$is
    f2$i1 <- cf2$i1
    f2$i2 <- cf2$i2
    f2$ds <- cf2$ds
    f2$d1 <- cf2$d1
    f2$d2 <- cf2$d2
    
    f2$si <- cf2$si
    f2$s1 <- cf2$s1
    ans2c <- f2$run()
      
    expect_equal(ans2c, cans2b)
})
    
  

## populating and using every kind of argument with a multi interface
test_that("populating member data in multi nimbleFunctionInterface", {
    f2G <- nimbleFunction(
        setup = function() {
            bs <- TRUE
            b1 <- c(FALSE, TRUE, TRUE)
            b2 <- matrix(c(FALSE, TRUE, TRUE, TRUE), nrow = 2)
            is <- as.integer(100)
            i1 <- as.integer(10107:10111)
            i2 <- matrix(as.integer(1234:1237), nrow = 2)
            ds <- 4315.3512
            d1 <- seq(15043.234, 15096, length = 3)
            d2 <- matrix(seq(492.1, 593.9, length = 4), nrow = 2)
            si <- 'hello'
            s1 <- c('hw','goodbye')
            setupOutputs(si, s1)
        },
        run = function() {
            returnType(double())
            return(bs + sum(b1) + sum(b2) + is + sum(i1) + sum(i2) + ds + sum(d1) + sum(d2))
        })
    
    f3G <- nimbleFunction(
        setup = function(f2G) {
            myf2 <- f2G() ## this will have a nested interface
        },
        run = function() {
            ans <- myf2$run()
            return(ans)
            returnType(double())
        })
    
    f3 <- f3G(f2G)
    cf3 <- compileNimble(f3)
    ans3 <- f3$run()
    cans3 <- cf3$run()

    expect_equal(ans3, cans3)

    myf2 <- f3$myf2
    cmyf2 <- cf3$myf2
    
    valueInCompiledNimbleFunction(cmyf2, 'bs', FALSE)
    myf2$bs <- valueInCompiledNimbleFunction(cmyf2, 'bs')
    
    valueInCompiledNimbleFunction(cmyf2, 'b1', c(FALSE, TRUE, FALSE, FALSE))
    myf2$b1 <- valueInCompiledNimbleFunction(cmyf2, 'b1')
    
    valueInCompiledNimbleFunction(cmyf2, 'b2', matrix(rep(c(FALSE, TRUE, FALSE), 3), nrow = 3))
    myf2$b2 <- valueInCompiledNimbleFunction(cmyf2, 'b2')
    
    valueInCompiledNimbleFunction(cmyf2, 'is', 200)
    myf2$is <- valueInCompiledNimbleFunction(cmyf2, 'is')
    
    valueInCompiledNimbleFunction(cmyf2, 'i1', as.integer(2003:2010))
    myf2$i1 <- valueInCompiledNimbleFunction(cmyf2, 'i1')
    
    valueInCompiledNimbleFunction(cmyf2, 'i2', matrix(as.integer(4321 + 1:9), nrow = 3))
    myf2$i2 <- valueInCompiledNimbleFunction(cmyf2, 'i2')
    
    valueInCompiledNimbleFunction(cmyf2, 'ds', 843.1)
    myf2$ds <- valueInCompiledNimbleFunction(cmyf2, 'ds')
    
    valueInCompiledNimbleFunction(cmyf2, 'd1', seq(5834.1, 60654.2, length = 4))
    myf2$d1 <- valueInCompiledNimbleFunction(cmyf2, 'd1')
    
    valueInCompiledNimbleFunction(cmyf2, 'd2', matrix(seq(4912.1, 5593.9, length = 9), nrow = 3))
    myf2$d2 <- valueInCompiledNimbleFunction(cmyf2, 'd2')
    
    valueInCompiledNimbleFunction(cmyf2, 'si', 'hello again')
    myf2$si <- valueInCompiledNimbleFunction(cmyf2, 'si')
    
    valueInCompiledNimbleFunction(cmyf2, 's1', c('hw again', 'goodbye again', 'and another thing'))
    myf2$s1 <-valueInCompiledNimbleFunction(cmyf2, 's1')
    
    ans3b <- f3$run()
    cans3b <- cf3$run()
    
    expect_equal(ans3b, cans3b)
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)

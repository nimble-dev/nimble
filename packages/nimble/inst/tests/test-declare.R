source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context("Testing of declare, setSize, numeric, integer, matrix, and array")

## I. 1D cases

## I.i 1D declare without sizes, then resize without c()

test_that("1D declare without sizes, then resize without c()", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function() {
            declare(z, double(1))
            setSize(z, 3)
            for(i in 1:3) z[i] = i
            return(z)
            returnType(double(1))
        })
    
    nf1 <- nf()
    expect_equivalent(nf1$run(), 1:3, info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(), 1:3, info = "compiled execution problem")
})

## I.ii 1D declare without sizes, then resize with c()

test_that("1D declare without sizes, then resize with c()", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function() {
            declare(z, double(1))
            setSize(z, c(3))
            for(i in 1:3) z[i] = i
            return(z)
            returnType(double(1))
        })
    
    nf1 <- nf()
    expect_equivalent(nf1$run(), 1:3, info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(), 1:3, info = "compiled execution problem")
})

## I.iii 1D declare with sizes, without c()
test_that("1D declare with sizes, without c()", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function() {
            declare(z, double(1, 3))
            for(i in 1:3) z[i] = i
            return(z)
            returnType(double(1))
        })

    nf1 <- nf()
    expect_equivalent(nf1$run(), array(1:3, dim=c(3)), info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(), 1:3, info = "compiled execution problem")
})
    
    ## I.iv 1D declare with sizes, with c()
test_that("1D declare with sizes, with c()", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function() {
            declare(z, double(1, c(3)))
            for(i in 1:3) z[i] = i
            return(z)
            returnType(double(1))
        })
    
    nf1 <- nf()
    expect_equivalent(nf1$run(), array(1:3, dim=c(3)), info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(), 1:3, info = "compiled execution problem")
})

## I.v 1D declare with sizes, with c() using a variable
test_that("1D declare with sizes, with c() using a variable", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function(s1 = integer()) {
            declare(z, double(1, c(s1)))
            for(i in 1:3) z[i] = i
            return(z)
            returnType(double(1))
        })
    
    nf1 <- nf()
    expect_equivalent(nf1$run(3), array(1:3, dim=c(3)), info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(3), 1:3, info = "compiled execution problem")
})


## I.vi 1D declare with sizes, followed by resize and size inference
test_that("1D declare with sizes, followed by resize and size inference", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function(m = double(2)) {
            declare(z, double(1, c(2)))
            setSize(z, 3)
            for(i in 1:3) z[i] = i
            ans <- m %*% z
            return(ans)
            returnType(double(2))
        })
    
    nf1 <- nf()
    m <- matrix(1:9, nrow = 3)
    expect_equivalent(nf1$run(m), m %*% 1:3, info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(m), m %*% 1:3, info = "compiled execution problem")
})

## II. 2D cases
## II.i 2D declare without sizes, then resize without c()

test_that("2D declare without sizes, then resize without c()", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function() {
            declare(z, double(2))
            setSize(z, 2, 4)
            for(i in 1:2) for(j in 1:4) z[i,j] = (j)*2 + (i) - 2
            return(z)
            returnType(double(2))
        })

    nf1 <- nf()
    expect_equivalent(nf1$run(), matrix(1:8, nrow = 2), info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(), matrix(1:8, nrow = 2), info = "compiled execution problem")
})

## II.ii 2D declare without sizes, then resize with c()

test_that("2D declare without sizes, then resize with c()", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function() {
            declare(z, double(2))
            setSize(z, c(2, 4))
            for(i in 1:2) for(j in 1:4) z[i,j] = (j)*2 + (i) - 2
            return(z)
            returnType(double(2))
        })

    nf1 <- nf()
    expect_equivalent(nf1$run(), matrix(1:8, nrow = 2), info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(), matrix(1:8, nrow = 2), info = "compiled execution problem")
})

## II.iii 2D declare with sizes, without c()
test_that("2D declare with sizes, without c()", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function() {
            declare(z, double(2, 2, 4))
            for(i in 1:2) for(j in 1:4) z[i,j] = (j)*2 + (i) - 2
            return(z)
            returnType(double(2))
        })
    
    nf1 <- nf()
    expect_equivalent(nf1$run(), matrix(1:8, nrow = 2), info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(), matrix(1:8, nrow = 2), info = "compiled execution problem")
})

## II.iv 2D declare with sizes, with c()
test_that("2D declare with sizes, with c()", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function() {
            declare(z, double(2, c(2, 4)))
            for(i in 1:2) for(j in 1:4) z[i,j] = (j)*2 + (i) - 2
            return(z)
            returnType(double(2))
        })
    
    nf1 <- nf()
    expect_equivalent(nf1$run(), matrix(1:8, nrow = 2), info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(), matrix(1:8, nrow = 2), info = "compiled execution problem")
})

## II.v 2D declare with sizes, with c() using a variable
test_that("2D declare with sizes, with c() using a variable", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function(s1 = integer()) {
            declare(z, double(2, c(2, s1)))
            for(i in 1:2) for(j in 1:4) z[i,j] = (j)*2 + (i) - 2
            return(z)
            returnType(double(2))
        })
    
    nf1 <- nf()
    expect_equivalent(nf1$run(4), matrix(1:8, nrow = 2), info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(4), matrix(1:8, nrow = 2), info = "compiled execution problem")
})


## II.vi 2D declare with sizes, followed by resize and size inference
test_that("2D declare with sizes, followed by resize and size inference", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function(s1 = integer(), m = double(2)) {
            declare(z, double(2, c(2, 2)))
            setSize(z, c(2, s1))
            for(i in 1:2) for(j in 1:4) z[i,j] = (j)*2 + (i) - 2
            ans <- m %*% z
            return(ans)
            returnType(double(2))
        })

    nf1 <- nf()
    m <- matrix(11:18, nrow = 4)
    z <- matrix(1:8, nrow = 2)
    expect_equivalent(nf1$run(4, m), m %*% z, info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(4, m), m %*% z, info = "compiled execution problem")
})
    
## III 3D cases
## III.i 3D declare without sizes, then resize without c()

test_that("3D declare without sizes, then resize without c()", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function() {
            declare(z, double(3))
            setSize(z, 2, 4, 3)
            for(i in 1:2) for(j in 1:4) for(k in 1:3) z[i, j, k] = (k)*2*4 + (j)*2 + (i) - 10
            return(z)
            returnType(double(3))
        })
    
    nf1 <- nf()
    expect_equivalent(nf1$run(), array(1:24, dim = c(2,4,3)), info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(), array(1:24, dim = c(2,4,3)), info = "compiled execution problem") 
})

## III.ii 3D declare without sizes, then resize with c()

test_that("3D declare without sizes, then resize with c()", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function() {
            declare(z, double(3))
            setSize(z, c(2, 4, 3))
            for(i in 1:2) for(j in 1:4) for(k in 1:3) z[i, j, k] = (k)*2*4 + (j)*2 + (i) - 10
            return(z)
            returnType(double(3))
        })
    
    nf1 <- nf()
    expect_equivalent(nf1$run(), array(1:24, dim = c(2,4,3)), info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(), array(1:24, dim = c(2,4,3)), info = "compiled execution problem") 
})

## III.iii 3D declare with sizes, without c()
test_that("3D declare with sizes, without c()", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function() {
            declare(z, double(3, 2, 4, 3))
            for(i in 1:2) for(j in 1:4) for(k in 1:3) z[i, j, k] = (k)*2*4 + (j)*2 + (i) - 10
            return(z)
            returnType(double(3))
        })
    
    nf1 <- nf()
    expect_equivalent(nf1$run(), array(1:24, dim = c(2,4,3)), info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(), array(1:24, dim = c(2,4,3)), info = "compiled execution problem") 
})

## III.iv 3D declare with sizes, with c()
test_that("3D declare with sizes, with c()" , {
    nf <- nimbleFunction(
        setup = function(){},
        run = function() {
            declare(z, double(3, c(2, 4, 3)))
            for(i in 1:2) for(j in 1:4) for(k in 1:3) z[i, j, k] = (k)*2*4 + (j)*2 + (i) - 10
            return(z)
            returnType(double(3))
        })
    
    nf1 <- nf()
    expect_equivalent(nf1$run(), array(1:24, dim = c(2,4,3)), info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(), array(1:24, dim = c(2,4,3)), info = "compiled execution problem") 
})

## III.v 3D declare with sizes, with c() using a variable
test_that("3D declare with sizes, with c() using a variable", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function(s1 = integer()) {
            declare(z, double(3, c(2, s1, 3)))
            for(i in 1:2) for(j in 1:4) for(k in 1:3) z[i, j, k] = (k)*2*4 + (j)*2 + (i) - 10
            return(z)
            returnType(double(3))
        })
    
    nf1 <- nf()
    expect_equivalent(nf1$run(4), array(1:24, dim = c(2,4,3)), info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(4), array(1:24, dim = c(2,4,3)), info = "compiled execution problem") 
})


## III.vi 3D declare with sizes, followed by resize and size inference
test_that("3D declare with sizes, followed by resize and size inferenc", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function(s1 = integer(), m = double(2)) {
            declare(z, double(3, c(1, 1, 2)))
            setSize(z, c(2, s1, 3))
            for(i in 1:2) for(j in 1:4) for(k in 1:3) z[i, j, k] = (k)*2*4 + (j)*2 + (i) - 10
            ans <- m %*% z[,,1]
            return(ans)
            returnType(double(2))
        })
    
    nf1 <- nf()
    m <- matrix(11:18, nrow = 4)
    z <- array(1:24, dim = c(2,4,3))
    expect_equivalent(nf1$run(4, m), m %*% z[,,1] , info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(4, m), m %*% z[,,1], info = "compiled execution problem") 
})

## IV. type declarations in arguments
## For this only the c() notation (or a single number for 1D case) works for sizes

## IV.i
test_that("type declarations in arguments, first case", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function(z = double(1)) {
            setSize(z, 3)
            for(i in 1:3) z[i] = i
            return(z)
            returnType(double(1))
        })
    
    nf1 <- nf()
    z <- 1:2
    expect_equivalent(nf1$run(z), 1:3, info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(z), 1:3, info = "compiled execution problem") 
})

## IV.ii ## this case arguably shouldn't work, but we allow it
test_that("type declarations in arguments, second case, allowed but arguable", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function(z = double(1, 4)) {
            setSize(z, 3)
            for(i in 1:3) z[i] = i
            return(z)
            returnType(double(1))
        })
    
    nf1 <- nf()
    z <- 1:2
    expect_equivalent(nf1$run(z), 1:3, info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(z), 1:3, info = "compiled execution problem") 
})


## IV.iii
test_that("type declarations in arguments, third case", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function(z = double(2, c(2, 4))) {
            for(i in 1:2) for(j in 1:4) z[i,j] = (j)*2 + (i) - 2
            return(z)
            returnType(double(2))
        })
    
    nf1 <- nf()
    z <- matrix(0, nrow = 2, ncol = 4)
    expect_equivalent(nf1$run(z),  matrix(1:8, nrow = 2), info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(z),  matrix(1:8, nrow = 2), info = "compiled execution problem") 
})

## IV.iv
test_that("type declarations in arguments, fourth case", {
    nf <- nimbleFunction(
        setup = function(){},
        run = function(z = double(0, default = 10)) {
            return(z/2)
            returnType(double())
        })
    
    nf1 <- nf()
    expect_equivalent(nf1$run(),  5, info = "R execution problem")
    cnf1 <- compileNimble(nf1)
    expect_equivalent(cnf1$run(),  5, info = "compiled execution problem") 
})

## testing of numeric(), integer(), matrix(), and array()

test_that("numeric cases 1-4", {
    expected <- numeric(10)
    Rfun <- nimbleFunction(run = function() {
        ans <- numeric(10)
        returnType(double(1))
        return(ans)
    })
    Cfun <- compileNimble(Rfun)
    expect_equal(Rfun(), expected, info = "case 1")
    expect_equal(Cfun(), expected, info = "case 2")
    expect_identical(class(Rfun()[1]), 'numeric', info = "case 3")
    expect_identical(class(Cfun()[1]), 'numeric', info = "case 4")
})

test_that("numeric cases 5-8", {
    expected <- rep(3, length = 2)
    Rfun <- nimbleFunction(run = function() {
        ans <- numeric(value = 3, length = 2)
        returnType(double(1))
        return(ans)
    })
    Cfun <- compileNimble(Rfun)
    expect_equal(Rfun(), expected, info = "case 5")
    expect_equal(Cfun(), expected, info = "case 6")
    expect_identical(class(Rfun()[1]), 'numeric', info = "case 7")
    expect_identical(class(Cfun()[1]), 'numeric', info = "case 8")
})

test_that("integer cases 9-12", {
    expected <- rep(9, 3)
    Rfun <- nimbleFunction(run = function() {
        x <- numeric(10, value = 3)
        ans <- integer(x[2], x[2]+x[3]*2)
        returnType(integer(1))
        return(ans)
    })
    Cfun <- compileNimble(Rfun)
    expect_equal(Rfun(), expected, info = "case 9")
    expect_equal(Cfun(), expected, info = "case 10")
    expect_identical(class(Rfun()[1]), 'integer', info = "case 11")
    expect_identical(class(Cfun()[1]), 'integer', info = "case 12")
})

test_that("numeric cases 13-16", {
    expected <- array(4, c(10,11))
    Rfun <- nimbleFunction(run = function() {
        ans <- array(4, c(10,11))
        returnType(double(2))
        return(ans)
    })
    Cfun <- compileNimble(Rfun)
    expect_equal(Rfun(), expected, info = "case 13")
    expect_equal(Cfun(), expected, info = "case 14")
    expect_identical(class(Rfun()[1]), 'numeric', info = "case 15")
    expect_identical(class(Cfun()[1]), 'numeric', info = "case 16")
})

test_that("integer cases 17-20", {
    expected <- matrix(as.integer(0), nrow=4, ncol=5)
    Rfun <- nimbleFunction(run = function() {
        x <- 4
        y <- 5
        ans <- matrix(init=FALSE, nrow=x, ncol=y, type='integer')
        returnType(integer(2))
        return(ans)
    })
    Cfun <- compileNimble(Rfun)
    expect_equal(Rfun(), expected, info = "case 17")
    expect_equal(dim(Cfun()), c(4, 5), info = "case 18")
    expect_identical(class(Rfun()[1,1]), 'integer', info = "case 19")
    expect_identical(class(Cfun()[1,1]), 'integer', info = "case 20")
})

test_that("numeric cases 21-24", {
    expected <- c(99, 10.5, 10.5, 10.5)
    Rfun <- nimbleFunction(run = function() {
        a <- numeric(4)
        a[1] <- 99
        b <- numeric(value=.5, 4)
        a[2:4] <- b[1:3] + numeric(3, 10)
        returnType(double(1))
        return(a)
    })
    Cfun <- compileNimble(Rfun)
    expect_equal(Rfun(), expected, info = "case 21")
    expect_equal(Cfun(), expected, info = "case 22")
    expect_identical(class(Rfun()[1]), 'numeric', info = "case 23")
    expect_identical(class(Cfun()[1]), 'numeric', info = "case 24")
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)


source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of declare, setSize, numeric, integer, matrix, and array")

caseCounter <- 1
makeCaseMsg <- function(i, xtra = character()) paste("declare and setSize testing case",i,xtra )
## I. 1D cases

## I.i 1D declare without sizes, then resize without c()

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
expect_equivalent(nf1$run(), 1:3, makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(), 1:3, makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1

## I.ii 1D declare without sizes, then resize with c()

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
expect_equivalent(nf1$run(), 1:3, makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(), 1:3, makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1

## I.iii 1D declare with sizes, without c()
nf <- nimbleFunction(
    setup = function(){},
    run = function() {
        declare(z, double(1, 3))
        for(i in 1:3) z[i] = i
        return(z)
        returnType(double(1))
    })

nf1 <- nf()
expect_equivalent(nf1$run(), array(1:3, dim=c(3)) , makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(), 1:3, makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1

## I.iv 1D declare with sizes, with c()
nf <- nimbleFunction(
    setup = function(){},
    run = function() {
        declare(z, double(1, c(3)))
        for(i in 1:3) z[i] = i
        return(z)
        returnType(double(1))
    })

nf1 <- nf()
expect_equivalent(nf1$run(), array(1:3, dim=c(3)) , makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(), 1:3, makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1

## I.v 1D declare with sizes, with c() using a variable
nf <- nimbleFunction(
    setup = function(){},
    run = function(s1 = integer()) {
        declare(z, double(1, c(s1)))
        for(i in 1:3) z[i] = i
        return(z)
        returnType(double(1))
    })

nf1 <- nf()
expect_equivalent(nf1$run(3), array(1:3, dim=c(3)) , makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(3), 1:3, makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1


## I.vi 1D declare with sizes, followed by resize and size inference
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
expect_equivalent(nf1$run(m), m %*% 1:3 , makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(m), m %*% 1:3, makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1

## II. 2D cases
## II.i 2D declare without sizes, then resize without c()

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
expect_equivalent(nf1$run(), matrix(1:8, nrow = 2), makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(), matrix(1:8, nrow = 2), makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1

## II.ii 2D declare without sizes, then resize with c()

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
expect_equivalent(nf1$run(), matrix(1:8, nrow = 2), makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(), matrix(1:8, nrow = 2), makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1

## II.iii 2D declare with sizes, without c()
nf <- nimbleFunction(
    setup = function(){},
    run = function() {
        declare(z, double(2, 2, 4))
        for(i in 1:2) for(j in 1:4) z[i,j] = (j)*2 + (i) - 2
        return(z)
        returnType(double(2))
    })

nf1 <- nf()
expect_equivalent(nf1$run(), matrix(1:8, nrow = 2) , makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(), matrix(1:8, nrow = 2), makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1

## II.iv 2D declare with sizes, with c()
nf <- nimbleFunction(
    setup = function(){},
    run = function() {
        declare(z, double(2, c(2, 4)))
        for(i in 1:2) for(j in 1:4) z[i,j] = (j)*2 + (i) - 2
        return(z)
        returnType(double(2))
    })

nf1 <- nf()
expect_equivalent(nf1$run(), matrix(1:8, nrow = 2) , makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(), matrix(1:8, nrow = 2), makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1

## II.v 2D declare with sizes, with c() using a variable
nf <- nimbleFunction(
    setup = function(){},
    run = function(s1 = integer()) {
        declare(z, double(2, c(2, s1)))
        for(i in 1:2) for(j in 1:4) z[i,j] = (j)*2 + (i) - 2
        return(z)
        returnType(double(2))
    })

nf1 <- nf()
expect_equivalent(nf1$run(4), matrix(1:8, nrow = 2) , makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(4), matrix(1:8, nrow = 2), makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1


## II.vi 2D declare with sizes, followed by resize and size inference
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
expect_equivalent(nf1$run(4, m), m %*% z , makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(4, m), m %*% z, makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1

## III 3D cases
## III.i 3D declare without sizes, then resize without c()

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
expect_equivalent(nf1$run(), array(1:24, dim = c(2,4,3)), makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(), array(1:24, dim = c(2,4,3)), makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1

## III.ii 3D declare without sizes, then resize with c()

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
expect_equivalent(nf1$run(), array(1:24, dim = c(2,4,3)), makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(), array(1:24, dim = c(2,4,3)), makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1

## III.iii 3D declare with sizes, without c()
nf <- nimbleFunction(
    setup = function(){},
    run = function() {
        declare(z, double(3, 2, 4, 3))
        for(i in 1:2) for(j in 1:4) for(k in 1:3) z[i, j, k] = (k)*2*4 + (j)*2 + (i) - 10
        return(z)
        returnType(double(3))
    })

nf1 <- nf()
expect_equivalent(nf1$run(), array(1:24, dim = c(2,4,3)), makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(), array(1:24, dim = c(2,4,3)), makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1

## III.iv 3D declare with sizes, with c()
nf <- nimbleFunction(
    setup = function(){},
    run = function() {
        declare(z, double(3, c(2, 4, 3)))
        for(i in 1:2) for(j in 1:4) for(k in 1:3) z[i, j, k] = (k)*2*4 + (j)*2 + (i) - 10
        return(z)
        returnType(double(3))
    })

nf1 <- nf()
expect_equivalent(nf1$run(), array(1:24, dim = c(2,4,3)), makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(), array(1:24, dim = c(2,4,3)), makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1

## III.v 3D declare with sizes, with c() using a variable
nf <- nimbleFunction(
    setup = function(){},
    run = function(s1 = integer()) {
        declare(z, double(3, c(2, s1, 3)))
        for(i in 1:2) for(j in 1:4) for(k in 1:3) z[i, j, k] = (k)*2*4 + (j)*2 + (i) - 10
        return(z)
        returnType(double(3))
    })

nf1 <- nf()
expect_equivalent(nf1$run(4), array(1:24, dim = c(2,4,3)), makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(4), array(1:24, dim = c(2,4,3)), makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1


## III.vi 3D declare with sizes, followed by resize and size inference
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
expect_equivalent(nf1$run(4, m), m %*% z[,,1] , makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(4, m), m %*% z[,,1], makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1

## IV. type declarations in arguments
## For this only the c() notation (or a single number for 1D case) works for sizes

## IV.i
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
expect_equivalent(nf1$run(z), 1:3, makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(z), 1:3, makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1

## IV.ii ## this case arguably shouldn't work, but we allow it
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
expect_equivalent(nf1$run(z), 1:3, makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(z), 1:3, makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1


## IV.iii
nf <- nimbleFunction(
    setup = function(){},
    run = function(z = double(2, c(2, 4))) {
        for(i in 1:2) for(j in 1:4) z[i,j] = (j)*2 + (i) - 2
        return(z)
        returnType(double(2))
    })

nf1 <- nf()
z <- matrix(0, nrow = 2, ncol = 4)
expect_equivalent(nf1$run(z),  matrix(1:8, nrow = 2), makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(z),  matrix(1:8, nrow = 2), makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1

## IV.iv
nf <- nimbleFunction(
    setup = function(){},
    run = function(z = double(0, default = 10)) {
        return(z/2)
        returnType(double())
    })

nf1 <- nf()
expect_equivalent(nf1$run(),  5, makeCaseMsg(caseCounter))
cnf1 <- compileNimble(nf1)
expect_equivalent(cnf1$run(),  5, makeCaseMsg(caseCounter, 'compiled')) 
caseCounter <- caseCounter + 1



## testing of numeric(), integer(), matrix(), and array()

expected <- numeric(10)
Rfun <- nimbleFunction(run = function() {
    ans <- numeric(10)
    returnType(double(1))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('numeric', expect_equal(Rfun(), expected))
test_that('numeric', expect_equal(Cfun(), expected))
test_that('numeric', expect_identical(class(Rfun()[1]), 'numeric'))
test_that('numeric', expect_identical(class(Cfun()[1]), 'numeric'))


expected <- rep(3, length = 2)
Rfun <- nimbleFunction(run = function() {
    ans <- numeric(value = 3, length = 2)
    returnType(double(1))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('numeric', expect_equal(Rfun(), expected))
test_that('numeric', expect_equal(Cfun(), expected))
test_that('numeric', expect_identical(class(Rfun()[1]), 'numeric'))
test_that('numeric', expect_identical(class(Cfun()[1]), 'numeric'))

expected <- rep(9, 3)
Rfun <- nimbleFunction(run = function() {
    x <- numeric(10, value = 3)
    ans <- integer(x[2], x[2]+x[3]*2)
    returnType(integer(1))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('integer', expect_equal(Rfun(), expected))
test_that('integer', expect_equal(Cfun(), expected))
test_that('integer', expect_identical(class(Rfun()[1]), 'integer'))
test_that('integer', expect_identical(class(Cfun()[1]), 'integer'))

expected <- array(4, c(10,11))
Rfun <- nimbleFunction(run = function() {
    ans <- array(4, c(10,11))
    returnType(double(2))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('integer', expect_equal(Rfun(), expected))
test_that('integer', expect_equal(Cfun(), expected))
test_that('integer', expect_identical(class(Rfun()[1]), 'numeric'))
test_that('integer', expect_identical(class(Cfun()[1]), 'numeric'))

expected <- matrix(as.integer(0), nrow=4, ncol=5)
Rfun <- nimbleFunction(run = function() {
    x <- 4
    y <- 5
    ans <- matrix(init=FALSE, nrow=x, ncol=y, type='integer')
    returnType(integer(2))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('integer', expect_equal(Rfun(), expected))
test_that('integer', expect_equal(Cfun(), expected))
test_that('integer', expect_identical(class(Rfun()[1,1]), 'integer'))
test_that('integer', expect_identical(class(Cfun()[1,1]), 'integer'))




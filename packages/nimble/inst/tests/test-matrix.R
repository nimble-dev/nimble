source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context('Testing Eigen matrix operations')

## To add a new test: simply add a character string in the 'tests' vector below.
## NOTE: if a character string in the 'tests' vector contains ' b ',
## that is, a lower-case 'b' surrounded by two spaces, then this testing framework assumes
## the return value of that test is 1-dimensional.  Otherwise, the return value is assumed
## to be 2-dimensional.  Violating this convention will break this testing framework.
tests <- c(
    'A %*% B', 'A %*% b', 't(b) %*% A',     ## note: this line all 2-dimensional
    'chol(A)',
    'forwardsolve(A, b )', 'forwardsolve(A, B)', 'forwardsolve(L, b )', 'forwardsolve(L, B)',
       'backsolve(A, b )',    'backsolve(A, B)',    'backsolve(R, b )',    'backsolve(R, B)',
           'solve(A, b )',        'solve(A, B)',        'solve(R, b )',        'solve(R, B)'
)

testNames <- as.character(sapply(tests, nimble:::Rname2CppName))

## (NOTE) that assumption is used on this next line:
testDims <- sapply(grepl(' b ', tests), function(x) if(x) 1 else 2)

methods <- mapply(function(test, dim) eval(substitute(function() { ans <- TEST; returnType(double(DIM)); return(ans) }, list(TEST = parse(text=test)[[1]], DIM = dim))), tests, testDims)
names(methods) <- testNames

n <- 4
set.seed(0)
tmp <- array(rnorm(n^2), c(n,n)) + n*diag(n)
A <- tmp + t(tmp)
R <- chol(A)
L <- t(R)
B <- array(as.numeric(1:n^2), c(n,n))
b <- as.numeric(1:n)

nfDef <- nimbleFunction(
    setup = function(A, R, L, B, b) {},
    run = function() {},
    methods = methods
)

Rnf <- nfDef(A, R, L, B, b)

Cnf <- compileNimble(Rnf)

testOneCase <- function(test, testName, Rnf, Cnf) {
    Rans <- eval(parse(text=test)[[1]])
    Rnfans <- eval(substitute(Rnf$TEST(), list(TEST=as.name(testName))))
    Cnfans <- eval(substitute(Cnf$TEST(), list(TEST=as.name(testName))))
    dif <- max(abs(Rans - Rnfans))
    if(dif > 1E-15) return(1)
    dif <- max(abs(Rans - Cnfans))
    if(dif > 1E-15) return(1)
    return(0)
}

for(i in seq_along(tests)) {
    test <- tests[i]
    testName <- testNames[i]
    testResult <- testOneCase(test, testName, Rnf, Cnf)
    ##if(testResult != 0) message('failed test: ', test)
    test_that(test, expect_equal(testResult, 0))
}









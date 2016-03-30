library(nimble)
source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context('Testing general numeric expressions in the DSL')

## To add a test, add a character string in the 'tests' vector below.
## The first character of each string will be stripped, and gives the dimension of the result.
## If there's an equal sign in the string, then the RHS gives the expected numeric result.

tests <- c(
    ## matrix-multiplication
    '2A %*% B', '2A %*% b', '2t(b) %*% A',
    ## chol(), including non-symmetric argument cases
    '2chol(A) = R',        '2chol(A + A)',        '2chol(A + L + t(L))', '2chol(3*A + L + t(L))',
    '2chol(A + L)',        '2chol(5*A - 2*L)',   ## non-symmetric
    ## forwardsolve(), backsolve(), solve()
    '1forwardsolve(A, b)', '2forwardsolve(A, B)', '1forwardsolve(L, b)', '2forwardsolve(L, B)',
       '1backsolve(A, b)',    '2backsolve(A, B)',    '1backsolve(R, b)',    '2backsolve(R, B)',
           '1solve(A, b)',        '2solve(A, B)',        '1solve(R, b)',        '2solve(R, B)',
    ## grab-bag
    '2forwardsolve(A, B)*2 - chol(A)',  '2inverse(A) + chol(inverse(A))'
)

testDims <- as.numeric(substring(tests, 1, 1))
testsWithoutDim <- substring(tests, 2)
testsSplitOnEqual <- strsplit(testsWithoutDim, '=')
testExprs <- sapply(testsSplitOnEqual, function(x) x[1])
testNames <- as.character(sapply(testExprs, nimble:::Rname2CppName))
ansExprs <- sapply(testsSplitOnEqual, function(x) if(length(x)>1) x[2] else NA)
methods <- mapply(function(test, dim) eval(substitute(function() { ans <- TEST; returnType(double(DIM)); return(ans) }, list(TEST = parse(text=test)[[1]], DIM = dim))), testExprs, testDims)
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

testOneCase <- function(test, testName, ans, Rnf, Cnf) {
    Rans <- eval(parse(text=test)[[1]])
    Rnfans <- eval(substitute(Rnf$TEST(), list(TEST=as.name(testName))))
    Cnfans <- eval(substitute(Cnf$TEST(), list(TEST=as.name(testName))))
    dif <- max(abs(Rans - Rnfans))
    if(dif > 1E-15) return(1)
    dif <- max(abs(Rans - Cnfans))
    if(dif > 1E-15) return(1)
    if(!is.na(ans)) {
        theans <- eval(parse(text=ans)[[1]])
        dif <- max(abs(Rans - theans))
        if(dif > 1E-15) return(1)
    }
    return(0)
}

for(i in seq_along(tests)) {
    test <- testExprs[i]
    testName <- testNames[i]
    ans <- ansExprs[i]
    testResult <- testOneCase(test, testName, ans, Rnf, Cnf)
    ##if(testResult != 0) message('failed test: ', test)
    test_that(test, expect_equal(testResult, 0))
}









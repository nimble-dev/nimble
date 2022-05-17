source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)

nimbleOptions(useADsolveAtomic = TRUE)
nimbleOptions(useADmatMultAtomic = TRUE) # used in meta-taping

## This returns
# y = forwardSolve(A, B)^2
# where A and B each have constant sections and
# variable sections.
# A is n1xn1. B is n1xn2.
# The variable sections are replaced with exp(-d * Ainput)
#    and exp(-2 * d * Binput)
#
# Derivatives are wrt d

argTypes <- list(d = "double()", Ainput = "double(2)", Binput = "double(2)")
op <- list(
  expr = quote({
    A <- matrix(nrow = n1, ncol = n1)
    B <- matrix(nrow = n1, ncol = n2)
    i <- j <- k <-  1L
    for(j in 1:n1) {
      for(i in 1:n1) A[i, j] <- Aconst[i,j]
      for(k in 1:n2) B[j, k] <- Bconst[j,k]
    }
    if(AupperLeft[1] != -1)
      A[ AupperLeft[1]:AlowerRight[1], AupperLeft[2]:AlowerRight[2] ] <- exp(-d * Ainput)
    if(BupperLeft[1] != -1)
      B[ BupperLeft[1]:BlowerRight[1], BupperLeft[2]:BlowerRight[2] ] <- exp(-2 * d * Binput)
    Y <- forwardsolve(A, B)
    out <- sum( Y * Y )
}),
  args = list(d = quote(double()),
              Ainput = quote(double(2)),
              Binput = quote(double(2))),
  outputType = quote(double())
)
forwardsolveTest_pieces <- make_AD_test2(op = op, wrt_args = "d", argTypes = argTypes, includeModelArgs = FALSE)

forwardsolveTest <- nimbleFunction(
  setup = function(Aconst, Bconst, AupperLeft, AlowerRight, BupperLeft, BlowerRight) { ## boundaries of non-constant region
    n1 <- nrow(Aconst)
    if(ncol(Aconst) != n1) stop("Aconst must be square")
    if(nrow(Bconst) != n1) stop("Bconst has wrong number of rows relative to Aconst")
    n2 <- ncol(Bconst)
  },
  run = forwardsolveTest_pieces$run,
  methods = forwardsolveTest_pieces$methods,
  buildDerivs = forwardsolveTest_pieces$buildDerivs
)


checkCase <- function(nf,
                     Aconst, Bconst, A_UL, A_LR,  B_UL, B_LR,
                     order = 0:2,
                     recordArgs, testArgs) {

    Rfxn <- nf(Aconst, Bconst, A_UL, A_LR,  B_UL, B_LR)
    Cfxn <- compileNimble(Rfxn)
    
    test_AD2_oneCall(Rfxn, Cfxn,
                     recordArgs = recordArgs, testArgs = testArgs,
                     order = order, wrt = 1)    
}

n1 <- 5
n2 <- 3

set.seed(1)
Aconst <- matrix(rnorm(n1*n1), nrow = n1)
Bconst <- matrix(rnorm(n1*n2), nrow = n1)
Bconst0 <- Bconst
Bconst0[1,] <-0


makeArgs = function(n1Ar, n1Ac, n1B, n2, d, Adiag = FALSE) {
  # These are replacement sections of A and B,
  # so n1 might differ for A and B, hence n1A and n1B
  # and rows and cols might differ in A even though actual A is square, hence n1Ar and n1Ac
  if(Adiag) Ain <- diag(runif(n1Ar))
  else Ain <- matrix(runif(n1Ar*n1Ac, min = 1, max = 3), nrow = n1Ar, ncol = n1Ac)
  list(
    Ain = Ain,
    Bin = matrix(runif(n1B*n2, min = 1, max = 3), nrow = n1B, ncol = n2),
    d = d
  )
}

## Case with all elements variable.  This works
recordArgs <- makeArgs(n1, n1, n1, n2, 1.2)
testArgs <- makeArgs(n1, n1, n1, n2, 1.4)

checkCase(forwardsolveTest, Aconst, Bconst, c(1, 1), c(n1, n1), c(1, 1), c(n1, n2),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with all of A constant.
recordArgs <- makeArgs(0, 0, n1, n2, 1.2)
testArgs <- makeArgs(0, 0, n1, n2, 1.4)

checkCase(forwardsolveTest, Aconst, Bconst, c(-1, 1), c(n1, n1), c(1, 1), c(n1, n2),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with all of A and B variable and A diagonal
recordArgs <- makeArgs(n1, n1, n1, n2, 1.2, Adiag=TRUE)
testArgs <- makeArgs(n1, n1, n1, n2, 1.4, Adiag=TRUE)

checkCase(forwardsolveTest, Aconst, Bconst, c(1, 1), c(n1, n1), c(1, 1), c(n1, n2),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with all of B constant.
recordArgs <- makeArgs(n1, n1, 0, 0, 1.2)
testArgs <- makeArgs(n1, n1, 0, 0, 1.4)

checkCase(forwardsolveTest, Aconst, Bconst, c(1, 1), c(n1, n1), c(-1, 1), c(n1, n2),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with all of B constant with some leading zeros.
recordArgs <- makeArgs(n1, n1, 0, 0, 1.2)
testArgs <- makeArgs(n1, n1, 0, 0, 1.4)

checkCase(forwardsolveTest, Aconst, Bconst0, c(1, 1), c(n1, n1), c(-1, 1), c(n1, n2),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with all of B variable except for some leading zeros
recordArgs <- makeArgs(n1, n1, n1-1, n2, 1.2)
testArgs <- makeArgs(n1, n1, n1-1, n2, 1.4)

checkCase(forwardsolveTest, Aconst, Bconst0, c(1, 1), c(n1, n1), c(2, 1), c(n1, n2),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with first rows of both A and B of constant elements. (I have a notes that this does not get special handling internally, but it might be out of date.).  My original testing code for this used Bconst0, but I'm not sure that was intended.
recordArgs <- makeArgs(n1-1, n1, n1-1, n2, 1.2)
testArgs <- makeArgs(n1-1, n1, n1-1, n2, 1.4)

checkCase(forwardsolveTest, Aconst, Bconst, c(2, 1), c(n1, n1), c(2, 1), c(n1, n2),
          recordArgs = recordArgs, testArgs = testArgs)

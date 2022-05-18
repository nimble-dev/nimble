source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)

nimbleOptions(useADmatMultAtomic = TRUE)

## This returns
# y = (A %*% B)^2
# where A and B each have constant sections and
# variable sections.
# The variable sections are replaced with exp(-d * Ainput)
#    and exp(-d * Binput)
#
# Derivatives are wrt d

argTypes <- list(d = "double()", Ainput = "double(2)", Binput = "double(2)")
op <- list(
  expr = quote({
    A <- matrix(nrow = n1, ncol = n2)
    B <- matrix(nrow = n2, ncol = n3)
    i <- j <- k <-  1L
    for(j in 1:n2) {
      for(i in 1:n1) A[i, j] <- Aconst[i,j]
      for(k in 1:n3) B[j, k] <- Bconst[j,k]
    }
    if(AupperLeft[1] != -1)
      A[ AupperLeft[1]:AlowerRight[1], AupperLeft[2]:AlowerRight[2] ] <- exp(-d * Ainput)
    if(BupperLeft[1] != -1)
      B[ BupperLeft[1]:BlowerRight[1], BupperLeft[2]:BlowerRight[2] ] <- exp(-d * Binput)
    Y <- A %*% B
    out <- sum( Y * Y )
  }),
  args = list(d = quote(double()),
              Ainput = quote(double(2)),
              Binput = quote(double(2))),
  outputType = quote(double())
)
matMultTestFxn_pieces <- make_AD_test2(op = op, wrt_args = "d", argTypes = argTypes, includeModelArgs = FALSE)

matMultTestFxn <- nimbleFunction(
  setup = function(Aconst, Bconst, AupperLeft, AlowerRight, BupperLeft, BlowerRight) { ## boundaries of non-constant region
    n1 <- nrow(Aconst)
    n2 <- ncol(Aconst)
    n3 <- ncol(Bconst)
  },
  run = matMultTestFxn_pieces$run,
  methods = matMultTestFxn_pieces$methods,
  buildDerivs = matMultTestFxn_pieces$buildDerivs
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
n2 <- 7
n3 <- 9

set.seed(1)
Aconst <- matrix(runif(n1*n2, min = 1, max = 3), nrow = n1)
Bconst <- matrix(runif(n2*n3, min = 1, max = 3), nrow = n2)
Aconst0 <- Aconst
Aconst0[1,] <- 0  # top row = 0
Aconst0[n1,] <- 0 # bottom row = 0

makeArgs = function(n1, n2A, n2B, n3, d) {
  # These are replacement sections of A and B,
  # so n2 might differ for A and B, hence n2A and n2B
  list(
    Ain = matrix(runif(n1*n2A, min = 1, max = 3), nrow = n1, ncol = n2A),
    Bin = matrix(runif(n2B*n3, min = 1, max = 3), nrow = n2B, ncol = n3),
    d = d
  )
}

## Warm-up case (test of the testing code) with some of A and B each constant and variable
recordArgs <- makeArgs(2, 2, 2, 2, 1.2)
testArgs <- makeArgs(2, 2, 2, 2, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(2, 2), c(3, 3), c(2, 2), c(3, 3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with all elements variable.  This works
set.seed(2)
recordArgs <- makeArgs(n1, n2, n2, n3, 1.2)
testArgs <- makeArgs(n1, n2, n2, n3, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(1, 1), c(n1, n2), c(1, 1), c(n2, n3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with all of A constant. 
set.seed(2)
recordArgs <- makeArgs(0, n2, n2, n3, 1.2)
testArgs <- makeArgs(0, n2, n2, n3, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(-1, 1), c(n1, n2), c(1, 1), c(n2, n3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with all of A constant, with some zero rows.
set.seed(2)
recordArgs <- makeArgs(0, n2, n2, n3, 1.2)
testArgs <- makeArgs(0, n2, n2, n3, 1.4)

checkCase(matMultTestFxn, Aconst0, Bconst, c(-1, 1), c(n1, n2), c(1, 1), c(n2, n3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with all of A constant and part of B constant.
set.seed(2)
recordArgs <- makeArgs(0, n2-2, n2-2, n3-2, 1.2)
testArgs <- makeArgs(0, n2-2, n2-2, n3-2, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(-1, 1), c(n1, n2), c(2, 2), c(n2-1, n3-1),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with all of B constant. This works
set.seed(2)
recordArgs <- makeArgs(n1, n2, n2, 0, 1.2)
testArgs <- makeArgs(n1, n2, n2, 0, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(1, 1), c(n1, n2), c(-1, 1), c(n2, n3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with all of B constant and part of A constant. 
set.seed(2)
recordArgs <- makeArgs(n1-2, n2-2, n2-2, 0, 1.2)
testArgs <- makeArgs(n1-2, n2-2, n2-2, 0, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(2, 2), c(n1-1, n2-1), c(-1, 1), c(n2, n3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with first row in A of constant elements.
set.seed(2)
recordArgs <- makeArgs(n1-1, n2, n2, n3, 1.2)
testArgs <- makeArgs(n1-1, n2, n2, n3, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(2, 1), c(n1, n2), c(1, 1), c(n2, n3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with first col in A of constant elements.
set.seed(2)
recordArgs <- makeArgs(n1, n2-1, n2, n3, 1.2)
testArgs <- makeArgs(n1, n2-1, n2, n3, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(1, 2), c(n1, n2), c(1, 1), c(n2, n3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with last row A of constant elements.
set.seed(2)
recordArgs <- makeArgs(n1-1, n2, n2, n3, 1.2)
testArgs <- makeArgs(n1-1, n2, n2, n3, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(1, 1), c(n1-1, n2), c(1, 1), c(n2, n3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with last col in A of constant elements.
set.seed(2)
recordArgs <- makeArgs(n1, n2-1, n2, n3, 1.2)
testArgs <- makeArgs(n1, n2-1, n2, n3, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(1, 1), c(n1, n2-1), c(1, 1), c(n2, n3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with both first and last row A of constant elements.
set.seed(2)
recordArgs <- makeArgs(n1-2, n2, n2, n3, 1.2)
testArgs <- makeArgs(n1-2, n2, n2, n3, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(2, 1), c(n1-1, n2), c(1, 1), c(n2, n3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with both first and last col of A of constant elements.
set.seed(2)
recordArgs <- makeArgs(n1, n2-2, n2, n3, 1.2)
testArgs <- makeArgs(n1, n2-2, n2, n3, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(1, 2), c(n1, n2-1), c(1, 1), c(n2, n3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with both first and last 2 COLs of A and first 2 and last ROWs of A all constant elements.
set.seed(2)
recordArgs <- makeArgs(n1-3, n2-3, n2, n3, 1.2)
testArgs <- makeArgs(n1-3, n2-3, n2, n3, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(3, 2), c(n1-1, n2-2), c(1, 1), c(n2, n3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with first col of B with constant elements.
set.seed(2)
recordArgs <- makeArgs(n1, n2, n2, n3-1, 1.2)
testArgs <- makeArgs(n1, n2, n2, n3-1, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(1, 1), c(n1, n2), c(1, 2), c(n2, n3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with last col of B with constant elements.
set.seed(2)
recordArgs <- makeArgs(n1, n2, n2, n3-1, 1.2)
testArgs <- makeArgs(n1, n2, n2, n3-1, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(1, 1), c(n1, n2), c(1, 1), c(n2, n3-1),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with first and last col of B with constant elements.
set.seed(2)
recordArgs <- makeArgs(n1, n2, n2, n3-2, 1.2)
testArgs <- makeArgs(n1, n2, n2, n3-2, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(1, 1), c(n1, n2), c(1, 2), c(n2, n3-1),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with first row of B with constant elements.
set.seed(2)
recordArgs <- makeArgs(n1, n2, n2-1, n3, 1.2)
testArgs <- makeArgs(n1, n2, n2-1, n3, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(1, 1), c(n1, n2), c(2, 1), c(n2, n3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with last row of B with constant elements.
set.seed(2)
recordArgs <- makeArgs(n1, n2, n2-1, n3, 1.2)
testArgs <- makeArgs(n1, n2, n2-1, n3, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(1, 1), c(n1, n2), c(1, 1), c(n2-1, n3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with first and last row of B with constant elements.
set.seed(2)
recordArgs <- makeArgs(n1, n2, n2-2, n3, 1.2)
testArgs <- makeArgs(n1, n2, n2-2, n3, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(1, 1), c(n1, n2), c(2, 1), c(n2-1, n3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with some rows and columns of B with constant elements.
set.seed(2)
recordArgs <- makeArgs(n1, n2, n2-3, n3-3, 1.2)
testArgs <- makeArgs(n1, n2, n2-3, n3-3, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(1, 1), c(n1, n2), c(3, 2), c(n2-1, n3-2),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with first row of A and first col of B with constant elements.
set.seed(2)
recordArgs <- makeArgs(n1-1, n2, n2, n3-1, 1.2)
testArgs <- makeArgs(n1-1, n2, n2, n3-1, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(2, 1), c(n1, n2), c(1, 2), c(n2, n3),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with last row of A and last col of B with constant elements. This works.
set.seed(2)
recordArgs <- makeArgs(n1-1, n2, n2, n3-1, 1.2)
testArgs <- makeArgs(n1-1, n2, n2, n3-1, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(1, 1), c(n1-1, n2), c(1, 1), c(n2, n3-1),
          recordArgs = recordArgs, testArgs = testArgs)

## Case with first row of A and last col of B with constant elements.
set.seed(2)
recordArgs <- makeArgs(n1-1, n2, n2, n3-1, 1.2)
testArgs <- makeArgs(n1-1, n2, n2, n3-1, 1.4)

checkCase(matMultTestFxn, Aconst, Bconst, c(2, 1), c(n1, n2), c(1, 1), c(n2, n3-1),
          recordArgs = recordArgs, testArgs = testArgs)

## STOPPED HERE: There are some random-configuration tests that couuld be ported over here too.

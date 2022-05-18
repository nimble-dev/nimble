source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)

nimbleOptions(useADmatMultAtomic = TRUE)
nimbleOptions(useADmatInverseAtomic = TRUE)

## This returns
# y = || exp(-d * A))^{-1} ||
# where A is n-x-n
argTypes <- list(d = "double()", Ainput = "double(2)")
op <- list(
  expr = quote({
    A <- matrix(nrow = n, ncol = n)
    i <- j <- 1L
    for(j in 1:n) {
      for(i in 1:n) A[i, j] <- Aconst[i,j]
    }
    if(AupperLeft[1] != -1)
      A[ AupperLeft[1]:AlowerRight[1], AupperLeft[2]:AlowerRight[2] ] <- exp(-d * Ainput)
    Y <- inverse(exp(-d * A))
    out <- sum( Y * Y )
  }),
  args = list(d = quote(double()),
              Ainput = quote(double(2))),
  outputType = quote(double())
)
matInverse_pieces <- make_AD_test2(op = op, wrt_args = "d", argTypes = argTypes, includeModelArgs = FALSE)

matInverse <- nimbleFunction(
  setup = function(Aconst, AupperLeft, AlowerRight) { ## boundaries of non-constant region
    n <- nrow(Aconst)
    if(ncol(Aconst) != n) stop("Aconst must be square")
  },
  run = matInverse_pieces$run,
  methods = matInverse_pieces$methods,
  buildDerivs = matInverse_pieces$buildDerivs
)

checkCase <- function(nf,
                     Aconst, A_UL, A_LR,
                     order = 0:2,
                     recordArgs, testArgs) {

    Rfxn <- nf(Aconst, A_UL, A_LR)
    Cfxn <- compileNimble(Rfxn)
    
    test_AD2_oneCall(Rfxn, Cfxn,
                     recordArgs = recordArgs, testArgs = testArgs,
                     order = order, wrt = 1,
                     RCrelTol = c(1e-12, 1e-05, 0.001))    
}

n <- 7

makeArgs = function(n1Ar, n1Ac, d, Adiag = FALSE) {
  # These are replacement sections of A and B,
  # so n1 might differ for A and B, hence n1A and n1B
  # and rows and cols might differ in A even though actual A is square, hence n1Ar and n1Ac
  if(Adiag) Ain <- diag(runif(n1Ar))
  else Ain <- matrix(runif(n1Ar*n1Ac, min = 1, max = 3), nrow = n1Ar, ncol = n1Ac)
  list(
    Ain = Ain,
    d = d
  )
}

set.seed(4)
Aconst <- matrix(runif(n*n, min = 1, max = 3), nrow = n)

## Case with all elements variable.
set.seed(3)
recordArgs <- makeArgs(n, n, 1.2)
testArgs <- makeArgs(n, n, 1.4)

checkCase(matInverse, Aconst, c(1, 1), c(n, n),
          recordArgs = recordArgs, testArgs = testArgs)

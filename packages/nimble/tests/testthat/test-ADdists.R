source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'testthat', 'AD_distribution_test_lists.R'), package = 'nimble'))
source(system.file(file.path('tests', 'testthat', 'AD_knownFailures.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)

## Some status notes:
## dcat: broken in double taping
## For many distributions, edge values of parameters (e.g. 0) may not
## work in derivatives
##
## Things to do or check:
## Add warning message if prec_param or scale_param is not CppAD::Constant
## Add test for invwishart
## Add test for lkj_chol
## Look for other gaps


test_that('Derivatives of dnorm function correctly.',
  {
    ADfun1 <- nimbleFunction(
      setup = function(){},
      run = function(y = double(1)) {
        outList <- derivs(testMethod(y), wrt = c('x'))
        returnType(ADNimbleList())
        return(outList)
      },
      methods = list(
        testMethod = function(x = double(1, 2)) {
          out <- dnorm(x[1],0,1)
          returnType(double())
          return(out)
        }
      ), buildDerivs = c('testMethod')
    )
    ADfunInst <- ADfun1()
    xRec <- matrix(c(1, -1))
    x <- matrix(c(2, -2))
    Rderivs <- ADfunInst$run(x)
    temporarilyAssignInGlobalEnv(ADfunInst)
    cADfunInst <- compileNimble(ADfunInst)
    cADfunInst$run(xRec)
    cderivs <- cADfunInst$run(x)
    expect_equal(cderivs$value, Rderivs$value)
    expect_equal(cderivs$jacobian, Rderivs$jacobian)
    expect_equal(cderivs$hessian, Rderivs$hessian)
  }
)


## Dirichlet with log argument
dirch_test_log <- make_AD_test2(
  op = list(
    name = "ddirch manual using additive log-ratio transformation",
    opParam = list(name = "ddirch manual"),
    # X = Xtrans = log(Xorig_i / Xorig_n), i = 1... n-1
    # Xorig[1:n-1] = exp(Xtrans)
    # Xorig[n] = 1/(1+sum(Xorig[1:n-1]))
    expr = quote({
      Xorig_1_nm1_over_Xn <- exp(x)
      Xorig_n <- 1/(1 + sum(Xorig_1_nm1_over_Xn))
      Xorig <- c(Xorig_1_nm1_over_Xn * Xorig_n, Xorig_n)
      out <- ddirch(x = Xorig, alpha=alpha, log = log)
    }),
    args = list(
      x = quote(double(1)),
      alpha = quote(double(1)),
      log = quote(double())
    ),
    outputType = quote(double())
  ),
  argTypes = c(x='double(1)', alpha='double(1)', log='double()'),
  wrt = c('x', 'alpha'),
  inputs = list(record = list(x = c(log(.2/.5), log(.3/.5)), alpha = c(2, 4, 4), log = 1),
                test   = list(x = c(log(.4/.45), log(.15/.45)), alpha = c(3, 2, 5), log = 0)))

dirch_test_out <- test_AD2(dirch_test_log)

## Dirichlet without log argument
dirch_test_fixedlog <- make_AD_test2(
  op = list(
    name = "ddirch manual using additive log-ratio transformation",
    opParam = list(name = "ddirch manual"),
    # X = Xtrans = log(Xorig_i / Xorig_n), i = 1... n-1
    # Xorig[1:n-1] = exp(Xtrans)
    # Xorig[n] = 1/(1+sum(Xorig[1:n-1]))
    expr = quote({
      Xorig_1_nm1_over_Xn <- exp(x)
      Xorig_n <- 1/(1 + sum(Xorig_1_nm1_over_Xn))
      Xorig <- c(Xorig_1_nm1_over_Xn * Xorig_n, Xorig_n)
      out <- ddirch(x = Xorig, alpha=alpha, log = 1)
    }),
    args = list(
      x = quote(double(1)),
      alpha = quote(double(1))
    ),
    outputType = quote(double())
  ),
  argTypes = c(x='double(1)', alpha='double(1)'),
  wrt = c('x', 'alpha'),
  inputs = list(record = list(x = c(log(.2/.5), log(.3/.5)), alpha = c(2, 4, 4)),
                test   = list(x = c(log(.4/.45), log(.15/.45)), alpha = c(3, 2, 5))))

dirch_test_out <- test_AD2(dirch_test_fixedlog)

## Wishart
## Testing Wishart is tricky because both x and cholesky are
## matrices and also because not all elements are indendent
makeARcov <- function(n, rho, sigma) {
  s2 <- sigma^2
  ans <- matrix(nrow = n, ncol = n)
  for(i in 1:n) {
    for(j in 1:i) {
      ans[i,j] <- ans[j,i] <- s2 * rho^(abs(i-j))
    }
  }
  ans
}
set.seed(123)
cholRec <- chol(makeARcov(4, .6, 2))
cholTest <- chol(makeARcov(4, .55, 3))

wRec <- rwish_chol(1, cholRec, df = 7, scale_param = FALSE)
wTest <- rwish_chol(1, cholTest, df = 7, scale_param = FALSE)

cholRecTri <- cholRec[ upper.tri(cholRec, TRUE)]
cholTestTri <-cholTest[ upper.tri(cholTest, TRUE)]

wRecTri <- wRec[ upper.tri(wRec, TRUE) ]
wTestTri <- wTest[ upper.tri(wTest, TRUE) ]

wish_test_log <- make_AD_test2(
  op = list(
    name = "dwish_chol manual",
    opParam = list(name = "dwish_chol manual"),
    expr = quote({
      # populate 2D matrices from the vectors
      # created from upper triangular values.
      x2D <- nimMatrix(nrow = 4, ncol = 4, init=FALSE)
      chol2D <- nimMatrix(nrow = 4, ncol = 4, init = FALSE)
      i <- 1L
      j <- 1L
      indOrig <- 1L
      for(j in 1:4) {
        for(i in 1:j) {
          x2D[i,j] <- x2D[j,i] <- x[indOrig]
          chol2D[i,j] <- chol[indOrig]  # chol2D does not need lower triangular entries
          indOrig <- indOrig + 1
        }
      }
      out <- dwish_chol(x = x2D, cholesky=chol2D, df = df, log = log)
    }),
    args = list(
      x = quote(double(1)),
      chol = quote(double(1)),
      df = quote(double()),
      log = quote(double())
    ),
    outputType = quote(double())
  ),
  argTypes = c(x='double(1)', chol='double(1)', df = 'double()', log='double()'),
  wrt = c('x', 'chol'),
  inputs = list(record = list(x = wRecTri, chol = cholRecTri, df = 7, log = 0),
                test   = list(x = wTestTri, chol = cholTestTri, df = 8, log = 1)),
)

wish_test_out <- test_AD2(wish_test_log)

## LKJ

p <- 5

lkj_test_log <- make_AD_test2(
  op = list(
    name = "dlkj_corr_cholesky manual",
    opParam = list(name = "dlkj_corr_cholesky manual"),
    expr = quote({
      # correlation matrix inverse transform 
      z <- nimMatrix(nrow = 5, ncol = 5, init=FALSE)
      u <- nimMatrix(nrow = 5, ncol = 5, init=FALSE)
        
      j <- 1L
      i <- 1L
      cnt <- 1L
      for(j in 2:5)
          for(i in 1:(j-1)) {
              z[i,j] <- tanh(x[cnt])
              cnt <- cnt + 1
          }
      u[1,2:5] <- z[1,2:5]
      u[2,3:5] <- z[2,3:5]*sqrt(1-u[1,3:5]^2)
      u[3,4:5] <- z[3,4:5]*sqrt(1-u[1,4:5]^2-u[2,4:5]^2)
      u[4,5] <- z[4,5]*sqrt(1-u[1,5]^2-u[2,5]^2-u[3,5]^2)
      for(j in 1:5)
          u[j,j] <- sqrt(1-sum(u[1:5,j]^2))
      out <- dlkj_corr_cholesky(x = u, eta = eta, p = 5, log = log)
    }),
    args = list(
      x = quote(double(1)),
      eta = quote(double()),
      p = quote(double()),
      log = quote(double())
    ),
    outputType = quote(double())
  ),
  argTypes = c(x='double(1)', eta='double()', p = 'double()', log='double()'),
  wrt = c('x', 'eta'),
  inputs = list(record = list(x = rnorm(choose(p,2)), eta = etaRec, p = 5, log = 0),
                test   = list(x = rnorm(choose(p,2)), eta = etaTest, p = 5, log = 1)),
)

lkj_test_out <- test_AD2(lkj_test_log)


## from AD_distribution_test_lists.R
resetTols()
testResults <- lapply(distn_tests2[1:140], test_AD2)
testResults <- lapply(distn_with_log_tests2[1:140], test_AD2)

ADtestEnv$RCrelTol[4] <- 1e-5 # set looser absolute tolerance for the dmnorm and dmvt tests
testResults <- lapply(distn_tests2[141:142], test_AD2)
testResults <- lapply(distn_with_log_tests2[141:142], test_AD2)
resetTols()

# test_AD_batch does not work with these.  Check that out
# make lapply run line for distn_tests2
# test_AD_batch(distn_with_log_tests2, knownFailures = AD_knownFailures, verbose = TRUE)

nimbleOptions(enableDerivs= EDopt)
nimbleOptions(buildModelDerivs = BMDopt)

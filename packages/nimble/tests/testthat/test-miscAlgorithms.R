source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

# Tests of Misc Algorithms which currently only includes expmAv
test_that("expmAv works", {
  
  ## Note that pracma was not accurate enough even for a small example
  ## expm::expm(A) %*% v
  A <- matrix(c(-0.120, 100, 0.9, -150), 2, 2)
  v <- c(120, 0)
  eAv <- c(192.793313160, 128.120518022) 
  
  ## Wrapper as derivs = TRUE with no setup.
  expAv_wrap <- nimbleFunction(
    run = function(A = double(2), v = double(1),
                 tol = double(0, default = 1e-8), 
                 rescaleFreq = double(0, default = 10),
                 Nmax = integer(0, default = 10000),
                 sparse = logical(0, default = FALSE)){
      ans <- expmAv(A, v, tol, rescaleFreq, Nmax, sparse)
      returnType(double(1))
      return(ans)
    }
  )
  expAvc <- compileNimble(expAv_wrap)

  eAvNimR <- expmAv(A, v)
  eAvNim <- expAvc(A, v)
  expect_equal(eAv, eAvNim, tol = 1e-6)
  expect_equal(eAv, eAvNimR, tol = 1e-6)

  ## Increase tolerance.
  eAvNim2 <- expAvc(A, v, tol = 1e-16)
  expect_equal(eAv, eAvNim2, tol = 1e-10)

  ## Compare sparse version with regular version:
  set.seed(1) 
  A <- Matrix::rsparsematrix(10, 10, .5)
  A <- A - diag(diag(A))
  v <- rnorm(10)
  d <- -rowSums(A) - 8
  diag(A) <- d
  eAvdense <- expAvc(as.matrix(A), v, sparse = FALSE)
  eAvsparse <- expAvc(as.matrix(A), v, sparse = TRUE)
  expect_equal(eAvdense, eAvsparse, tol = 1e-16)

  ## Increase tolerance:
  eAvsparse2 <- expAvc(as.matrix(A), v, sparse = TRUE, tol = 1e-16)
  expect_equal(eAvsparse2, eAvsparse, tol = 1e-8) ## Comparing 1e-8 and 1e-16
  eAvsparse3 <- expAvc(as.matrix(A), v, sparse = TRUE, tol = 1e-5)
  expect_equal(eAvsparse3, eAvsparse2, tol = 1e-5) ## Comparing 1e-8 and 1e-16

  ## Check rescaling:
  eAvsparse_scale1 <- expAvc(as.matrix(A), v, sparse = TRUE, tol = 1e-8, rescaleFreq = 1)
  eAvsparse_scale10 <- expAvc(as.matrix(A), v, sparse = TRUE, tol = 1e-8, rescaleFreq = 10)
  eAvsparse_scale50 <- expAvc(as.matrix(A), v, sparse = TRUE, tol = 1e-8, rescaleFreq = 50)
  expect_equal(eAvsparse_scale1, eAvsparse, tol = 1e-16) ##  Approximation is best with frequent rescaling.
  expect_equal(eAvsparse_scale10, eAvsparse, tol = 1e-12)
  expect_equal(eAvsparse_scale50, eAvsparse, tol = 1e-10)
})

# Test full matrix exponential using scaling and squaring.
test_that("Matrix Exponential Works", {
  
  ## expm::expm(A) %*% v
  A <- matrix(c(-0.120, 100, 0.9, -150), 2, 2)
  v <- c(120, 0)
  eAv <- c(192.793313160, 128.120518022) 
  
  ## Wrapper as derivs = TRUE with no setup.
  expmAc <- compileNimble(expmA)

  eAR <- expmA(A)
  eAC <- expmAc(A)
  eAvR <- (eAR %*% v)[,1]
  eAvC <- (eAC %*% v)[,1]
  expect_equal(eAv, eAvR, tol = 1e-6)
  expect_equal(eAv, eAvC, tol = 1e-6)

  ## Increase tolerance.
  eAvNim2 <- (expmAc(A, tol = 1e-16)  %*% v)[,1]
  expect_equal(eAv, eAvNim2, tol = 1e-10)

  ## Bigger Matrix:
  set.seed(101)
  n <- 10
  qd <- rgamma(n, 1, 1)
  A <- matrix(0, n, n)
  A[1,2] <- qd[1]
  for( i in 2:(n-1) ){  
    A[i, i+1] <- qd[i]/2
    A[i, i-1] <- qd[i]/2
  }
  A[n, n-1] <- qd[n]
  diag(A) <- -qd
  v <- rgamma(n, 5, 1)
  eAv <- expmAv(A, v, tol = 1e-16)
  eAv2 <- (expmAc(A, tol = 1e-16)  %*% v)[,1]
  expect_equal(eAv, eAv2, tol = 1e-14)
})


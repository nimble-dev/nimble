source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

# Tests of Misc Algorithms which currently includes expAv and expm
test_that("expAv works.", {
    
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
      ans <- expAv(A, v, tol, rescaleFreq, Nmax, sparse)
      returnType(double(1))
      return(ans)
    }
  )
  expAvc <- compileNimble(expAv_wrap)

  eAvNimR <- expAv(A, v)
  eAvNim <- expAvc(A, v)
  expect_equal(eAv, eAvNim, tol = 1e-8)
  expect_equal(eAv, eAvNimR, tol = 1e-8)

  ## Increase tolerance.
  eAvNim2 <- expAvc(A, v, tol = 1e-16)
  expect_equal(eAv, eAvNim2, tol = 1e-10)

  ## Compare sparse version with regular version:
  set.seed(1) 
  A <- as.matrix(Matrix::rsparsematrix(10, 10, .5))
  A <- A - diag(diag(A))
  v <- rnorm(10)
  d <- -rowSums(A) - 8
  diag(A) <- d
  eAvdense <- expAvc(A, v, sparse = FALSE)
  eAvsparse <- expAvc(A, v, sparse = TRUE)
  expect_equal(eAvdense, eAvsparse, tol = 1e-15)

  ## Increase tolerance:
  eAvsparse2 <- expAvc(as.matrix(A), v, sparse = TRUE, tol = 1e-16)
  expect_equal(eAvsparse2, eAvsparse, tol = 1e-8) ## Comparing 1e-8 and 1e-16
  eAvsparse3 <- expAvc(as.matrix(A), v, sparse = TRUE, tol = 1e-5)
  expect_equal(eAvsparse3, eAvsparse2, tol = 1e-5) ## Comparing 1e-8 and 1e-16

  ## Check rescaling:
  eAvsparse_scale1 <- expAvc(as.matrix(A), v, sparse = TRUE, tol = 1e-8, rescaleFreq = 1)
  eAvsparse_scale10 <- expAvc(as.matrix(A), v, sparse = TRUE, tol = 1e-8, rescaleFreq = 10)
  eAvsparse_scale50 <- expAvc(as.matrix(A), v, sparse = TRUE, tol = 1e-8, rescaleFreq = 50)
  expect_equal(eAvsparse_scale1, eAvsparse, tol = 1e-15) ##  Approximation is best with frequent rescaling.
  expect_equal(eAvsparse_scale10, eAvsparse, tol = 1e-12)
  expect_equal(eAvsparse_scale50, eAvsparse, tol = 1e-10)
  
  ## Test larger matrix that has positive diagonals, does not require negative.
  set.seed(2) 
  A <- as.matrix(Matrix::rsparsematrix(50, 50, 0.25))
  v <- rnorm(50)
  eAvdense <- expAvc(A, v, sparse = FALSE)
  eAvsparse <- expAvc(A, v, sparse = TRUE)
  expect_equal(eAvdense, eAvsparse, tol = 1e-15)
  
  # expAv_rtmb <- RTMB::expAv(A, v)
  ## Indices 1,5, 10 Check sprintf("%.16f",expAv_rtmb[c(1,5,10)])
  expAv_rtmb <- c(-10.4849004612874221, 5.7463833745583406, -4.3508466614873917)
  expect_equal(eAvsparse[c(1,5,10)], expAv_rtmb, tol = 1e-8)
  
  ## Larger numerical test with bigger time period, leads to overflow for both expm package and ours, with same result.
  # expAv_rtmb <- RTMB::expAv(A*10, v)
  # Indices 1,5, 10 Check sprintf("%.16f",expAv_rtmb[c(1,5,10)])
  vals <- c(8591151029674.1933593750000000, -2283297500254.1206054687500000, -16239759412070.7285156250000000)
  eAv <- expAvc(A*10, v, tol = 1e-8)
  expect_equal(eAv[c(1,5,10)], vals, tol = 1e-8)   
})

# Test full matrix exponential using scaling and squaring.
test_that("Matrix Exponential Works", {
  
  ## expm::expm(A) %*% v
  A <- matrix(c(-0.120, 100, 0.9, -150), 2, 2)
  v <- c(120, 0)
  eAv <- c(192.793313160, 128.120518022) 
  
  ## Wrapper as derivs = TRUE with no setup.
  expm_wrap <- nimbleFunction(
    run = function(A = double(2), tol = double(0, default = 1e-8)){
      returnType(double(2))
      ans <- expm(A, tol)
      return(ans)
    }
  )
  expmc <- compileNimble(expm_wrap)

  eAR <- expm(A)
  eAC <- expmc(A)
  eAvR <- (eAR %*% v)[,1]
  eAvC <- (eAC %*% v)[,1]
  expect_equal(eAv, eAvR, tol = 1e-6)
  expect_equal(eAv, eAvC, tol = 1e-6)

  ## Increase tolerance.
  eAvNim2 <- (expmc(A, tol = 1e-16)  %*% v)[,1]
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
  eAv <- expAv(A, v, tol = 1e-16)
  eAv2 <- (expmc(A, tol = 1e-16)  %*% v)[,1]
  expect_equal(eAv, eAv2, tol = 1e-14)
  
  ## Now a much larger matrix and compare with expm::expm...
  set.seed(2) 
  A <- as.matrix(Matrix::rsparsematrix(50, 50, 0.25))

  # test <- expm::expm(A)
  # sprintf("%.16f",test[cbind(c(1,10,20,32),c(1,15,25,32))])
  vals <- c(1.4785161123142350, 2.9154566237879371, -0.9980758066172498, 1.1608012543866268)
  eAv <- expm(A, tol = 1e-14)
  expect_equal(eAv[cbind(c(1,10,20,32),c(1,15,25,32))], vals, tol = 1e-14) 
  
  ## Larger numerical test with bigger time period, leads to overflow for both expm package and ours, with same result.
  # test <- expm::expm(A*10)
  # sprintf("%.16f",test[cbind(c(1,10,20,32),c(1,15,25,32))])
  vals <- c(2247213364689.7807617187500000, 1321534261508.6760253906250000, 81388538896.9135742187500000, -666846013929.7093505859375000)
  eAv <- expm(A*10, tol = 1e-14)
  expect_equal(eAv[cbind(c(1,10,20,32),c(1,15,25,32))], vals, tol = 1e-12) 
})

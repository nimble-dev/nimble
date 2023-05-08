source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)

nimbleOptions(useADcholAtomic = TRUE)
nimbleOptions(useADsolveAtomic  = TRUE)                              
nimbleOptions(useADmatMultAtomic = TRUE)                             
nimbleOptions(useADmatInverseAtomic  = TRUE)                       

relTol <- eval(formals(test_ADModelCalculate)$relTol)
relTol[3] <- 1e-6
relTol[4] <- 1e-4

verbose <- FALSE

context("Testing of derivatives for calculate() for nimbleModel with various mv distributions")

code <- nimbleCode({
    Sigma1[1:n,1:n] <- exp(-dist[1:n,1:n]/rho)

    Q[1:n,1:n] <- inverse(Sigma1[1:n, 1:n])
    y[4, 1:n] ~ dmnorm(mu4[1:n], Q[1:n,1:n])

    Uprec[1:n, 1:n] <- chol(Q[1:n,1:n])
    Ucov[1:n, 1:n] <- chol(Sigma1[1:n,1:n])
    y[5, 1:n] ~ dmnorm(mu5[1:n], cholesky = Uprec[1:n,1:n], prec_param = 1)
    y[6, 1:n] ~ dmnorm(mu6[1:n], cholesky = Ucov[1:n,1:n], prec_param = 0)

    W1[1:n, 1:n] ~ dinvwish(R = R[1:n,1:n], df = nu)

    UR[1:n, 1:n] <- chol(R[1:n,1:n])
    US[1:n, 1:n] <- chol(inverse(R[1:n,1:n]))
    W2[1:n, 1:n] ~ dinvwish(cholesky = UR[1:n,1:n], df = nu, scale_param = 0)
    W3[1:n, 1:n] ~ dinvwish(cholesky = US[1:n,1:n], df = nu, scale_param = 1)

    W4[1:n, 1:5] ~ dwish(R[1:n, 1:n], df = nu)
    W5[1:n, 1:5] ~ dwish(cholesky = UR[1:n, 1:n], df = nu, scale_param = 0)
    W6[1:n, 1:5] ~ dwish(cholesky = US[1:n, 1:n], df = nu, scale_param = 1)
    
    mu4[1:n] ~ dmnorm(z[1:n], W4[1:n,1:n])
    mu5[1:n] ~ dmnorm(z[1:n], W5[1:n,1:n])
    mu6[1:n] ~ dmnorm(z[1:n], W6[1:n,1:n])
    rho ~ dgamma(2, 3)
    nu ~ dunif(0, 100)
})

set.seed(1)
n <- 5
locs <- runif(n)
dd <- fields::rdist(locs)
R <- crossprod(matrix(rnorm(n^2), n, n))
model <- nimbleModel(code, constants = list(n = n),
                     inits = list(dist = dd, R = R, nu = 8, rho = rgamma(1, 1, 1),
                                                                 z = rep(1, n)))
model$simulate()
model$calculate()
model$setData('y')

newDist <- as.matrix(dist(runif(n)))
newR <- crossprod(matrix(rnorm(n*n), n))
newW4 <- crossprod(matrix(rnorm(n*n), n))
newW5 <- crossprod(matrix(rnorm(n*n), n))
newW6 <- crossprod(matrix(rnorm(n*n), n))
relTolTmp <- relTol
relTolTmp[1] <- 1e-14
relTolTmp[2] <- 1e-6
relTolTmp[3] <- 1e-1
relTolTmp[4] <- 1e-1
relTolTmp[5] <- 1e-11



## rOutput2d11 result can be wildly out of tolerance, so not checking it.

test_ADModelCalculate(model, newUpdateNodes = list(nu = 12.1, dist = newDist, R = newR, W4 = newW4, W5 = newW5, W6 = newW6),
                      x = 'prior', absTolThreshold = 1e-12, checkCompiledValuesIdentical = FALSE,
                      useParamTransform = TRUE, useFasterRderivs = TRUE, checkDoubleUncHessian = FALSE,
                      relTol = relTolTmp, verbose = verbose,
                      name = 'various multivariate dists')
## 1310 seconds.

## This segfaults as of 2023-05-08 (and did as well 2023-03-25), with libnimble.a (or with libnimble.so).
## There is no call to `clearCompiled` in the models testing, so that shouldn't be the issue.
## Not sure how able to get the timing above.

nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)






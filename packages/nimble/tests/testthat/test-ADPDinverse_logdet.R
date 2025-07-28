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

# Standalone test
cov <- crossprod(matrix(rnorm(25), 5))

PDinverse_logdet_test <- make_AD_test2(
  op = list(
    name = "PDinverse_logdet with prec_param FALSE: positive definite inverse and log determinant test",
    opParam = list(name = "PDinverse_logdet"),
    expr = quote({
      pdl <- PDinverse_logdet(mat, 0)
      out <- pdl
    }),
    args = list(
      mat = quote(double(2))
    ),
    outputType = quote(double(1))
  ),
  argTypes = c(mat='double(2)'),
  wrt = c('mat'),
  inputs = list(record = list(mat = cov),
                test   = list(mat = cov+0.1))
)
PDinverse_logdet_test_out <- test_AD2(PDinverse_logdet_test)

# Test in a model with dmnormAD

if(FALSE) { # This is seg-faulting with: "corrupted size vs. prev_size while consolidating"
set.seed(1)
n <- 3
locs <- cbind(runif(n),runif(n))
dd <- matrix(nrow = n, ncol = n)
for(i in 1:n) for(j in 1:n) dd[i,j] <- sqrt(sum((locs[i,]-locs[j,])^2))
code <- nimbleCode({
  sigma ~ dunif(0, 20)
  rho ~ dunif(0,3)
  cov[1:3,1:3] <- sigma*sigma*exp(-dd[1:3,1:3]/rho)
  prec_ldet[1:10] <- PDinverse_logdet(cov[1:3, 1:3], 0)
  x[1:3] ~ dmnormAD(mean = mu[1:3], prec_ldet = prec_ldet[1:10])
})
inits <- list(sigma = 0.8, rho = 1.2, mu = c(.2, .3, .1))
constants = list(dd=dd)
data = list(x = c(.1, .2, .3))
model <- nimbleModel(code, data = data, constants = constants, inits = inits, buildDerivs=TRUE)

relTolTmp <- relTol
relTolTmp[2] <- 1e-6
relTolTmp[3] <- 1e-2
relTolTmp[4] <- 1e-2
relTolTmp[5] <- 1e-13
test_ADModelCalculate(model, useParamTransform = TRUE,
                      newConstantNodes = list(dd = dd*1.1, x=c(.15, .25, .35)),
                      newUpdateNodes = list(mu = c(.22, .32, .12)),
                      checkCompiledValuesIdentical = FALSE, checkDoubleUncHessian = TRUE,
                      useFasterRderivs = TRUE, verbose = TRUE, name = 'dmnormAD with prec_ldet')

source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
nimbleOptions(experimentalEnableDerivs = TRUE)
context("Testing of derivatives for nimbleFunctions.")

test_that('Derivatives of dnorm function correctly.',
  {
    ADfun1 <- nimbleFunction(
      setup = function(){},
      run = function(x = double(0)) {
        outList <- derivs(testMethod(x))
        returnType(ADNimbleList())
        return(outList)
      },
      methods = list(
        testMethod = function(x = double(0)) {
          out <- dnorm(x,0,1)
          returnType(double())
          return(out)
        }
      ), enableDerivs = list('testMethod')
    )
    
    ADfunInst <- ADfun1()
    temporarilyAssignInGlobalEnv(ADfunInst)  
    cADfunInst <- compileNimble(ADfunInst)
    
    Rderiv <- D(expression((1/(sqrt(2*pi)))*exp(-(x^2)/2)), 'x') ## Can be replaced by NIMBLE's R version of derivs in the future.
    x <- 1.4
    expect_equal(cADfunInst$run(x)$gradient[1], eval(Rderiv)) ## Temporary simple test to make sure compilation and gradient calculation work.
  }
)

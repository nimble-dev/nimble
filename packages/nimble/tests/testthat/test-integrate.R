
library(nimble)

# integrand should take a vector and return a vector
integrand <- nimbleFunction(
  run = function(x = double(1)) {
    return(3.1415927 * x)
    returnType(double(1))
  }
)

foo <- nimbleFunction(
  run = function() {
    return(integrate(integrand, 0, 1))
    returnType(double())
  }
)
nimbleOptions(showCompilerOutput = TRUE)

comp <- compileNimble(foo, integrand)

comp$integrand(c( .5, 1)) # check that integrand works
foo() # should return 0.5 * pi
comp$foo()

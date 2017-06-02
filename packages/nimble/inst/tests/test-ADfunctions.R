nf <- nimbleFunction(
  setup = function(){},
  run = function(x = double(0)) {
    outList <- derivs(testMethod(x))
    returnType(ADNimbleList())
    return(outList)
  },
  methods = list(
    testMethod = function(x = double(0)) {
      out <- dnorm(x,0.1,1.1)
      returnType(double())
      return(out)
    }
  ), enableDerivs = list('testMethod')
)
nf1 <- nf()
cnf1 <- compileNimble(nf1)
x <- 2.0
cnf1$run(x)



nf <- nimbleFunction(
  setup = function(){},
  run = function(x = double(0)) {
    out <- testMethod(x)
    returnType(double())
    return(out)
  },
  methods = list(
    testMethod = function(x = double(0)) {
      out <- dnorm(x,0.1,1.1)
      returnType(double())
      return(out)
    }
  ), enableDerivs = list('testMethod')
)
nf1 <- nf()
cnf1 <- compileNimble(nf1)
x <- 2.0
cnf1$run(x)

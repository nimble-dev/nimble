test_that('Addition works', {
      fun <- nimbleFunction(run = function(x = double(0)) {
          returnType(double(0))
          return(x + x)
      })
      expect_equal(fun(0), compileNimble(fun)(0))
})

test_that('One equals zero', {
      expect_failure(
          expect_equal(1, 0),
          info = 'KNOWN ISSUE https://github.com/nimble-dev/nimble/issues/10'
      )
})


  test_that('Function works with vectors', {
      sizes <- c(1,2,3,10,100)
      for (size in sizes) {
          x <- rep(1, size)
          y <- rep(2, size)
          expect_equal(x + x, y, info = paste(' where size =', size))
      })
  }

test_that('Compiler fails when no returnType is provided', {
      nf <- nimbleFunction(run = function(){ return(0) })
      expect_error(compileNimble(nf))
  })

code = nimbleCode({
    for(i in 1:4) {
        y[i] ~ dnorm(mu[k[i]], sd = 1)
    }
})

dims = list(mu = 5)
inits = list(mu = rnorm(5), k = rep(1,4))
data = list(y = rnorm(4))


    # test case where necessary dims not provided

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of default MCMC")

RwarnLevel <- options('warn')$warn
options(warn = -1)
nimbleVerboseSetting <- nimbleOptions('verbose')

testsDynIndex <- list(
    list(
        case = 'basic dynamic index',
        code = nimbleCode({
            for(i in 1:4) {
                y[i] ~ dnorm(mu[k[i]], sd = 1)
            }
            for(j in 1:5) {
                mu[j] ~ dnorm(0, 1)
            }
        }), 
        dims = list(mu = 5), inits = list(mu = rnorm(5), k = rep(1,4)), 
        data = list(y = rnorm(4)),
        expectedDeps = list(list(parent = 'k[1]', result = 'y[1]'),
                            list(parent = 'mu[1]', result = c('mu[1]', paste0('y[',1:4,']'))),
                            list(parent = 'mu', result = c(paste0('mu[',1:5,']'),
                                                           paste0('y[',1:4,']')))),
        validIndexes =list(list(var = 'k[1]', value = 1),
                             list(var = 'k[1]', value = 5)),
        invalidIndexes =list(list(var = 'k[1]', value = 0),
                             list(var = 'k[1]', value = NA),
                             list(var = 'k[1]', value = 6))
    ),
    list( 
        case = 'dynamic index of multivariate node',
        code = nimbleCode({
            for(i in 1:4) {
                y[i] ~ dnorm(mu[k[i]], sd = 1)
            }
            for(j in 1:5) {
                mu[1:5] ~ dmnorm(z[1:5], pr[1:5,1:5])
            }
        }), 
        dims = list(mu = 5), inits = list(mu = rnorm(5), k = rep(1,4),
                                          z = rep(0,5), pr = diag(5)), 
        data = list(y = rnorm(4)),
        expectedDeps = list(list(parent = 'k[1]', result = 'y[1]'),
                            list(parent = 'mu[1]', result = c('mu[1:5]', paste0('y[',1:4,']'))),
                            list(parent = 'mu', result = c(paste0('mu[1:5]'),
                                                           paste0('y[',1:4,']')))),
        validIndexes =list(list(var = 'k[1]', value = 1),
                           list(var = 'k[1]', value = 5)),
        invalidIndexes =list(list(var = 'k[1]', value = 0),
                             list(var = 'k[1]', value = 6))
    ),
    list(
        case = 'dynamic index functional',
        code = nimbleCode({
            for(i in 1:4) {
                y[i] ~ dnorm(mu[k[i]+1], sd = 1)
            }
            for(j in 1:5) {
                mu[j] ~ dnorm(0, 1)
            }
        }), 
        dims = list(mu = 5), inits = list(mu = rnorm(5), k = rep(1,4)), 
        data = list(y = rnorm(4)),
        expectedDeps = list(list(parent = 'k[1]', result = 'y[1]'),
                            list(parent = 'mu[1]', result = c('mu[1]', paste0('y[',1:4,']'))),
                            list(parent = 'mu', result = c(paste0('mu[', 1:5, ']'),
                                                           paste0('y[',1:4,']')))),
        validIndexes =list(list(var = 'k[1]', value = 0),
                           list(var = 'k[1]', value = 4)),
        invalidIndexes =list(list(var = 'k[1]', value = -1),
                             list(var = 'k[1]', value = 5))

    ),
    list(
        case = 'dynamic index multiple input functional',
        code = nimbleCode({
            for(i in 1:4) {
                y[i] ~ dnorm(mu[k[i]+j[i]], sd = 1)
            }
            for(j in 1:5) {
                mu[j] ~ dnorm(0, 1)
            }
        }), 
        dims = list(mu = 5), inits = list(mu = rnorm(5), k = rep(1,4), j = rep(1,4)), 
        data = list(y = rnorm(4)),
        expectedDeps = list(list(parent = 'k[1]', result = 'y[1]'),
                            list(parent = 'j[2]', result = 'y[2]'),
                            list(parent = 'mu[1]', result = c('mu[1]', paste0('y[',1:4,']'))),
                            list(parent = 'mu', result = c(paste0('mu[', 1:5, ']'),
                                                           paste0('y[',1:4,']')))),
        validIndexes =list(list(var = c('k[1]','j[1]'), value = c(0,2))),
        invalidIndexes =list(list(var = c('k[1]', 'j[1]'), value = c(0,6)))
    ),
    list(  
        case = 'dynamic index multiple input functional, multivariate indexing',
        code = nimbleCode({
            for(i in 1:4) {
                y[i] ~ dnorm(mu[k[i], 2, k[i]+j[i]], sd = 1)
            }
            for(j in 1:5) 
                for(l in 1:3) 
                    for(m in 1:8)
                mu[j,l,m] ~ dnorm(0, 1)
        }), 
        inits = list(mu = array(rnorm(120), c(5,3,8)), k = rep(1,4), j = rep(1:4)), 
        data = list(y = rnorm(4)),
        expectedDeps = list(list(parent = 'k[1]', result = 'y[1]'),
                            list(parent = 'j[2]', result = 'y[2]'),
                            list(parent = 'mu[1,1,1]', result = c('mu[1, 1, 1]')),
                            list(parent = 'mu[1,2,1]', result = c('mu[1, 2, 1]',
                                 expandNames('y', 1:4))),
                            list(parent = 'mu', result = c(expandNames('mu', 1:5, 1:3, 1:8),
                                                           expandNames('y',1:4)))),
        validIndexes =list(list(var = c('k[1]','j[1]'), value = c(3,5))),
        invalidIndexes =list(list(var = c('k[1]', 'j[1]'), value = c(3,6)),
                             list(var = c('k[1]', 'j[1]'), value = c(6,2)))
    )
    list( 
        case = 'dynamic index, multivariate node indexing',
        code = nimbleCode({
            for(i in 1:4) 
                y[i, 1:3] ~ dmnorm(mu[k[i], 1:3], pr[1:3, 1:3])
            for(i in 1:3)
                mu[i, 1:3] ~ dmnorm(z[1:3], pr[1:3, 1:3])
        }), 
        inits = list(mu = matrix(rnorm(9), 3), k = rep(1, 4), z = rep(0, 3),
                     pr = diag(3)),
        data = list(y = matrix(rnorm(12), 4, 3)),
        expectedDeps = list(list(parent = 'k[1]', result = 'y[1, 1:3]'),

                            list(parent = 'mu[1,2]', result = c('mu[1, 1:3]',
                                                                expandNames('y', 1:4, "1:3"))),
                            list(parent = 'mu', result = c(expandNames('mu', 1:3, "1:3"),
                                                           expandNames('y', 1:4, "1:3")))),
        validIndexes =list(list(var = c('k[1]'), value = 3)),
        invalidIndexes =list(list(var = c('k[1]'), value = 0),
                             list(var = c('k[1]'), value = 4))
    ),                             
    list( 
        case = 'dynamic index, crossed multivariate node indexing',
        code = nimbleCode({
            for(i in 1:4) {
                y[i, 1:3] ~ dmnorm(mu[k[i], 1:3], pr[1:3, 1:3])
            }
            for(i in 1:3)
                mu[1:3, i] ~ dmnorm(z[1:3], pr[1:3, 1:3])
        }), 
        inits = list(mu = matrix(rnorm(9), 3), k = rep(1, 4), z = rep(0, 3),
                     pr = diag(3)),
        data = list(y = matrix(rnorm(12), 4, 3)),
        expectedDeps = list(list(parent = 'k[1]', result = 'y[1, 1:3]'),
                            list(parent = 'mu[2,2]', result = c('mu[1:3, 2]',
                                                                expandNames('y', 1:4, "1:3"))),
                            list(parent = 'mu', result = c(expandNames('mu', "1:3", 1:3),
                                                           expandNames('y', 1:4, "1:3")))),
        validIndexes =list(list(var = c('k[1]'), value = 3)),
        invalidIndexes =list(list(var = c('k[1]'), value = 0),
                             list(var = c('k[1]'), value = 4))
    ),     
    list(  # This should pass, but only if we don't test bounds for nested index, d[i].
        case = 'dynamic index multiple input functional, multiple indexing, with incorrect nested index',
        code = nimbleCode({
            for(i in 1:4) {
                y[i] ~ dnorm(mu[k[d[i]], 2, k[i]+j[i]], sd = 1)
            }
            ## k[i]+j[i] \in 1:8; d[i] \in 1:2; k[i] \in 1:5
            for(ii in 1:5) 
                for(jj in 1:3) 
                    for(kk in 1:8)
                        mu[ii,jj,kk] ~ dnorm(0, 1)
            for(i in 1:4)
                k[i] ~ dcat(p[1:5])
        }), 
        inits = list(mu = array(rnorm(120), c(5,3,8)), k = rep(1,4), j = rep(1, 4), d = rep(1,4)),
        data = list(y = rnorm(4)),
        expectedDeps = list(list(parent = 'd[2]', result = 'y[2]'),
                            list(parent = 'k[2]', result = c('k[2]', expandNames('y', 1:4))),
                            list(parent = 'j[2]', result = 'y[2]'),
                            list(parent = 'mu[1,1,1]', result = c('mu[1, 1, 1]')),
                            list(parent = 'mu[1,2,1]', result = c('mu[1, 2, 1]',
                                                                  expandNames('y', 1:4))),
                            list(parent = 'mu', result = c(expandNames('mu', 1:5, 1:3, 1:8),
                                                           expandNames('y', 1:4)))),
        validIndexes =list(list(var = c('d[1]', 'k[1]','j[1]'), value = c(2,5,2))),
        invalidIndexes =list(list(var = c('d[1]', 'k[1]','j[1]'), value = c(1,5,4)),
                             list(var = c('d[1]', 'k[1]','j[1]'), value = c(1,6,2)))
    )
    
)


testsInvalidDynIndex <- list(
    list(
        case = 'vector dynamic index',
        code = nimbleCode({
            y[1:3] ~ dmnorm(mu[k[1:3]], pr[1:3,1:3])
            mu[1:5] ~ dmnorm(z[1:5], pr[1:5,1:5])
        }), 
        inits = list(mu = rnorm(5)), 
        data = list(y = rnorm(3)),
        expectError = TRUE
    )
)


## These are cases NIMBLE fails on that we are aware of. XFAIL should be reported.
testsInvalidDynIndexExpectedFailures <- list(
    list(  # This should pass with correct error trapping but it does not because we don't test bounds for nested index, d[i], so failure occurs in R execution (C execution fails gracefully because of how C++ handles k[d[0]].
        case = 'dynamic index multiple input functional, multiple indexing, with nested index',
        code = nimbleCode({
            for(i in 1:4) {
                y[i] ~ dnorm(mu[k[d[i]], 2, k[i]+j[i]], sd = 1)
            }
            ## k[i]+j[i] \in 1:8; d[i] \in 1:2; k[i] \in 1:5
            for(ii in 1:5) 
                for(jj in 1:3) 
                    for(kk in 1:8)
                        mu[ii,jj,kk] ~ dnorm(0, 1)
            for(i in 1:4)
                k[i] ~ dcat(p[1:5])
        }), 
        inits = list(mu = array(rnorm(120), c(5,3,8)), k = rep(1,4), j = rep(1, 4), d = rep(1,4)),
        data = list(y = rnorm(4)),
        invalidIndexes =list(list(var = c('d[1]', 'k[1]','j[1]'), value = c(5,4,4))),
        expectFailure = TRUE
    )


)


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



    # test case where necessary dims not provided

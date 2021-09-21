## These tests check building, compilation and calculation with dynamically-indexed models. Testing includes that errors are correctly emitted when invalid indexes are used via the 'invalidIndexes' element. In a few cases with nested indexing, R calculation fails with error in condition of an if() and we don't catch that error, but compiled calculation fails gracefully.

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
        case = 'basic dynamic index, no context',
        code = nimbleCode({
            y ~ dnorm(mu[k+1], 1)
            k ~ dcat(p[1:3])
            for(i in 1:3)
                mu[i] ~ dnorm(0,1)
        }), 
        dims = list(mu = 3), inits = list(mu = rnorm(3), k = 1), 
        data = list(y = rnorm(1)),
        expectedDeps = list(list(parent = 'k', result = c('k','y')),
                            list(parent = 'mu[1]', result = c('mu[1]', 'y')),
                            list(parent = 'mu', result = c(expandNames('mu', 1:3), 'y'))),
        validIndexes =list(list(var = 'k', value = 0)),
        invalidIndexes =list(list(var = 'k[1]', value = -1),
                             list(var = 'k[1]', value = NA),
                             list(var = 'k[1]', value = 4))
    ),
    list(
        case = 'basic dynamic index in deterministic expression',
        code = nimbleCode({
            for(i in 1:4) {
                y[i] <- exp(mu[k[i]+1])
                z[i] ~ dnorm(y[i], 1)
            }
            for(j in 1:5) {
                mu[j] ~ dnorm(0, 1)
            }
        }), 
        inits = list(mu = rnorm(5), k = rep(1,4)),
        data = list(z = rnorm(4)),
        expectedDeps = list(list(parent = 'k[1]', result = 'z[1]'),
                            list(parent = 'mu[1]', result = c('mu[1]', paste0('z[',1:4,']'))),
                            list(parent = 'mu', result = c(paste0('mu[',1:5,']'),
                                                           paste0('z[',1:4,']')))),
        validIndexes =list(list(var = 'k[1]', value = 0),
                             list(var = 'k[1]', value = 4)),
        invalidIndexes =list(list(var = 'k[1]', value = -1),
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
        case = 'basic dynamic index with calculated index',
        code = nimbleCode({
            for(i in 1:4) {
                k[i] <- floor(j[i] + 1)
                y[i] ~ dnorm(mu[k[i]], sd = 1)
                j[i] ~ dcat(p[1:4])
            }
            for(i in 1:5) {
                mu[i] ~ dnorm(0, 1)
            }
        }), 
        inits = list(mu = rnorm(5), j = rep(1,4)), 
        data = list(y = rnorm(4)),
        expectedDeps = list(list(parent = 'k[1]', result = 'y[1]'),
                            list(parent = 'j[1]', result = c('j[1]','y[1]')),
                            list(parent = 'mu[1]', result = c('mu[1]', paste0('y[',1:4,']'))),
                            list(parent = 'mu', result = c(paste0('mu[',1:5,']'),
                                                           paste0('y[',1:4,']')))),
        validIndexes =list(list(var = 'j[1]', value = 0),
                             list(var = 'j[1]', value = 4)),
        invalidIndexes =list(list(var = 'j[1]', value = -1),
                             list(var = 'j[1]', value = NA),
                             list(var = 'j[1]', value = 5))
    ),
    list(
        case = 'basic dynamic index with non-simple index',
        code = nimbleCode({
            for(i in 1:4) {
                y[i] ~ dnorm(mu[k[9-i]], 1)
            }
            for(i in 1:5) {
                mu[i] ~ dnorm(0, 1)
            }
            for(i in 2:9)
                k[i] ~ dcat(p[1:5])
        }), 
        inits = list(mu = rnorm(5), k = rep(1, 9), p = rep(1/5,5)),
        data = list(y = rnorm(4)),
        expectedDeps = list(list(parent = 'k[2]', result = 'k[2]'),
                            list(parent = 'k[6]', result = c('k[6]','y[3]')),
                            list(parent = 'mu[1]', result = c('mu[1]', paste0('y[',4:1,']'))),
                            list(parent = 'mu', result = c(paste0('mu[',1:5,']'),
                                                           paste0('y[',4:1,']')))),
        validIndexes =list(list(var = 'k[1]', value = 0),
                             list(var = 'k[6]', value = 5)),
        invalidIndexes =list(list(var = 'k[6]', value = -1),
                             list(var = 'k[6]', value = NA),
                             list(var = 'k[6]', value = 6))
    ),
    list(
        case = 'basic dynamic index with lifted node',
        code = nimbleCode({
            for(i in 1:4) {
                y[i] ~ dnorm(mu, tau[k[i]])
            }
            for(i in 1:6)
                tau[i] ~ dgamma(1, 1)
        }), 
        inits = list(mu = rnorm(1), tau = rep(1,6), k = rep(1,4)), 
        data = list(y = rnorm(4)),
        expectedDeps = list(list(parent = 'k[1]', result = 'y[1]'),
                            list(parent = 'mu', result = expandNames('y', 1:4)),
                            list(parent = 'tau[1]', result = c('tau[1]',
                                                               expandNames('y', 1:4))),
                            list(parent = 'tau', result = c(expandNames('tau', 1:6),
                                                            expandNames('y', 1:4)))),
        validIndexes =list(list(var = 'k[1]', value = 6)),
        invalidIndexes =list(list(var = 'k[1]', value = -1),
                             list(var = 'k[1]', value = NA),
                             list(var = 'k[1]', value = 7))
    ),
    list(
        case = 'dynamic index nested indexing',
        code = nimbleCode({
            for(i in 1:4) {
                y[i] ~ dnorm(mu[k[j[i]+1]], sd = 1)
            }
            for(i in 1:5) {
                mu[i] ~ dnorm(0, 1)
            }
            for(i in 1:8)
                k[i] ~ dcat(p[1:5])
        }), 
        inits = list(mu = rnorm(5), k = rep(1,8), j = rep(1,4), p = rep(1/5,5)), 
        data = list(y = rnorm(4)),
        expectedDeps = list(list(parent = 'j[2]', result = 'y[2]'),
                            list(parent = 'k[2]', result = c('k[2]', expandNames('y', 1:4))),
                            list(parent = 'mu[1]', result = c('mu[1]', paste0('y[',1:4,']'))),
                            list(parent = 'mu', result = c(paste0('mu[', 1:5, ']'),
                                                           paste0('y[',1:4,']')))),
        validIndexes =list(list(var = c('k[8]','j[1]'), value = c(5,7))),
        invalidIndexes = list(list(var = c('k[1]', 'j[1]'), value = c(1,8)),  ## mu[k[j[1]+1]] = mu[k[9]] = mu[NA] 
                              list(var = c('k[2]', 'j[1]'), value = c(6,1)))
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
    ),
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
    ),
    list(
        case = 'multiple dynamic indexes on one node',
        code = nimbleCode({
            for(i in 1:4) {
                y[i] ~ dnorm(mu[k[i]+1], sd = exp(mu[k[i]]))
            }
            for(j in 1:5) {
                mu[j] ~ dnorm(0, 1)
            }
        }), 
        dims = list(mu = 5), inits = list(mu = rnorm(5), k = rep(1,4)), 
        data = list(y = rnorm(4)),
        expectedDeps = list(list(parent = 'k[1]', result = 'y[1]'),
                            list(parent = 'mu[1]', result = c('mu[1]', expandNames('y', 1:4))),
                            list(parent = 'mu', result = c(expandNames('mu', 1:5),
                                                           expandNames('y', 1:4)))),
        validIndexes =list(list(var = 'k[1]', value = 1)),
        invalidIndexes =list(list(var = 'k[1]', value = 5))
    ),
    list(
        case = 'multiple dynamic indexes on two nodes',
        code = nimbleCode({
            for(i in 1:4) {
                y[i] ~ dnorm(mu[k[i]], sd = sigma[k[i]])
                k[i] ~ dcat(p[1:5])
            }
            for(j in 1:5) {
                mu[j] ~ dnorm(0, 1)
            }
            for(j in 1:6) {
                sigma[j] ~ dgamma(1, 1)
            }
        }), 
        inits = list(mu = rnorm(5), sigma = rep(1,6), k = rep(1,4)), 
        data = list(y = rnorm(4)),
        expectedDeps = list(list(parent = 'k[1]', result = c('k[1]', 'y[1]')),
                            list(parent = 'sigma[1]', result = c('sigma[1]', expandNames('y', 1:4))),
                            list(parent = 'mu[2]', result = c('mu[2]', expandNames('y', 1:4))),
                            list(parent = 'sigma', result = c(expandNames('sigma', 1:6),
                                                           expandNames('y', 1:4)))),
        validIndexes =list(list(var = 'k[4]', value = 5)),
        invalidIndexes =list(list(var = 'k[4]', value = 6))
    ),
    list(  
        case = 'dynamic index multiple indexing of various types',
        code = nimbleCode({
            for(i in 1:4) {
                y[i] ~ dnorm(mu[i, 2, k[i]], sd = 1)
                k[i] ~ dcat(p[1:5])
            }
            for(ii in 1:5) 
                for(jj in 1:3) 
                    for(kk in 1:8)
                        mu[ii,jj,kk] ~ dnorm(0, 1)
        }), 
        inits = list(mu = array(rnorm(120), c(5,3,8)), k = rep(1,4)),
        data = list(y = rnorm(4)),
        expectedDeps = list(list(parent = 'k[2]', result = c('k[2]', 'y[2]')),
                            list(parent = 'mu[1,1,1]', result = c('mu[1, 1, 1]')),
                            list(parent = 'mu[1,2,1]', result = c('mu[1, 2, 1]', 'y[1]')),
                            list(parent = 'mu', result = c(expandNames('mu', 1:5, 1:3, 1:8),
                                                           expandNames('y', 1:4)))),
        validIndexes =list(list(var = 'k[4]', value = 8)),
        invalidIndexes =list(list(var = 'k[4]', value = 9))
    ),
    list(  
        case = 'dynamic index multivariate nested indexing',
        code = nimbleCode({
            for(i in 1:4) {
                y[i] ~ dnorm(mu[1, k[2, j[i]]], 1)
                j[i] ~ dcat(p[1:5])
            }
            for(ii in 1:2)
                for(jj in 1:3) 
                    mu[ii, jj] ~ dnorm(0,1)
            p[1:5] ~ ddirch(alpha[1:5])
            for(i in 1:5)
                k[2, i] ~ dcat(q[1:3])
        }),
        inits = list(mu = matrix(rnorm(6),2,3), p = rep(1/5,5), q = rep(1/3,3),
                     j = rep(1, 4), k = matrix(1,2,5)),
        data = list(y = rnorm(4)),
        expectedDeps = list(list(parent = 'j[2]', result = c('j[2]', 'y[2]')),
                            list(parent = 'k[2,2]', result = c('k[2, 2]',
                                                               expandNames('y', 1:4))),
                            list(parent = 'mu[2,1]', result = c('mu[2, 1]')),
                            list(parent = 'mu[1,1]', result = c('mu[1, 1]',
                                                                expandNames('y', 1:4))),
                            list(parent = 'mu', result = c(expandNames('mu', 1:2, 1:3),
                                                           expandNames('y', 1:4)))),
        validIndexes =list(list(var = c('j[1]','k[2,4]'), value = c(4, 3))),
        invalidIndexes =list(list(var = c('j[1]','k[2,4]'), value = c(0,1)),  ## mu[1, k[2, j[1]]] = mu[1, k[2, 0]] = mu[1, NA]  
                             list(var = c('j[1]','k[2,1]'), value = c(1,4)))
    ),
    list(  
        case = 'dynamic index bivariate index',
        code = nimbleCode({
            for(i in 1:2) {
                for(j in 1:4) {
                    y[i, j] ~ dnorm(mu[i, k[i,j], 3], sd = sigma[k[i,j]])
                    k[i,j] ~ dcat(p[i,1:3])
                }
                p[i,1:3] ~ ddirch(alpha[1:3])
            }
            for(i in 1:2) {
                for(j in 1:3) {
                    mu[i,j,3] ~ dnorm(0,1)
                }}
            for(j in 1:3)
                sigma[j] ~ dhalfflat()
        }), 
        inits = list(mu = array(rnorm(12), c(2,3,3)), k = matrix(1, 2, 4)),
        data = list(y = matrix(rnorm(8),2,4)),
        expectedDeps = list(list(parent = 'k[2,4]', result = c('k[2, 4]', 'y[2, 4]')),
                            list(parent = 'sigma[2]', result = c('sigma[2]',
                                                                 expandNames('y', 1:2, 1:4))),
                            list(parent = 'mu[1, 1, 3]', result = c('mu[1, 1, 3]',
                                                                    expandNames('y', 1, 1:4))),
                            list(parent = 'mu', result = c(expandNames('mu', 1:2, 1:3, 3),
                                                           expandNames('y', 1:2, 1:4)))),
        validIndexes =list(list(var = 'k[2,4]', value = 3)),
        invalidIndexes =list(list(var = 'k[2,4]', value = 4))
    ),
    list(  
        case = 'dynamic index with indexes in separate statements same indexed var',
        code = nimbleCode({
            for(i in 1:4) {
                y[i] ~ dnorm(mu[k[i], 1], sd = 1)
                k[i] ~ dcat(p[1:5])
            }
            for(i in 1:8) {
                z[i] ~ dnorm(mu[1, j[i]], 1)
                j[i] ~ dcat(q[1:3])
            }
            for(i in 1:5)
                for(j in 1:3)
                    mu[i,j] ~ dnorm(0,1)
        }),
        data = list(y = rnorm(4), z = rnorm(8)),
        inits = list(k = rep(1,4), j = rep(1, 8), mu = matrix(rnorm(15), 5, 3)),
        expectedDeps = list(list(parent = 'k[2]', result = c('k[2]', 'y[2]')),
                            list(parent = 'j[8]', result = c('j[8]', 'z[8]')),
                            list(parent = 'mu[1,1]', result = c('mu[1, 1]',
                                                                expandNames('y', 1:4),
                                                                expandNames('z', 1:8))),
                            list(parent = 'mu[2,1]', result = c('mu[2, 1]',
                                                                expandNames('y', 1:4))),
                            list(parent = 'mu[1,3]', result = c('mu[1, 3]',
                                                                expandNames('z', 1:8))),
                            list(parent = 'mu', result = c(expandNames('mu', 1:5, 1:3),
                                                           expandNames('y', 1:4),
                                                           expandNames('z', 1:8)))),
        validIndexes =list(list(var = c('k[1]','j[1]'), value = c(4,3))),
        invalidIndexes =list(list(var = c('k[1]','j[1]'), value = c(0,1)),
                             list(var = c('k[1]','j[1]'), value = c(1,0)),
                             list(var = c('k[1]','j[1]'), value = c(NA, 1)),
                             list(var = c('k[1]','j[1]'), value = c(1, NA)))
    ),
    list(
        case = 'dynamic index with indexes with chaining',
        code = nimbleCode({
            for(i in 1:5)
                mu[i] ~ dnorm(0, sd = 1)
            for(i in 1:9) {
                k[i] ~ dcat(p[1:5])
                y[i] ~ dnorm( mu[k[i]], sd = 1)
            }
            for(i in 1:6) {
                j[i] ~ dcat(q[1:9])
                z[i] ~ dnorm( y[j[i]], sd = 1)
            }
        }),
        data = list(z = rnorm(6)),
        inits = list(y = rnorm(9), k = rep(1,9), j = rep(1, 6), mu = rnorm(5),
                     p = rep(1/5,5), q = rep(1/9,9)),
        expectedDeps = list(list(parent = 'k[2]', result = c('k[2]', 'y[2]')),
                            list(parent = 'j[6]', result = c('j[6]', 'z[6]')),
                            list(parent = 'y[1]', result = c('y[1]',
                                                             expandNames('z', 1:6))),
                            list(parent = 'mu[2]', result = c('mu[2]',
                                                              expandNames('y', 1:9)))),
        validIndexes =list(list(var = c('k[1]','j[1]'), value = c(5,9))),
        invalidIndexes =list(list(var = c('k[1]','j[1]'), value = c(0,1)),
                             list(var = c('k[1]','j[1]'), value = c(1,0)),
                             list(var = c('k[1]','j[1]'), value = c(NA, 1)),
                             list(var = c('k[1]','j[1]'), value = c(1, NA)))
    ),
    list(
        case = 'dynamic index with truncation',
        code = nimbleCode({
            for(i in 1:5) {
                y[i] ~ T(dnorm(mu[k[i]],1), a[i], 10)
                a[i] ~ dunif(0,1)
            }
            for(i in 1:7)
                mu[i] ~ T(dnorm(0, 1), 0, 1)
        }),
        inits = list(mu = rnorm(7), k = rep(1,5), a = rep(0.5, 5)), 
        data = list(y = rnorm(5)),
        expectedDeps = list(list(parent = 'k[1]', result = c('y[1]')),
                            list(parent = 'a[2]', result = c('a[2]', 'y[2]')),
                            list(parent = 'mu[2]', result = c('mu[2]', expandNames('y', 1:5)))),
        validIndexes =list(list(var = 'k[4]', value = 7)),
        invalidIndexes =list(list(var = 'k[4]', value = 8))
    ),
    list(
        case = 'dynamic index plus index specified by constant',
        code = nimbleCode({
            for(i in 1:5) {
                y[i] ~ dnorm(mu[k[i], block[i]], 1)
                k[i] ~ dcat(p[1:4])
            }
        }),
        inits = list(mu = matrix(rnorm(4*2), 4, 2), k = rep(1,5)),
        data = list(y = rnorm(5)),
        constants = list(block = c(1,1,2,1,2)),
        expectedDeps = list(list(parent = 'mu[1,1]', result = c('y[1]', 'y[2]', 'y[4]')),
                            list(parent = 'mu[2,1]', result = c('y[1]', 'y[2]', 'y[4]')),
                            list(parent = 'mu[2,2]', result = c('y[3]', 'y[5]')),
                            list(parent = 'k[2]', result = c('k[2]', 'y[2]'))),
        validIndexes =list(list(var = 'k[4]', value = 4)),
        invalidIndexes =list(list(var = 'k[4]', value = 5))
    ),
    list(  
        case = 'dynamic index multiple input functional, multiple indexing, with invalid nested index',
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
        inits = list(mu = array(rnorm(120), c(5,3,8)), k = rep(1,4), j = rep(1, 4), d = rep(1,4), p = rep(1/5, 5)),
        data = list(y = rnorm(4)),
        invalidIndexes =list(list(var = c('d[1]', 'k[1]','j[1]'), value = c(5,4,4)))  ## mu[k[d[1]]] = mu[k[5]] = mu[NA]  
    )
)

## These are cases that NIMBLE should error out of.
testsInvalidDynIndex <- list(
    list(
        case = 'vector dynamic index (invalid)',
        code = nimbleCode({
            y[1:3] ~ dmnorm(mu[k[1:3]], pr[1:3,1:3])
            mu[1:5] ~ dmnorm(z[1:5], pr[1:5,1:5])
        }), 
        inits = list(mu = rnorm(5)), 
        data = list(y = rnorm(3)),
        expectError = TRUE,
        expectErrorMsg = "only scalar indices are allowed"
    ),
    list(
        case = 'vector dynamic index in multidimensional (invalid)',
        code = nimbleCode({
            for(i in 1:2) 
                y[i, 1:3] ~ dmnorm(mu[k[i, 1:3]], pr[1:3,1:3])
        }), 
        inits = list(mu = rnorm(10)), 
        data = list(y = matrix(rnorm(6), 2, 3)),
        expectError = TRUE,
        expectErrorMsg = "only scalar indices are allowed"
    ),
    list(
        case = 'dynamic index sequence (invalid)',
        code = nimbleCode({
            for(i in 1:2) {
                y[i,1:2] ~ dnorm(mu[k[i]:(k[i]+1)], sd = sigma)
            }
        }), 
        inits = list(mu = rnorm(5)),
        expectError = TRUE,
        expectErrorMsg = "Dynamic indexing found in a vector"
    ),
    list(
        case = 'dynamic index sequence, multiple dynamic indexes (invalid)',
        code = nimbleCode({
            for(i in 1:2) {
                y[i,1:2] ~ dnorm(mu[k[i]:j[i]], sd = sigma)
            }
        }), 
        inits = list(mu = rnorm(5)),
        expectError = TRUE,
        expectErrorMsg = "Dynamic indexing found in a vector"
    ),
    list(
        case = 'dynamic indexing of constants',
        code = nimbleCode({
            for(i in 1:2) {
                y[i] ~ dnorm(mu[k[i]], sd = sigma)
            }
        }), 
        constants = list(mu = rnorm(5)),
        expectError = TRUE,
        expectErrorMsg = "dynamic indexing of constants"
    )
)





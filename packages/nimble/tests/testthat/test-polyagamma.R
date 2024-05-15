source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

context("Testing of Polya-gamma sampler")

RwarnLevel <- options('warn')$warn
options(warn = 1)

## verbose: set to FALSE
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

## MCMC progress bar: set to FALSE
nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

test_that('polyagamma validity checks', {
    ## Various valid basic model structures
    code <- nimbleCode({
        for(i in 1:n) {
            p0[i] <- expit(b0+b1*x[i])
            y0[i] ~ dbern(p0[i])
            
            logit(p1[i]) <- a0 + a1*x[i]
            y1[i]~dbern(p1[i])
            
            logit(p2[i]) <- inprod(x2[i,1:p], b2[1:p])
            y2[i]~dbin(p2[i], size = m)
            
            logit(p3[i])  <- x2[i,1:p] %*% b3[1:p]
            y3[i]~dbin(p3[i], 1)            
        }
        m ~ dpois(5)
        
        b0~dnorm(0, sd=10)
        b1~dnorm(0, sd=10)
        a0~dnorm(0, sd=10)
        a1~dnorm(0, sd=10)
        
        b2[1:p] ~ dmnorm(mu[1:p], Q[1:p,1:p])
        for(i in 1:p)
            b3[i] ~ dnorm(0, sd = 10) 
    })

    n <- 10
    p <- 3
    constants <- list(n=n, p=p, x = runif(n), x2 = matrix(runif(n*p),n))
    ys <- rep(1,n)
    data <- list(y0 = ys, y1 = ys, y2 = ys, y3 = ys)
    
    m <- nimbleModel(code, constants = constants, data = data)
    conf <- configureMCMC(m, nodes = NULL)
    expect_silent(conf$addSampler(type='polyagamma', target=c('b0','b1')))
    expect_silent(conf$addSampler(type='polyagamma', target=c('a0','a1')))
    expect_silent(conf$addSampler(type='polyagamma', target=c('b2')))
    expect_silent(conf$addSampler(type='polyagamma', target=c('b3')))
    expect_silent(mcmc <- buildMCMC(conf))

    ## random effects check design matrix
    code <- nimbleCode({
        for(i in 1:n) {
            p0[i] <- expit(b0+b1*x[i]+u[k[i]])
            y0[i] ~ dbern(p0[i])
        }
                
        b0~dnorm(0, sd=10)
        b1~dnorm(0, sd=10)

        for(i in 1:p)
            u[i] ~ dnorm(0,1)
    })

    
    n <- 9
    p <- 3
    constants <- list(n=n, p=p, x = runif(n), k = rep(1:3, each = 3))
    ys <- rep(1,n)
    data <- list(y0 = ys)
    
    m <- nimbleModel(code, constants = constants, data = data)
    conf <- configureMCMC(m, nodes = NULL)
    expect_silent(conf$addSampler(type='polyagamma', target=c('b0','u','b1')))
    expect_silent(mcmc <- buildMCMC(conf))
    mcmc$run(1)
    expect_equal(mcmc$samplerFunctions[[1]]$X,
                     cbind(rep(1,n), c(rep(1,3),rep(0,6)),
                           c(rep(0,3),rep(1,3),rep(0,3)),
                           c(rep(0,6),rep(1,3)), constants$x))

    ## dnorm + dmnorm
    code <- nimbleCode({
        for(i in 1:n) {
            p0[i] <- expit(b[1]+u[k[i]]+b[2]*x[i])
            y0[i] ~ dbern(p0[i])
        }

        b[1:2] ~ dmnorm(z[1:2], pr[1:2,1:2])
        
        for(i in 1:p)
            u[i] ~ dnorm(0,1)
    })

    
    n <- 9
    p <- 3
    constants <- list(n=n, p=p, x = runif(n), k = rep(1:3, each = 3), z=rep(0,2), pr=diag(2))
    ys <- rep(1,n)
    data <- list(y0 = ys)
    
    m <- nimbleModel(code, constants = constants, data = data)
    conf <- configureMCMC(m, nodes = NULL)
    expect_silent(conf$addSampler(type='polyagamma', target=c('b','u')))
    expect_silent(mcmc <- buildMCMC(conf))
    mcmc$run(1)
    expect_equal(mcmc$samplerFunctions[[1]]$X,
                     cbind(rep(1,n), constants$x, c(rep(1,3),rep(0,6)),
                           c(rep(0,3),rep(1,3),rep(0,3)),
                           c(rep(0,6),rep(1,3))))
    lens <- c(2,1,1,1); names(lens) <- c('b[1:2]','u[1]','u[2]','u[3]')
    expect_equal(mcmc$samplerFunctions[[1]]$nodeLengths, lens)

    
    ## zero-inflated
    ## despite presence of initial structural zeroes.
    
    code <- nimbleCode({
        for(i in 1:n) {
            p0[i] <- expit(b0+b1*x[i])
            y0[i] ~ dbern(z1[i]*z2[i]*p0[i])
            z1[i] ~ dbin(q, 1)
            z2[i] ~ dbern(q)
        }
        
        b0~dnorm(0, sd=10)
        b1~dnorm(0, sd=10)
    })

    n <- 10
    constants <- list(n=n, x = runif(n), q = .75)
    data <- list(y0 = c(1,0,1,0,1,0,1,0,1,0))
    inits <- list(z1 = c(1,1,1,1,1,0,1,0,1,0), z2 = c(1,0,1,0,1,0,1,0,1,1))

    set.seed(1)
    m <- nimbleModel(code, constants = constants, data = data, inits = inits)
    conf <- configureMCMC(m, nodes = c('z1','z2'))
    expect_silent(conf$addSampler(type='polyagamma', target=c('b0','b1'), control = list(fixedDesignColumns = TRUE)))
    expect_silent(mcmc <- buildMCMC(conf))
    mcmc$run(1)
    ## Check that designMatrix properly filled in initially, despite presence of initial structural zeroes.
    expect_equal(mcmc$samplerFunctions[[21]]$X,
                 cbind(rep(1,n), constants$x))
    expect_equal(mcmc$samplerFunctions[[21]]$n, sum(m$z1*m$z2))
    expect_equal(mcmc$samplerFunctions[[21]]$probNonZero[1:mcmc$samplerFunctions[[21]]$n],
                 which(m$z1*m$z2 == 1))
    mcmc$samplerFunctions[[21]]$setProbParam()
    expect_equal(mcmc$samplerFunctions[[21]]$psi[1:mcmc$samplerFunctions[[21]]$n],
                 m$b0 + m$b1*constants$x[mcmc$samplerFunctions[[21]]$probNonZero[1:mcmc$samplerFunctions[[21]]$n]])


    code <- nimbleCode({
        for(i in 1:n) {
            p0[i] <- expit(b0+b1*x[i]+u[k[i]])
            y0[i] ~ dbern(p0[i])
            k[i] ~ dcat(q[1:3])
        }
        for(i in 1:p)
            u[i] ~ dnorm(0, 1)
        b0~dnorm(0, sd=10)
        b1~dnorm(0, sd=10)
    })

    n <- 10
    p <- 3
    constants <- list(n=n, x = runif(n), p = p, q=c(.2,.4,.4))
    data <- list(y0 = c(1,0,1,0,1,0,1,0,1,0))
    inits <- list(k = c(1,2,3,1,2,3,1,2,3,3))

    set.seed(1)
    m <- nimbleModel(code, constants = constants, data = data, inits = inits)
    conf <- configureMCMC(m, nodes = c('k'))
    expect_silent(conf$addSampler(type='polyagamma', target=c('b0','b1','u'), control = list(nonTargetNodes = 'k')))
    expect_silent(mcmc <- buildMCMC(conf))
    mcmc$run(3)
    ## Check designMatrix filled correctly based on stochastic indexing.
    result <- matrix(0, n, p)
    result[which(m$k==1),1] <- 1
    result[which(m$k==2),2] <- 1
    result[which(m$k==3),3] <- 1
    expect_equal(mcmc$samplerFunctions[[11]]$X, cbind(1, constants$x, result))
    
    
    ## Various invalid models.
    code <- nimbleCode({
        for(i in 1:n) {
            p0[i] <- z2[i]*expit(b0+b1*x[i])
            y0[i] ~ dbern(z1[i]*p0[i])
            z1[i] ~ dbin(q, 1)
            z2[i] ~ dbern(q)
        }
        
        b0~dnorm(0, sd=10)
        b1~dnorm(0, sd=10)
    })

    n <- 10
    constants <- list(n=n, x = runif(n), q = .75)
    data <- list(y0 = c(1,0,1,0,1,0,1,0,1,0))
    inits <- list(z1 = c(1,1,1,1,1,0,1,0,1,0), z2 = c(1,0,1,0,1,0,1,0,1,1))

    set.seed(1)
    m <- nimbleModel(code, constants = constants, data = data, inits = inits)
    conf <- configureMCMC(m, nodes = c('z1','z2'))
    expect_silent(conf$addSampler(type='polyagamma', target=c('b0','b1'), control = list(fixedDesignColumns = TRUE)))
    expect_error(mcmc <- buildMCMC(conf), "zero inflation cannot be specified directly")

    code <- nimbleCode({
        for(i in 1:n) {
            p0[i] <- expit(b0+exp(b1)*x[i])
            y0[i] ~ dbern(z1[i]*z2[i]*p0[i])
            z1[i] ~ dbin(q, 1)
            z2[i] ~ dbern(q)
        }
        
        b0~dnorm(0, sd=10)
        b1~dnorm(0, sd=10)
    })

    n <- 10
    constants <- list(n=n, x = runif(n), q = .75)
    data <- list(y0 = c(1,0,1,0,1,0,1,0,1,0))
    inits <- list(z1 = c(1,1,1,1,1,0,1,0,1,0), z2 = c(1,0,1,0,1,0,1,0,1,1))

    set.seed(1)
    m <- nimbleModel(code, constants = constants, data = data, inits = inits)
    conf <- configureMCMC(m, nodes = c('z1','z2'))
    expect_silent(conf$addSampler(type='polyagamma', target=c('b0','b1'), control = list(fixedDesignColumns = TRUE)))
    expect_error(mcmc <- buildMCMC(conf), "as a linear function")

    code <- nimbleCode({
        for(i in 1:n) {
            p0[i] <- expit(b0+b1*x[i])
            y0[i] ~ dbern(z1[i]*z2[i]*p0[i])
            z1[i] ~ dbin(q, 3)
            z2[i] ~ dbern(q)
        }
        
        b0~dnorm(0, sd=10)
        b1~dnorm(0, sd=10)
    })

    n <- 10
    constants <- list(n=n, x = runif(n), q = .75)
    data <- list(y0 = c(1,0,1,0,1,0,1,0,1,0))
    inits <- list(z1 = c(1,1,1,1,1,0,1,0,1,0), z2 = c(1,0,1,0,1,0,1,0,1,1))

    set.seed(1)
    m <- nimbleModel(code, constants = constants, data = data, inits = inits)
    conf <- configureMCMC(m, nodes = c('z1','z2'))
    expect_silent(conf$addSampler(type='polyagamma', target=c('b0','b1'), control = list(fixedDesignColumns = TRUE)))
    expect_error(mcmc <- buildMCMC(conf), "Zero inflation nodes must be `dbern`")
    
    code <- nimbleCode({
        for(i in 1:n) {
            p0[i] <- probit(b0+b1*x[i])
            y0[i] ~ dbern(z1[i]*z2[i]*p0[i])
            z1[i] ~ dbin(q, 1)
            z2[i] ~ dbern(q)
        }
        
        b0~dnorm(0, sd=10)
        b1~dnorm(0, sd=10)
    })

    n <- 10
    constants <- list(n=n, x = runif(n), q = .75)
    data <- list(y0 = c(1,0,1,0,1,0,1,0,1,0))
    inits <- list(z1 = c(1,1,1,1,1,0,1,0,1,0), z2 = c(1,0,1,0,1,0,1,0,1,1))

    set.seed(1)
    m <- nimbleModel(code, constants = constants, data = data, inits = inits)
    conf <- configureMCMC(m, nodes = c('z1','z2'))
    expect_silent(conf$addSampler(type='polyagamma', target=c('b0','b1'), control = list(fixedDesignColumns = TRUE)))
    expect_error(mcmc <- buildMCMC(conf), "via logit link")
    
     code <- nimbleCode({
        for(i in 1:n) {
            p0[i] <- probit(b0+b1*x[i])
            y0[i] ~ dbern(z1[i]*z2[i]*p0[i])
            z1[i] ~ dbin(q, 1)
            z2[i] ~ dbern(q)
        }
        
        b0~dflat()
        b1~dnorm(0, sd=10)
    })

    n <- 10
    constants <- list(n=n, x = runif(n), q = .75)
    data <- list(y0 = c(1,0,1,0,1,0,1,0,1,0))
    inits <- list(z1 = c(1,1,1,1,1,0,1,0,1,0), z2 = c(1,0,1,0,1,0,1,0,1,1))

    set.seed(1)
    m <- nimbleModel(code, constants = constants, data = data, inits = inits)
    conf <- configureMCMC(m, nodes = c('z1','z2'))
    expect_silent(conf$addSampler(type='polyagamma', target=c('b0','b1'), control = list(fixedDesignColumns = TRUE)))
    expect_error(mcmc <- buildMCMC(conf), "must have `dnorm`")   
                    
    code <- nimbleCode({
        for(i in 1:n) {
            p0[i] <- expit(b0+b1*x[i])
            y0[i] ~ dbern(z1[i]*z2[i]*p0[i])
            z1[i] ~ dbin(q, 1)
            z2[i] ~ dbern(q)
        }
        w ~ dnorm(b0, 3)
        b0~dnorm(0, sd=10)
        b1~dnorm(0, sd=10)
    })

    n <- 10
    constants <- list(n=n, x = runif(n), q = .75)
    data <- list(y0 = c(1,0,1,0,1,0,1,0,1,0), w = 1)
    inits <- list(z1 = c(1,1,1,1,1,0,1,0,1,0), z2 = c(1,0,1,0,1,0,1,0,1,1))

    set.seed(1)
    m <- nimbleModel(code, constants = constants, data = data, inits = inits)
    conf <- configureMCMC(m, nodes = c('z1','z2'))
    expect_silent(conf$addSampler(type='polyagamma', target=c('b0','b1'), control = list(fixedDesignColumns = TRUE)))
    expect_error(mcmc <- buildMCMC(conf), "must be distributed `dbern`")

    code <- nimbleCode({
        for(i in 1:n) {
            p0[i] <- expit(b0+b1*x[i])
            y0[i] ~ dbern(z1[i]*z2[i]*p0[i])
            z1[i] ~ dbin(q, 1)
            z2[i] ~ dbern(q)
        }
        w ~ dbern(expit(b0))
        b0~dnorm(0, sd=10)
        b1~dnorm(0, sd=10)
    })

    n <- 10
    constants <- list(n=n, x = runif(n), q = .75)
    data <- list(y0 = c(1,0,1,0,1,0,1,0,1,0), w = 1)
    inits <- list(z1 = c(1,1,1,1,1,0,1,0,1,0), z2 = c(1,0,1,0,1,0,1,0,1,1))

    set.seed(1)
    m <- nimbleModel(code, constants = constants, data = data, inits = inits)
    conf <- configureMCMC(m, nodes = c('z1','z2'))
    expect_silent(conf$addSampler(type='polyagamma', target=c('b0','b1'), control = list(fixedDesignColumns = TRUE)))
    expect_error(mcmc <- buildMCMC(conf), "should all be part of the same declaration")
})


test_that('polyagamma MCMC results', {
    ## Compare to non-PG.
    code <- nimbleCode({
        for(i in 1:n) {
            p0[i] <- expit(b0+b1*x[i,1]+b2*x[i,2])
            y0[i] ~ dbern(p0[i])
        }
        x[1,1] ~ dnorm(0, sd = 3)
        x[3,1] ~ dnorm(0, sd = 3)
        b0~dnorm(0, sd=10)
        b1~dnorm(0, sd=10)
        b2~dnorm(0, sd=10)
    })

    set.seed(1)
    n <- 1000
    b0 <- 0.3
    b1 <- 0.5
    b2 <- -0.2
    x = matrix(runif(n*2), n)
    inits <- list(b0 = 0, b1 = 0, b2 = 0)
    constants <- list(n=n)
    linpred <- b0 + b1*x[,1]+x[,2]
    x[1,1] <- NA
    x[3,1] <- NA
    data <- list(y0 = rbinom(n, size = 1, prob = expit(linpred)), x = x)

    set.seed(1)
    m <- nimbleModel(code, constants = constants, data = data, inits = inits)
    conf <- configureMCMC(m, nodes = c('x'))
    conf$addSampler(type='polyagamma', target=c('b0','b1','b2'),
                    control = list(fixedDesignColumns = c(TRUE,FALSE,TRUE),
                                   nonTargetNodes = 'x'))
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc, project = m)
    outPG <- runMCMC(cmcmc, niter=10000)

    m <- nimbleModel(code, constants = constants, data = data, inits = inits)
    conf <- configureMCMC(m, sliceOnly = TRUE)
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc, project = m)
    outDefault <- runMCMC(cmcmc, niter=10000)

    betas <- 1:3
    expect_equal(colMeans(outPG[1001:10000,betas]), colMeans(outDefault[1001:10000,betas]),
                 tolerance = 0.003)
    expect_equal(apply(outPG[1001:10000,betas], 2, sd), apply(outDefault[1001:10000,betas], 2, sd),
                 tolerance = 0.003)
    xs <- c(4,6)
    expect_equal(colMeans(outPG[1001:10000,xs]), colMeans(outDefault[1001:10000,xs]),
                 tolerance = .2)
    expect_equal(apply(outPG[1001:10000,xs], 2, sd), apply(outDefault[1001:10000,xs], 2, sd),
                 tolerance = .2)
    
})

# check Paul examples
# dbin case?

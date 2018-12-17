source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of dynamic indexing")

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

source(system.file(file.path('tests', 'dynamicIndexingTestLists.R'), package = 'nimble'))

nimbleAllowDynamicIndexingSetting <- nimbleOptions('allowDynamicIndexing')
nimbleOptions(allowDynamicIndexing = TRUE)

goldFileName <- 'dynamicIndexingTestLog_Correct.Rout'
tempFileName <- 'dynamicIndexingTestLog.Rout'
generatingGoldFile <- !is.null(nimbleOptions('generateGoldFileForDynamicIndexingTesting'))
outputFile <- if(generatingGoldFile) file.path(nimbleOptions('generateGoldFileForDynamicIndexingTesting'), goldFileName) else tempFileName

## capture warnings
sink_with_messages(outputFile)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

## check variations on use of dynamic indexing in BUGS code including building and compiling,
## dependencies, and valid and invalid dynamic index values

ans1 <- sapply(testsDynIndex, test_dynamic_indexing_model)
ans2 <- sapply(testsInvalidDynIndex, test_dynamic_indexing_model)

## check conjugacy detection

test_that("Testing normal-normal conjugacy detection with dynamic indexing", { 
    code = nimbleCode({
        for(i in 1:4) 
            y[i] ~ dnorm(mu[k[i]], sd = 1)
        for(j in 1:5) 
            mu[j] ~ dnorm(0, 1)
    })
    m = nimbleModel(code, data = list(y = rnorm(4)),
                    inits = list(k = rep(1,4)))
    conf <- configureMCMC(m)
    expect_match(conf$getSamplers()[[1]]$name, "conjugate_dnorm_dnorm",
                 info = 'failed to detect conjugacy')
})

    
test_that("Testing Poisson-normal non-conjugacy detection with dynamic indexing", { 
    code = nimbleCode({
        for(i in 1:4) 
            y[i] ~ dpois(mu[k[i]])
        for(j in 1:5)
            mu[j] ~ dnorm(0, 1)
    })
    m = nimbleModel(code, data = list(y = rpois(4, 1)),
                    inits = list(k = rep(1,4)))
    conf <- configureMCMC(m)
    expect_equal(length(grep("conjugate", conf$getSamplers()[[1]]$name)), 0,
                 info = 'incorrectly detected conjugacy')
})
    
test_that("Testing truncated normal-normal non-conjugacy detection with dynamic indexing", { 
    code = nimbleCode({
        for(i in 1:4) 
            y[i] ~ T(dnorm(mu[k[i]], sd = 1), 0, 1)
        for(j in 1:5) 
            mu[j] ~ dnorm(0, 1)
    })
    m = nimbleModel(code, data = list(y = rnorm(4)),
                    inits = list(k = rep(1,4)))
    conf <- configureMCMC(m)
    expect_equal(length(grep("conjugate", conf$getSamplers()[[1]]$name)), 0,
                 info = "incorrectly detected conjugacy")
})

test_that("Testing normal-truncated normal non-conjugacy detection with dynamic indexing", {     
    code = nimbleCode({
        for(i in 1:4) 
            y[i] ~ dnorm(mu[k[i]], sd = 1)
        for(j in 1:5) 
            mu[j] ~ T(dnorm(0, 1), 0 ,1)
    })
    m = nimbleModel(code, data = list(y = rnorm(4)),
                    inits = list(k = rep(1,4)))
    conf <- configureMCMC(m)
    expect_equal(length(grep("conjugate", conf$getSamplers()[[1]]$name)), 0,
                 info = "incorrectly detected conjugacy")
})

test_that("Testing normal-normal-IG multi-conjugacy detection with dynamic indexing", {     
    code = nimbleCode({
        for(i in 1:2) {
            y[i] ~ dnorm(cc*mu[xi[i]+offset], var = s0)
            xi[i] ~ dcat(p[1:3])
        }
        for(j in 1:4)
            mu[j] ~ dnorm(0,1)
        s0 ~ dinvgamma(1,1)
    })
    m <- nimbleModel(code, data = list(y = rep(1,2)), inits = list(xi = 1:2, offset = 1))
    conf <- configureMCMC(m)
    expect_match(conf$getSamplers()[[1]]$name, "conjugate_dnorm_dnorm",
                 info = "failed to detect normal-normal conjugacy")
    expect_match(conf$getSamplers()[[5]]$name, "conjugate_dinvgamma_dnorm",
                 info = "failed to detect invgamma-normal conjugacy")
})

test_that("Testing normal-normal-IG mixed-conjugacy detection with dynamic indexing", {     
    code = nimbleCode({
        for(i in 1:2) {
            y[i] ~ dnorm(s0 + mu[xi[i]], var = s0)
            xi[i] ~ dcat(p[1:3])
        }
        for(j in 1:4)
            mu[j] ~ dnorm(0,1)
        s0 ~ dinvgamma(1,1)
    })
    m <- nimbleModel(code, data = list(y = rep(1,2)), inits = list(xi = 1:2))
    conf <- configureMCMC(m)
    expect_match(conf$getSamplers()[[1]]$name, "conjugate_dnorm_dnorm",
                 info = "failed to detect normal-normal conjugacy when var param in mean and var")
    expect_equal(length(grep("conjugate", conf$getSamplers()[[5]]$name)), 0,
                 info = "incorrectly detected conjugacy when var param in dnorm mean and var")
})

test_that("Testing normal-normal-IG mixed-conjugacy detection (case 2) with dynamic indexing", {     
    code = nimbleCode({
        for(i in 1:2) {
            y[i] ~ dnorm(mu[round(xi[i]+s0)], var = s0)
            xi[i] ~ dcat(p[1:3])
        }
        for(j in 1:4)
            mu[j] ~ dnorm(0,1)
        s0 ~ dinvgamma(1,1)
    })
    m <- nimbleModel(code, data = list(y = rep(1,2)), inits = list(s0 = 1, xi = 1:2))
    conf <- configureMCMC(m)
    expect_match(conf$getSamplers()[[1]]$name, "conjugate_dnorm_dnorm",
                 info = "failed to detect normal-normal conjugacy when var param in mean index")
    expect_equal(length(grep("conjugate", conf$getSamplers()[[5]]$name)), 0,
                 info = "incorrectly detected conjugacy when var param in mean index")
})

test_that("Testing normal-normal-IG complicated dependency non-conjugacy detection with dynamic indexing", {     
    code = nimbleCode({
        for(i in 1:2) {
            y[i] ~ dnorm(sqrt(s0 + mu[xi[i]]), var = s0)
            xi[i] ~ dcat(p[1:3])
        }
        for(j in 1:4)
            mu[j] ~ dnorm(0,1)
        s0 ~ dinvgamma(1,1)
    })
    m <- nimbleModel(code, data = list(y = rep(1,2)), inits = list(xi = 1:2))
    conf <- configureMCMC(m)
    expect_equal(length(grep("conjugate", conf$getSamplers()[[1]]$name)), 0,
                 info = "incorrectly detected conjugacy wth nonlinear functional of mean param in dnorm")
    expect_equal(length(grep("conjugate", conf$getSamplers()[[5]]$name)), 0,
                 info = "incorrectly detected conjugacy when var param in dnorm mean and var, case 2")
})

test_that("Testing normal-normal-IG complicated dependency mixed-conjugacy detection with dynamic indexing", {    
    code = nimbleCode({
        for(i in 1:2) {
            y[i] ~ dnorm(mu[round(xi[i]+mu[i])], var = s0)
            xi[i] ~ dcat(p[1:3])
        }
        for(j in 1:4)
            mu[j] ~ dnorm(0,1)
        s0 ~ dinvgamma(1,1)
    })
    m <- nimbleModel(code, data = list(y = rep(1,2)), inits = list(mu = rep(1, 4), xi = 1:2))
    conf <- configureMCMC(m)
    expect_equal(length(grep("conjugate", conf$getSamplers()[[1]]$name)), 0,
                 info = "incorrectly detected conjugacy when mean param in own index")
    expect_match(conf$getSamplers()[[5]]$name, "conjugate_dinvgamma_dnorm",
                 info = "failed to detect invgamma-normal conjugacy when mean param in own index")
})

test_that("Testing multivariate normal-Wishart dependency conjugacy detection with dynamic indexing", {    
    code = nimbleCode({
        y[1:3] ~ dmnorm(mu[k, 1:3], pr[1:3,1:3])
        pr[1:3,1:3] ~ dwish(pr0[1:3,1:3], 5)
        for(i in 1:3)
            mu[i, 1:3] ~ dmnorm(z[1:3], pr0[1:3, 1:3])
        k ~ dcat(p[1:3])
    })
    m = nimbleModel(code, inits = list(k = 1, pr0 = diag(3), pr = diag(3),
                                       z = rep(1,3)))
    conf <- configureMCMC(m)
    expect_match(conf$getSamplers()[[2]]$name, "conjugate_dwish_dmnorm",
                 info = "failed to detect wishart-dmnorm conjugacy in multivariate multi-conjugacy setting")
    expect_match(conf$getSamplers()[[3]]$name, "conjugate_dmnorm_dmnorm",
                 info = "failed to detect dmnorm-dmnorm conjugacy in multivariate multi-conjugacy setting")
})

test_that("Expected failure: Testing multivariate normal-Wishart dependency conjugacy detection with dynamic indexing", {    
    code = nimbleCode({
        y[1:3] ~ dmnorm(mu[k, 1:3], pr[1:3,1:3])
        pr[1:3,1:3] ~ dwish(pr0[1:3,1:3], 5)
        for(i in 1:3)
            mu[1:3, i] ~ dmnorm(z[1:3], pr0[1:3, 1:3])
        k ~ dcat(p[1:3])
    })
    m = nimbleModel(code, inits = list(k = 1, pr0 = diag(3), pr = diag(3),
                                       z = rep(1,3)))
    conf <- configureMCMC(m)
    expect_failure(expect_match(conf$getSamplers()[[3]]$name, "conjugate_dmnorm_dmnorm",
            info = "failed to detect dmnorm-dmnorm conjugacy  in multivariate crossed multi-conjugacy setting"))
})

test_that("Testing normal-normal-IG multi-index multi-conjugacy detection with dynamic indexing", {    
    code = nimbleCode({
        for(i in 1:2) {
            k2[i] <- 1 + k1[i]
            k1[i] ~ dbern(p)
            for(j in 1:3)
                for(l in 1:2) {
                    y[i,j,l] ~ dnorm(mu[1+k1[i], 2, k2[i]], var = s0)
                    mu[i, j, l] ~ dnorm(0, 1)
                }
        }
        s0 ~ dinvgamma(1,1)
    })
    m = nimbleModel(code, data = list(y = array(rnorm(2*2*3), c(2,3,2))),
                    inits = list(k1 = rep(0,2)))
    conf <- configureMCMC(m)
    expect_match(conf$getSamplers()[[3]]$name, "conjugate_dnorm_dnorm",
                 info = "failed to detect normal-normal conjugacy")
    expect_match(conf$getSamplers()[[13]]$name, "conjugate_dinvgamma_dnorm",
                 info = "failed to detect invgamma-normal conjugacy")
})

test_that("Testing normal-normal-IG multi-index mixed-conjugacy detection with dynamic indexing", {    
    code = nimbleCode({
        for(i in 1:2) {
            k2[i] <- s0 + k1[i]
            k1[i] ~ dbern(p)
            for(j in 1:3)
                for(l in 1:2) {
                    y[i,j,l] ~ dnorm(mu[1+k1[i], 2, k2[i]], var = s0)
                    mu[i, j, l] ~ dnorm(0, 1)
                }
        }
        s0 ~ dinvgamma(1,1)
    })
    m = nimbleModel(code, data = list(y = array(rnorm(2*2*3), c(2,3,2))),
                    inits = list(s0 = 1, k1 = rep(0,2)))
    conf <- configureMCMC(m)
    expect_match(conf$getSamplers()[[3]]$name, "conjugate_dnorm_dnorm",
                 info = "failed to detect normal-normal conjugacy, var not conjugate")
    expect_equal(length(grep("conjugate", conf$getSamplers()[[13]]$name)), 0,
                 info = "incorrectly detected conjugacy for variance, var not conjugate") 

})

## Testing that uninitialized dynamic indexes are initialized at start of MCMC and MCMC runs

test_that("Testing initialization of uninitialized dynamic indexes", {
    mix_normal <- nimbleCode({
        for(i in 1:n) {
            z[i] ~ dcat(lambda[1:K])
            y[i] ~ dnorm(mu[z[i]],1)       
        }
        lambda[1:K] ~ ddirch(alpha[1:K])
        for(k in 1:K) {
            mu[k] ~ dnorm(0.0,0.0001);
        }   
    })
    
    n <- 200; K <- 2;
    y1 <- c(rnorm(100, mean = 3, sd = 1), rnorm(100, mean = -3, sd = 1))
    consts <- list(n = n, K = K)
    data <- list(y = y1, alpha = rep(1,K))

    if(nimbleOptions('verbose')) {
        expect_message(model1 <- nimbleModel(code = mix_normal, constants = consts, data = data), "missing value where TRUE/FALSE needed")
    } else expect_silent(model1 <- nimbleModel(code = mix_normal, constants = consts, data = data))
    
    conf1 <- configureMCMC(model1)
    Rmcmc <- buildMCMC(model1)
    Cmodel <- compileNimble(model1)
    Cmcmc <- compileNimble(Rmcmc, project = model1)
    expect_equal(Cmcmc$run(100), NULL)
})
    
## MCMC testing

test_that("MCMC with invalid indexes produce warning, but runs", {
    code <- nimbleCode({
        y ~ dnorm(x[k-1], 1)
        k ~ dcat(p[1:3])
    })
    m <- nimbleModel(code, data = list(y = 0), inits = list(k = 2, x = c(0,1), p = rep(1/3,3)))
    cm <- compileNimble(m)
    mcmc = buildMCMC(m)
    cmcmc = compileNimble(mcmc ,project=m)
    set.seed(1)
    expect_output(cmcmc$run(10000), "dynamic index out of bounds")
    out <- as.matrix(cmcmc$mvSamples)
    expect_equal(sum(out == 2) / sum(out == 3), dnorm(0)/dnorm(1), tolerance = .02)
})


models <- c('hearts')

### test BUGS examples - models and MCMC
sapply(models, testBUGSmodel, useInits = TRUE)
sapply(models, test_mcmc, numItsC = 1000)

## beta0C is beta0 in inits file
system.in.dir(paste("sed 's/beta0/beta0C/g' cervix-inits.R > ", file.path(tempdir(), "cervix-inits.R")), dir = system.file('classic-bugs','vol2','cervix', package = 'nimble'))
## need missing x's to have initial values or initial model$calculate fails
system.in.dir(paste("echo 'x <- c(rep(as.numeric(NA), 115), rep(1, 1929))' >> ", file.path(tempdir(), "cervix-inits.R")), dir = system.file('classic-bugs','vol2','cervix', package = 'nimble'))

test_that('cervix model and MCMC test', {
    testBUGSmodel('cervix', dir = "", model = system.file('classic-bugs','vol2','cervix','cervix.bug', package = 'nimble'), data = system.file('classic-bugs','vol2','cervix','cervix-data.R', package = 'nimble'),  inits = file.path(tempdir(), "cervix-inits.R"),  useInits = TRUE)
    test_mcmc(model = system.file('classic-bugs','vol2','cervix','cervix.bug', package = 'nimble'), name = 'cervix', inits = file.path(tempdir(), "cervix-inits.R"), data = system.file('classic-bugs', 'vol2', 'cervix','cervix-data.R', package = 'nimble'), numItsC = 1000, avoidNestedTest = TRUE)
})

## There are some issues with using the model as provided in the BUGS example because of use of zeros
## to exclude categories in dirichlet-distributed vector, therefore recreate BUGS code here.
test_that('biops model and MCMC test', {
    system.in.dir(paste("echo 'var\n
    biopsies[ns,4], #  grades observed in ith session (multinomial)\n
    nbiops[ns],     # total number of biopsies in ith session\n
    truex[ns],       # true state in ith session\n
    error[4,4],     # error matrix in taking biopsies\n
    prior[4,4],     # prior parameters for rows of error[,]\n
    p[4];           # underlying   incidence of true  states\n
model {\n
   for (i in 1:ns){\n
      truex[i]       ~ dcat(p[]);\n
      biopsies[i,]  ~ dmulti(error[truex[i],],nbiops[i]); \n
   }\n
   error[1, 1:4] <- c(1, 0, 0, 0)\n
    error[2, 1:2] ~ ddirch(prior[2, 1:2])\n
    error[3, 1:3] ~ ddirch(prior[3, 1:3])\n
    error[4, 1:4] ~ ddirch(prior[4, 1:4])\n
   p[]       ~ ddirch(prior[4,]);     # prior for p\n
   }' >> ", file.path(tempdir(), "biops.bug")), dir = system.file('classic-bugs','vol2','biops', package = 'nimble'))
    system.in.dir(paste("sed 's/true/truex/g' biops-inits.R > ", file.path(tempdir(), "biops-inits.R")), dir = system.file('classic-bugs','vol2','biops', package = 'nimble'))
system.in.dir(paste("echo 'error <- matrix(c(1,0,0,0, .5, .5, 0, 0, 1/3,1/3,1/3,0,1/4,1/4,1/4,1/4), 4,4, byrow=T)'  >> ", file.path(tempdir(), "biops-inits.R")), dir = system.file('classic-bugs','vol2','cervix', package = 'nimble'))
    testBUGSmodel('biops', dir = "", model = file.path(tempdir(), "biops.bug"), data = system.file('classic-bugs','vol2','biops','biops-data.R', package = 'nimble'),  inits = file.path(tempdir(), "biops-inits.R"),  useInits = TRUE)
    test_mcmc(model = file.path(tempdir(), "biops.bug"), name = 'biops', inits = file.path(tempdir(), "biops-inits.R"), data = system.file('classic-bugs', 'vol2', 'biops','biops-data.R', package = 'nimble'), numItsC = 1000, avoidNestedTest = TRUE)
})


test_that('basic mixture model with conjugacy', {
    n <- 1000
    d <- 4
    set.seed(2)
    mns <- c(-.9, .2, 1.6, -1.1)
    p <- c(.4, .14, .1, .36)
    sds <- c(.2, .4, .2, .1)
    k <- sample(1:d, n, replace = TRUE, prob = p)
    y <- rnorm(n, mns[k], sds[k])
    code <- nimbleCode({
        for(i in 1:n) {
            y[i] ~ dnorm(mu[k[i]], sd = sigma[k[i]])
            k[i] ~ dcat(p[1:d])
        }
        for(j in 1:d) {
            mu[j] ~ dnorm(mu0, sd = tau)
            sigma[j] ~ dgamma(a, b)
        }
        mu0 ~ dnorm(0, sd = 100)
        tau ~ dhalfflat()
        a ~ dhalfflat()
        b ~ dhalfflat()
        p[1:d] ~ ddirch(alpha[1:d])
    })
        
    test_mcmc(name = 'basic mixture model with conjugacy',
              model = code, seed = 1, numItsC_results = 20000,
              data = list(d=d,n=n,y=y),
              inits = list(sigma = rep(1,4), mu = c(-1, 0, 1, 2),
                           p = rep(.25, 4), k = sample(1:4, n, replace = TRUE),
                           mu0 = 0, a = 1, b = 1, tau = 1, alpha = rep(1, 4)),
              results = list(mean = list("mu[1]" = mns[4],
                                         "mu[2]" = mns[3],
                                         "mu[3]" = mns[1],
                                         "mu[4]" = mns[2],
                                         "p[1]" = p[4],
                                         "p[2]" = p[3],
                                         "p[3]" = p[1],
                                         "p[4]" = p[2])),
              resultsTolerance = list(mean = list("mu[1]" = .05,
                                                  "mu[2]" = .02,
                                                  "mu[3]" = .05,
                                                  "mu[4]" = .05,
                                                  "p[1]" = .05,
                                                  "p[2]" = .02,
                                                  "p[3]" = .1,
                                                  "p[4]" = .02)),
              avoidNestedTest = TRUE)

})

test_that('basic mixture model without conjugacy', {
    n <- 1000; d <- 4
    set.seed(2)
    mns <- c(8, 15, 0.5, 4)
    p <- c(.45, .14, .05, .36)
    k <- sample(1:d, n, replace = TRUE, prob = p)
    y <- rpois(n, mns[k])
    code <- nimbleCode({
        for(i in 1:n) {
            y[i] ~ dpois(mu[k[i]])
            k[i] ~ dcat(p[1:d])
        }
        for(j in 1:d) {
            mu[j] ~ dnorm(mu0, sd = tau)
        }
        mu0 ~ dnorm(0, sd = 100)
        tau ~ dhalfflat()
        p[1:d] ~ ddirch(alpha[1:d])
    })
    test_mcmc(name = 'basic mixture model without conjugacy',
              model = code, seed = 1, numItsC_results = 20000,
              data = list(d=d,n=n,y=y),
              inits = list(mu = rep(3, 4),
                           p = rep(.25, 4), k = sample(1:4, n, replace = TRUE),
                           mu0 = 4, tau = 1, alpha = rep(1, 4)),
              results = list(mean = list("mu[1]" = mns[3],
                                         "mu[2]" = mns[4],
                                         "mu[3]" = mns[1],
                                         "mu[4]" = mns[2],
                                         "p[1]" = p[3],
                                         "p[2]" = p[4],
                                         "p[3]" = p[1],
                                         "p[4]" = p[2])),
              resultsTolerance = list(mean = list("mu[1]" = .5,
                                                  "mu[2]" = .8,
                                                  "mu[3]" = 1.5,
                                                  "mu[4]" = 1,
                                                  "p[1]" = .03,
                                                  "p[2]" = .12,
                                                  "p[3]" = .12,
                                                  "p[4]" = .05)),
              avoidNestedTest=TRUE)
})


test_that('range checking with dynamic indexing', {
    code <- nimbleCode({
        for(i in 3:4) {
            z[i-1] ~ dcat(alpha[z[i-2], 1:2])
        }
    }
    )
    
    m <- nimbleModel(code, inits = list(alpha = matrix(c(0.9,.1,0.1,.9), 2),
                                        z = c(1,2,1)), calculate = FALSE)
    ## will give error if issue #790 is not fixed:
    expect_silent(output <- m$calculate())

    code <- nimbleCode({
        for(i in 1:2) {
            z[i+1] ~ dcat(alpha[z[i], 1:2])
        }
    }
    )
    
    m <- nimbleModel(code, inits = list(alpha = matrix(c(0.9,.1,0.1,.9), 2),
                                        z = c(1,2,1)), calculate = FALSE)
    expect_silent(output <- m$calculate())
    
    code <- nimbleCode({
        for(i in 2:3) {
            z[i] ~ dcat(alpha[z[i-1], 1:2])
        }
    }
    )
    
    
    m <- nimbleModel(code, inits = list(alpha = matrix(c(0.9,.1,0.1,.9), 2),
                                        z = c(1,2,1)), calculate = FALSE)
    expect_silent(output <- m$calculate())
    
})


if(FALSE) {
    ## Heisenbug here - running manually vs. via testthat causes different ordering of mixture components even though we are setting seed before each use of RNG; this test _does_ pass when run manually.
test_that('basic multivariate mixture model with conjugacy', {
    n <- 1000; d <- 4
    set.seed(2)
    mns <- cbind(c(1.5, 1.5, -1.5, -1.5), c(1.5, -1.5, 1.5, -1.5))
    p <- c(.45, .14, .05, .36)
    k <- sample(1:d, n, replace = TRUE, prob = p)
    y <- cbind(rnorm(n), rnorm(n))
    y <- y + mns[k, ]
    code <- nimbleCode({
        for(i in 1:n) {
            y[i, 1:2] ~ dmnorm(mu[k[i], 1:2], pr[1:2, 1:2])
            k[i] ~ dcat(p[1:d])
        }
        for(i in 1:d) {
                mu[i, 1:2] ~ dmnorm(z[1:2], pr0[1:2, 1:2])
        }
        ## if put a positive prior on diag elements directly can get negative samples and
        ## failure in R MCMC because of Cholesky
        myvars[1] ~ dunif(-2, 4)
        myvars[2] ~ dunif(-2, 4)
        pr[1, 1] <- exp(myvars[1])
        pr[2, 2] <- exp(myvars[2])
        pr[1, 2] <- 0
        pr[2, 1] <- 0
        p[1:d] ~ ddirch(alpha[1:d])
    })
    test_mcmc(name = 'basic multivariate mixture model with conjugacy',
              model = code, seed = 1, numItsC_results = 20000,
              data = list(y = y, z = rep(0, 2), pr0 = diag(rep(1e-4, 2)), n = n, d = d),
              inits = list(mu = cbind(rnorm(4), rnorm(4)), k = sample(1:4, n, replace = TRUE),
                           p = rep(.25, 4), alpha = rep(1,4), myvars = rep(1,2)),
              results = list(mean = list("mu[1, 1]" = mns[2,1], "mu[1, 2]" = mns[2,2], 
                                         "mu[2, 1]" = mns[3,1], "mu[2, 2]" = mns[3,2],
                                         "mu[3, 1]" = mns[1,1], "mu[3, 2]" = mns[1,2],
                                         "mu[4, 1]" = mns[4,1], "mu[4, 2]" = mns[4,2],
                                         "p[1]" = p[2],
                                         "p[2]" = p[3],
                                         "p[3]" = p[1],
                                         "p[4]" = p[4],
                                         "myvars[1]" = 0,
                                         "myvars[2]" = 0)),
              resultsTolerance = list(mean = list("mu[1, 1]" = .2, "mu[1, 2]" = .1, 
                                         "mu[2, 1]" = .6, "mu[2, 2]" = .05,
                                         "mu[3, 1]" = .05, "mu[3, 2]" = .1,
                                         "mu[4, 1]" = .05, "mu[4, 2]" = .2,
                                         "p[1]" = .04,
                                         "p[2]" = .01,
                                         "p[3]" = .01,
                                         "p[4]" = .03,
                                         "myvars[1]" = .1,
                                         "myvars[2]" = .1)), avoidNestedTest = TRUE)
})
}



sink(NULL)

if(!generatingGoldFile) {
    trialResults <- readLines(tempFileName)
    correctResults <- readLines(system.file(file.path('tests', goldFileName), package = 'nimble'))
    compareFilesByLine(trialResults, correctResults)
}

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)
nimbleOptions(allowDynamicIndexing = nimbleAllowDynamicIndexingSetting)

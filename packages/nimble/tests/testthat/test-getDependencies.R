source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context("Testing of getDependencies")

## note this testing is not intended to blur into general model processing testing.
## It assumes the model processing is ok and really tests traversal of the graph to get dependencies

## also self=FALSE is not a deep testing need because this is processed in R after the deeper processing (graph traversal) in C++

## first model, with no criss-crossing dependencies.
test_that("getDependencies in first model with no criss-crossing dependencies", {
    c1 <- nimbleCode({
        a1 ~ dnorm(0,1)
        a2 ~ dnorm(0,1)
        b1 ~ dnorm(a1, 1) ## only a1 
        b2 ~ dnorm(a1, sd = a2) ## a1 and a2
        c1 <- a1 + 1 
        d1 ~ dnorm(c1, 1) ## only a1 via c1
        c2 <- a1 + 2
        c3 <- c2^2
        e1 ~ dnorm(c3, 1) ## only a1 via c2 and d1
        c4 <- a2 + 3
        f1 ~ dnorm(c3, sd = c4) ## a1 via c1 and d1; a2 via c3
        g1 ~ dnorm(f1, 1)
    })
    
    m1 <- nimbleModel(c1)
    ans1 <- m1$topologicallySortNodes(c('a1','b1','b2','c1','d1','c2','c3','e1','f1'))
    
    ## basic cases from a1
    expect_identical(m1$getDependencies('a1'), ans1)
    expect_identical(m1$getDependencies(c('a1','a1')), ans1)
    expect_identical(m1$getDependencies(c('a1','c1')), ans1)
    expect_identical(m1$getDependencies(c('c1','a1')), ans1)
    expect_identical(m1$getDependencies(c('a1', 'f1')), c(ans1, 'g1'))
    expect_identical(m1$getDependencies(c('f1', 'a1')), c(ans1, 'g1'))
    expect_identical(m1$getDependencies(c('a1'), downstream=TRUE), c(ans1, 'g1'))
    
    ## basic cases from a1 and a2
    ans2 <- m1$topologicallySortNodes(c('a1','a2','b1','b2','c4','c1','d1','c2','c3','e1','f1'))
    expect_identical(m1$getDependencies(c('a1','a2')), ans2)
    expect_identical(m1$getDependencies(c('a2','a1')), ans2)
    
    ## some omit cases
    ans3 <- m1$topologicallySortNodes(c('a1','b1','b2','c1','d1','c2'))
    expect_identical(m1$getDependencies('a1', omit = 'c3'), ans3)
    
    ans4 <- m1$topologicallySortNodes(c('a1','a2','b1','b2','c1','d1','c2','c3','e1','f1'))
    expect_identical(m1$getDependencies(c('a2','a1'), omit = 'c4', downstream = TRUE), c(ans4, 'g1'))
    
    ## some cases starting from deterministic
    expect_identical(m1$getDependencies('c1'), c('c1','d1'))
    expect_identical(m1$getDependencies('c2'), c('c2','c3','e1','f1'))
    expect_identical(m1$getDependencies(c('c1','c2')), c('c1','c2','d1','c3','e1','f1'))
})
    
################
    ## second model, with some criss-crossing dependencies
test_that("getDependencies in second model with some criss-crossing dependencies", {
    c2 <- nimbleCode({
        a1 <- 1
        b1 ~ dnorm(a1, 1)
        c1 <- b1 + 1
        c2 <- b1 + 2
        d1 ~ dnorm(c1, sd = c2)
        e1 ~ dnorm(d1, 1)
        f1 ~ dnorm(d1, sd = c1)
        g1 ~ dnorm(e1, sd = c1)
        c3 <- c1 + c2
        h1 ~ dnorm(c3, 1)
        h2 ~ dnorm(c3, sd = c1)
    })

    m2 <- nimbleModel(c2)

    ans1 <- m2$topologicallySortNodes(c('a1','b1'))
    expect_identical(m2$getDependencies('a1'), ans1)

    ans2 <- m2$topologicallySortNodes(c('b1','c1','c2','d1','f1','g1','c3','h1','h2'))
    expect_identical(m2$getDependencies('b1'), ans2)
    expect_identical(m2$getDependencies(c('b1','c1')), ans2)
    expect_identical(m2$getDependencies(c('c1','b1')), ans2)
    expect_identical(m2$getDependencies(c('b1','c1','f1')), ans2)
    expect_identical(m2$getDependencies(c('f1','c1','b1')), ans2)
    expect_identical(m2$getDependencies(c('f1','b1','c1')), ans2)

    ans3 <- m2$topologicallySortNodes(c('b1','c1','c2','d1','e1','f1','g1','c3','h1','h2'))
    expect_identical(m2$getDependencies(c('b1','c1','d1')), ans3)
    expect_identical(m2$getDependencies(rev(ans2)), ans3)
    expect_identical(m2$getDependencies(c('b1','c1','e1')), ans3)
    expect_identical(m2$getDependencies(c('e1','c1','b1')), ans3)

    ans4 <- m2$topologicallySortNodes(c('c2','d1','c3','h1','h2'))
    expect_identical(m2$getDependencies(c('c2')), ans4)
    expect_identical(m2$getDependencies(c('c3','c2')), ans4)
    expect_identical(m2$getDependencies(c('c2','h1','c3')), ans4)
})

test_that("getParentNodes works", {
    code=nimbleCode({
        for(j in 1:J) {
            for(i in 1:I)
                y[j,i] ~ dnorm(theta[j], sigma)
            theta[j] ~ dnorm(mu, sd = tau)
        }
        mu ~ dnorm(mu0,1)
        mu0 <- mu00
        sigma ~ dunif(sigma0,1)
        sigma1 ~ dunif(0,1)
        tau ~ dunif(0,1)
    })
    I <- 4
    J <- 3
    constants <- list(I = I, J = J)
    m <- nimbleModel(code,
                     data = list(y = matrix(rnorm(I*J), J, I)),
                     constants = constants)
    expect_identical(getParentNodes('mu', m, stochOnly =  TRUE), character(0))
    expect_identical(getParentNodes('mu', m), c('mu0', 'mu00'))
    expect_identical(getParentNodes('y', m, stochOnly = TRUE),
                     c('theta[1]','theta[2]','theta[3]', 'sigma'))
    expect_identical(getParentNodes('y', m),
                     c('lifted_d1_over_sqrt_oPsigma_cP', 'theta[1]','theta[2]','theta[3]', 'sigma'))
    expect_identical(getParentNodes('y[2, 1:3]', m, stochOnly = TRUE),
                     c('theta[2]', 'sigma'))
    expect_identical(getParentNodes('theta', m),
                     c('tau', 'mu'))
    
})


options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)

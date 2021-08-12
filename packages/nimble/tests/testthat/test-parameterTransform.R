source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context('Testing of parameterTransform nimbleFunction')


##
## parameterTransform testing code below
##

test_that('exact values from uncompiled and compiled parameterTransform', {
    code <- nimbleCode({
        a ~ dnorm(0, 1)
        b ~ dgamma(1, 1)
        c ~ dunif(2, 10)
        d[1:3] ~ dmnorm(mu[1:3], cov = C[1:3,1:3])
        e[1:3,1:3] ~ dwish(R = C[1:3,1:3], df = 5)
    })
    constants <- list(mu=rep(0,3), C=diag(3))
    data <- list()
    U <- matrix(c(.2,2,4,0,1,1,0,0,4), nrow=3, byrow=TRUE)
    eInit <- t(U) %*% U
    inits <- list(a=0, b=1, c=5, d=rep(0,3), e=eInit)
    Rmodel <- nimbleModel(code, constants, data, inits)
    ##
    expect_equal(Rmodel$calculate(), -33.077938542)
    ##
    nodes <- letters[1:5]
    pt <- parameterTransform(Rmodel, nodes)
    ## 
    Cmodel <- compileNimble(Rmodel)
    Cpt <- compileNimble(pt, project = Rmodel)
    ##
    vals <- values(Rmodel, nodes)
    ##
    ## test uncompiled
    theNF <- pt
    tVals <- theNF$transform(vals)
    vals2 <- theNF$inverseTransform(tVals)
    ##
    expect_true(theNF$getOriginalLength() == 15)
    expect_true(theNF$getTransformedLength() == 12)
    expect_true(all(round(vals - c(0.00, 1.00, 5.00, 0.00, 0.00, 0.00, 0.04, 0.40, 0.80, 0.40, 5.00, 9.00, 0.80, 9.00, 33.00), 16) == 0))
    expect_true(all(round(tVals - c(0, 0, -0.510825623765991, 0, 0, 0, -1.6094379124341, 2, 0, 4, 1, 1.38629436111989), 14) == 0))
    expect_true(all(vals - vals2 == 0))
    expect_true(round(theNF$logDetJacobian(tVals), 7) == -0.9571127)
    ##
    ## test compiled
    theNF <- Cpt
    tVals <- theNF$transform(vals)
    vals2 <- theNF$inverseTransform(tVals)
    ##
    expect_true(theNF$getOriginalLength() == 15)
    expect_true(theNF$getTransformedLength() == 12)
    expect_true(all(round(vals - c(0.00, 1.00, 5.00, 0.00, 0.00, 0.00, 0.04, 0.40, 0.80, 0.40, 5.00, 9.00, 0.80, 9.00, 33.00), 16) == 0))
    expect_true(all(round(tVals - c(0, 0, -0.510825623765991, 0, 0, 0, -1.6094379124341, 2, 0, 4, 1, 1.38629436111989), 14) == 0))
    expect_true(all(vals - vals2 == 0))
    expect_true(round(theNF$logDetJacobian(tVals), 7) == -0.9571127)
})


test_that('exact values parameterTransform of multivariate dmnorm and dwish nodes', {
    code <- nimbleCode({
        a ~ dnorm(0, 1)
        b ~ dgamma(1, 1)
        c ~ dunif(2, 10)
        d[1:3] ~ dmnorm(mu[1:3], cov = C[1:3,1:3])
        e[1:3,1:3] ~ dwish(R = C[1:3,1:3], df = 5)
        f ~ dunif(0, 1)
        g ~ dunif(0, 5)
        h ~ dt(2, 2, 4)
        ii[1:5,1:5] ~ dwish(R = Ci[1:5,1:5], df = 10)
        j ~ dunif(-5, 5)
    })
    ##
    Ci <- diag(5)
    Ci[1,2] <- Ci[2,1] <- 0.2
    Ci[1,3] <- Ci[3,1] <- 0.1
    Ci[4,5] <- Ci[5,4] <- 0.3
    ##
    constants <- list(mu=rep(0,3), C=diag(3), Ci=Ci)
    data <- list()
    U <- matrix(c(.4,3,5,0,1,-1,0,0,4), nrow=3, byrow=TRUE)
    eInit <- t(U) %*% U
    inits <- list(a=0, b=1, c=5, d=rep(0,3), e=eInit, f=0.5, g=4, h=1, ii=diag(5), j=-1)
    ##
    Rmodel <- nimbleModel(code, constants, data, inits)
    ##
    expect_equal(Rmodel$calculate(), -80.6027522654)
    ##
    nodes <- Rmodel$getNodeNames(stochOnly = TRUE)
    pt <- parameterTransform(Rmodel, nodes)
    ## 
    Cmodel <- compileNimble(Rmodel)
    Cpt <- compileNimble(pt, project = Rmodel)
    ##
    vals <- values(Rmodel, nodes)
    ##
    ## test uncompiled
    theNF <- pt
    tVals <- theNF$transform(vals)
    vals2 <- theNF$inverseTransform(tVals)
    ##
    expect_true(theNF$getOriginalLength() == 44)
    expect_true(theNF$getTransformedLength() == 31)
    expect_true(all(round(vals,15) - c(0, 1, 5, 0.5, 4, 1, -1, 0, 0, 0, 0.16, 1.2, 2, 1.2, 10, 14, 2, 14, 42, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1) == 0))
    expect_true(all(round(tVals - c(0, 0, -0.510825623765991, 0, 1.38629436111989, 1, -0.405465108108164, 0, 0, 0, -0.916290731874155, 3, -0.00000000000000177635683940025, 5, -1, 1.38629436111989, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 14) == 0))
    expect_true(all(vals - vals2 == 0))
    expect_true(round(theNF$logDetJacobian(tVals) - 4.54724272356, 11) == 0)
    ##
    ## test compiled
    theNF <- Cpt
    tVals <- theNF$transform(vals)
    vals2 <- theNF$inverseTransform(tVals)
    ##
    expect_true(theNF$getOriginalLength() == 44)
    expect_true(theNF$getTransformedLength() == 31)
    expect_true(all(round(vals,15) - c(0, 1, 5, 0.5, 4, 1, -1, 0, 0, 0, 0.16, 1.2, 2, 1.2, 10, 14, 2, 14, 42, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1) == 0))
    expect_true(all(round(tVals - c(0, 0, -0.510825623765991, 0, 1.38629436111989, 1, -0.405465108108164, 0, 0, 0, -0.916290731874155, 3, -0.00000000000000177635683940025, 5, -1, 1.38629436111989, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 14) == 0))
    expect_true(all(vals - vals2 == 0))
    expect_true(round(theNF$logDetJacobian(tVals) - 4.54724272356, 11) == 0)
})



test_that('parameterTransform for dirichlet distribution', {
    code <- nimbleCode({
        x2[1:2] ~ ddirich(a2[1:2])
        x3[1:3] ~ ddirich(a3[1:3])
        x5[1:5] ~ ddirich(a5[1:5])
    })
    constants <- list(a2 = rep(1, 2),
                      a3 = rep(1, 3),
                      a5 = rep(1, 5))
    data <- list()
    inits <- list(x2 = c(.5, .5),
                  x3 = c(.2, .3, .5),
                  x5 = c(.1, .2, .3, .05, .35))
    ##
    Rmodel <- nimbleModel(code, constants, data, inits)
    ##
    expect_equal(Rmodel$calculate(), 3.871201)
    ##
    nodes <- Rmodel$getNodeNames(stochOnly = TRUE)
    pt <- parameterTransform(Rmodel, nodes)
    ##
    Cmodel <- compileNimble(Rmodel)
    Cpt <- compileNimble(pt, project = Rmodel)
    ##
    vals <- values(Rmodel, nodes)
    ##
    ## test uncompiled
    theNF <- pt
    tVals <- theNF$transform(vals)
    vals2 <- theNF$inverseTransform(tVals)
    ##
    expect_true(theNF$getOriginalLength() == 10)
    expect_true(theNF$getTransformedLength() == 7)
    expect_true(all(round(vals,15) - c(.5, .5, .2, .3, .5, .1, .2, .3, .05, .35) == 0))
    expect_true(all(round(tVals - c(0.0000000, -1.3862944, -0.5108256, -2.1972246, -1.2527630, -0.2876821, -1.9459101)) == 0))
    expect_true(all(round(vals - vals2, 16) == 0))
    expect_true(theNF$logDetJacobian(tVals) + 14.0544024662466205 == 0)
    ##
    ## test compiled
    theNF <- Cpt
    tVals <- theNF$transform(vals)
    vals2 <- theNF$inverseTransform(tVals)
    ##
    expect_true(theNF$getOriginalLength() == 10)
    expect_true(theNF$getTransformedLength() == 7)
    expect_true(all(round(vals,15) - c(.5, .5, .2, .3, .5, .1, .2, .3, .05, .35) == 0))
    expect_true(all(round(tVals - c(0.0000000, -1.3862944, -0.5108256, -2.1972246, -1.2527630, -0.2876821, -1.9459101)) == 0))
    expect_true(all(round(vals - vals2, 16) == 0))
    expect_true(theNF$logDetJacobian(tVals) + 14.0544024662466205 == 0)
})


test_that('parameterTransform for lkj correlation matrices', {
    ## Correlation matrix
    R <- matrix(c(
        1, 0.9, .3, -.5, .1,
        0.9, 1, .15, -.3, .1,
        .3, .15, 1, .3, .1,
        -.5, -.3, .3, 1, .1,
        .1,.1,.1,.1, 1)
      , 5, 5)

    U <- chol(R)

    ## Partial sums of columns of U
    PS <- matrix(0, 5, 5)
    PS[1,] <- 1
    PS[2, 2:5] <- 1-U[1, 2:5]^2
    PS[3, 3] <- 1 - sum(U[1:2, 3]^2)
    PS[3, 4] <- 1-U[1,4]^2-U[2,4]^2
    PS[4,4] <- 1 - sum(U[1:3, 4]^2)
    PS[3, 5] <- 1-U[1,5]^2-U[2,5]^2
    PS[4, 5]<- 1-U[1,5]^2-U[2,5]^2-U[3,5]^2
    PS[5, 5]<- 1 - sum(U[1:4, 5]^2)
    ## Canonical partial correlations
    Z <- diag(5)
    Z[1,2:5] <- U[1, 2:5]
    Z[2,3] <- U[2,3]/sqrt(PS[2,3])
    Z[2,4] <- U[2,4]/sqrt(PS[2,4])
    Z[3,4] <- U[3,4]/sqrt(PS[3,4])
    Z[2,5] <- U[2,5]/sqrt(PS[2,5])
    Z[3,5] <- U[3,5]/sqrt(PS[3,5])
    Z[4,5] <- U[4,5]/sqrt(PS[4,5])

    ## Unconstrained transformed parameters
    yt <- atanh(Z)
    diag(yt) <- 0
    yt <- yt[yt!=0]

    ## log of determinant of Jacobian per Stan Ref Manual Section 10.12 of version 2.27
    logDetJac <- 0.5*(log(PS[2,3])+log(PS[2,4])+log(PS[3,4])+log(PS[2,5]) + log(PS[3,5]) + log(PS[4,5])) -
        2*sum(log(cosh(yt)))

    set.seed(1)
    n <- 100
    p <- 5
    y <- t(t(U)%*%matrix(rnorm(p*n),p,n))
    
    code <- nimbleCode({
        for(i in 1:n)
            y[i, 1:p] ~ dmnorm(mu[1:p], cholesky = U[1:p, 1:p], prec_param = 0)
        U[1:p,1:p] ~ dlkj_corr_cholesky(eta = 1.3, p)
    })

    ## Can't use dlkj until lkj merged into devel and then into ADoak
    ## so expect failure until then.
    expect_failure({
        m <- nimbleModel(code, constants = list(n = n, p = p, mu = rep(0, p)),
                     data = list(y = y), inits = list(U = U))

        pt <- parameterTransform(m, nodes = 'Ustar')
        yt_from_pt <- pt$transform(m$Ustar)
        U_from_pt  <- pt$inverseTransform(yt)
        logDetJac_from_pt  <- pt$logDetJacobian(yt)
        
        cm <- compileNimble(m)
        cpt <- compileNimble(pt, project = m)
        cyt_from_pt <- cpt$transform(m$Ustar)
        cU_from_pt  <- cpt$inverseTransform(yt)
        clogDetJac_from_pt  <- cpt$logDetJacobian(yt)
        
        expect_identical(yt, yt_from_pt)
        expect_identical(U, matrix(U_from_pt, p, p))
        expect_identical(logDetJac, logDetJac_from_pt)
        
        expect_identical(yt, cyt_from_pt)
        expect_identical(U, matrix(cU_from_pt, p, p))
        expect_identical(logDetJac, clogDetJac_from_pt)
    })
})
    
options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)



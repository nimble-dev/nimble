source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

context("Testing of CAR distributions")

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)
nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)




## testing dcar_normal sampling
test_that('dcar_normal sampling', {
    cat('===== Starting MCMC test dcar_normal sampling. =====')
    
    code <- nimbleCode({
        alpha0 ~ dflat()
        S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], 3)
        for(i in 1:N)
            mu[i] <- alpha0 + S[i]
        for(i in 1:2) {
            log(lambda[i]) <- mu[i]
            Y[i] ~ dpois(lambda[i])
        }
        Y[3] ~ dnorm(mu[3], 3)
        ymean4 <- 5*mu[4]
        Y[4] ~ dnorm(ymean4, 7)
        ymean5 <- 2*mu[5]
        Y[5] ~ dnorm(ymean5, 1)
    })
    
    constants <- list(N = 6,
                      num = c(3,4,4,3,2,2),
                      adj = c(2,3,4,   1,3,5,6,   1,2,4,5,   1,3,6,  2,3,   2,4),
                      weights = rep(1, 18),
                      L = 18)
    data <- list(Y = c(10,12,15,20,24))
    inits <- list(alpha0 = 0,
                  S = c(0,0,0,0,0,0))
    
    Rmodel <- nimbleModel(code, constants, data, inits)
    conf <- configureMCMC(Rmodel)
    Rmcmc <- buildMCMC(conf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    
    niter <- 20
    set.seed(0); Rmcmc$run(niter)
    set.seed(0); Cmcmc$run(niter)
    
    Rsamples <- as.matrix(Rmcmc$mvSamples)
    Csamples <- as.matrix(Cmcmc$mvSamples)
    
    sampleNames <- colnames(Rsamples)
    
    expect_true(all(Rsamples[, sampleNames] - Csamples[, sampleNames] == 0),
                info = 'agreement between R and C sampling of dcar_normal')
    
    expect_lt(max(abs(as.numeric(Csamples[20, sampleNames]) - 
                 c(1.639127, 1.815422, 1.676655, 5.099797, 2.345276, 7.018026, 2.696936))),  1e-6,
                 label = 'exact sample values for dcar_normal')
})



## testing dcar_proper density evaluation,
## and generation of default values of C and M
test_that('dcar_proper density evaluation', {
    cat('===== Starting test dcar_proper density evaluations. =====')

    x <- c(1, 3, 3, 4)
    mu <- rep(3, 4)
    adj <- c(2, 1,3, 2,4, 3)
    num <- c(1, 2, 2, 1)
    lp <- 0.004158475
    expect_equal(dcar_proper(x, mu, adj=adj, num=num, tau=1, gamma=0), lp,
                 info = 'C density evaluation for dcar_proper() omitting C and M')

    weights <- rep(1, 6)
    CM <- as.carCM(adj, weights, num)
    C <- CM$C
    M <- CM$M
    expect_equal(dcar_proper(x, mu, C, adj, num, M, tau=1, gamma=0), lp,
                 info = 'C density evaluation for dcar_proper() weights all one')

    weights <- c(2, 2, 3, 3, 4, 4)
    CM2 <- as.carCM(adj, weights, num)
    C2 <- CM2$C
    M2 <- CM2$M
    lp2 <- 0.001050636
    expect_equal(dcar_proper(x, mu, C2, adj, num, M2, tau=1, gamma=0), lp2,
                 info = 'C density evaluation for dcar_proper() weights different')
})



## testing dcar_proper sampling
test_that('dcar_proper sampling', {
    cat('===== Starting MCMC test dcar_proper sampling. =====')

    code <- nimbleCode({
        tau ~ dgamma(0.001, 0.001)
        gamma ~ dunif(-1, 1)
        x[1:N] ~ dcar_proper(mu[1:N], adj=adj[1:L], num=num[1:N], tau=tau, gamma=gamma)
        y[1] ~ dnorm(x[1], 1)
        y[2] ~ dnorm(3*x[2] + 5, 10)
        y[3] ~ dnorm(x[3]^2, 1)
        y[4] ~ dnorm(x[4]^2, 10)
    })

    mu <- 1:4
    adj <- c(2, 1, 3, 2, 4, 3)
    num <- c(1, 2, 2, 1)
    tau <- 1
    gamma <- 0
    y <- c(3, 6, 8, 10)
    x <- rep(0, 4)
    constants <- list(mu = mu, adj = adj, num = num, N = 4, L = 6)
    data <- list(y = y)
    inits <- list(tau = tau, gamma = gamma, x = x)

    Rmodel <- nimbleModel(code, constants, data, inits)
    Cmodel <- compileNimble(Rmodel)
    lp <- -574.964
    
    expect_lt(abs(calculate(Rmodel) - lp), 1E-5,
                 label = 'calculate for dcar_proper()')
    
    expect_lt(abs(calculate(Cmodel) - lp), 1E-5,
                 label = 'calculate for dcar_proper(), compiled')
    
    weights <- rep(1, 6)
    CM <- as.carCM(adj, weights, num)
    C <- CM$C
    M <- CM$M
    Q <- tau * diag(1/M) %*% (diag(4) - gamma*CAR_calcCmatrix(C, adj, num))
    lp <- dmnorm_chol(x, mu, chol = chol(Q), prec_param = TRUE, log = TRUE)
    
    expect_equal(calculate(Rmodel, 'x[1:4]'), lp,
                 info = 'R density evaluation for dcar_proper()')
    
    expect_equal(calculate(Cmodel, 'x[1:4]'), lp,
                 info = 'C density evaluation for dcar_proper()')

    set.seed(0); xnew <- rmnorm_chol(n = 1, mu, chol = chol(Q), prec_param = TRUE)
    set.seed(0); simulate(Rmodel, 'x[1:4]')
    set.seed(0); simulate(Cmodel, 'x[1:4]')
    
    expect_equal(xnew, Rmodel$x, info = 'R dcar_proper() simulate function')
    expect_equal(xnew, Cmodel$x, info = 'R dcar_proper() simulate function')
    
    Rmodel$x <- x
    Cmodel$x <- x
    
    conf <- configureMCMC(Rmodel)
    conf$addMonitors('x')
    Rmcmc <- buildMCMC(conf)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

    niter <- 20
    set.seed(0); Rmcmc$run(niter)
    set.seed(0); Cmcmc$run(niter)

    Rsamples <- as.matrix(Rmcmc$mvSamples)
    Csamples <- as.matrix(Cmcmc$mvSamples)

    sampleNames <- colnames(Rsamples)

    expect_true(all(Rsamples[, sampleNames] - Csamples[, sampleNames] == 0),
                info = 'agreement between R and C sampling of dcar_proper')

    expect_lt(max(abs(as.numeric(Csamples[20, sampleNames]) -
                 c(-0.86201288, 0.07689823, 2.04074467, 0.24380342, 3.00405982, -3.25336913))), 1e-6,
                 label = 'exact sample values for dcar_proper')
})


## testing dcar_proper distribution gives correct
## likelihood evaluation, when Cmatrix is singular
test_that('dcar_proper gives correct likelihood with singular Cmatrix', {
    nHabRows <- nHabCols <- 5
    adj <- NULL
    numadj <- NULL
    numHabWindows <- nHabRows * nHabCols 
    for(i in 1:numHabWindows){
        {
            lenadj <- length(adj)
            if(i == 1) adj <- c(adj, c(2, nHabCols + 1))
            else if((2 <= i) & (i <= nHabCols - 1)) adj <- c(adj, c(i-1, i+1, i+nHabCols))
            else if(i == nHabCols) adj <- c(adj, c(i-1, i+nHabCols))
            else if(i %in% c(1:(nHabRows-2)*nHabCols+1)) adj <- c(adj, c(i-nHabCols, i+1, i+nHabCols))
            else if(i %in% c(2:(nHabRows-1)*nHabCols)) adj <- c(adj, c(i-nHabCols, i-1, i+nHabCols))
            else if(i == (nHabRows-1)*nHabCols+1) adj <- c(adj, c(i-nHabCols, i+1))
            else if(((nHabRows-1)*nHabCols+2 <= i) & (i <= nHabCols*nHabRows-1)) adj <- c(adj, c(i-nHabCols, i-1, i+1))
            else if(i == nHabCols*nHabRows) adj <- c(adj, c(i-nHabCols, i-1))
            else adj <- c(adj, c(i-nHabCols, i-1, i+1, i+nHabCols))
            numadj[i] <- length(adj) - lenadj
        }
    }
    num <- numadj
    ##
    tau <- 1
    gamma <- 0.5
    mu <- rep(0, numHabWindows)
    x <- rep(0, numHabWindows)
    ##
    C <- CAR_calcC(adj, num)
    M <- CAR_calcM(num)
    Cmatrix <- CAR_calcCmatrix(C, adj, num)
    ##
    ## dmvnorm
    ## function (x, mean = rep(0, p), sigma = diag(p), log = FALSE)
    ## cov = (I-gamma*Cmatrix)^-1* %*% Mmatrix / tau )
    CovMatrix <- solve(diag(numHabWindows) - gamma*Cmatrix) %*% diag(M) / tau
    ## dmvnorm(x, mean = mu, sigma = CovMatrix)   ## 0.00009134347
    ## (2pi)^(-N/2) * |Cov|^(-1/2) * exp(-1/2 * (x-mu)' %*% Cov^-1 %*% (x-mu))
    lp_true <- (2*pi)^(-numHabWindows/2) * det(CovMatrix)^(-1/2)   ## 0.00009134347
    ##
    lp_dcar <- dcar_proper(x, mu = mu, adj = adj, num = num, tau = tau, gamma = gamma)
    ##
    expect_equal(lp_true, lp_dcar)
    ##
    set.seed(0)
    xnew <- rcar_proper(n = 1, mu = mu, adj = adj, num = num, tau = tau, gamma = gamma)
    expect_true(all(round(xnew, 8) == c(0.60896921, -0.18071413, 0.86723430, 0.87426463, 0.66326831, -0.95558320, -0.56113360, -0.22185377, 0.14990019, 1.37682998, 0.32056427, -0.48639310, -0.65728773, -0.16961897, -0.30264915, -0.23875949, 0.10629314, -0.40382663, 0.18898157, -0.72969397, -0.09361655, 0.25556273, 0.16314071, 0.46935282, -0.04233136)))
})


test_that('CAR conjugacy checking new skipExpansionsNode system', {
    ##
    code <- nimbleCode({
        S[1:N] ~ dcar_normal(adj[1:L], weights[1:L], numneighbours[1:N], 1)
        for(i in 1:K) {
            beta[i] ~ dnorm(0, 1)
        }
        for(i in 1:N){
            eta[i] <- inprod(beta[1:K], x[1:K])
            mu[i] <- S[i] + eta[i]
            y[i] ~ dnorm(mu[i], 1)
        }
    })
    ##
    N <- 3
    L <- 4
    K <- 7
    ##
    constants <- list(N=N, L=L, K=K, adj=c(2,1,3,2), weights=rep(1,L), numneighbours=c(1,2,1))
    data <- list(y = rep(0,N))
    inits <- list(S = rep(0,N), beta = rep(0,K), x=1:K)
    ##
    Rmodel <- nimbleModel(code, constants, data, inits)
    conf <- configureMCMC(Rmodel)
    Rmcmc <- buildMCMC(conf)
    ##
    expect_true(class(Rmcmc) == 'MCMC')
    expect_true(conf$samplerConfs[[8]]$name == 'CAR_normal')
    expect_true(class(Rmcmc$samplerFunctions$contentsList[[8]]$componentSamplerFunctions$contentsList[[1]]) == 'CAR_scalar_conjugate')
    expect_true(class(Rmcmc$samplerFunctions$contentsList[[8]]$componentSamplerFunctions$contentsList[[2]]) == 'CAR_scalar_conjugate')
    expect_true(class(Rmcmc$samplerFunctions$contentsList[[8]]$componentSamplerFunctions$contentsList[[3]]) == 'CAR_scalar_conjugate')
    ##
    compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
    Cmodel <- compiledList$model; Cmcmc <- compiledList$mcmc
    ##
    set.seed(0); Rsamples <- runMCMC(Rmcmc, 10)
    set.seed(0); Csamples <- runMCMC(Cmcmc, 10)
    ##
    expectedSamples <- c(0.97357331, 0.07601302, 0.10439196, -0.37719856, 0.15912985, 0.03509085, -0.01162275, 0.17958068, -0.34811805, 0.10319592)
    Rcolnames <- colnames(Rsamples)
    ##
    expect_true(all(round(as.numeric(Rsamples[10,Rcolnames]),8) == expectedSamples))
    expect_true(all(round(as.numeric(Csamples[10,Rcolnames]),8) == expectedSamples))
})




options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)

## File for testing distributions provided by NIMBLE

source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
context('Testing NIMBLE distributions')

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

goldFileName <- 'distributionsTestLog_Correct.Rout'
tempFileName <- 'distributionsTestLog.Rout'
generatingGoldFile <- !is.null(nimbleOptions('generateGoldFileForDistributionsTesting'))
outputFile <- if(generatingGoldFile) file.path(nimbleOptions('generateGoldFileForDistributionsTesting'), goldFileName) else tempFileName

sink(outputFile)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)


## mvt
x <- c(1, 1, 2)
mn <- c(1, 2, 3)
sc <- diag(c(1, 2, 3))
sc[1, 3] <- sc[3, 1] <- 0.1
df <- 5

## test use directly from R

truth <- mvtnorm::dmvt(x, delta = mn, sigma = sc, df = df, log = FALSE)

test_that("dmvt_chol calculates density correctly in R",

          expect_equal(dmvt_chol(x, mn, chol(sc), df, prec_param = FALSE),
                       truth,
                       info = paste0("incorrect dmvt calculation in R")))

## test use through nimble function


test_that("Test that dmvt_chol works correctly in nimbleFunction", {
          nf <- nimbleFunction(
              run = function(x = double(1), mn = double(1),
                             scale = double(2), df = double(0)) {
                  returnType(double(0))
                  ch <- chol(scale)
                  out <- dmvt_chol(x = x, mu = mn, cholesky = ch,
                                   df = df, prec_param = FALSE, log = FALSE)
                  return(out)
              }
          )
          expect_equal(nf(x, mn, sc, df), (truth), 
                       info = paste0("incorrect dmvt value in R nimble function"))
          cnf <- compileNimble(nf)
          expect_equal(cnf(x, mn, sc, df), (truth), 
                       info = paste0("incorrect dmvt value in compiled nimble function"))
})

## test use in model

test_that("Test that mvt calculate and simulate are correct in nodeFunctions", {
    mvt_code <- nimbleCode({
        x[1:3] ~ dmvt(mn[1:3], scale = sc[1:3,1:3], df = df)
    })
    mvt_model <- nimbleModel(mvt_code, constants = list(mn = mn, sc = sc, prec = FALSE, df = df))
    mvt_model$x <- x
    expect_equal(exp(mvt_model$calculate()), (truth),
                 info = paste0("incorrect likelihood value for uncompiled dmvt"))
    c_mvt_model <- compileNimble(mvt_model)
    c_mvt_model$x
    expect_equal(exp(c_mvt_model$calculate()), (truth),
                 info = paste0("incorrect likelihood value for compiled dmvt"))

    ## random sampling
    set.seed(0)
    r_samps <- t(replicate(10000, rmvt_chol(n = 1, mn, chol(sc), df, prec_param = FALSE)))
    true_cov <- sc*df/(df-2)

    expect_lt(max(abs(colMeans(r_samps)- mn)) , 0.03,
                 label = "Difference in random sample means in R exceeds tolerance")
    expect_lt(max(abs(cov(r_samps) - true_cov)) , 0.1,
                 label = "Difference in random sample covs in R exceeds tolerance")

    nf_sampling <- nimbleFunction(
        run = function(mn = double(1), scale = double(2), df = double(0)) {
            returnType(double(1))
            ch <- chol(scale)
            out <- rmvt_chol(n = 1, mu = mn, cholesky = ch,
                             df = df, prec_param = FALSE)
            return(out)
        }
    )
    nf_samps <- t(replicate(10000, nf_sampling(mn, sc, df)))
    expect_lt(abs(colMeans(nf_samps)-mn) , 0.03,
                 label = "Difference in means in nf exceeds tolerance")
    expect_lt(abs(cov(nf_samps)- true_cov) , 0.1,
                 label = "Difference in covs in nf exceeds tolerance")

    ## sampling via `simulate`
    set.seed(0)
    simul_samp <- function(model) {
        model$simulate()
        return(model$x)
    }
    
    simul_samps <- t(replicate(10000, simul_samp(c_mvt_model)))
    
    expect_lt(abs(colMeans(simul_samps)-mn) , 0.03,
                 label = "Difference in means in random samples (simulate) exceeds tolerance")
    
    expect_lt(abs(cov(simul_samps)-true_cov) , 0.1,
                 label = "Difference in covs in random samples (simulate) exceeds tolerance")
})

## dlkj

## test use directly from R

eta <- 3.3
k <- 5

test_that("dlkj_corr_cholesky calculates density correctly", {
    set.seed(1)
    U <- rlkj_corr_cholesky(1, eta, k)
    truth <- sum(log(diag(U)) * (k-(1:k)+2*eta-2))
    expect_equal(dlkj_corr_cholesky(U, eta, k, log = TRUE),
                 truth,
                 info = paste0("incorrect dlkj calculation in R"))
    
    nf <- nimbleFunction(
        run = function(eta = double(0), p = double(0)) {
            returnType(double(0))
            x <- rlkj_corr_cholesky(n = 1, eta, p)
            out <- dlkj_corr_cholesky(x = x, eta, p, log = TRUE)
            return(out)
        }
    )

    set.seed(1)
    expect_equal(nf(eta, k), truth,
                 info = paste0("incorrect dlkj value in R nimble function"))
    cnf <- compileNimble(nf)
    set.seed(1)
    expect_equal(cnf(eta, k), truth,
                 info = paste0("incorrect dlkj value in compiled nimble function"))
})

## test use in model

## code from rethinking package for comparison
rlkjcorr_rethinking <- function ( n , K , eta = 1 ) {

    stopifnot(is.numeric(K), K >= 2, K == as.integer(K))
    stopifnot(eta > 0)
    #if (K == 1) return(matrix(1, 1, 1))
    
    f <- function() {
        alpha <- eta + (K - 2)/2
        r12 <- 2 * rbeta(1, alpha, alpha) - 1
        R <- matrix(0, K, K) # upper triangular Cholesky factor until return()
        R[1,1] <- 1
        R[1,2] <- r12
        R[2,2] <- sqrt(1 - r12^2)
        if(K > 2) for (m in 2:(K - 1)) {
            alpha <- alpha - 0.5
            y <- rbeta(1, m / 2, alpha)
            # Draw uniformally on a hypersphere
            z <- rnorm(m, 0, 1)
            z <- z / sqrt(crossprod(z)[1])
            
            R[1:m,m+1] <- sqrt(y) * z
            R[m+1,m+1] <- sqrt(1 - y)
        }
        return(crossprod(R))
    }
    R <- replicate( n , f() )
    if ( dim(R)[3]==1 ) {
        R <- R[,,1]
    } else {
        # need to move 3rd dimension to front, so conforms to array structure that Stan uses
        R <- aperm(R,c(3,1,2))
    }
    return(R)
}


test_that("Test that dlkj calculate and simulate are correct in nodeFunctions and nimbleFunctions", {
    set.seed(1)
    U <- rlkj_corr_cholesky(1, eta, k)
    truth <- sum(log(diag(U)) * (k-(1:k)+2*eta-2))
    lkj_code <- nimbleCode({
        U[1:k,1:k] ~ dlkj_corr_cholesky(eta = eta, p = k)
    })
    lkj_model <- nimbleModel(lkj_code, constants = list(eta = eta, k = k))
    lkj_model$U <- U
    expect_equal(lkj_model$calculate(), truth,
                 info = paste0("incorrect log-likelihood value for uncompiled dlkj"))
    c_lkj_model <- compileNimble(lkj_model)
    expect_equal(c_lkj_model$calculate(), truth,
                 info = paste0("incorrect log-likelihood value for compiled dlkj"))

    ## random sampling - compare to code from rethinking package
    
    nf_sampling <- nimbleFunction(
        run = function(eta = double(0), p = double(0)) {
            returnType(double(2))
            out <- rlkj_corr_cholesky(n = 1, eta, p)
            return(out)
        }
    )

    set.seed(1)
    R <- rlkjcorr_rethinking(1, K = k, eta = eta)
    set.seed(1)
    testR <- crossprod(rlkj_corr_cholesky(1, eta, k))
    set.seed(1)
    testnf <- crossprod(nf_sampling(eta, k))
    expect_equal(R, testR, info = "nimble-based and rethinking-based rlkj simulations differ")
    expect_equal(R, testnf, info = "nimble-based and rethinking-based rlkj simulations differ")
})

test_that("rlkj_corr_cholesky size processing works", {
    nf <- nimbleFunction(
        run = function(eta = double(0), p = double(0), A = double(2)) {
            returnType(double(3))
            tmp <- rlkj_corr_cholesky(n = 1, eta, p)
            x <- nimArray(0, c(p+3,p+3,2))
            x[2:(p+1),3:(p+2), 2] <- t(tmp) %*% tmp %*% A
            return(x)
        }
    )
    cnf <- compileNimble(nf)

    set.seed(1)
    A <- matrix(rnorm(k*k), k)
    set.seed(1)
    x <- crossprod(rlkj_corr_cholesky(1, eta, k)) %*% A
    set.seed(1)
    x2 <- cnf(eta, k, A)[2:(k+1),3:(k+2),2]
    expect_equal(x, x2)
})

## dmulti and dcat

set.seed(0)
normGen <- rmulti(1, 1000, prob = c(.1, .1, .8))
set.seed(0)
unnormGen <- rmulti(1, 1000, prob = c(10, 10, 80))
normResult <- dmulti(normGen, prob = c(.1, .1, .8), log = TRUE)
unnormResult <- dmulti(normGen, prob = c(10, 10, 80), log = TRUE)

test_that("rmulti handles 'probs' that do not sum to one: ",
              expect_identical(normGen, unnormGen,
                               info = "normalized and unnormalized probabilities give different results"))
test_that("dmulti handles 'probs' that do not sum to one: ",
              expect_equal(normResult, unnormResult,
                           info = "normalized and unnormalized probabilities give different results"))

set.seed(0)
normGen <- rcat(1, prob = c(.1, .1, .8))
set.seed(0)
unnormGen <- rcat(1, prob = c(10, 10, 80))
normResult <- dcat(normGen, prob = c(.1, .1, .8), log = TRUE)
unnormResult <- dcat(normGen, prob = c(10, 10, 80), log = TRUE)

test_that("rcat handles 'probs' that do not sum to one: ",
              expect_identical(normGen, unnormGen,
                               info = "normalized and unnormalized probabilities give different results"))
test_that("dcat handles 'probs' that do not sum to one: ",
              expect_equal(normResult, unnormResult,
                           info = "normalized and unnormalized probabilities give different results"))

## dinvgamma
set.seed(0)
y <- 1.1; a <- 1; c <- 2; alpha <- 3; beta <- 2; theta <- 1

manDens <- alpha*log(beta) - lgamma(alpha) - (alpha+1)*log(y) - beta/y
test_that("Test that dinvgamma gets correct result: ",
              expect_lt(abs(manDens - dinvgamma(y, alpha, scale = beta, log = TRUE)) , 1e-15))

set.seed(0)
smp1 <- rinvgamma(100000, shape = alpha, scale = beta)
set.seed(0)
smp2 <- rinvgamma(100000, shape = alpha, rate = 1/beta)

test_that("Test that rinvgamma with scale gets correct result: ",
              expect_lt(abs(beta / (alpha-1)-mean(smp1)) , 0.01,
                           label = "Difference in mean exceeds tolerance"))
test_that("Test that rinvgamma with rate gets correct result: ",
              expect_lt(abs(beta / (alpha-1)-mean(smp2)) , 0.01,
                           label = "Difference in mean exceeds tolerance"))
test_that("Test that rinvgamma with scale gets correct result: ",
              expect_lt(abs(beta/((alpha-1)*sqrt(alpha-2))-sd(smp1)) ,  0.1,
                           label = "Difference in sd exceeds tolerance"))
test_that("Test that rinvgamma with rate gets correct result: ",
              expect_lt(abs(beta/((alpha-1)*sqrt(alpha-2))-sd(smp2)) , 0.1,
                           label = "Difference in sd exceeds tolerance"))

quantile <- quantile(smp1, .15)
attributes(quantile) <- NULL
test_that("Test that pinvgamma gets correct result: ",
              expect_lt(abs(qinvgamma(.15, alpha, scale = beta) - quantile) , 0.005,
                           label = "Difference in quantile exceeds tolerance"))
p <- mean(smp1 < .5)
test_that("Test that qinvgamma gets correct result: ",
              expect_lt(abs(pinvgamma(.5, alpha, scale = beta)-p) , 0.005,
                           label = "Difference in probability exceeds tolerance"))

test_that("dinvgamma-dinvgamma conjugacy with dependency using scale", {
    code <- nimbleCode({
        y ~ dinvgamma(a, scale = c*theta)
        theta ~ dinvgamma(alpha, beta)
    })
    m = nimbleModel(code, data = list(y = y),
                    inits = list(theta = theta, a = a, c = c, alpha = alpha, beta = beta))
    conf <- configureMCMC(m)
    samplers <- conf$getSamplers()
    expect_identical(samplers[[1]]$name, 'RW',
                     info = "conjugacy improperly detected")
})

test_that("dinvgamma-dinvgamma conjugacy with linear dependency", {
    code <- nimbleCode({
        y ~ dinvgamma(a, rate = 2 + c*theta)
        theta ~ dinvgamma(alpha, beta)
    })
    m = nimbleModel(code, data = list(y = y),
                    inits = list(theta = theta, a = a, c = c, alpha = alpha, beta = beta))
    conf <- configureMCMC(m)
    samplers <- conf$getSamplers()
    expect_identical(samplers[[1]]$name, 'RW',
                     info = "conjugacy improperly detected")
})

test_that("dinvgamma-dinvgamma conjugacy with dependency using rate", {
    code <- nimbleCode({
        y ~ dinvgamma(a, rate = c*theta)
        theta ~ dinvgamma(alpha, beta)
    })
    m = nimbleModel(code, data = list(y=y),
                    inits = list(theta = theta, a = a, c = c, alpha = alpha, beta = beta))
    conf <- configureMCMC(m)
    samplers <- conf$getSamplers()
    expect_identical(samplers[[1]]$name, 'conjugate_dinvgamma_dinvgamma_multiplicative',
                     info = "conjugacy not detected")

    mcmc <- buildMCMC(conf)
    comp <- compileNimble(m, mcmc)
    set.seed(0)
    comp$mcmc$run(100)
    smp <- as.matrix(comp$mcmc$mvSamples)
    
    manualSampler <- function(n, y, a, c, alpha, beta) {
        out <- rep(0, n)
        shape = a + alpha
        scale = beta + 1/(c*y)
        set.seed(0)
        out1 <- 1/rgamma(n, shape, rate = scale)
        set.seed(0)
        out2 <- rinvgamma(n, shape, scale = scale)
        return(list(out1, out2))
    }
    smpMan <- manualSampler(100, y, a, c, alpha, beta) 
    expect_identical(smp[,1], smpMan[[1]],
                               info = "NIMBLE conjugate sampler and manual sampler results differ (1)")
    expect_identical(smp[,1], smpMan[[2]],
                     info = "NIMBLE conjugate sampler and manual sampler results differ (2)")
})


test_that("dgamma-dinvgamma conjugacy with dependency using rate", {
    code <- nimbleCode({
        y ~ dinvgamma(a, scale = c*theta)
        theta ~ dgamma(alpha, beta)
    })
    m = nimbleModel(code, data = list(y=y),
                    inits = list(theta = theta, a = a, c = c, alpha = alpha, beta = beta))
    conf <- configureMCMC(m)
    samplers <- conf$getSamplers()
    expect_identical(samplers[[1]]$name, 'conjugate_dgamma_dinvgamma_multiplicative',
                     info = "conjugacy not detected")

    mcmc <- buildMCMC(conf)
    comp <- compileNimble(m, mcmc)
    set.seed(0)
    comp$mcmc$run(10)
    smp <- as.matrix(comp$mcmc$mvSamples)
    
    manualSampler <- function(n, y, a, c, alpha, beta) {
        out <- rep(0, n)
        shape = a + alpha
        rate = beta + c/y
        set.seed(0)
        out <- rgamma(n, shape, rate = rate)
        return(out)
    }
    smpMan <- manualSampler(10, y, a, c, alpha, beta)
    expect_identical(smp[,1], smpMan,
                     info = "NIMBLE gamma conjugate sampler and manual sampler results differ")
})

# dinvwish_chol

test_that("Test that rinvwish_chol, dinvwish_chol, dwish_chol, and rwish_chol give correct results", {
    set.seed(1)
    df <- 20.3
    d <- 5
    C <- crossprod(matrix(rnorm(d^2, 0), d))
    pm <- C / (df - d - 1)
    
    n <- 100
    draws1 <- draws2 <- draws3 <- array(0, c(d, d, n))
    set.seed(1)
    for(i in 1:n)
        draws1[,,i] <- solve(rwish_chol(1, chol(C), df, scale_param = FALSE))
    pmean1 <- apply(draws1, c(1,2), mean)
    
    set.seed(1)
    for(i in 1:n)
        draws2[,,i] <- rinvwish_chol(1, chol(C), df, scale_param = TRUE)
    pmean2 <- apply(draws2, c(1,2), mean)
    
    set.seed(1)
    for(i in 1:n)
        draws3[,,i] <- rinvwish_chol(1, chol(solve(C)), df, scale_param = FALSE)
    pmean3 <- apply(draws3, c(1,2), mean)
    expect_lt(abs(max(abs(pmean1 - pm))) , 0.03,
                 label = "mean of inverse of rwish draws differs from truth")
    expect_lt(abs(max(abs(pmean2 - pm))) , 0.03,
                 label = "mean of rinvwish with scale draws differs from truth")
    expect_lt(abs(max(abs(pmean3 - pm))) , 0.03,
                 label = "mean of rinvwish with rate differs from truth")

    # draws1 is from rwish not rinvwish, but for comparing density values it doesn't matter their origin.
    dens1 <- dinvwish_chol(draws1[,,1], chol(C), df, scale_param = TRUE, log = TRUE)
    dens2 <- dinvwish_chol(draws1[,,1], chol(solve(C)), df, scale_param = FALSE, log = TRUE)
    
    dfun <- function(W, S, nu) {
        k <- nrow(W)
        U = chol(S)
        return(-(log(2)*nu*k/2+(k*(k-1)/4)*log(pi) +sum(lgamma((nu + 1 - 1:k)/2))) + nu*sum(log(diag(U))) -
               (nu+k+1)*sum(log(diag(chol(W)))) -0.5*sum(diag(S %*% solve(W))))
    }
    
    dens3 <-  dfun(draws1[,,1], C, df)
    expect_lt(abs(dens1-dens3) ,  0.000001,
                 label = "dinvwish with scale differs from truth")
    expect_lt(abs(dens2-dens3) ,  0.000001,
                 label = "dinvwish with rate differs from truth")

    dens1 <- dwish_chol(draws1[,,1], chol(C), df, scale_param = TRUE, log = TRUE)
    dens2 <- dwish_chol(draws1[,,1], chol(solve(C)), df, scale_param = FALSE, log = TRUE)

    dfun <- function(W, S, nu) {
        k <- nrow(W)
        U = chol(S)
        return(-(log(2)*nu*k/2+(k*(k-1)/4)*log(pi) +sum(lgamma((nu + 1 - 1:k)/2))) - nu*sum(log(diag(U))) +
               (nu-k-1)*sum(log(diag(chol(W)))) -0.5*sum(diag(solve(S, W))))               
    }

    dens3 <-  dfun(draws1[,,1], C, df)
    expect_lt(abs(dens1-dens3) , 0.000001,
                 label = "dinvwish with scale differs from truth")
    expect_lt(abs(dens2-dens3) , 0.000001,
                 label = "dinvwish with rate differs from truth")

})
                                                                        
set.seed(1)

trueCor <- matrix(c(1, .3, .7, .3, 1, -0.2, .7, -0.2, 1), 3)
covs <- c(3, 2, .5)

trueCov = diag(sqrt(covs)) %*% trueCor %*% diag(sqrt(covs))

M = 3
n = 20
S = crossprod(matrix(rnorm(M*M), M)) # diag(rep(1,3))
nu = 4                                    
mu = 1:3
Y = mu + t(chol(trueCov)) %*% matrix(rnorm(3*n), ncol = n)
data <- list(Y = t(Y), n = n, M = M)

code <- nimbleCode( {
  for(i in 1:n) {
    Y[i, 1:M] ~ dmnorm(mu[1:M], cov = Omega[1:M,1:M])
  }
  Omega[1:M,1:M] ~ dinvwish(S[1:M,1:M], nu)
})

codeR <- nimbleCode( {
  for(i in 1:n) {
    Y[i, 1:M] ~ dmnorm(mu[1:M], cov = Omega[1:M,1:M])
  }
  Omega[1:M,1:M] ~ dinvwish(R = R[1:M,1:M], df = nu)
})

newDf = nu + n
newS = S + tcrossprod(Y- mu)
OmegaTrueMean = newS / (newDf - M - 1)
                   
invwishRV <- array(0, c(M, M, 10000))
for(i in 1:10000) {
  invwishRV[,,i] <- rinvwish_chol(1, chol(newS), df = newDf, scale_param = TRUE)
}
OmegaSimTrueSDs = apply(invwishRV, c(1,2), sd)

test_mcmc(model = code, name = 'conjugate inverse Wishart', data = data,
          seed = 0, numItsC = 1000, inits = list(Omega = trueCov, mu = mu, S = S, nu = nu),
          results = list(mean = list(Omega = OmegaTrueMean),
            sd = list(Omega = OmegaSimTrueSDs)),
          resultsTolerance = list(mean = list(Omega = matrix(.05, M,M)),
            sd = list(Omega = matrix(0.1, M, M))))

test_mcmc(model = codeR, name = 'conjugate inverse Wishart', data = data,
          seed = 0, numItsC = 1000, inits = list(Omega = trueCov, mu = mu, R = solve(S), nu = nu),
          results = list(mean = list(Omega = OmegaTrueMean),
            sd = list(Omega = OmegaSimTrueSDs)),
          resultsTolerance = list(mean = list(Omega = matrix(.05, M,M)),
            sd = list(Omega = matrix(0.1, M, M))))

test_that("dinvwish-dmnorm conjugacy", {
    data = list(Y=t(Y))
    constants = list(n = n, M = M)
    m = nimbleModel(code, data = data, inits = list(nu = nu, Omega = trueCov, mu = mu, S = S),
                    constants = constants)
    conf <- configureMCMC(m)
    expect_identical(conf$getSamplers()[[1]]$name, 'conjugate_dinvwish_dmnorm_identity',
                     info = "conjugacy not detected")
})


## dflat
test_that("dflat and dhalfflat usage", {
    set.seed(0)
    code <- nimbleCode({
        for(i in 1:n) 
            y[i] ~ dnorm(mu, sd = sigma)
        mu ~ dflat()
        sigma ~ dhalfflat()
    })
    
    n <- 100
    inits <- list(mu = 0, sigma = 1)
    m <- nimbleModel(code, constants = list(n = n), data = list(y = rnorm(n)), inits = inits)
    cm <- compileNimble(m)
    
    
    expect_identical(m$calculate('mu'), 0, "incorrect R calculate for dflat")
    expect_identical(m$calculate('mu'), cm$calculate('mu'), "incorrect compiled calculate for dflat")
    expect_identical(m$calculate('sigma'), 0, "incorrect R calculate for dhalfflat")
    expect_identical(m$calculate('sigma'), cm$calculate('sigma'), "incorrect compiled calculate for dhalfflat")
    m$simulate('mu'); m$simulate('sigma');  cm$simulate('mu'); cm$simulate('sigma')
    expect_identical(m$mu, NaN, "incorrect R simulate for dflat")
    expect_identical(m$mu, cm$mu, "incorrect compiled simulate for dflat")
    expect_identical(m$sigma, NaN, "incorrect R simulate for dhalfflat")
    expect_identical(m$sigma, cm$sigma, "incorrect compiled simulate for dhalfflat")

    m$setInits(inits)
    cm$setInits(inits)
    
    conf <- configureMCMC(m, nodes = 'mu')
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project = m)
    set.seed(1)
    mcmc$run(10)
    set.seed(1)
    cmcmc$run(10000)
    rsmp <- c(as.matrix(mcmc$mvSamples)[,'mu'])
    csmp <- c(as.matrix(cmcmc$mvSamples)[,'mu'])
    
    
    expect_identical(conf$getSamplers()[[1]]$name, 'conjugate_dflat_dnorm_identity', info = "did not detect conjugacy")
    expect_identical(rsmp, csmp[1:10], "R and compiled samples don't match")
    expect_lt(abs(mean(csmp) - mean(m$y)) , 0.001, label = "posterior mean for mu not correct")
    expect_lt(abs(sd(csmp) - 1/sqrt(n)) , 0.002, label = "posterior sd for mu not correct")
    
    m$setInits(inits)
    cm$setInits(inits)
    conf <- configureMCMC(m, nodes = 'sigma')
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project = m, resetFunctions = TRUE)
    set.seed(1)
    mcmc$run(10)
    set.seed(1)
    cmcmc$run(10000)
    rsmp <- c(as.matrix(mcmc$mvSamples)[,'sigma'])
    csmp <- c(as.matrix(cmcmc$mvSamples)[,'sigma'])
    
    expect_identical(conf$getSamplers()[[1]]$name, 'conjugate_dhalfflat_dnorm_identity', info = "did not detect conjugacy")
    expect_identical(rsmp, csmp[1:10], info = "R and compiled samples don't match")
    expect_lt(abs(mean(csmp^2) - var(m$y)) ,  0.03, label = "posterior mean for sigma not correct")
    
    m$setInits(inits)
    cm$setInits(inits)
    conf <- configureMCMC(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project = m, resetFunctions = TRUE)
    set.seed(1)
    mcmc$run(10)
    set.seed(1)
    cmcmc$run(10000)
    rsmp <- as.matrix(mcmc$mvSamples)
    csmp <- as.matrix(cmcmc$mvSamples)
    
    expect_equivalent(rsmp, csmp[1:10, ], info = "R and compiled samples don't match")
# for some reason, these are not identical - after 8-10 iterations, R and C values diverge in 8th decimal place or so
    expect_lt(abs(mean(csmp[ , 'mu']) - mean(m$y)) ,  0.001, label = "posterior mean for mu not correct")
    expect_lt(abs(sd(csmp[ , 'mu']) - sd(m$y)/sqrt(n)) ,  0.003, label = "posterior sd for mu not correct")

    expect_lt(abs(mean(csmp[ , 'sigma']^2) - var(m$y)) ,  0.05, label = "posterior mean for sigma not correct")
})

test_that("ddexp usage", {
    n <- 100000
    mean <- 1
    scale <- 3
    p_exp <- 0.8
    value <- mean + qexp(p_exp, 1/scale)
    symm_value <- 2*mean - value
    p <- 1-(1-p_exp)/2
    dens <- 0.5 * dexp(value-mean, 1/scale) 
    expect_identical(ddexp(value, mean, scale), dens, "basic use of ddexp")
    expect_equal(ddexp(value, mean, scale, log = TRUE), log(dens), info = "basic use of ddexp, log scale")
    expect_identical(ddexp(symm_value, mean, scale), dens, "basic use of ddexp, left half")
    expect_equal(ddexp(symm_value, mean, scale, log = TRUE), log(dens), info = "basic use of ddexp, left half, log scale")

    set.seed(1)
    smp <- rdexp(n, mean, scale)
    expect_lt(abs(mean(smp) - mean) ,  0.02, label = "mean of random sample")
    expect_lt(abs(mean(smp > value) - (1-p)) ,  0.01, label = "upper quantile of random sample")
    expect_lt(abs(mean(smp < symm_value) - (1-p)) ,  0.01, label = "lower quantile of random sample")

    expect_identical(pdexp(value, mean, scale), p, "basic use of pdexp")
    expect_identical(pdexp(value, mean, scale, log.p = TRUE), log(p), "basic use of pdexp, log scale")
    expect_lt(abs(pdexp(value, mean, scale, lower.tail = FALSE) -  (1-p)) ,  0.001, label = "basic use of pdexp, upper tail")
    expect_identical(pdexp(value, mean, scale, lower.tail = FALSE, log.p = TRUE), log(1-p), "basic use of pdexp, upper tail, log scale")

    expect_lt(abs(pdexp(symm_value, mean, scale)- (1-p)) ,  0.001, label = "basic use of pdexp, left half")
    expect_lt(abs(pdexp(symm_value, mean, scale, log.p = TRUE) - log(1-p)) ,  0.001, label = "basic use of pdexp, left half, log scale")
    expect_identical(pdexp(symm_value, mean, scale, lower.tail = FALSE), p, "basic use of pdexp, left half, upper tail")
    expect_identical(pdexp(symm_value, mean, scale, lower.tail = FALSE, log.p = TRUE), log(p),
                 "basic use of pdexp, left half, upper tail, log scale")

    expect_identical(qdexp(p, mean, scale), value, "basic use of qdexp")
    expect_identical(qdexp(log(p), mean, scale, log.p = TRUE), value, "basic use of qdexp, log scale")
    expect_identical(qdexp(p, mean, scale, lower.tail = FALSE), symm_value, "basic use of qdexp, upper tail")
    expect_identical(qdexp(log(p), mean, scale, lower.tail = FALSE, log.p = TRUE), symm_value, "basic use of qdexp, upper tail, log scale")
    
    expect_identical(qdexp(1-p, mean, scale), symm_value, "basic use of qdexp, left half")
    expect_identical(qdexp(log(1-p), mean, scale, log.p = TRUE), symm_value, "basic use of qdexp, left half, log scale")
    expect_identical(qdexp(1-p, mean, scale, lower.tail = FALSE), value, "basic use of qdexp, left half, upper tail")
    expect_identical(qdexp(log(1-p), mean, scale, lower.tail = FALSE, log.p = TRUE), value, "basic use of qdexp, left half, upper tail, log scale")

    a <- nimbleFunction(
        run = function(x = double(0), p = double(0), mu = double(0), scale = double(0)) {
            returnType(double(1))
            dens <- ddexp(x, mu, scale, log = FALSE)
            smp <- rdexp(10, mu, scale)
            prob <- pdexp(x, mu, scale)
            quantile <- qdexp(p, mu, scale)
            return(c(dens, smp, prob, quantile))
        })
    ca <- compileNimble(a)
    set.seed(1)
    output <- a(value, p, mean, scale)
    set.seed(1)
    coutput <- ca(value, p, mean, scale)
    expect_identical(output, coutput, "compiled and uncompiled use of ddexp don't match")
    expect_identical(output, c(dens, smp[1:10], p, value), "use of ddexp in nimbleFunction does not match correct values") 
    
    code <- nimbleCode({
        for(i in 1:n) 
            y[i] ~ dnorm(mu, 1)
        mu ~ ddexp(mu0, scale = sigma)
    })
    
    n <- 10
    inits <- list(mu0 = mean, sigma = scale)
    m <- nimbleModel(code, constants = list(n = n), data = list(y = rnorm(n)), inits = inits)
    cm <- compileNimble(m)
    m$mu <- cm$mu <- value
    expect_equal(m$calculate('mu'), log(dens), info = "incorrect R calculate for ddexp")
    expect_identical(m$calculate('mu'), cm$calculate('mu'), "incorrect compiled calculate for ddexp")
    set.seed(1)
    cm$simulate('mu')
    expect_identical(smp[1], cm$mu, "incorrect compiled simulate for ddexp")

    ## check use of rate
    set.seed(1)
    x <- rdexp(1, 0.5, scale = 0.25)
    d1 <- ddexp(x, 0.5, 0.25)
    d2 <- ddexp(x, 0.5, rate = 4)
    expect_identical(d1, d2)

    ## nimbleFunction
    f <- nimbleFunction(
        run = function(x = double(), location = double(), par2 = double()) {
            returnType(double(1))
            out1 <- ddexp(x, location, par2)
            out2 <- ddexp(x, location, rate = 1/par2)
            return(c(out1, out2))
        })
    cf <- compileNimble(f)
    out <- cf(x, 0.5, 0.25)
    expect_identical(out, c(d1, d2))

    ## in a model
    code <- nimbleCode({
        x1 ~ ddexp(mu, scale = scale)
        x2 ~ ddexp(mu, rate)
    })
    m <- nimbleModel(code, data = list(x1 = x, x2 = x), inits = list(mu = 0.5, scale = 0.25, rate = 4))
    cm <- compileNimble(m)
    out1 <- cm$calculate('x1')
    out2 <- cm$calculate('x2')
    expect_identical(exp(c(out1, out2)), c(d1, d2))

})


test_that("recycling behavior from R and within nimbleFunctions for non-R-native distributions", {

    ## dt_nonstandard
    set.seed(1)
    param <- runif(3)
    x <- rt_nonstandard(6, 3, param, 3.5)
    expect_equal(length(x), 6)
    param <- rep(param, 2)
    d <- dt_nonstandard(x[1:3], 3, param, 3.5)
    expect_identical(d[1:3], d[4:6])
    p <- pt_nonstandard(x[1:3], 3, param, 3.5)
    q <- qt_nonstandard(p, 3, param, 3.5)
    expect_equal(rep(x[1:3], 2), q)
    expect_identical(p[1:3], p[4:6])
    expect_identical(q[1:3], q[4:6])

    ## Use of recycling
    f <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            dd <- dt_nonstandard(x, 3, theta, 3.5)
            pp <- pt_nonstandard(x, 3, theta, 3.5)
            qq <- qt_nonstandard(pp, 3, theta, 3.5)
            returnType(double(1))
            return(c(dd, pp, qq))
        })
    cf <- compileNimble(f)
    out <- cf(x[1:3], param) 
    expect_identical(out, c(d, p, q), info = 'dt_nonstandard nf')

    ## dexp_nimble
    set.seed(1)
    param <- runif(3)
    x <- rexp_nimble(6, param)
    expect_equal(length(x), 6)
    param <- rep(param, 2)
    d <- dexp_nimble(x[1:3], param)
    expect_identical(d[1:3], d[4:6])
    p <- pexp_nimble(x[1:3], param)
    q <- qexp_nimble(p, param)
    expect_equal(rep(x[1:3], 2), q)
    expect_identical(p[1:3], p[4:6])
    expect_identical(q[1:3], q[4:6])

    f <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            d <- dexp_nimble(x, theta)
            p <- pexp_nimble(x, theta)
            q <- qexp_nimble(p, theta)
            returnType(double(1))
            return(c(d, p, q))
        })
    cf <- compileNimble(f)
    out <- cf(x[1:3], param) 
    expect_identical(out, c(d, p, q), info = 'dexp_nimble nf')

    ## ddexp
    set.seed(1)
    param <- runif(3)
    x <- rdexp(6, 0.5, param)
    expect_equal(length(x), 6)
    param <- rep(param, 2)
    d <- ddexp(x[1:3], 0.5, param)
    expect_identical(d[1:3], d[4:6])
    p <- pdexp(x[1:3], 0.5, param)
    q <- qdexp(p, 0.5, param)
    expect_equal(rep(x[1:3], 2), q)
    expect_identical(p[1:3], p[4:6])
    expect_identical(q[1:3], q[4:6])

    f <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            d <- ddexp(x, 0.5, theta)
            p <- pdexp(x, 0.5, theta)
            q <- qdexp(p, 0.5, theta)
            returnType(double(1))
            return(c(d, p, q))
        })
    cf <- compileNimble(f)
    out <- cf(x[1:3], param) 
    expect_identical(out, c(d, p, q), info = 'ddexp nf')

    ## dsqrt_invgamma (no p or q functions; not available in nimbleFunction)
    set.seed(1)
    param <- runif(3)
    x <- rsqrtinvgamma(6, 0.5, param)
    expect_equal(length(x), 6)
    param <- rep(param, 2)
    d <- dsqrtinvgamma(x[1:3], 0.5, param)
    expect_identical(d[1:3], d[4:6])
    
    ## dinvgamma
    set.seed(1)
    param <- runif(3)
    x <- rinvgamma(6, 0.5, param)
    expect_equal(length(x), 6)
    param <- rep(param, 2)
    d <- dinvgamma(x[1:3], 0.5, param)
    expect_identical(d[1:3], d[4:6])
    p <- pinvgamma(x[1:3], 0.5, param)
    q <- qinvgamma(p, 0.5, param)
    expect_equal(rep(x[1:3], 2), q)
    expect_identical(p[1:3], p[4:6])
    expect_identical(q[1:3], q[4:6])

    f <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            d <- dinvgamma(x, 0.5, theta)
            p <- pinvgamma(x, 0.5, theta)
            q <- qinvgamma(p, 0.5, theta)
            returnType(double(1))
            return(c(d, p, q))
        })
    cf <- compileNimble(f)
    out <- cf(x[1:3], param) 
    expect_identical(out, c(d, p, q), info = 'dinvgamma nf')
})
    
   
sink(NULL)

if(!generatingGoldFile) {
    test_that("Log file matches gold file", {
        trialResults <- readLines(tempFileName)
        correctResults <- readLines(system.file(file.path('tests', 'testthat', goldFileName), package = 'nimble'))
        compareFilesByLine(trialResults, correctResults)
    })
}

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)

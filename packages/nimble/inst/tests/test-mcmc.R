source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of default MCMC")

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

## If you do *not* want to write to results files
##    comment out the sink() call below.  And consider setting verbose = FALSE 
## To record a new gold file, nimbleOptions('generateGoldFileForMCMCtesting') should contain the path to the directory where you want to put it
## e.g. nimbleOptions(generateGoldFileForMCMCtesting = getwd())
## Comparison to the gold file won't work until it is installed with the package.

goldFileName <- 'mcmcTestLog_Correct.Rout'
tempFileName <- 'mcmcTestLog.Rout'
generatingGoldFile <- !is.null(nimbleOptions('generateGoldFileForMCMCtesting'))
outputFile <- if(generatingGoldFile) file.path(nimbleOptions('generateGoldFileForMCMCtesting'), goldFileName) else tempFileName

## capture warnings
sink_with_messages(outputFile)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)
## tests of classic BUGS examples

test_mcmc('blocker', numItsC = 1000, resampleData = TRUE)
# 100% coverage; looks fine

test_mcmc('bones', numItsC = 10000, resampleData = TRUE)
# 100% coverage; looks fine

test_mcmc('dyes', numItsC = 1000, resampleData = TRUE)
# 100% coverage; looks fine

test_mcmc('equiv', numItsC = 1000, resampleData = TRUE)
# looks good
# testing: tau[2]=97.95, 198.8 ; tau[1]=102.2,55
# phi = -.008,.052; pi = -.1805,.052

test_mcmc('line', numItsC = 1000, resampleData = TRUE)
# 100% coverage; looks fine

test_mcmc('oxford', numItsC = 1000, resampleData = TRUE)
# probably ok; seems to overcover for 'b', but 'b' in this
# parameteriz'n is a top-level node and the multiplic'n
# by sigma seems to lead to frequentist overcoverage
# similar results in JAGS

test_mcmc('pump', numItsC = 1000, resampleData = TRUE)
# 100% coverage; looks fine

test_mcmc('rats', numItsC = 1000, resampleData = TRUE)
# 93.8% coverage; looks fine and compares well to JAGS
# however in resampleData, one of the taus wildly misses

test_mcmc('seeds', numItsC = 1000, resampleData = TRUE)
# fine

test_mcmc('dugongs', numItsC = 1000, resampleData = TRUE)
# 100% coverage; looks fine


test_mcmc('epil', model = 'epil2.bug', inits = 'epil-inits.R',
              data = 'epil-data.R', numItsC = 1000, resampleData = TRUE)
# looks ok

test_mcmc('epil', model = 'epil3.bug', inits = 'epil-inits.R',
              data = 'epil-data.R', numItsC = 1000, resampleData = TRUE)
# looks ok

test_mcmc('seeds', model = 'seedsuni.bug', inits = 'seeds-init.R',
              data = 'seeds-data.R', numItsC = 1000, resampleData = TRUE)
# looks fine - intervals for b's seem a bit large but probably ok
# particularly since default seeds.bug seems fine
# results compared to JAGS look fine
test_mcmc('seeds', model = 'seedssig.bug', inits = 'seeds-init.R',
              data = 'seeds-data.R', numItsC = 1000, resampleData = TRUE)
# looks fine - intervals for b's seem a bit large but probably ok

test_mcmc('birats', model = 'birats1.bug', inits = 'birats-inits.R',
              data = 'birats-data.R', numItsC = 1000, resampleData = TRUE)
# seems fine

test_mcmc('birats', model = 'birats3.bug', inits = 'birats-inits.R',
              data = 'birats-data.R', numItsC = 1000, resampleData = TRUE)
# seems fine

test_mcmc('birats', model = 'birats2.bug', inits = 'birats-inits.R',
            data = 'birats-data.R', numItsC = 1000, resampleData = TRUE)
# looks fine now that values() returns in order
# result changes as of v0.4 because in v0.3-1 'omega.beta' was found
# as both topNode and nontopNode and was being simulated into
# incorrectly in resampleData - this affected values further downstream

test_mcmc('ice', model = 'icear.bug', inits = 'ice-inits.R',
          data = 'ice-data.R', numItsC = 1000, resampleData = TRUE,
          knownFailures = list('coverage' = 'KNOWN ISSUE: coverage is low'))
# resampleData gives very large magnitude betas because beta[1],beta[2] are not
# actually topNodes because of (weak) dependence on tau, and
# are simulated from their priors to have large magnitude values

test_that('ice example reworked', {
                                        # rework ice example so that beta[1] and beta[2] will be top nodes
    system.in.dir(paste("sed 's/tau\\*1.0E-6/1.0E-6/g' icear.bug > ", file.path(tempdir(), "icear.bug")), dir = system.file('classic-bugs','vol2','ice', package = 'nimble'))
    test_mcmc(model = file.path(tempdir(), "icear.bug"), inits = system.file('classic-bugs', 'vol2', 'ice','ice-inits.R', package = 'nimble'), data = system.file('classic-bugs', 'vol2', 'ice','ice-data.R', package = 'nimble'), numItsC = 1000, resampleData = TRUE, avoidNestedTest = TRUE)
                                        # looks fine, but alpha and beta values shifted a bit (systematically) relative to JAGS results - on further inspection this is because mixing for this model is poor in both NIMBLE and JAGS - with longer runs they seem to agree (as best as one can tell given the mixing without doing a super long run)
})

test_mcmc('beetles', model = 'beetles-logit.bug', inits = 'beetles-inits.R',
          data = 'beetles-data.R', numItsC = 1000, resampleData = TRUE)
                                        # getting warning; deterministic model node is NA or NaN in model initialization
                                        # weirdness with llike.sat[8] being NaN on init (actually that makes sense), and with weird lifting of RHS of llike.sat

test_that('leuk example setup', {
    writeLines(c("var","Y[N,T],","dN[N,T];"), con = file.path(tempdir(), "leuk.bug")) ## echo doesn't seem to work on Windows
                                        # need nimStep in data block as we no longer have step
    system.in.dir(paste("cat leuk.bug >> ", file.path(tempdir(), "leuk.bug")), dir = system.file('classic-bugs','vol1','leuk',package = 'nimble'))
                                        # need nimStep in data block as we no longer have step
    system.in.dir(paste("sed -i -e 's/step/nimStep/g'", file.path(tempdir(), "leuk.bug")))
    
    test_mcmc(model = file.path(tempdir(), "leuk.bug"), name = 'leuk', inits = system.file('classic-bugs', 'vol1', 'leuk','leuk-init.R', package = 'nimble'), data = system.file('classic-bugs', 'vol1', 'leuk','leuk-data.R', package = 'nimble'), numItsC = 1000,
              results = list(mean = list(beta = 1.58), sd = list(beta = 0.43)),
              resultsTolerance = list(mean = list(beta = 0.02), sd = list(beta = 0.02)), avoidNestedTest = TRUE)
})

test_that('salm example setup', {
    writeLines(paste("var","logx[doses];"), con = file.path(tempdir(), "salm.bug"))
    system.in.dir(paste("cat salm.bug >>", file.path(tempdir(), "salm.bug")), dir = system.file('classic-bugs','vol1','salm', package = 'nimble'))
    test_mcmc(model = file.path(tempdir(), "salm.bug"), name = 'salm', inits = system.file('classic-bugs', 'vol1', 'salm','salm-init.R', package = 'nimble'), data = system.file('classic-bugs', 'vol1', 'salm','salm-data.R', package = 'nimble'), numItsC = 1000, avoidNestedTest = TRUE)
                                        # looks good compared to JAGS
})

test_that('air example setup', {
    file.copy(system.file('classic-bugs','vol2','air','air.bug', package = 'nimble'), file.path(tempdir(), "air.bug"), overwrite=TRUE)
    system.in.dir(paste("sed -i -e 's/mean(X)/mean(X\\[\\])/g'", file.path(tempdir(), "air.bug")))
    test_mcmc(model = file.path(tempdir(), "air.bug"), name = 'air', inits = system.file('classic-bugs', 'vol2', 'air','air-inits.R', package = 'nimble'), data = system.file('classic-bugs', 'vol2', 'air','air-data.R', package = 'nimble'), numItsC = 1000, avoidNestedTest = TRUE)
                                        # theta[2] posterior is a bit off from JAGS - would be worth more investigation
})

test_that('jaw-linear setup', {
    system.in.dir(paste("sed 's/mean(age)/mean(age\\[1:M\\])/g' jaw-linear.bug > ", file.path(tempdir(), "jaw-linear.bug")), dir = system.file('classic-bugs','vol2','jaw', package = 'nimble')) # alternative way to get size info in there
    test_mcmc(model = file.path(tempdir(), "jaw-linear.bug"), name = 'jaw-linear', inits = system.file('classic-bugs', 'vol2', 'jaw','jaw-inits.R', package = 'nimble'), data = system.file('classic-bugs', 'vol2', 'jaw','jaw-data.R', package = 'nimble'), numItsC = 1000, avoidNestedTest = TRUE) # , knownFailures = list('R MCMC' = 'Cholesky of NA matrix fails in R 3.4.2 in calculate(model) of initializeModel() but not in R 3.4.1'))
})
## note R MCMC used to fail when tried to do Cholesky of 0 matrix in 2-point method, but no longer doing multiplicative link for Wishart targets
                                      
test_mcmc('pump',
          resampleData = TRUE,
          results = list(mean = list(
                             "theta[1]" = 0.06,
                             "theta[2]" = 0.10,
                             "theta[9]" = 1.58,
                             "theta[10]" = 1.97,
                             alpha = 0.73,
                             beta = 0.98)),
          resultsTolerance = list(mean = list(
                                      "theta[1]" = 0.01,
                                      "theta[2]" = 0.01,
                                      "theta[9]" = 0.05,
                                      "theta[10]" = 0.05,
                                      alpha = 0.1,
                                      beta = 0.1)))


test_that('gap setup', {
    ## LogProb gap: bug fixed in after v0.3
    ## Problem that occurred in v0.3: because of gap in logProb_a (i.e. logProb_a[2]
    ## is defined but logProb_a[1] is not)
    ## Because logProbs get scrambled, the random walk sampler would always accept,
    ## meaning the sd of proposal steps approaches Inf
    gapCode <- nimbleCode({
	a[1] <- 1
	a[2] ~ dnorm(0,1)
    })
    
    test_mcmc(model = gapCode, seed = 0, numItsC = 100000,
              results = list(mean = list(`a[2]` = 0) ),
              resultsTolerance = list(mean = list(`a[2]` = 0.1)),
              samplers = list(list(type = 'RW', target = 'a[2]')),
              avoidNestedTest = TRUE
              )
})

if(.Platform$OS.type == 'windows') {
    message("Stopping tests now in Windows to avoid crashing until we can unload compiled projects")
    message("To continue testing use 'mcmc2' tests")
    q("no")
}

### Daniel's world's simplest MCMC demo

test_that('very simple example setup', {
    code <- nimbleCode({
        x ~ dnorm(0, 2)
        y ~ dnorm(x+1, 3)
        z ~ dnorm(y+2, 4)
    })
    data = list(y = 3)
    
    test_mcmc(model = code,
              name = 'very simple example',
              data = data,
              resampleData = FALSE,
              results = list(
                  mean = list(x = 6/5, z = 5),
                  sd = list(x = 1/sqrt(5), z = 1/2)),
              resultsTolerance = list(mean = list(x = .1, z = .1),
                                      sd = list(x = .05, z = .05)),
              avoidNestedTest = TRUE)
})

### basic block sampler example

test_that('basic no-block sampler setup', {
    code <- nimbleCode({
        for(i in 1:3) {
            x[i] ~ dnorm(0, 1)
            y[i] ~ dnorm(x[i], 2)
        }
    })
    data = list(y = -1:1)
    
    test_mcmc(model = code,
              name = 'basic no-block sampler',
              data = data,
              resampleData = FALSE,
              results = list(
                  mean = list(x = c(-2/3,0,2/3)),
                  var = list(x = rep(1/3,3))),
              resultsTolerance = list(mean = list(x = rep(.1,3)),
                                      var = list(x = rep(.05,3))),
              avoidNestedTest = TRUE)
    
    
    test_mcmc(model = code,
              name = 'basic block sampler on scalars',
              data = data, resampleData = FALSE,
              results = list(
                  mean = list(x = c(-2/3,0,2/3)),
                  var = list(x = rep(1/3,3))),
              resultsTolerance = list(mean = list(x = rep(.1,3)),
                                      var = list(x = rep(.05,3))),
              samplers = list(
                  list(type = 'RW_block', target = 'x[1]'),
                  list(type = 'RW_block', target = 'x[2]'),
                  list(type = 'RW_block', target = 'x[3]')
              ), removeAllDefaultSamplers = TRUE, numItsC = 10000,
              avoidNestedTest = TRUE)
    
    test_mcmc(model = code,
              name = 'basic block sampler on vector',
              data = data,
              resampleData = FALSE,
              results = list(
                  mean = list(x = c(-2/3,0,2/3)),
                  var = list(x = rep(1/3,3))),
              resultsTolerance = list(mean = list(x = rep(.1,3)),
                                      var = list(x = rep(.05,3))),
              samplers = list(
                  list(type = 'RW_block', target = 'x', control = list(adaptInterval = 500))
              ), numItsC = 10000, avoidNestedTest = TRUE)  
})

### slice sampler example

test_that('slice sampler example setup', {
    code <- nimbleCode({
        z ~ dnorm(0, 1)
        normal5_10 ~ dnorm(5, sd = 10)
        beta1_1 ~ dbeta(1, 1)
        beta3_5 ~ dbeta(3, 5)
        binom10_p5 ~ dbin(size=10, prob=0.5)
        binom20_p3 ~ dbin(size=20, prob=0.3)
    })
    
    test_mcmc(model = code,
              name = "slice sampler example",
              resampleData = FALSE,
              results = list(
                  mean = list(z = 0, "beta1_1" = 0.5, "beta3_5" = 3/(3+5),
                              "binom10_p5" = 10*.5, "binom20_p3" = 20*.3),
                  sd = list(z = 1, "beta1_1" = sqrt(1/12),
                            "beta3_5" = sqrt(3*5/((3+5)^2*(3+5+1))),
                            "binom10_p5" = sqrt(10*.5*.5),
                            "binom20_p3" = sqrt(20*.3*.7))),
              resultsTolerance = list(
                  mean = list(z = 0.1, "beta1_1" = 0.5, "beta3_5" = .2,
                              "binom10_p5" = .25, "binom20_p3" = .25),
                  sd = list(z = .1, "beta1_1" = .05, "beta3_5" = .03,
                            "binom10_p5" = .2, "binom20_p3" = .25)),
              samplers = list(list(type = 'slice', target = 'z', control = list(adaptInterval = 10)),
                              list(type = 'slice', target = 'normal5_10', control = list(adaptInterval = 10)),
                              list(type = 'slice', target = 'beta1_1', control = list(adaptInterval = 10)),
                              list(type = 'slice', target = 'beta3_5', control = list(adaptInterval = 10)),
                              list(type = 'slice', target = 'binom10_p5', control = list(adaptInterval = 10)),
                              list(type = 'slice', target = 'binom20_p3', control = list(adaptInterval = 10))),
              avoidNestedTest = TRUE)
    })


### elliptical slice sampler 'ess'

test_that('elliptical slice sampler setup', {
    set.seed(0)
    ESScode <- quote({
        x[1:d] ~ dmnorm(mu_x[1:d], prec = prec_x[1:d, 1:d])
        y[1:d] ~ dmnorm(x[1:d], prec = prec_y[1:d, 1:d])
    })
    d <- 3
    mu_x <- rnorm(d)
    temp <- array(rnorm(d^2), c(d,d))
    prec_x <- solve(temp %*% t(temp))
    temp <- array(rnorm(d^2), c(d,d))
    prec_y <- solve(temp %*% t(temp))
    y <- rnorm(d)
    ESSconstants <- list(d = d, mu_x = mu_x, prec_x = prec_x, prec_y = prec_y)
    ESSdata <- list(y = y)
    ESSinits <- list(x = rep(0, d))
    
    test_mcmc(model = ESScode, data = c(ESSconstants, ESSdata), inits = ESSinits,
              name = 'exact values of elliptical slice sampler',
              seed = 0,
              exactSample = list(`x[1]` = c(-0.492880566939352, -0.214539223107114, 1.79345037297218, 1.17324496091208, 2.14095077672555, 1.60417482445964, 1.94196916651627, 2.66737323347255, 2.66744178776022, 0.253966883192744), `x[2]` = c(-0.161210109217102, -0.0726534676226932, 0.338308532423757, -0.823652445515156, -0.344130712698579, -0.132642244861469, -0.0253168895009594, 0.0701624130921676, 0.0796842215444978, -0.66369112443311), `x[3]` = c(0.278627475932455, 0.0661336950029345, 0.407055002920732, 1.98761228946318, 1.0839897275519, 1.00262648370199, 0.459841485268785, 2.59229443025387, 1.83769567435409, 1.92954706515119)),
              samplers = list(list(type = 'ess', target = 'x')))
    
    test_mcmc(model = ESScode, data = c(ESSconstants, ESSdata), inits = ESSinits,
              name = 'results to tolerance of elliptical slice sampler',
              results = list(mean = list(x = c(1.0216463, -0.4007247, 1.1416904))),
              resultsTolerance = list(mean = list(x = c(0.01, 0.01, 0.01))),
              numItsC = 100000,
              samplers = list(list(type = 'ess', target = 'x')), avoidNestedTest = TRUE)
    })


### demo2 of check conjugacy

test_that('beta-binom conjugacy setup', {
    code <- nimbleCode({
        x ~ dbeta(3, 13)
        y[1] ~ dbin(x, 10)
        y[2] ~ dbin(x, 20)
    })
    data = list(y = c(3,4))
    
    test_mcmc(model = code, name = 'check of beta-binom conjugacy', data = data, exactSample = list(x = c(0.195510839527966, 0.332847482503424,0.247768152764931, 0.121748195439553, 0.157842271774841, 0.197566496350904, 0.216991517500577, 0.276609942874852, 0.165733872345582, 0.144695512780252)), seed = 0, avoidNestedTest = TRUE)
})
### checkConjugacy_demo3_run.R - various conjugacies

test_that('various conjugacies setup', {
    code <- nimbleCode({
        x ~ dgamma(1, 1)       # should satisfy 'gamma' conjugacy class
        a  ~ dnorm(0, x)     # should satisfy 'norm' conjugacy class
        a2 ~ dnorm(0, tau = 3*x+0)
        b  ~ dpois(0+5*x)
        b2 ~ dpois(1*x*1)
        c ~ dgamma(1, 7*x*5)
        for(i in 2:3) {
            jTau[i] <- 1
            jNorm[i] ~ dnorm(c * (a+3) - i, var = jTau[i])
            kTauSd[i] <- 2
            kLogNorm[i] ~ dlnorm(0 - a - 6*i, kTauSd[i])
        }
        jNorm[1] <- 0
        kLogNorm[1] <- 0
    })

    sampleVals = list(x = c(3.950556165467749, 1.556947815895538, 1.371834934033851, 2.036442813764752, 2.247416118159410, 2.537131924778210, 2.382184991769738, 2.653737836857812, 2.934255734970981, 3.007873553270551),
                      c = c(0.010341199485849559, 0.010341199485849559, 0.003846483017887228, 0.003846483017887228, 0.003846483017887228, 0.006269117826484087, 0.009183580181658716, 0.009183580181658716, 0.006361841408434201, 0.006361841408434201))
    
    test_mcmc(model = code, name = 'check various conjugacies', exactSample = sampleVals, seed = 0, mcmcControl = list(scale=0.01), avoidNestedTest = TRUE)
    ## with fixing of jNorm[1] and kLogNorm[1] we no longer have: knownFailures = list('R C samples match' = "KNOWN ISSUE: R and C posterior samples are not equal for 'various conjugacies'"))
})

### Weibull-gamma conjugacy
test_that('Weibull-gamma conjugacy setup', {
    y <- 3; depShape <- 2; c <- 2; shape <- 1; rate <- 2
    code <- nimbleCode({
        y ~ dweib(shape = depShape, lambda = c*theta)
        theta ~ dgamma(shape, rate = rate)
    })
    m <- nimbleModel(code, data = list(y = y), inits = list(c = c, theta = 1),
                     constants = list(depShape = depShape, shape = shape, rate = rate))
    conf <- configureMCMC(m)
    samplers <- conf$getSamplers()
    expect_identical(samplers[[1]]$name, 'conjugate_dgamma_dweib',
                                   info = "dweibull-dgamma conjugacy with dependency using lambda not detected")
    mcmc <- buildMCMC(conf)
    comp <- compileNimble(m, mcmc)
    set.seed(0)
    comp$mcmc$run(10)
    smp <- as.matrix(comp$mcmc$mvSamples)

    manualSampler <- function(n, y, depShape, c, shape, rate) {
        out <- rep(0, n)
        shape = shape + 1
        rate = rate + c*y^depShape
        set.seed(0)
        out <- rgamma(n, shape, rate = rate)
        return(out)
    }
    smpMan <- manualSampler(10, y, depShape, c, shape, rate)

    expect_identical(smp[,1], smpMan,
                                   info = "NIMBLE gamma-Weibull conjugate sampler and manual sampler results differ")
})    

test_that('Dirichlet-multinomial conjugacy setup', {
### Dirichlet-multinomial conjugacy

# as of v0.4, exact numerical results here have changed because
# ddirch now sometimes returns NaN rather than -Inf (when an
# alpha is proposed to be negative) -- this changes the RNG
# sequence because NaN values result in no runif() call in decide()

# single multinomial
    set.seed(0)
    n <- 100
    alpha <- c(10, 30, 15, 60, 1)
    K <- length(alpha)
    p <- c(.12, .24, .09, .54, .01)
    y <- rmulti(1, n, p)
    
    code <- function() {
        y[1:K] ~ dmulti(p[1:K], n);
        p[1:K] ~ ddirch(alpha[1:K]);
        for(i in 1:K) {
            alpha[i] ~ dgamma(.001, .001);
        }
    }
    
    inits <- list(p = rep(1/K, K), alpha = rep(K, K))
    data <- list(n = n, K = K, y = y)
    
    test_mcmc(model = code, name = 'Dirichlet-multinomial example', data= data, seed = 0, numItsC = 10000,
              inits = inits,
              results = list(mean = list(p = p)),
              resultsTolerance = list(mean = list(p = rep(.06, K))), avoidNestedTest = TRUE)
})
## bad mixing for alphas; probably explains why posterior estimates for alphas changed so much as of v 0.4
  
## with replication
test_that('Dirichlet-multinomial with replication setup', {
    set.seed(0)
    n <- 100
    m <- 20
    alpha <- c(10, 30, 15, 60, 1)
    K <- length(alpha)
    y <- p <- matrix(0, m, K)
    for(i in 1:m) {
        p[i, ] <- rdirch(1, alpha)
        y[i, ] <- rmulti(1, n, p[i, ])
    }
    
    code <- function() {
        for(i in 1:m) {
            y[i, 1:K] ~ dmulti(p[i, 1:K], n);
            p[i, 1:K] ~ ddirch(alpha[1:K]);
        }
        for(i in 1:K) {
            alpha[i] ~ dgamma(.001, .001);
        }
    }
    
    inits <- list(p = matrix(1/K, m, K), alpha = rep(1/K, K))
    data <- list(n = n, K = K, m = m, y = y)

    ## two tolerance failures are known, for p[39] and p[76]
    test_mcmc(model = code, name = 'Dirichlet-multinomial with replication', data= data,
              seed = 0, numItsC = 1000,
              inits = inits, numItsC_results = 100000,
              results = list(mean = list(p = p, alpha = alpha)),
              resultsTolerance = list(mean = list(p = matrix(.05, m, K),
                                                  alpha = c(5,10,10,20,.5))),
              knownFailures = list('MCMC match to known posterior: p mean 39' = 'KNOWN ISSUE: two samples outside resultsTolerance',
                                  'MCMC match to known posterior: p mean 76' = 'KNOWN ISSUE: two samples outside resultsTolerance'), avoidNestedTest = TRUE)
})
# note alphas mix poorly (and are highly correlated),
# presumably because of cross-level dependence between
# p's and alphas.  cross-level sampler would probably work well here,
# or, of course, integrating over the p's

test_that('Dirichlet-categorical conjugacy setup', {
### Dirichlet-categorical conjugacy

# single multinomial represented as categorical
    set.seed(0)
    n <- 100
    alpha <- c(10, 30, 15, 60, 1)
    K <- length(alpha)
    p <- c(.12, .24, .09, .54, .01)
    y <- rmulti(1, n, p)
    y <- rep(seq_along(y), times = y)
    
    code <- function() {
        for(i in 1:n)
          y[i] ~ dcat(p[1:K])
        p[1:K] ~ ddirch(alpha[1:K])
        for(i in 1:K) {
            alpha[i] ~ dgamma(.001, .001);
        }
    }
    
    inits <- list(p = rep(1/K, K), alpha = rep(K, K))
    data <- list(n = n, K = K, y = y)
    
    test_mcmc(model = code, name = 'Dirichlet-categorical example', data= data, seed = 0, numItsC = 10000,
              inits = inits,
              results = list(mean = list(p = p)),
              resultsTolerance = list(mean = list(p = rep(.06, K))), avoidNestedTest = TRUE)
})

## also note that MCMC results here should be identical to those from
## Dirichlet-multinomial case two tests up from this


### block sampler on MVN node
test_that('block sampler on MVN node setup', {
    code <- nimbleCode({
        mu[1] <- 10
        mu[2] <- 20
        mu[3] <- 30
        x[1:3] ~ dmnorm(mu[1:3], prec = Q[1:3,1:3])
    })
    
    Q = matrix(c(1.0,0.2,-1.0,0.2,4.04,1.6,-1.0,1.6,10.81), nrow=3)
    data = list(Q = Q)
    inits = list(x = c(10, 20, 30))
    
    test_mcmc(model = code, name = 'block sampler on multivariate node', data = data, seed = 0, numItsC = 10000,
              results = list(mean = list(x = c(10,20,30)),
                             var = list(x = diag(solve(Q)))),
              resultsTolerance = list(mean = list(x = rep(1,3)),
                                      var = list(x = c(.1, .03, .01))),
              samplers = list(
                  list(type = 'RW_block', target = 'x[1:3]')), avoidNestedTest = TRUE)
                                        # caution: setting targetNodes='x' works but the initial end sampler is not removed because x[1:3] in targetNode in default sampler != 'x' in targetNodes passed in
    if(FALSE) {
        Rmodel <- nimbleModel(code, constants = list(Q=Q))
        mcmcspec <- MCMCspec(Rmodel, nodes = NULL)
        mcmcspec$addSampler(type = 'RW_block', target = 'x', control = list(adaptInterval=500))
        mcmcspec$getMonitors()
        Rmcmc <- buildMCMC(mcmcspec)
        Cmodel <- compileNimble(Rmodel)
        Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
        Cmcmc(200000)    ## this runs nearly instantaneously on my computer -DT
        samples <- as.matrix(nfVar(Cmcmc, 'mvSamples'))
        samples <- samples[50001:200000,]
        dim(samples)
        apply(samples, 2, mean)
        solve(Q)
        cov(samples)
        propCov <- nfVar(Cmcmc, 'samplerFunctions')[[1]]$propCov
        scale <- nfVar(Cmcmc, 'samplerFunctions')[[1]]$scale
        propCov * scale^2
        
        nfVar(Cmcmc, 'samplerFunctions')[[1]]$scaleHistory
        nfVar(Cmcmc, 'samplerFunctions')[[1]]$acceptanceRateHistory
        nfVar(Cmcmc, 'samplerFunctions')[[1]]$scale
        nfVar(Cmcmc, 'samplerFunctions')[[1]]$propCov
        ## why is the proposal cov w/ .99 cross-corrs?
        ## also MCMC in C takes a surprisingly long time - this might be threaded lin alg behaving badly on small matrices
    }
})

test_that('second block sampler on multivariate node', {
### DT's model
    mu <- c(1,2,3)
    corr <- matrix(c(1,.8,0.3,.8,1,0,0.3,0,1), nrow=3)
    varr <- c(1,2,3)
    Sig <- diag(sqrt(varr))
    Q <- Sig %*% corr %*% Sig
    P <- solve(Q)
    
    code <- nimbleCode({
        x[1:3] ~ dmnorm(mu[1:3], prec = P[1:3,1:3])
    })
    data = list(P = P, mu = mu)
    
    test_mcmc(model = code, name = 'second block sampler on multivariate node', data = data, seed = 0, numItsC = 100000,
              results = list(mean = list(x = mu),
                             var = list(x = varr)),
              resultsTolerance = list(mean = list(x = rep(.1,3)),
                                      var = list(x = c(.1,.1,.1))),
              samplers = list(
                  list(type = 'RW_block', target = 'x[1:3]')), avoidNestedTest = TRUE)
    })

test_that('MVN conjugate setup', {
### MVN conjugate update
    
    set.seed(0)
    mu0 = 1:3
    Q0 = matrix(c(1, .2, .8, .2, 2, 1, .8, 1, 2), nrow = 3)
    Q = solve(matrix(c(3, 1.7, .9, 1.7, 2, .6, .9, .6, 1), nrow = 3))
    a = c(-2, .5, 1)
    B = matrix(rnorm(9), 3)
    
##### not currently working - see Perry's email of ~ 10/6/14
    ## code <- nimbleCode({
    ##   mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
    ##   y[1:3] ~ dmnorm(asCol(a[1:3]) + B[1:3, 1:3] %*% asCol(mu[1:3]), Q[1:3, 1:3])
    ## })
    
    code <- nimbleCode({
        mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
        y_mean[1:3] <- asCol(a[1:3]) + B[1:3, 1:3] %*% asCol(mu[1:3])
        y[1:3] ~ dmnorm(y_mean[1:3], Q[1:3, 1:3])
    })
    
    
    mu <- mu0 + chol(solve(Q0)) %*% rnorm(3)
                                        # make sure y is a vec not a 1-col matrix or get a dimensionality error
    y <- c(a + B%*%mu + chol(solve(Q)) %*% rnorm(3))
    data = list(mu0 = mu0, Q0 = Q0, Q = Q, a = a, B = B, y = y)
    
    muQtrue = t(B) %*% Q%*%B + Q0
    muMeanTrue = c(solve(muQtrue, crossprod(B, Q%*%(y-a)) + Q0%*%mu0))
    
    test_mcmc(model = code, name = 'two-level multivariate normal', data = data, seed = 0, numItsC = 10000,
              results = list(mean = list(mu = muMeanTrue),
                             cov = list(mu = solve(muQtrue))),
              resultsTolerance = list(mean = list(mu = rep(.02,3)),
                                      cov = list(mu = matrix(.01, 3, 3))), avoidNestedTest = TRUE)
    

### scalar RW updates in place of conjugate mv update

    test_mcmc(model = code, name = 'two-level multivariate normal with scalar updaters', data = data, seed = 0, numItsC = 100000,
              results = list(mean = list(mu = muMeanTrue),
                             cov = list(mu = solve(muQtrue))),
              resultsTolerance = list(mean = list(mu = rep(.03,3)),
                                      cov = list(mu = matrix(.03, 3, 3))),
              samplers = list(list(type = 'RW', target = 'mu[1]'),
                              list(type = 'RW', target = 'mu[2]'),
                              list(type = 'RW', target = 'mu[3]')),
              removeAllDefaultSamplers = TRUE, avoidNestedTest = TRUE)
    
})


## another example of MVN conjugate sampler, for test-mcmc.R
## using both cov and prec parametrizaions of MVN,
## and various linear links
test_that('another MVN conjugate sampler setup', {
    set.seed(0)
    prior_mean <- rep(0,5)
    tmp <- array(rnorm(25), c(5,5))
    tmp <- tmp + t(tmp) + 5*diag(5)
    prior_cov <- tmp
    a <- array(rnorm(20), c(4,5))
    B <- array(NA, c(4,5,5))
    for(i in c(2,4))   B[i,,] <- array(rnorm(25), c(5,5))
    B[1,,] <- diag(5)
    B[3,,] <- diag(5)
    M_y <- array(NA, c(4,5,5))
    for(i in 1:4) {
        tmp <- array(rnorm(25,i), c(5,5))
        tmp <- tmp + t(tmp) + 5*i*diag(5)
        M_y[i,,] <- tmp
    }
    x <- rep(0, 5)
    y <- array(rnorm(20), c(4,5))
    
    code <- nimbleCode({
        x[1:5] ~ dmnorm(mean = prior_mean[1:5], cov = prior_cov[1:5,1:5])
        for(i in 1:4)
            mu_y[i,1:5] <- asCol(a[i,1:5]) + B[i,1:5,1:5] %*% asCol(x[1:5])
        y[1,1:5] ~ dmnorm(mu_y[1,1:5], prec = M_y[1,1:5,1:5])
        y[2,1:5] ~ dmnorm(mu_y[2,1:5], cov  = M_y[2,1:5,1:5])
        y[3,1:5] ~ dmnorm(mu_y[3,1:5], prec = M_y[3,1:5,1:5])
        y[4,1:5] ~ dmnorm(mu_y[4,1:5], cov  = M_y[4,1:5,1:5])
    })
    constants <- list(prior_mean=prior_mean, prior_cov=prior_cov, a=a, B=B, M_y=M_y)
    data <- list(y=y)
    inits <- list(x=x)
    Rmodel <- nimbleModel(code, constants, data, inits)
    spec <- configureMCMC(Rmodel)
    Rmcmc <- buildMCMC(spec)
    
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    
    set.seed(0)
    Rmcmc$run(10)
    Rsamples <- as.matrix(Rmcmc$mvSamples)
    set.seed(0)
    Cmcmc$run(10)
    Csamples <- as.matrix(Cmcmc$mvSamples)

    expect_equal(round(as.numeric(Rsamples), 8),
                 c(0.97473128, 0.50438666, 1.1251132, 0.83830666, 0.74077066, 0.92935482, 0.83758372, 0.98708273, 1.24199937, 0.67348127, -0.54387714, -0.60713969, -0.51392796, -0.3176801, -0.34416529, -0.08530564, -0.47160157, -0.21996584, -0.20504917, -0.77287122, 0.78462584, 0.46103509, 0.43862813, 0.49343096, 0.61020864, 0.55088287, 0.53887202, 0.49863894, 0.62691318, 0.80142839, 0.34941152, 0.06623608, 0.05624477, 0.21369178, 0.26585415, -0.1439989, -0.03133488, 0.3544062, -0.03518959, 0.27415746, 0.40977, 0.8351078, 0.25719293, 0.05663917, 0.30894028, 0.33113315, 0.47647909, 0.26143962, 0.07180759, 0.27255767),
                 info = 'R sample not correct compared to known result')
    
    dif <- as.numeric(Rsamples - Csamples)
    expect_true(max(abs(dif)) < 1E-15, info = 'R and C equiv')
    
    y_prec <- array(NA, c(4,5,5))
    y_prec[1,,] <-       M_y[1,,]
    y_prec[2,,] <- solve(M_y[2,,])
    y_prec[3,,] <-       M_y[3,,]
    y_prec[4,,] <- solve(M_y[4,,])
    contribution_mean <- array(NA, c(4,5))
    for(i in 1:4)   contribution_mean[i,] <- t(B[i,,]) %*% y_prec[i,,] %*% (y[i,] - a[i,])
    contribution_prec <- array(NA, c(4,5,5))
    for(i in 1:4)   contribution_prec[i,,] <- t(B[i,,]) %*% y_prec[i,,] %*% B[i,,]
    prior_prec <- solve(prior_cov)
    post_prec <- prior_prec + apply(contribution_prec, c(2,3), sum)
    post_cov <- solve(post_prec)
    post_mean <- (post_cov %*% (prior_prec %*% prior_mean + apply(contribution_mean, 2, sum)))[,1]
    
    Cmcmc$run(100000)
    Csamples <- as.matrix(Cmcmc$mvSamples)
    
    dif_mean <- as.numeric(apply(Csamples, 2, mean)) - post_mean
    expect_true(max(abs(dif_mean)) < 0.001, info = 'posterior mean')
    
    dif_cov <- as.numeric(cov(Csamples) - post_cov)
    expect_true(max(abs(dif_cov)) < 0.001, info = 'posterior cov')
})


### test of conjugate Wishart
test_that('conjugate Wishart setup', {
    set.seed(0)
    
    trueCor <- matrix(c(1, .3, .7, .3, 1, -0.2, .7, -0.2, 1), 3)
    covs <- c(3, 2, .5)
    
    trueCov = diag(sqrt(covs)) %*% trueCor %*% diag(sqrt(covs))
    Omega = solve(trueCov)
    
    n = 20
    R = diag(rep(1,3))
    mu = 1:3
    Y = mu + t(chol(trueCov)) %*% matrix(rnorm(3*n), ncol = n)
    M = 3
    data <- list(Y = t(Y), n = n, M = M, mu = mu, R = R)
    
    code <- nimbleCode( {
        for(i in 1:n) {
            Y[i, 1:M] ~ dmnorm(mu[1:M], Omega[1:M,1:M]);
        }
        Omega[1:M,1:M] ~ dwish(R[1:M,1:M], 4);
    })
    
    newDf = 4 + n
    newR = R + tcrossprod(Y- mu)
    OmegaTrueMean = newDf * solve(newR)
    
    wishRV <- array(0, c(M, M, 10000))
    for(i in 1:10000) {
        z <- t(chol(solve(newR))) %*% matrix(rnorm(3*newDf), ncol = newDf)
        wishRV[ , , i] <- tcrossprod(z)
    }
    OmegaSimTrueSDs = apply(wishRV, c(1,2), sd)
    
    test_mcmc(model = code, name = 'conjugate Wishart', data = data, seed = 0, numItsC = 1000, inits = list(Omega = OmegaTrueMean),
              results = list(mean = list(Omega = OmegaTrueMean ),
                             sd = list(Omega = OmegaSimTrueSDs)),
              resultsTolerance = list(mean = list(Omega = matrix(.05, M,M)),
                                      sd = list(Omega = matrix(0.06, M, M))), avoidNestedTest = TRUE)
                                        # issue with Chol in R MCMC - probably same issue as in jaw-linear
    
})

test_that('conjugate Wishart setup with scaling', {
    set.seed(0)
    
    trueCor <- matrix(c(1, .3, .7, .3, 1, -0.2, .7, -0.2, 1), 3)
    covs <- c(3, 2, .5)
    tau <- 4
    
    trueCov = diag(sqrt(covs)) %*% trueCor %*% diag(sqrt(covs))
    Omega = solve(trueCov) / tau
    
    n = 20
    R = diag(rep(1,3))
    mu = 1:3
    Y = mu + t(chol(trueCov)) %*% matrix(rnorm(3*n), ncol = n)
    M = 3
    data <- list(Y = t(Y), n = n, M = M, mu = mu, R = R)
    
    code <- nimbleCode( {
        for(i in 1:n) {
            Y[i, 1:M] ~ dmnorm(mu[1:M], tmp[1:M,1:M])
        }
        tmp[1:M,1:M] <- tau * Omega[1:M,1:M]
        Omega[1:M,1:M] ~ dwish(R[1:M,1:M], 4);
    })
    
    newDf = 4 + n
    newR = R + tcrossprod(Y - mu)*tau
    OmegaTrueMean = newDf * solve(newR)
    
    wishRV <- array(0, c(M, M, 10000))
    for(i in 1:10000) {
        z <- t(chol(solve(newR))) %*% matrix(rnorm(3*newDf), ncol = newDf)
        wishRV[ , , i] <- tcrossprod(z)
    }
    OmegaSimTrueSDs = apply(wishRV, c(1,2), sd)

    m <- nimbleModel(code, data = data[1],constants=data[2:5],inits = list(Omega = OmegaTrueMean, tau = tau))
    conf <- configureMCMC(m)
    expect_equal(conf$getSamplers()[[1]]$name, "conjugate_dwish_dmnorm", info = "conjugate dmnorm-dwish with scaling not detected")
    
    test_mcmc(model = code, name = 'conjugate Wishart, scaled', data = data, seed = 0, numItsC = 1000, inits = list(Omega = OmegaTrueMean, tau = tau),
              results = list(mean = list(Omega = OmegaTrueMean ),
                             sd = list(Omega = OmegaSimTrueSDs)),
              resultsTolerance = list(mean = list(Omega = matrix(.05, M,M)),
                                      sd = list(Omega = matrix(0.06, M, M))), avoidNestedTest = TRUE)
                                        # issue with Chol in R MCMC - probably same issue as in jaw-linear
    
})

test_that('using RW_wishart sampler on non-conjugate Wishart node', {
    set.seed(0)
    trueCor <- matrix(c(1, .3, .7, .3, 1, -0.2, .7, -0.2, 1), 3)
    covs <- c(3, 2, .5)
    trueCov = diag(sqrt(covs)) %*% trueCor %*% diag(sqrt(covs))
    Omega = solve(trueCov)
    n = 20
    R = diag(rep(1,3))
    mu = 1:3
    Y = mu + t(chol(trueCov)) %*% matrix(rnorm(3*n), ncol = n)
    M = 3
    data <- list(Y = t(Y), n = n, M = M, mu = mu, R = R)
    code <- nimbleCode( {
        for(i in 1:n) {
            Y[i, 1:M] ~ dmnorm(mu[1:M], Omega[1:M,1:M])
        }
        Omega[1:M,1:M] ~ dwish(R[1:M,1:M], 4)
    })
    newDf = 4 + n
    newR = R + tcrossprod(Y- mu)
    OmegaTrueMean = newDf * solve(newR)
    wishRV <- array(0, c(M, M, 10000))
    for(i in 1:10000) {
        z <- t(chol(solve(newR))) %*% matrix(rnorm(3*newDf), ncol = newDf)
        wishRV[ , , i] <- tcrossprod(z)
    }
    OmegaSimTrueSDs = apply(wishRV, c(1,2), sd)
    allData <- data
    constants <- list(n = allData$n, M = allData$M, mu = allData$mu, R = allData$R)
    data <- list(Y = allData$Y)
    inits <- list(Omega = OmegaTrueMean)

    Rmodel <- nimbleModel(code, constants, data, inits)
    Rmodel$calculate()
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('Omega', 'RW_wishart')
    Rmcmc <- buildMCMC(conf)
    compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
    Cmodel <- compiledList$model; Cmcmc <- compiledList$mcmc
    set.seed(0)
    samples <- runMCMC(Cmcmc, 10000)
    d1 <- as.numeric(apply(samples, 2, mean)) - as.numeric(OmegaTrueMean)
    difference <- sum(round(d1,9) - c(0.024469145, -0.011872571, -0.045035297, -0.011872571, -0.003443918, 0.009363410, -0.045035297,  0.009363410,  0.049971420))
    expect_true(difference == 0)
})


test_that('using RW_wishart sampler on inverse-Wishart distribution', {
    code <- nimbleCode( {
        for(i in 1:n) {
            Y[i, 1:M] ~ dmnorm(mu[1:M], cov = C[1:M,1:M])
        }
        C[1:M,1:M] ~ dinvwish(R[1:M,1:M], 4)
    })
    
    set.seed(0)
    trueCor <- matrix(c(1, .3, .7, .3, 1, -0.2, .7, -0.2, 1), 3)
    covs <- c(3, 2, .5)
    trueCov = diag(sqrt(covs)) %*% trueCor %*% diag(sqrt(covs))
    n = 20
    R = diag(rep(1,3))
    mu = 1:3
    Y = mu + t(chol(trueCov)) %*% matrix(rnorm(3*n), ncol = n)
    M = 3
    data <- list(Y = t(Y), n = n, M = M, mu = mu, R = R)

    newDf = 4 + n
    newR = R + tcrossprod(Y- mu)
    CTrueMean <- newR / newDf
    CTrueMeanVec <- as.numeric(CTrueMean)
    allData <- data
    constants <- list(n = allData$n, M = allData$M, mu = allData$mu, R = allData$R)
    data <- list(Y = allData$Y)
    inits <- list(C = CTrueMean)
    niter <- 50000

    Rmodel <- nimbleModel(code, constants, data, inits)
    Rmodel$calculate()
    conf <- configureMCMC(Rmodel)
    Rmcmc <- buildMCMC(conf)
    compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
    Cmodel <- compiledList$model; Cmcmc <- compiledList$mcmc
    set.seed(0)
    samples <- runMCMC(Cmcmc, niter)
    conjMean <- as.numeric(apply(samples, 2, mean))
    conjSD <- as.numeric(apply(samples, 2, sd))

    Rmodel <- nimbleModel(code, constants, data, inits)
    Rmodel$calculate()
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('C', 'RW_wishart')
    Rmcmc <- buildMCMC(conf)
    compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
    Cmodel <- compiledList$model; Cmcmc <- compiledList$mcmc
    set.seed(0)
    samples <- runMCMC(Cmcmc, niter)
    RWMean <- as.numeric(apply(samples, 2, mean))
    RWSD <- as.numeric(apply(samples, 2, sd))

    expect_true(all(round(as.numeric(RWMean - conjMean), 9) == c(-0.001651758, -0.009675571, 0.004894809, -0.009675571, 0.015533882, -0.008256095, 0.004894809, -0.008256095, 0.002119615)))
    expect_true(all(round(as.numeric(RWSD - conjSD), 9) == c(0.022803503, -0.010107015, 0.012342044, -0.010107015, 0.006191412, -0.000091101, 0.012342044, -0.000091101, 0.001340032)))
})

test_that('detect conjugacy when scaling Wishart, inverse Wishart cases', {
    code <- nimbleCode({
        mycov[1:p, 1:p]  <- lambda * Sigma[1:p,1:p] / eta
        y[1:p] ~ dmnorm(z[1:p], cov = mycov[1:p, 1:p])
        Sigma[1:p,1:p] ~ dinvwish(S[1:p,1:p], nu)
    })
    m  <- nimbleModel(code, constants = list(p = 3))
    expect_identical(length(m$checkConjugacy('Sigma')), 1L, 'inverse-Wishart case')

    code <- nimbleCode({
        mycov[1:p, 1:p]  <- lambda * Sigma[1:p,1:p] / eta
        y[1:p] ~ dmnorm(z[1:p], prec = mycov[1:p, 1:p])
        Sigma[1:p,1:p] ~ dwish(S[1:p,1:p], nu)
    })
    m  <- nimbleModel(code, constants = list(p = 3))
    expect_identical(length(m$checkConjugacy('Sigma')), 1L, 'Wishart case')

    code <- nimbleCode({
        mycov[1:p, 1:p]  <- lambda * Sigma[1:p,1:p] / eta
        y[1:p] ~ dmnorm(z[1:p], cov = mycov[1:p, 1:p])
        Sigma[1:p,1:p] ~ dwish(S[1:p,1:p], nu)
    })
    m  <- nimbleModel(code, constants = list(p = 3))
    expect_identical(length(m$checkConjugacy('Sigma')), 0L, 'inverse-Wishart not-conj case')

    code <- nimbleCode({
        mycov[1:p, 1:p]  <- lambda * Sigma[1:p,1:p] / eta
        y[1:p] ~ dmnorm(z[1:p], prec = mycov[1:p, 1:p])
        Sigma[1:p,1:p] ~ dinvwish(S[1:p,1:p], nu)
    })
    m  <- nimbleModel(code, constants = list(p = 3))
    expect_identical(length(m$checkConjugacy('Sigma')), 0L, 'Wishart not-conj case')

    code <- nimbleCode({
        mycov[1:p, 1:p]  <- lambda[1:p,1:p] * Sigma[1:p,1:p]
        y[1:p] ~ dmnorm(z[1:p], cov = mycov[1:p, 1:p])
        Sigma[1:p,1:p] ~ dinvwish(S[1:p,1:p], nu)
    })
    m  <- nimbleModel(code, constants = list(p = 3))
    expect_identical(length(m$checkConjugacy('Sigma')), 0L, 'Wishart case')

    code <- nimbleCode({
        mycov[1:p, 1:p]  <- lambda[1:p,1:p] %*% Sigma[1:p,1:p]
        y[1:p] ~ dmnorm(z[1:p], cov = mycov[1:p, 1:p])
        Sigma[1:p,1:p] ~ dinvwish(S[1:p,1:p], nu)
    })
    m  <- nimbleModel(code, constants = list(p = 3))
    expect_identical(length(m$checkConjugacy('Sigma')), 0L, 'Wishart case')
})


## testing conjugate MVN updating with ragged dependencies;
## that is, dmnorm dependents of different lengths from the target node
test_that('conjugate MVN with ragged dependencies', {
    cat('===== Starting MCMC test for conjugate MVN with ragged dependencies. =====')
    code <- nimbleCode({
        x[1:3] ~ dmnorm(mu0[1:3], prec = ident[1:3,1:3])
        mu_y2[1:2] <- asCol(a[1:2]) + B[1:2,1:3] %*% asCol(x[1:3])
        mu_y3[1:3] <- asCol(a[1:3]) + B[1:3,1:3] %*% asCol(x[1:3])
        mu_y5[1:5] <- asCol(a[1:5]) + B[1:5,1:3] %*% asCol(x[1:3])
        y2[1:2] ~ dmnorm(mu_y2[1:2], prec = prec_y[1:2,1:2])
        y3[1:3] ~ dmnorm(mu_y3[1:3], prec = prec_y[1:3,1:3])
        y5[1:5] ~ dmnorm(mu_y5[1:5], prec = prec_y[1:5,1:5])
    })
    
    mu0 <- rep(0,3)
    ident <- diag(3)
    a <- 11:15
    B <- matrix(1:15, nrow=5, ncol=3, byrow=TRUE)
    prec_y <- diag(1:5)
    
    constants <- list(mu0=mu0, ident=ident, a=a, B=B, prec_y=prec_y)
    data <- list(y2=1:2, y3=1:3, y5=1:5)
    inits <- list(x=rep(0,3))
    
    Rmodel <- nimbleModel(code, constants, data, inits)
    
    spec <- configureMCMC(Rmodel)
    Rmcmc <- buildMCMC(spec)
    
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    
    set.seed(0)
    Rmcmc$run(10)
    
    set.seed(0)
    Cmcmc$run(10)
    
    Rsamples <- as.matrix(Rmcmc$mvSamples)
    Csamples <- as.matrix(Cmcmc$mvSamples)

    expect_true(all(abs(as.numeric(Rsamples[,]) - c(4.96686874, 3.94112676, 4.55975130, 4.01930176, 4.47744412, 4.12927167, 4.91242131, 4.62837537, 4.54227859, 4.97237602, -1.12524733, 1.24545265, -0.13454814, 0.82755276, 0.08252775, 0.71187071, -0.31322184, -0.57462284, -0.64800963, -0.52885823, -3.92276916, -5.23904995, -4.53535941, -4.89919931, -4.66995650, -4.94181562, -4.63558011, -4.16385294, -4.03469945, -4.51128205)) < 1E-8),
                info = 'correct samples for ragged dmnorm conjugate update')

    dif <- Rsamples - Csamples

    expect_true(all(abs(dif) < 2E-13), info = 'R and C samples same for ragged dmnorm conjugate update')
    
    set.seed(0)
    Cmcmc$run(200000)
    Csamples <- as.matrix(Cmcmc$mvSamples)
    
    obsmean <- apply(Csamples, 2, mean)
    
    obsprec <- inverse(cov(Csamples))
    
    pprec <- ident +
        t(B[1:2,1:3]) %*% prec_y[1:2,1:2] %*% B[1:2,1:3] +
        t(B[1:3,1:3]) %*% prec_y[1:3,1:3] %*% B[1:3,1:3] +
        t(B[1:5,1:3]) %*% prec_y[1:5,1:5] %*% B[1:5,1:3]
    
    
    pmean <- inverse(pprec) %*% (ident %*% mu0 +
                                 t(B[1:2,1:3]) %*% prec_y[1:2,1:2] %*% (1:2 - a[1:2]) +
                                 t(B[1:3,1:3]) %*% prec_y[1:3,1:3] %*% (1:3 - a[1:3]) +
                                 t(B[1:5,1:3]) %*% prec_y[1:5,1:5] %*% (1:5 - a[1:5])   )
    

    expect_true(all(abs(pmean - obsmean) / pmean < 0.01), info = 'ragged dmnorm conjugate posterior mean')
    expect_true(all(abs(pprec - obsprec) / pprec < 0.005), info = 'ragged dmnorm conjugate posterior precision')    
})
    
## testing binary sampler
test_that('binary sampler setup', {
    cat('===== Starting MCMC test for binary sampler. =====')
    code <- nimbleCode({
        a ~ dbern(0.5)
        b ~ dbern(0.6)
        c ~ dbern(0.05)
        d ~ dbin(prob=0.2, size=1)
        e ~ dbinom(prob=0.9, size=1)
        f ~ dbern(0.5)
        g ~ dbern(0.5)
        h ~ dbern(0.5)
        for(i in 1:10)
            yf[i] ~ dnorm(f, sd = 1)
        for(i in 1:10)
            yg[i] ~ dnorm(g, sd = 1)
        for(i in 1:10)
            yh[i] ~ dnorm(h, sd = 1)
    })
    constants <- list()
    data <- list(yf = c(rep(0,2), rep(1,8)), yg = c(rep(0,8), rep(1,2)), yh = c(rep(0,5), rep(1,5)))
    inits <- list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, h=0)
    
    Rmodel <- nimbleModel(code, constants, data, inits)
    
    expect_true(Rmodel$isBinary('a'), info = 'model$isBinary')
    expect_true(Rmodel$isBinary('b'), info = 'model$isBinary')
    expect_true(Rmodel$isBinary('c'), info = 'model$isBinary')
    expect_true(Rmodel$isBinary('d'), info = 'model$isBinary')
    expect_true(Rmodel$isBinary('e'), info = 'model$isBinary')
    expect_true(Rmodel$isBinary('f'), info = 'model$isBinary')
    expect_true(Rmodel$isBinary('g'), info = 'model$isBinary')
    expect_true(Rmodel$isBinary('h'), info = 'model$isBinary')
    
    spec <- configureMCMC(Rmodel, nodes = NULL)
    spec$addSampler('a', 'binary', print=FALSE)
    spec$addSampler('b', 'binary', print=FALSE)
    spec$addSampler('c', 'binary', print=FALSE)
    spec$addSampler('d', 'binary', print=FALSE)
    spec$addSampler('e', 'binary', print=FALSE)
    spec$addSampler('f', 'binary', print=FALSE)
    spec$addSampler('g', 'binary', print=FALSE)
    spec$addSampler('h', 'binary', print=FALSE)
    ##spec$printSamplers()
    
    Rmcmc <- buildMCMC(spec)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    
    set.seed(0)
    Cmcmc$run(100000)
    samples <- as.matrix(Cmcmc$mvSamples)
    means <- apply(samples, 2, mean)
    ##means
    
    tol <- 0.0025
    test_that('binary sampler posterior', expect_lt(abs(means[['a']] - 0.5), tol))
    test_that('binary sampler posterior', expect_lt(abs(means[['b']] - 0.6), tol))
    test_that('binary sampler posterior', expect_lt(abs(means[['c']] - 0.05), tol))
    test_that('binary sampler posterior', expect_lt(abs(means[['d']] - 0.2), tol))
    test_that('binary sampler posterior', expect_lt(abs(means[['e']] - 0.9), tol))
    test_that('binary sampler posterior', expect_lt(abs(means[['f']] - 0.9525), tol))
    test_that('binary sampler posterior', expect_lt(abs(means[['g']] - 0.0475), tol))
    test_that('binary sampler posterior', expect_lt(abs(means[['h']] - 0.5), tol))
    
})
    
    ## testing the binary sampler handles 'out of bounds' ok
test_that('binary sampler handles out of bounds', {
    cat('===== Starting MCMC test for binary sampler handles out of bounds. =====')
    code <- nimbleCode({
        px ~ dbern(0.5)
        py ~ dbern(0.5)
        x ~ dnorm(0, sd = px - 0.5)
        y ~ dnorm(0, tau = py)
    })
    constants <- list()
    data <- list(x = 0, y = 0)
    inits <- list(px = 1, py = 1)
    Rmodel <- nimbleModel(code, constants, data, inits)
    
    spec <- configureMCMC(Rmodel)
    if(nimbleOptions('verbose')) spec$printSamplers()
    Rmcmc <- buildMCMC(spec)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    
    set.seed(0)
    Rmcmc$run(100)
    Rsamples <- as.matrix(Rmcmc$mvSamples)
    test_that('binary sampler out-of-bounds', expect_true(all(as.numeric(Rsamples) == 1)))
    
    set.seed(0)
    Cmcmc$run(100)
    Csamples <- as.matrix(Cmcmc$mvSamples)
    test_that('binary sampler out-of-bounds', expect_true(all(as.numeric(Csamples) == 1)))
})
    
    ## testing the RW_multinomial sampler
test_that('RW_multinomial sampler', {
    cat('===== Starting MCMC test for RW_multinomial sampler. =====')
    codeTest <- nimbleCode ({
        X[1:nGroups] ~ dmultinom(size=N, prob=pVecX[1:nGroups])
        Y[1:nGroups] ~ dmultinom(size=N, prob=pVecY[1:nGroups])
        for (ii in 1:nGroups) {
            Z[ii] ~ dbeta(1 + X[ii], 1 + Y[ii])
        }
    })
    
    set.seed(0)
    nGroups   <- 5
    N         <- 1E6
    pVecX     <- rdirch(1, rep(1, nGroups))
    pVecY     <- rdirch(1, rep(1, nGroups))
    X         <- rmultinom(1, N, pVecX)[,1]
    Y         <- rmultinom(1, N, pVecY)[,1]
    Z         <- rbeta(nGroups, 1+X, 1+Y)
    ## Hard code in the results of sample() since output from sample
    ## changed as of R 3.6.0 to fix a long-standing bug in R.
    smpX <- pVecX[c(2,1,4,3,5)]
    smpY <- pVecY[c(1,4,2,3,5)]
    fakeSample <- sample(pVecX)  # to keep random number stream as before
    Xini      <- rmultinom(1, N, smpX)[,1]
    fakeSample <- sample(pVecY)  # to keep random number stream as before
    Yini      <- rmultinom(1, N, smpY)[,1]
    Constants <- list(nGroups=nGroups)
    Inits     <- list(X=Xini, Y=Yini, pVecX=pVecX, pVecY=pVecY, N=N)
    Data      <- list(Z=Z)
    modelTest <- nimbleModel(codeTest, constants=Constants, inits=Inits, data=Data, check=TRUE)
    cModelTest <- compileNimble(modelTest)
    
    mcmcTestConfig <- configureMCMC(cModelTest, print = nimbleOptions('verbose'))
    samplers <- mcmcTestConfig$getSamplers()
    test_that('assign RW_multinomial sampler', expect_equal(samplers[[1]]$name, 'RW_multinomial'))
    test_that('assign RW_multinomial sampler', expect_equal(samplers[[2]]$name, 'RW_multinomial'))
    mcmcTest  <- buildMCMC(mcmcTestConfig)
    cMcmcTest <- compileNimble(mcmcTest, project=modelTest)
    
    ## Optionally resample data
    cModelTest$N      <- N <- 1E3
    (cModelTest$pVecX <- sort(rdirch(1, rep(1, nGroups))))
    (cModelTest$pVecY <- sort(rdirch(1, rep(1, nGroups))))
    simulate(cModelTest, "X", includeData=TRUE); (X <- cModelTest$X)
    simulate(cModelTest, "Y", includeData=TRUE); (Y <- cModelTest$Y)
    simulate(cModelTest, "Z", includeData=TRUE); (Z <- cModelTest$Z)
    
    niter  <- 1E4
    cMcmcTest$run(niter)
    samples <- as.matrix(cMcmcTest$mvSamples)
    
    expect_identical(as.numeric(samples[10000,]), c(8, 25, 31, 115, 821, 25,19, 84, 510, 362),
                     info = 'exact results of RW_multinomial sampler')
})

## testing the RW_multinomial sampler on distribution of size 2
test_that('RW_multinomial sampler on distribution of size 2', {
    cat('===== Starting MCMC test RW_multinomial sampler on distribution of size 2. =====')
    code <- nimbleCode({
        prob[1] <- p
        prob[2] <- 1-p
        x[1:2] ~ dmultinom(size = N, prob = prob[1:2])
        y      ~ dbinom(   size = N, prob = p)
    })
    
    set.seed(0)
    N <- 100
    p <- 0.3
    x1 <- rbinom(1, size=N, prob=p)
    x2 <- N - x1
    inits <- list(N = N, p = p, x = c(x1, x2), y = x1)
    Rmodel <- nimbleModel(code, constants=list(), data=list(), inits=inits)
    Cmodel <- compileNimble(Rmodel)
    
    conf <- configureMCMC(Rmodel)
    if(nimbleOptions('verbose')) conf$printSamplers()
    conf$removeSamplers()
    if(nimbleOptions('verbose')) conf$printSamplers()
    conf$addSampler(target = 'x', type = 'RW_multinomial', print = nimbleOptions('verbose'))
    conf$addSampler(target = 'y', type = 'slice', print = nimbleOptions('verbose'))
    if(nimbleOptions('verbose')) conf$printSamplers()
    Rmcmc  <- buildMCMC(conf)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    
    Cmcmc$run(100000)
    
    samples <- as.matrix(Cmcmc$mvSamples)
    fracs <- apply(samples, 2, mean) / N
    expect_true(all(abs(as.numeric(fracs[c(1,3)]) - p) < 0.01),
                info = 'RW_multinomial sampler results within tolerance')
})

## testing RW_dirichlet sampler
## consistent with conjugate multinomial sampler
test_that('RW_dirichlet sampler consistent with conjugate multinomial sampler', {
    cat('===== Starting MCMC test if RW_dirichlet sampler consistent with conjugate multinomial sampler. =====')
    n <- 100
    alpha <- c(10, 30, 15)
    K <- length(alpha)
    p <- c(.12, .24, .09)
    y <- c(23, 67, 10)
    code <- quote({
        p[1:K]  ~ ddirch(alpha[1:K])
        p2[1:K] ~ ddirch(alpha[1:K])
        y[1:K]  ~ dmulti(p[1:K], n)
        y2[1:K] ~ dmulti(p2[1:K], n)
    })
    inits <- list(p = rep(1/K, K), alpha = alpha, p2 = rep(1/K,K))
    constants <- list(n=n, K=K)
    data <- list(y = y, y2 = y)
    Rmodel <- nimbleModel(code, constants, data, inits)
    Cmodel <- compileNimble(Rmodel)
    conf <- configureMCMC(Rmodel, nodes=NULL)
    conf$addSampler('p[1:3]',  'RW_dirichlet')
    conf$addSampler('p2[1:3]', 'conjugate')
    Rmcmc <- buildMCMC(conf)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    
    niter <- 30
    set.seed(0)
    Rmcmc$run(niter)
    Rsamples <- as.matrix(Rmcmc$mvSamples)
    set.seed(0)
    Cmcmc$run(niter)
    Csamples <- as.matrix(Cmcmc$mvSamples)
    
    nodes <- c('p[1]','p[2]','p[3]')
    ans <- c(0.12812261, 0.6728109, 0.19906652)
    tol <- 1e-6
    expect_equal(as.numeric(Rsamples[30, nodes]), ans, tolerance=tol, info = 'correct R RW_dirichlet samples')
    expect_equal(as.numeric(Csamples[30, nodes]), ans, tolerance=tol, info = 'correct C RW_dirichlet samples')
    
    Cmcmc$run(100000)
    Csamples <- as.matrix(Cmcmc$mvSamples)
    means <- apply(Csamples, 2, mean)

    expect_true(all(abs(means[c('p[1]','p[2]','p[3]')] - means[c('p2[1]','p2[2]','p2[3]')]) < 0.001),
                info = 'agreement between RW_dirichlet and conjugate dirichlet sampling' )
})
    
## testing RW_dirichlet sampler
## more complicated -- intermediate deterministic nodes, and non-conjugate
## test agreement between RW_dirichlet sampler, and writing model with component gammas
test_that('RW_dirichlet sampler more complicated', {
    cat('===== Starting MCMC test RW_dirichlet sampler more complicated. =====')
    code <- nimbleCode({
        alpha[1] <- 4
        alpha[2] ~ dgamma(1,1)
        alpha[3] ~ dgamma(0.1, 0.1)
        alpha[4] <- alpha[2] + alpha[3]
        p1[1:4] ~ ddirch(alpha[1:4])
        q1[1] <- p1[1]
        q1[2] <- p1[2]
        q1[3] <- p1[3]/2
        q1[4] <- p1[3]/2
        q1[5] <- p1[4]
        y1[1:5] ~ dmulti(q1[1:5], n)
        for(i in 1:4) {
            theta[i] ~ dgamma(alpha[i], 1)
            p2[i] <- theta[i] / V
        }
        V <- sum(theta[1:4])
        q2[1] <- p2[1]
        q2[2] <- p2[2]
        q2[3] <- p2[3]/2
        q2[4] <- p2[3]/2
        q2[5] <- p2[4]
        y2[1:5] ~ dmulti(q2[1:5], n)
    })
    y <- c(50, 10, 5, 5, 30)
    constants <- list(n=100)
    data <- list(y1=y, y2=y)
    inits <- list(alpha=c(4,1,1,2), p1=rep(0.25,4), p2=rep(0.25,4), theta=rep(1,4))
    Rmodel <- nimbleModel(code, constants, data, inits)
    conf <- configureMCMC(Rmodel)
    if(nimbleOptions('verbose')) conf$printSamplers()
    conf$addMonitors('p1','p2')
    Rmcmc <- buildMCMC(conf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    
    Cmcmc$run(100000)
    samples <- as.matrix(Cmcmc$mvSamples)
    means <- apply(samples, 2, mean)
    sds <- apply(samples, 2, sd)

    expect_true(all(abs(means[c('p1[1]','p1[2]','p1[3]','p1[4]')] - means[c('p2[1]','p2[2]','p2[3]','p2[4]')]) < 0.01),
                info = 'non-conjugate agreement between RW_dirichlet and component gamma sampling: mean')
    expect_true(all(abs(sds[c('p1[1]','p1[2]','p1[3]','p1[4]')] - sds[c('p2[1]','p2[2]','p2[3]','p2[4]')]) < 0.01),
                info = 'non-conjugate agreement between RW_dirichlet and component gamma sampling: sd')
})


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
    
    expect_equal(as.numeric(Csamples[20, sampleNames]),
                 c(1.639127, 1.815422, 1.676655, 5.099797, 2.345276, 7.018026, 2.696936),
                 tolerance = 1e-6,
                 info = 'exact sample values for dcar_normal')
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
    
    expect_equal(calculate(Rmodel), lp, tol = 1E-5,
                 info = 'calculate for dcar_proper()')
    
    expect_equal(calculate(Cmodel), lp, tol = 1E-5,
                 info = 'calculate for dcar_proper(), compiled')
    
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

    expect_equal(as.numeric(Csamples[20, sampleNames]),
                 c(0.0978025, -0.6643286, 1.9510954, 0.2413084, 2.6684426, -3.2533691),
                 tolerance = 1e-6,
                 info = 'exact sample values for dcar_proper')
})


## testing dmnorm-dnorm conjugacies we don't detect

test_that('dnorm-dmnorm conjugacies NIMBLE fails to detect', {
    code = nimbleCode({
        y[1:3] ~ dmnorm(mu[1:3], pr[1:3,1:3])
        for(i in 1:3)
            mu[i] ~ dnorm(0,1)
        pr[1:3,1:3] ~ dwish(pr0[1:3,1:3], 5)
    })
    m = nimbleModel(code, inits = list(pr0 = diag(3), pr = diag(3)))
    conf <- configureMCMC(m)
    expect_failure(expect_match(conf$getSamplers()[[1]]$name, "conjugate_dmnorm_dnorm",
         info = "failed to detect dmnorm-dnorm conjugacy"),
         info = "EXPECTED FAILURE NOT FAILING: this known failure should occur because of limitations in conjugacy detection with dmnorm dependents of dnorm target")

    code = nimbleCode({
        for(i in 1:3)
            y[i] ~ dnorm(mu[i], var = s0)
        mu[1:3] ~ dmnorm(z[1:3],pr[1:3,1:3])
        s0 ~ dinvgamma(1,1)
    })
    m = nimbleModel(code, inits = list(z = rep(0,3), pr = diag(3)))
    conf <- configureMCMC(m)
    expect_failure(expect_match(conf$getSamplers()[[2]]$name, "conjugate_dnorm_dmnorm",
         info = "failed to detect dmnorm-dnorm conjugacy"),
         info = "EXPECTED FAILURE NOT FAILING: this known failure should occur because of limitations in conjugacy detection with dnorm dependents of dmnorm target")
})

## dnorm prior in vectorized regression mean (inprod, matrix multiplication)

test_that('NIMBLE detects dnorm-dnorm conjugacy via inprod() or %*%', {

    code <- nimbleCode({
        for(i in 1:n) 
            y[i] ~ dnorm(b0 + inprod(beta[1:p], X[i, 1:p]), 1)
        for(i in 1:p) 
            beta[i] ~ dnorm(0, 1)
        b0 ~ dnorm(0, 1)
    })
    constants <- list(n = 5, p = 3)
    data <- list(y = rnorm(constants$n),
                 X = matrix(rnorm(constants$n * constants$p), constants$n))
    inits <- list(b0 = 1, beta = rnorm(constants$p))
    m <- nimbleModel(code, data = data, constants = constants)
    conf <- configureMCMC(m)
    expect_identical(conf$getSamplers()[[1]]$name, 'conjugate_dnorm_dnorm',
                                   info = "conjugacy with inprod not detected")

    code <- nimbleCode({
        for(i in 1:n) 
            y[i] ~ dnorm(b0 + sum(beta[1:p]*X[i, 1:p]), 1)
        for(i in 1:p) 
            beta[i] ~ dnorm(0, 1)
        b0 ~ dnorm(0, 1)
    })
    constants <- list(n = 5, p = 3)
    data <- list(y = rnorm(constants$n),
                 X = matrix(rnorm(constants$n * constants$p), constants$n))
    m <- nimbleModel(code, data = data, constants = constants)
    conf <- configureMCMC(m)
    expect_identical(conf$getSamplers()[[1]]$name, 'conjugate_dnorm_dnorm',
                                   info = "conjugacy with inprod not detected")


    ## compare to conjugate sampler using summed contributions
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc, project = m)
    smp1 <- runMCMC(cmcmc, 50, setSeed = 1)

    code <- nimbleCode({
        for(i in 1:n) 
            y[i] ~ dnorm(b0 + sum(beta[1:p]*X[i, 1:p]), 1)
        for(i in 1:p) 
            beta[i] ~ dnorm(0, 1)
        b0 ~ dnorm(0, 1)
    })
    m <- nimbleModel(code, data = data, constants = constants)
    conf <- configureMCMC(m)
    expect_identical(conf$getSamplers()[[1]]$name, 'conjugate_dnorm_dnorm',
                                   info = "conjugacy with sum not detected")


    ## compare to conjugate sampler using summed contributions
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc, project = m)
    smp2 <- runMCMC(cmcmc, 50, setSeed = 1)

    code <- nimbleCode({
        for(i in 1:n) 
            y[i] ~ dnorm(b0 + beta[1]*X[i,1] + beta[2]*X[i,2] + beta[3]*X[i,3], 1)
        for(i in 1:p) 
            beta[i] ~ dnorm(0, 1)
        b0 ~ dnorm(0, 1)
    })
    m <- nimbleModel(code, data = data, constants = constants)
    conf <- configureMCMC(m)
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc, project = m)
    smp3 <- runMCMC(cmcmc, 50, setSeed = 1)
    expect_equal(smp1, smp3, info = 'conjugate sampler with inprod does not match summation')
    expect_equal(smp2, smp3, info = 'conjugate sampler with sum does not match summation')
    
    code <- nimbleCode({
        for(i in 1:n) 
            y[i] ~ dnorm(b0 + inprod(exp(beta[1:p]), X[i, 1:p]), 1)
        for(i in 1:p) 
            beta[i] ~ dnorm(0, 1)
        b0 ~ dnorm(0, 1)
    })
    m <- nimbleModel(code, data = data, constants = constants)
    conf <- configureMCMC(m)
    expect_identical(conf$getSamplers()[[1]]$name, 'RW',
                                   info = "conjugacy with inprod improperly detected")

    code <- nimbleCode({
        for(i in 1:n) 
            y[i] ~ dnorm(b0 + exp(inprod(beta[1:p], X[i, 1:p])), 1)
        for(i in 1:p) 
            beta[i] ~ dnorm(0, 1)
        b0 ~ dnorm(0, 1)
    })
    m <- nimbleModel(code, data = data, constants = constants)
    conf <- configureMCMC(m)
    expect_identical(conf$getSamplers()[[1]]$name, 'RW',
                                   info = "conjugacy with inprod improperly detected")

    code <- nimbleCode({
        for(i in 1:n) 
            y[i] ~ dnorm(b0 + sum(exp(beta[1:p])*X[i, 1:p]), 1)
        for(i in 1:p) 
            beta[i] ~ dnorm(0, 1)
        b0 ~ dnorm(0, 1)
    })
    m <- nimbleModel(code, data = data, constants = constants)
    conf <- configureMCMC(m)
    expect_identical(conf$getSamplers()[[1]]$name, 'RW',
                                   info = "conjugacy with sum improperly detected")

    code <- nimbleCode({
        for(i in 1:n) 
            y[i] ~ dnorm(b0 + (X[i, 1:p] %*% beta[1:p])[1,1], 1)
        for(i in 1:p) 
            beta[i] ~ dnorm(0, 1)
        b0 ~ dnorm(0, 1)
    })
    m <- nimbleModel(code, data = data, constants = constants)
    conf <- configureMCMC(m)
    expect_identical(conf$getSamplers()[[1]]$name, 'conjugate_dnorm_dnorm',
                                   info = "conjugacy with matrix multiplication not detected")
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc, project = m)
    smp3 <- runMCMC(cmcmc, 50, setSeed = 1)
    expect_equal(smp1, smp3, info = 'conjugate sampler with matrix mult. does not match summation')

    check <- nimble:::cc_checkLinearity(quote(b0 + (X[1, 1:3] %*% structureExpr(beta[1], beta[2], beta[3]))[1,1]), 'beta[1]')
    expect_identical(check, list(offset = quote(b0 + X[1, 1:3] * structureExpr(beta[1], beta[2], beta[3])),
                                 scale = quote(X[1, 1:3])))

    check <- nimble:::cc_checkLinearity(quote(b0 + inprod(structureExpr(beta[1], beta[2], beta[3]), X[1, 1:3])), 'beta[1]')
    expect_identical(check, list(offset = quote(b0 + structureExpr(beta[1], beta[2], beta[3]) * X[1, 1:3]),
                                 scale = quote(X[1, 1:3])))

    ## check nested specifications
    
    code <- nimbleCode({
        for(i in 1:n) 
            y[i] ~ dnorm(b0 + inprod(zbeta[1:p], X[i, 1:p]), 1)
        for(i in 1:p) {
            beta[i] ~ dnorm(0, 1)
            zbeta[i] <- z[i] * beta[i]
        }
        b0 ~ dnorm(0, 1)
    })
    constants <- list(n = 5, p = 3)
    data <- list(y = rnorm(constants$n),
                 X = matrix(rnorm(constants$n * constants$p), constants$n))
    inits <- list(b0 = 1, beta = rnorm(constants$p))
    m <- nimbleModel(code, data = data, constants = constants)
    conf <- configureMCMC(m)
    expect_identical(conf$getSamplers()[[1]]$name, 'conjugate_dnorm_dnorm',
                                   info = "conjugacy with inprod not detected")
   
    code <- nimbleCode({
        for(i in 1:n) 
            y[i] ~ dnorm(b0 + inprod(zbeta[1:p], X[i, 1:p]), 1)
        for(i in 1:p) {
            beta[i] ~ dnorm(0, 1)
            wbeta[i] <- a + w * beta[i]
            zbeta[i] <- z[i] * wbeta[i]
        }
        b0 ~ dnorm(0, 1)
    })
    constants <- list(n = 5, p = 3)
    data <- list(y = rnorm(constants$n),
                 X = matrix(rnorm(constants$n * constants$p), constants$n))
    inits <- list(b0 = 1, beta = rnorm(constants$p))
    m <- nimbleModel(code, data = data, constants = constants)
    conf <- configureMCMC(m)
    expect_identical(conf$getSamplers()[[1]]$name, 'conjugate_dnorm_dnorm',
                                   info = "conjugacy with inprod not detected")

    code <- nimbleCode({
        for(i in 1:n) 
            y[i] ~ dnorm((X[i, 1:p] %*% zbeta[1:p])[1], 1)
        for(i in 1:p) {
            beta[i] ~ dnorm(0, 1)
            zbeta[i] <- z[i] * beta[i]
        }
    })
    constants <- list(n = 5, p = 3)
    data <- list(y = rnorm(constants$n),
                 X = matrix(rnorm(constants$n * constants$p), constants$n))
    inits <- list(b0 = 1, beta = rnorm(constants$p))
    m <- nimbleModel(code, data = data, constants = constants)
    conf <- configureMCMC(m)
    expect_identical(conf$getSamplers()[[1]]$name, 'conjugate_dnorm_dnorm',
                     info = "conjugacy with inprod not detected")

   code <- nimbleCode({
        for(i in 1:n) 
            y[i] ~ dnorm(b0 + inprod(zbeta[1:p], X[i, 1:p]), 1)
        for(i in 1:p) {
            beta[i] ~ dnorm(0, 1)
            wbeta[i] <- exp(w * beta[i])
            zbeta[i] <- z[i] * wbeta[i]
        }
        b0 ~ dnorm(0, 1)
    })
    constants <- list(n = 5, p = 3)
    data <- list(y = rnorm(constants$n),
                 X = matrix(rnorm(constants$n * constants$p), constants$n))
    inits <- list(b0 = 1, beta = rnorm(constants$p))
    m <- nimbleModel(code, data = data, constants = constants)
    conf <- configureMCMC(m)
    expect_identical(conf$getSamplers()[[1]]$name, 'RW',
                     info = "conjugacy with inprod mistakenly detected")

    expect_identical(nimble:::cc_checkLinearity(
        quote(structureExpr(z[1] * exp(w * beta[1]), z[2] * exp(w * beta[2]))),
        'beta[2]'), NULL)
    output <- nimble:::cc_checkLinearity(
        quote(structureExpr(z[1] * exp(w * beta[1]), a + z[2] * (d + w * beta[2]))),
        'beta[2]')
    expect_identical(is.list(output), TRUE)  ## should be a list with scale/offset
    
})


test_that('MCMC with logProb variable being monitored builds and compiles.', {
  cat('===== Starting MCMC test of logProb variable monitoring =====')
  code <- nimbleCode({
    prob[1] <- p
    prob[2] <- 1-p
    x[1:2] ~ dmultinom(size = N, prob = prob[1:2])
    y      ~ dbinom(   size = N, prob = p)
  })
  set.seed(0)
  N <- 100
  p <- 0.3
  x1 <- rbinom(1, size=N, prob=p)
  x2 <- N - x1
  inits <- list(N = N, p = p, x = c(x1, x2), y = x1)
  Rmodel <- nimbleModel(code, constants=list(), data=list(), inits=inits)
  Cmodel <- compileNimble(Rmodel)
  expect_silent(conf <- configureMCMC(Rmodel, monitors = 'logProb_y'))
  expect_silent(Rmcmc  <- buildMCMC(conf))
  if(nimbleOptions('verbose')) {
      expect_message(Cmcmc <- compileNimble(Rmcmc, project = Rmodel), "compilation finished")
  } else expect_silent(Cmcmc <- compileNimble(Rmcmc, project = Rmodel))
  Cmcmc$run(10)
})

test_that('slice sampler bails out of loop', {
    code <- nimbleCode({
        y ~ dnorm(0, sd = sigma)
        sigma ~ dunif(0.5-1e-15, 0.5+1e-15)
    })
    Rmodel <- nimbleModel(code, data = list(y = 0), inits = list(sigma = 0.5))
    Cmodel <- compileNimble(Rmodel)
    conf <- configureMCMC(Rmodel, monitors = 'logProb_y', onlySlice = TRUE)
    Rmcmc  <- buildMCMC(conf)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    expect_silent(Cmcmc$run(10))

    conf$removeSamplers('sigma')
    conf$addSampler('sigma','slice', control = list(maxContractions = 10))
    Rmcmc  <- buildMCMC(conf)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel, resetFunctions = TRUE)
    expect_output(Cmcmc$run(1), "slice sampler reached maximum number of contractions")
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

test_that('checkConjugacy corner case when linear scale is identically zero', {
    targetNode <- 'beta[4]'
    linearityCheckExpr <- quote(beta[4] * 0 * alpha.smrcent[3])
    conjugacyCheck <- nimble:::cc_checkLinearity(linearityCheckExpr, targetNode)
    expect_identical(conjugacyCheck, list(offset = 0, scale = 0))
    
    targetNode <- 'beta[4]'
    linearityCheckExpr <- quote(beta[1] + beta[2] * 0 + beta[3] * alpha.smrcent[3] + beta[4] * 0 * alpha.smrcent[3] + alpha.stream[1] + alpha.family[3, 1])
    conjugacyCheck <- nimble:::cc_checkLinearity(linearityCheckExpr, targetNode)
    expect_identical(conjugacyCheck,
                     list(offset = quote(beta[1] + beta[2] * 0 + beta[3] * alpha.smrcent[3] + alpha.stream[1] + alpha.family[3, 1]),
                          scale = 0))
})

test_that('cc_checkScalar operates correctly', {
    expect_true(cc_checkScalar(quote(lambda)))
    expect_true(cc_checkScalar(quote(lambda*eta)))
    expect_true(cc_checkScalar(quote(exp(lambda))))
    expect_true(cc_checkScalar(quote(5+exp(lambda[2]))))

    expect_false(cc_checkScalar(quote(exp(lambda[2:3]))))
    expect_false(cc_checkScalar(quote(lambda[1:2,1:2])))
    expect_false(cc_checkScalar(quote(lambda[1:2,1:2]/eta)))
    expect_false(cc_checkScalar(quote(eta*(theta*lambda[1:2,1:2]))))
    expect_false(cc_checkScalar(quote(lambda[1:2,1:2,1:5])))

    expect_true(cc_checkScalar(quote(lambda[xi[i]])))
    expect_true(cc_checkScalar(quote(lambda[xi[i],xi[j]])))
    expect_false(cc_checkScalar(quote(lambda[xi[i]:3])))
    expect_false(cc_checkScalar(quote(lambda[xi[i]:3,2])))

    ## Ideally this case would evaluate to TRUE, but we would
    ## have to handle knowing output dims of user-defined fxns.
    expect_false(cc_checkScalar(quote(sum(lambda[1:5]))))
    expect_false(cc_checkScalar(quote(foo(lambda))))

})

test_that('cc_stripExpr operates correctly', {
    expr <- 'coeff * (log(value)-offset) * taulog'
    expect_identical(cc_stripExpr(parse(text = 'coeff^2 * tau')[[1]], TRUE, TRUE),
                 'tau')
    expect_identical(cc_stripExpr(parse(text = expr)[[1]], TRUE, TRUE),
                 '(log(value)) * taulog')
    expect_identical(cc_stripExpr(parse(text = expr)[[1]], TRUE, FALSE),
                 'coeff * (log(value)) * taulog')
    expect_identical(cc_stripExpr(parse(text = expr)[[1]], FALSE, TRUE),
                 '(log(value)-offset) * taulog')
    expect_identical(cc_stripExpr(parse(text = expr)[[1]], FALSE, FALSE),
                 expr)
})

sink(NULL)

if(!generatingGoldFile) {
    trialResults <- readLines(tempFileName)
    trialResults <- trialResults[grep('Error in x$.self$finalize() : attempt to apply non-function', trialResults, invert = TRUE, fixed = TRUE)]
    correctResults <- readLines(system.file(file.path('tests', goldFileName), package = 'nimble'))
    compareFilesByLine(trialResults, correctResults)
}

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

## Tests for copy() [implemented as nimCopy], values(), and values()<-
## These use some of the same internals (accessors), so they are in the same testing file.
## These tests use lists of nimbleFunctions, initialization code, and testing code.
## They include copying to and from models and/or modelValues, using default arguments, using logProb = [TRUE|FALSE], and using the same or different blocks of variables.
## Checks are made internally for uncompiled and compiled cases.  Then uncompiled and compiled outcomes are compared to check that they behaved identically.

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context('Testing of nimCopy and values')




#############
## Here is a model with deterministic and stochastic variables of dimensions 0-4, including a multivariate node
copyTestModelCode <- nimbleCode({
    x0 ~ dnorm(0,1); d0 <- x0 + 10000
    for(i in 1:4) {x1[i] ~ dnorm(0,1); d1[i] <- x1[i]+10000}
    for(i in 1:4) for(j in 1:4) {x2[i,j] ~ dnorm(0,1); d2[i,j] <- x2[i,j]+10000}
    for(i in 1:4) for(j in 1:4) for(k in 1:4) {x3[i,j,k] ~ dnorm(0, 1); d3[i,j,k] <- x3[i,j,k]+10000}
    for(i in 1:4) for(j in 1:4) for(k in 1:4) for(l in 1:4) {x4[i,j,k,l] ~ dnorm(0, 1); d4[i,j,k,l] <- x4[i,j,k,l]+10000}
    for(i in 1:4) for(j in 1:4) for(k in 1:4) for(l in 1:4) for(m in 1:4) {x5[i,j,k,l,m] ~ dnorm(0, 1); d5[i,j,k,l,m] <- x5[i,j,k,l,m]+10000}
    for(i in 1:4) for(j in 1:4) for(k in 1:4) for(l in 1:4) for(m in 1:4) for(n in 1:4) {x6[i,j,k,l,m,n] ~ dnorm(0, 1); d6[i,j,k,l,m,n] <- x6[i,j,k,l,m,n]+10000}

    v1[1:4] ~ dmnorm(v1mu[1:4], v1sigma[1:4, 1:4]) ## not testing indexing here!
    for(i in 1:5) w1[ 2:5, i] ~ dmnorm(v1mu[1:4], v1sigma[1:4, 1:4])
    for(i in 1:4) foo[i] <- w1[i+1, 1]
})

copyTestConstants <- list()
copyTestData <- list(v1mu = rep(0,4), v1sigma = diag(4))

## A list of nimbleFunctions used for various cases
copyTestNFcodeList <- list(
    ## First 3 are ModelToMV
    nfModelToMVall = quote({
        nimbleFunction(
            setup = function(from, to, rowTo = 3, logProb = TRUE) {},
            run = function() {
                nimCopy(from = from, to = to, rowTo = rowTo, logProb = logProb)
            }
        )
    }),
    nfModelToMVsomeSame = quote({
        nimbleFunction(
            setup = function(from, nodes, to, rowTo = 3, logProb = TRUE) {},
            run = function() {
                nimCopy(from = from, nodes = nodes, to = to, rowTo = rowTo, logProb = logProb)
            }
        )
    }),
    nfModelToMVsomeDiff = quote({
        nimbleFunction(
            setup = function(from, nodes, to, nodesTo, rowTo = 3, logProb = TRUE) {},
            run = function() {
                nimCopy(from = from, nodes = nodes, to = to, rowTo = rowTo, nodesTo = nodesTo, logProb = logProb)
            }
        )
    }),
    ## Next 3 are MVToModel
    nfMVToModelAll = quote({
        nimbleFunction(
            setup = function(from, to, row = 3, logProb = TRUE) {},
            run = function() {
                nimCopy(from = from, to = to, row = row, logProb = logProb)
            }
        )
    }),
    nfMVToModelSomeSame = quote({
        nimbleFunction(
            setup = function(from, nodes, to, row = 3, logProb = TRUE) {},
            run = function() {
                nimCopy(from = from, nodes = nodes, to = to, row = row, logProb = logProb)
            }
        )
    }),
    nfMVToModelSomeDiff = quote({
        nimbleFunction(
            setup = function(from, nodes, to, nodesTo, row = 3, logProb = TRUE) {},
            run = function() {
                nimCopy(from = from, nodes = nodes, to = to, nodesTo = nodesTo, row = row, logProb = logProb)
            }
        )
    })
)

## A list of inputs for runOneCopyTest
copyTestCaseList <- list(
## First 2 are Model2MV checking default of all nodes if nodes arg is not provided
    Model2MVall_logProbTRUE = list(
        label = 'Model2MVall_logProbTRUE', ## A label
        nfName = 'nfModelToMVall',         ## Name of a nimbleFunction from above list to define as nf
        seed = round(runif(1, 1, 10000)),  ## A seed generated once then the same for compiled and uncompiled cases
        initCode = quote({simulate(m); calculate(m)}),  ## Code to initiate the model and/or modelValues
        compile = c(FALSE, TRUE),          ## sequence of yes/no compile
        compareRtoCpp = TRUE,              ## Should uncompiled and compiled be compared? When TRUE, only works if compile == c(FALSE, TRUE), in that order
        nfMcode = quote({nf(m, mv, logProb = TRUE)}), ## code for specializing the nf
        testThatLines = quote({            ## code for testing
            for(oneName in m$getVarNames(includeLogProb = TRUE)) {
                test_that('single', expect_identical(m[[oneName]], mv[[oneName]][[3]],
                                                     info = 'Model2MVall_logProbTRUE'))
            }
        })
    ),
    Model2MVall_logProbFALSE = list(
        label = 'Model2MVall_logProbFALSE',
        nfName = 'nfModelToMVall',
        seed = round(runif(1, 1, 10000)),
        initCode = quote({simulate(m); calculate(m)}),
        compile = c(FALSE, TRUE),
        compareRtoCpp = TRUE,
        nfMcode = quote({nf(m, mv, logProb = FALSE)}),
        testThatLines = quote({
            for(oneName in m$getVarNames(includeLogProb = FALSE)) {
                test_that('single', expect_identical(m[[oneName]], mv[[oneName]][[3]],
                                                     info = 'Model2MVall_logProbFALSE'))
            }
            lpNames <- m$getVarNames(includeLogProb = TRUE)
            lpNames <- lpNames[grep('logProb_', lpNames)]
            for(oneName in lpNames) {
                test_that('single', expect_false(identical(m[[oneName]], mv[[oneName]][[3]])))
            }
        })
    ),
## Next 2 copy some blocks of some variables, using the same for from and to (nodes provided but nodesTo not provided)
    Model2MVsomeSame_logProbFALSE = list( ## some = some blocks of nodes, ## same = same nodes from and to 
        label = 'Model2MVsomeSame_logProbFALSE',
        nfName = 'nfModelToMVsomeSame',
        seed = round(runif(1, 1, 10000)),
        initCode = quote({simulate(m); calculate(m)}),
        compile = c(FALSE, TRUE),
        compareRtoCpp = TRUE,
        nfMcode = quote({nf(m, mv, nodes = c('x0','d0','x1[2:3]','d1[2:3]','x2[2:3,2:3]','d2[2:3, 2:3]',
                                             'x3[2:3,2:3,2:3]','d3[2:3,2:3,2:3]', 'x4[2:3,2:3,2:3,2:3]','d4[2:3,2:3,2:3,2:3]',
                                             'x5[2:3, 2:3, 2:3, 2:3, 2:3]','d5[2:3, 2:3, 2:3, 2:3, 2:3]',
                                             'x6[2:3, 2:3, 2:3, 2:3, 2:3, 2:3]','d6[2:3, 2:3, 2:3, 2:3, 2:3, 2:3]',
                                             'v1[2:3]', 'w1[2:3,2:3]'), logProb = FALSE)}),
        testThatLines = quote({
            info <- 'Model2MVsomeSame_logProbFALSE'
            test_that('single', expect_identical(as.numeric(m[['x0']]), as.numeric(mv[['x0']][[3]]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d0']]), as.numeric(mv[['d0']][[3]]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['x1']])[2:3], as.numeric(mv[['x1']][[3]])[2:3], info = info ))
            test_that('single', expect_identical(as.numeric(m[['d1']])[2:3], as.numeric(mv[['d1']][[3]])[2:3], info = info ))
            test_that('single', expect_identical(as.numeric(m[['x2']][2:3, 2:3]), as.numeric(mv[['x2']][[3]][2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d2']][2:3, 2:3]), as.numeric(mv[['d2']][[3]][2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['x3']][2:3, 2:3, 2:3]), as.numeric(mv[['x3']][[3]][2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d3']][2:3, 2:3, 2:3]), as.numeric(mv[['d3']][[3]][2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['x4']][2:3, 2:3, 2:3, 2:3]), as.numeric(mv[['x4']][[3]][2:3, 2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d4']][2:3, 2:3, 2:3, 2:3]), as.numeric(mv[['d4']][[3]][2:3, 2:3, 2:3, 2:3]), info = info ))

            test_that('single', expect_identical(as.numeric(m[['x5']][2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(mv[['x5']][[3]][2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d5']][2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(mv[['d5']][[3]][2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))

            test_that('single', expect_identical(as.numeric(m[['x6']][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(mv[['x6']][[3]][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d6']][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(mv[['d6']][[3]][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))

            test_that('single', expect_identical(as.numeric(m[['v1']][2:3]), as.numeric(mv[['v1']][[3]][2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['w1']][2:3, 2:3]), as.numeric(mv[['w1']][[3]][2:3, 2:3]), info = info ))
        })
       ),
        Model2MVsomeSame_logProbTRUE = list(
        label = 'Model2MVsomeSame_logProbTRUE',
        nfName = 'nfModelToMVsomeSame',
        seed = round(runif(1, 1, 10000)),
        initCode = quote({simulate(m); calculate(m)}),
        compile = c(FALSE, TRUE),
        compareRtoCpp = TRUE,
        nfMcode = quote({nf(m, mv, nodes = c('x0','d0','x1[2:3]','d1[2:3]','x2[2:3,2:3]','d2[2:3, 2:3]',
                                             'x3[2:3,2:3,2:3]','d3[2:3,2:3,2:3]', 'x4[2:3,2:3,2:3,2:3]','d4[2:3,2:3,2:3,2:3]',
                                             'x5[2:3, 2:3, 2:3, 2:3, 2:3]','d5[2:3, 2:3, 2:3, 2:3, 2:3]',
                                             'x6[2:3, 2:3, 2:3, 2:3, 2:3, 2:3]','d6[2:3, 2:3, 2:3, 2:3, 2:3, 2:3]',
                                             'v1[2:3]', 'w1[2:3,2:3]'), logProb = TRUE)}),
        testThatLines = quote({
            info = 'Model2MVsomeSame_logProbTRUE'
            test_that('single', expect_identical(as.numeric(m[['x0']]), as.numeric(mv[['x0']][[3]]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d0']]), as.numeric(mv[['d0']][[3]]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['x1']])[2:3], as.numeric(mv[['x1']][[3]])[2:3], info = info ))
            test_that('single', expect_identical(as.numeric(m[['d1']])[2:3], as.numeric(mv[['d1']][[3]])[2:3], info = info ))
            test_that('single', expect_identical(as.numeric(m[['x2']][2:3, 2:3]), as.numeric(mv[['x2']][[3]][2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d2']][2:3, 2:3]), as.numeric(mv[['d2']][[3]][2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['x3']][2:3, 2:3, 2:3]), as.numeric(mv[['x3']][[3]][2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d3']][2:3, 2:3, 2:3]), as.numeric(mv[['d3']][[3]][2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['x4']][2:3, 2:3, 2:3, 2:3]), as.numeric(mv[['x4']][[3]][2:3, 2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d4']][2:3, 2:3, 2:3, 2:3]), as.numeric(mv[['d4']][[3]][2:3, 2:3, 2:3, 2:3]), info = info ))

            test_that('single', expect_identical(as.numeric(m[['x5']][2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(mv[['x5']][[3]][2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d5']][2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(mv[['d5']][[3]][2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))

            test_that('single', expect_identical(as.numeric(m[['x6']][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(mv[['x6']][[3]][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d6']][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(mv[['d6']][[3]][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))

            test_that('single', expect_identical(as.numeric(m[['v1']][2:3]), as.numeric(mv[['v1']][[3]][2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['w1']][2:3, 2:3]), as.numeric(mv[['w1']][[3]][2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['logProb_x0']]), as.numeric(mv[['logProb_x0']][[3]]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['logProb_x1']])[2:3], as.numeric(mv[['logProb_x1']][[3]])[2:3], info = info ))
            test_that('single', expect_identical(as.numeric(m[['logProb_x2']][2:3, 2:3]), as.numeric(mv[['logProb_x2']][[3]][2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['logProb_x3']][2:3, 2:3, 2:3]), as.numeric(mv[['logProb_x3']][[3]][2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['logProb_x4']][2:3, 2:3, 2:3, 2:3]), as.numeric(mv[['logProb_x4']][[3]][2:3, 2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['logProb_x5']][2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(mv[['logProb_x5']][[3]][2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['logProb_x6']][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(mv[['logProb_x6']][[3]][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))

            test_that('single', expect_identical(as.numeric(m[['logProb_v1']][1]), as.numeric(mv[['logProb_v1']][[3]][1]), info = info )) ## Note these last two use the collapsing of logProb vars for multivariate nodes
            test_that('single', expect_identical(as.numeric(m[['logProb_w1']][2, 2:3]), as.numeric(mv[['logProb_w1']][[3]][2, 2:3]), info = info ))
        })
        ),
## Next 2 copy some node blocks, with different nodes and nodesTo
    Model2MVsomeDiff_logProbFALSE = list( ## diff = different nodesTo from nodes
        label = 'Model2MVsomeDiff_logProbFALSE',
        nfName = 'nfModelToMVsomeDiff',
        seed = round(runif(1, 1, 10000)),
        initCode = quote({simulate(m); calculate(m)}),
        compile = c(FALSE, TRUE),
        compareRtoCpp = TRUE,
        nfMcode = quote({nf(m, mv, nodes = c('x0','d0','x1[1:2]','d1[1:2]','x2[1:2,1:2]','d2[1:2, 1:2]',
                                             'x3[1:2,1:2,1:2]','d3[1:2,1:2,1:2]', 'x4[1:2,1:2,1:2,1:2]','d4[1:2,1:2,1:2,1:2]',
                                             'x5[1:2, 1:2, 1:2, 1:2, 1:2]','d5[1:2, 1:2, 1:2, 1:2, 1:2]',
                                             'x6[1:2, 1:2, 1:2, 1:2, 1:2, 1:2]','d6[1:2, 1:2, 1:2, 1:2, 1:2, 1:2]',
                                             'v1[1:2]', 'w1[1:2,1:2]'),
                            nodesTo = c('x0','d0','x1[3:4]','d1[3:4]','x2[3:4,3:4]','d2[3:4, 3:4]',
                                        'x3[3:4,3:4,3:4]','d3[3:4,3:4,3:4]', 'x4[3:4,3:4,3:4,3:4]','d4[3:4,3:4,3:4,3:4]',
                                        'x5[3:4, 3:4, 3:4, 3:4, 3:4]','d5[3:4, 3:4, 3:4, 3:4, 3:4]',
                                        'x6[3:4, 3:4, 3:4, 3:4, 3:4, 3:4]','d6[3:4, 3:4, 3:4, 3:4, 3:4, 3:4]',
                                        'v1[3:4]', 'w1[3:4,3:4]'),
                            logProb = FALSE)}),
        testThatLines = quote({
            info = 'Model2MVsomeDiff_logProbFALSE'
            test_that('single', expect_identical(as.numeric(m[['x0']]), as.numeric(mv[['x0']][[3]]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d0']]), as.numeric(mv[['d0']][[3]]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['x1']])[1:2], as.numeric(mv[['x1']][[3]])[3:4], info = info ))
            test_that('single', expect_identical(as.numeric(m[['d1']])[1:2], as.numeric(mv[['d1']][[3]])[3:4], info = info ))
            test_that('single', expect_identical(as.numeric(m[['x2']][1:2, 1:2]), as.numeric(mv[['x2']][[3]][3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d2']][1:2, 1:2]), as.numeric(mv[['d2']][[3]][3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['x3']][1:2, 1:2, 1:2]), as.numeric(mv[['x3']][[3]][3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d3']][1:2, 1:2, 1:2]), as.numeric(mv[['d3']][[3]][3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['x4']][1:2, 1:2, 1:2, 1:2]), as.numeric(mv[['x4']][[3]][3:4, 3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d4']][1:2, 1:2, 1:2, 1:2]), as.numeric(mv[['d4']][[3]][3:4, 3:4, 3:4, 3:4]), info = info ))

            test_that('single', expect_identical(as.numeric(m[['x5']][1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(mv[['x5']][[3]][3:4, 3:4, 3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d5']][1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(mv[['d5']][[3]][3:4, 3:4, 3:4, 3:4, 3:4]), info = info ))

            test_that('single', expect_identical(as.numeric(m[['x6']][1:2, 1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(mv[['x6']][[3]][3:4, 3:4, 3:4, 3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d6']][1:2, 1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(mv[['d6']][[3]][3:4, 3:4, 3:4, 3:4, 3:4, 3:4]), info = info ))

            test_that('single', expect_identical(as.numeric(m[['v1']][1:2]), as.numeric(mv[['v1']][[3]][3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['w1']][1:2, 1:2]), as.numeric(mv[['w1']][[3]][3:4, 3:4]), info = info ))
        })
    ),
    Model2MVsomeDiff_logProbTRUE = list( ## diff = different nodesTo from nodes
        label = 'Model2MVsomeDiff_logProbTRUE',
        nfName = 'nfModelToMVsomeDiff',
        seed = round(runif(1, 1, 10000)),
        compile = c(FALSE, TRUE),
        compareRtoCpp = TRUE,
        nfMcode = quote({nf(m, mv, nodes = c('x0','d0','x1[1:2]','d1[1:2]','x2[1:2,1:2]','d2[1:2, 1:2]',
                                             'x3[1:2,1:2,1:2]','d3[1:2,1:2,1:2]', 'x4[1:2,1:2,1:2,1:2]','d4[1:2,1:2,1:2,1:2]',
                                             'x5[1:2, 1:2, 1:2, 1:2, 1:2]','d5[1:2, 1:2, 1:2, 1:2, 1:2]',
                                             'x6[1:2, 1:2, 1:2, 1:2, 1:2, 1:2]','d6[1:2, 1:2, 1:2, 1:2, 1:2, 1:2]',
                                             'v1[1:2]', 'w1[1:2,1:2]'),
                            nodesTo = c('x0','d0','x1[3:4]','d1[3:4]','x2[3:4,3:4]','d2[3:4, 3:4]',
                                        'x3[3:4,3:4,3:4]','d3[3:4,3:4,3:4]', 'x4[3:4,3:4,3:4,3:4]','d4[3:4,3:4,3:4,3:4]',
                                        'x5[3:4, 3:4, 3:4, 3:4, 3:4]','d5[3:4, 3:4, 3:4, 3:4, 3:4]',
                                        'x6[3:4, 3:4, 3:4, 3:4, 3:4, 3:4]','d6[3:4, 3:4, 3:4, 3:4, 3:4, 3:4]',
                                        'v1[3:4]', 'w1[3:4,3:4]'),
                            logProb = TRUE)}),
        testThatLines = quote({
            info = 'Model2MVsomeDiff_logProbTRUE'
            test_that('single', expect_identical(as.numeric(m[['x0']]), as.numeric(mv[['x0']][[3]]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d0']]), as.numeric(mv[['d0']][[3]]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['x1']])[1:2], as.numeric(mv[['x1']][[3]])[3:4], info = info ))
            test_that('single', expect_identical(as.numeric(m[['d1']])[1:2], as.numeric(mv[['d1']][[3]])[3:4], info = info ))
            test_that('single', expect_identical(as.numeric(m[['x2']][1:2, 1:2]), as.numeric(mv[['x2']][[3]][3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d2']][1:2, 1:2]), as.numeric(mv[['d2']][[3]][3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['x3']][1:2, 1:2, 1:2]), as.numeric(mv[['x3']][[3]][3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d3']][1:2, 1:2, 1:2]), as.numeric(mv[['d3']][[3]][3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['x4']][1:2, 1:2, 1:2, 1:2]), as.numeric(mv[['x4']][[3]][3:4, 3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d4']][1:2, 1:2, 1:2, 1:2]), as.numeric(mv[['d4']][[3]][3:4, 3:4, 3:4, 3:4]), info = info ))

            test_that('single', expect_identical(as.numeric(m[['x5']][1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(mv[['x5']][[3]][3:4, 3:4, 3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d5']][1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(mv[['d5']][[3]][3:4, 3:4, 3:4, 3:4, 3:4]), info = info ))

            test_that('single', expect_identical(as.numeric(m[['x6']][1:2, 1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(mv[['x6']][[3]][3:4, 3:4, 3:4, 3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['d6']][1:2, 1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(mv[['d6']][[3]][3:4, 3:4, 3:4, 3:4, 3:4, 3:4]), info = info ))

            test_that('single', expect_identical(as.numeric(m[['v1']][1:2]), as.numeric(mv[['v1']][[3]][3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['w1']][1:2, 1:2]), as.numeric(mv[['w1']][[3]][3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['logProb_x0']]), as.numeric(mv[['logProb_x0']][[3]]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['logProb_x1']])[1:2], as.numeric(mv[['logProb_x1']][[3]])[3:4], info = info ))
            test_that('single', expect_identical(as.numeric(m[['logProb_x2']][1:2, 1:2]), as.numeric(mv[['logProb_x2']][[3]][3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['logProb_x3']][1:2, 1:2, 1:2]), as.numeric(mv[['logProb_x3']][[3]][3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['logProb_x4']][1:2, 1:2, 1:2, 1:2]), as.numeric(mv[['logProb_x4']][[3]][3:4, 3:4, 3:4, 3:4]) ))
            test_that('single', expect_identical(as.numeric(m[['logProb_x5']][1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(mv[['logProb_x5']][[3]][3:4, 3:4, 3:4, 3:4, 3:4]) ))
            test_that('single', expect_identical(as.numeric(m[['logProb_x6']][1:2, 1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(mv[['logProb_x6']][[3]][3:4, 3:4, 3:4, 3:4, 3:4, 3:4]) ))
            test_that('single', expect_identical(as.numeric(m[['logProb_v1']][1]), as.numeric(mv[['logProb_v1']][[3]][1]), info = info ))
            test_that('single', expect_identical(as.numeric(m[['logProb_w1']][2, 1:2]), as.numeric(mv[['logProb_w1']][[3]][2, 3:4]), info = info )) ## curious case, logProbs go in [2, i]
        })
    )
)

## This list of cases is for copying from a modelValues to a model
copyTestCaseListMVtoModel <- list(
## First 2 are MV2Model with all nodes (default blank nodes arg)
    MV2ModelAll_logProbTRUE = list(
        label = 'MV2ModelAll_logProbTRUE',
        nfName = 'nfMVToModelAll',
        seed = round(runif(1, 1, 10000)),
        initCode = quote({simulate(m); calculate(m); nimCopy(from = m, to = mv, rowTo = 3, logProb = TRUE); simulate(m); calculate(m)}),
        compile = c(FALSE, TRUE),
        compareRtoCpp = TRUE, 
        nfMcode = quote({nf(mv, m, logProb = TRUE)}),
        testThatLines = quote({
            for(oneName in mv$getVarNames(includeLogProb = TRUE)) {
                test_that('single', expect_identical(m[[oneName]], mv[[oneName]][[3]], info = 'MV2ModelAll_logProbTRUE' ))
            }
        })
    ),
    MV2ModelAll_logProbFALSE = list(
        label = 'MV2ModelAll_logProbFALSE',
        nfName = 'nfMVToModelAll',
        seed = round(runif(1, 1, 10000)),
        initCode = quote({simulate(m); calculate(m); nimCopy(from = m, to = mv, rowTo = 3, logProb = TRUE); simulate(m); calculate(m)}),
        compile = c(FALSE, TRUE),
        compareRtoCpp = TRUE,
        nfMcode = quote({nf(mv, m, logProb = FALSE)}),
        testThatLines = quote({
            for(oneName in mv$getVarNames(includeLogProb = FALSE)) {
                test_that('single', expect_identical(m[[oneName]], mv[[oneName]][[3]] ))
            }
            lpNames <- mv$getVarNames(includeLogProb = TRUE)
            lpNames <- lpNames[grep('logProb_', lpNames)]
            for(oneName in lpNames) {
                test_that('single', expect_false(identical(m[[oneName]], mv[[oneName]][[3]])))
            }
        })
    ),
## Next 2 are MV2Model with nodes but not nodesTo
    MV2ModelSomeSame_logProbFALSE = list( ## some = some blocks of nodes, ## same = same nodes from and to 
        label = 'MV2ModelSomeSame_logProbFALSE',
        nfName = 'nfMVToModelSomeSame',
        seed = round(runif(1, 1, 10000)),
        initCode = quote({simulate(m); calculate(m); nimCopy(from = m, to = mv, rowTo = 3, logProb = TRUE); simulate(m); calculate(m)}),
        compile = c(FALSE, TRUE),
        compareRtoCpp = TRUE,
        nfMcode = quote({nf(mv, m, nodes = c('x0','d0','x1[2:3]','d1[2:3]','x2[2:3,2:3]','d2[2:3, 2:3]',
                                             'x3[2:3,2:3,2:3]','d3[2:3,2:3,2:3]', 'x4[2:3,2:3,2:3,2:3]','d4[2:3,2:3,2:3,2:3]',
                                             'x5[2:3, 2:3, 2:3, 2:3, 2:3]','d5[2:3, 2:3, 2:3, 2:3, 2:3]',
                                             'x6[2:3, 2:3, 2:3, 2:3, 2:3, 2:3]','d6[2:3, 2:3, 2:3, 2:3, 2:3, 2:3]',
                                             'v1[2:3]', 'w1[2:3,2:3]'), logProb = FALSE)}),
        testThatLines = quote({
            info <- 'MV2ModelSomeSame_logProbFALSE'
            test_that('single', expect_identical(as.numeric(mv[['x0']][[3]]), as.numeric(m[['x0']]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d0']][[3]]), as.numeric(m[['d0']]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['x1']][[3]])[2:3], as.numeric(m[['x1']])[2:3], info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d1']][[3]])[2:3], as.numeric(m[['d1']])[2:3], info = info ))
            test_that('single', expect_identical(as.numeric(mv[['x2']][[3]][2:3, 2:3]), as.numeric(m[['x2']][2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d2']][[3]][2:3, 2:3]), as.numeric(m[['d2']][2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['x3']][[3]][2:3, 2:3, 2:3]), as.numeric(m[['x3']][2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d3']][[3]][2:3, 2:3, 2:3]), as.numeric(m[['d3']][2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['x4']][[3]][2:3, 2:3, 2:3, 2:3]), as.numeric(m[['x4']][2:3, 2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d4']][[3]][2:3, 2:3, 2:3, 2:3]), as.numeric(m[['d4']][2:3, 2:3, 2:3, 2:3]), info = info ))

            test_that('single', expect_identical(as.numeric(mv[['x5']][[3]][2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(m[['x5']][2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d5']][[3]][2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(m[['d5']][2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))

            test_that('single', expect_identical(as.numeric(mv[['x6']][[3]][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(m[['x6']][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d6']][[3]][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(m[['d6']][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))

            test_that('single', expect_identical(as.numeric(mv[['v1']][[3]][2:3]), as.numeric(m[['v1']][2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['w1']][[3]][2:3, 2:3]), as.numeric(m[['w1']][2:3, 2:3]), info = info ))
        })
       ),
    MV2ModelsomeSame_logProbTRUE = list(
        label = 'MV2ModelSomeSame_logProbTRUE',
        nfName = 'nfMVToModelSomeSame',
        seed = round(runif(1, 1, 10000)),
        initCode = quote({simulate(m); calculate(m); nimCopy(from = m, to = mv, rowTo = 3, logProb = TRUE); simulate(m); calculate(m)}),
        compile = c(FALSE, TRUE),
        compareRtoCpp = TRUE,
        nfMcode = quote({nf(mv, m, nodes = c('x0','d0','x1[2:3]','d1[2:3]','x2[2:3,2:3]','d2[2:3, 2:3]',
                                             'x3[2:3,2:3,2:3]','d3[2:3,2:3,2:3]', 'x4[2:3,2:3,2:3,2:3]','d4[2:3,2:3,2:3,2:3]',
                                             'x5[2:3, 2:3, 2:3, 2:3, 2:3]','d5[2:3, 2:3, 2:3, 2:3, 2:3]',
                                             'x6[2:3, 2:3, 2:3, 2:3, 2:3, 2:3]','d6[2:3, 2:3, 2:3, 2:3, 2:3, 2:3]',
                                             'v1[2:3]', 'w1[2:3,2:3]'), logProb = TRUE)}),
        testThatLines = quote({
            info <- 'MV2ModelSomeSame_logProbTRUE'
            test_that('single', expect_identical(as.numeric(mv[['x0']][[3]]), as.numeric(m[['x0']]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d0']][[3]]), as.numeric(m[['d0']]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['x1']][[3]])[2:3], as.numeric(m[['x1']])[2:3], info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d1']][[3]])[2:3], as.numeric(m[['d1']])[2:3], info = info ))
            test_that('single', expect_identical(as.numeric(mv[['x2']][[3]][2:3, 2:3]), as.numeric(m[['x2']][2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d2']][[3]][2:3, 2:3]), as.numeric(m[['d2']][2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['x3']][[3]][2:3, 2:3, 2:3]), as.numeric(m[['x3']][2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d3']][[3]][2:3, 2:3, 2:3]), as.numeric(m[['d3']][2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['x4']][[3]][2:3, 2:3, 2:3, 2:3]), as.numeric(m[['x4']][2:3, 2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d4']][[3]][2:3, 2:3, 2:3, 2:3]), as.numeric(m[['d4']][2:3, 2:3, 2:3, 2:3]), info = info ))

            test_that('single', expect_identical(as.numeric(mv[['x5']][[3]][2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(m[['x5']][2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d5']][[3]][2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(m[['d5']][2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))

            test_that('single', expect_identical(as.numeric(mv[['x6']][[3]][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(m[['x6']][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d6']][[3]][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(m[['d6']][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))

            test_that('single', expect_identical(as.numeric(mv[['v1']][[3]][2:3]), as.numeric(m[['v1']][2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['w1']][[3]][2:3, 2:3]), as.numeric(m[['w1']][2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_x0']][[3]]), as.numeric(m[['logProb_x0']]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_x1']][[3]])[2:3], as.numeric(m[['logProb_x1']])[2:3], info = info ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_x2']][[3]][2:3, 2:3]), as.numeric(m[['logProb_x2']][2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_x3']][[3]][2:3, 2:3, 2:3]), as.numeric(m[['logProb_x3']][2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_x4']][[3]][2:3, 2:3, 2:3, 2:3]), as.numeric(m[['logProb_x4']][2:3, 2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_x5']][[3]][2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(m[['logProb_x5']][2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_x6']][[3]][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), as.numeric(m[['logProb_x6']][2:3, 2:3, 2:3, 2:3, 2:3, 2:3]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_v1']][[3]][1]), as.numeric(m[['logProb_v1']][1]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_w1']][[3]][2, 2:3]), as.numeric(m[['logProb_w1']][2, 2:3]), info = info ))
        })
    ),
## Next 2 are MV2Model with nodes and nodesTo provided
    MV2ModelSomeDiff_logProbFALSE = list( ## diff = different nodesTo from nodes
        label = 'MV2ModelsomeDiff_logProbFALSE',
        nfName = 'nfMVToModelSomeDiff',
        seed = round(runif(1, 1, 10000)),
        initCode = quote({simulate(m); calculate(m); nimCopy(from = m, to = mv, rowTo = 3, logProb = TRUE); simulate(m); calculate(m)}),
        compile = c(FALSE, TRUE),
        compareRtoCpp = TRUE,
        nfMcode = quote({nf(mv, m, nodes = c('x0','d0','x1[1:2]','d1[1:2]','x2[1:2,1:2]','d2[1:2, 1:2]',
                                             'x3[1:2,1:2,1:2]','d3[1:2,1:2,1:2]', 'x4[1:2,1:2,1:2,1:2]','d4[1:2,1:2,1:2,1:2]',
                                             'x5[1:2, 1:2, 1:2, 1:2, 1:2]','d5[1:2, 1:2, 1:2, 1:2, 1:2]',
                                             'x6[1:2, 1:2, 1:2, 1:2, 1:2, 1:2]','d6[1:2, 1:2, 1:2, 1:2, 1:2, 1:2]',
                                             'v1[1:2]', 'w1[1:2,1:2]'),
                            nodesTo = c('x0','d0','x1[3:4]','d1[3:4]','x2[3:4,3:4]','d2[3:4, 3:4]',
                                        'x3[3:4,3:4,3:4]','d3[3:4,3:4,3:4]', 'x4[3:4,3:4,3:4,3:4]','d4[3:4,3:4,3:4,3:4]',
                                        'x5[3:4, 3:4, 3:4, 3:4, 3:4]','d5[3:4, 3:4, 3:4, 3:4, 3:4]',
                                        'x6[3:4, 3:4, 3:4, 3:4, 3:4, 3:4]','d6[3:4, 3:4, 3:4, 3:4, 3:4, 3:4]',
                                        'v1[3:4]', 'w1[3:4,3:4]'),
                            logProb = FALSE)}),
        testThatLines = quote({
            info <- 'MV2ModelsomeDiff_logProbFALSE'
            test_that('single', expect_identical(as.numeric(mv[['x0']][[3]]), as.numeric(m[['x0']]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d0']][[3]]), as.numeric(m[['d0']]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['x1']][[3]])[1:2], as.numeric(m[['x1']])[3:4], info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d1']][[3]])[1:2], as.numeric(m[['d1']])[3:4], info = info ))
            test_that('single', expect_identical(as.numeric(mv[['x2']][[3]][1:2, 1:2]), as.numeric(m[['x2']][3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d2']][[3]][1:2, 1:2]), as.numeric(m[['d2']][3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['x3']][[3]][1:2, 1:2, 1:2]), as.numeric(m[['x3']][3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d3']][[3]][1:2, 1:2, 1:2]), as.numeric(m[['d3']][3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['x4']][[3]][1:2, 1:2, 1:2, 1:2]), as.numeric(m[['x4']][3:4, 3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d4']][[3]][1:2, 1:2, 1:2, 1:2]), as.numeric(m[['d4']][3:4, 3:4, 3:4, 3:4]), info = info ))

            test_that('single', expect_identical(as.numeric(mv[['x5']][[3]][1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(m[['x5']][3:4, 3:4, 3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d5']][[3]][1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(m[['d5']][3:4, 3:4, 3:4, 3:4, 3:4]), info = info ))

            test_that('single', expect_identical(as.numeric(mv[['x6']][[3]][1:2, 1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(m[['x6']][3:4, 3:4, 3:4, 3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d6']][[3]][1:2, 1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(m[['d6']][3:4, 3:4, 3:4, 3:4, 3:4, 3:4]), info = info ))

            test_that('single', expect_identical(as.numeric(mv[['v1']][[3]][1:2]), as.numeric(m[['v1']][3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['w1']][[3]][1:2, 1:2]), as.numeric(m[['w1']][3:4, 3:4]), info = info ))
        })
    ),
    MV2ModelSomeDiff_logProbTRUE = list( ## diff = different nodesTo from nodes
        label = 'MV2ModelsomeDiff_logProbTRUE',
        nfName = 'nfMVToModelSomeDiff',
        seed = round(runif(1, 1, 10000)),
        initCode = quote({simulate(m); calculate(m); nimCopy(from = m, to = mv, rowTo = 3, logProb = TRUE); simulate(m); calculate(m)}),
        compile = c(FALSE, TRUE),
        compareRtoCpp = TRUE,
        nfMcode = quote({nf(mv, m, nodes = c('x0','d0','x1[1:2]','d1[1:2]','x2[1:2,1:2]','d2[1:2, 1:2]',
                                             'x3[1:2,1:2,1:2]','d3[1:2,1:2,1:2]', 'x4[1:2,1:2,1:2,1:2]','d4[1:2,1:2,1:2,1:2]',
                                             'x5[1:2, 1:2, 1:2, 1:2, 1:2]','d5[1:2, 1:2, 1:2, 1:2, 1:2]',
                                             'x6[1:2, 1:2, 1:2, 1:2, 1:2, 1:2]','d6[1:2, 1:2, 1:2, 1:2, 1:2, 1:2]',
                                             'v1[1:2]', 'w1[1:2,1:2]'),
                            nodesTo = c('x0','d0','x1[3:4]','d1[3:4]','x2[3:4,3:4]','d2[3:4, 3:4]',
                                        'x3[3:4,3:4,3:4]','d3[3:4,3:4,3:4]', 'x4[3:4,3:4,3:4,3:4]','d4[3:4,3:4,3:4,3:4]',
                                        'x5[3:4, 3:4, 3:4, 3:4, 3:4]','d5[3:4, 3:4, 3:4, 3:4, 3:4]',
                                        'x6[3:4, 3:4, 3:4, 3:4, 3:4, 3:4]','d6[3:4, 3:4, 3:4, 3:4, 3:4, 3:4]',
                                        'v1[3:4]', 'w1[3:4,3:4]'),
            logProb = TRUE)}),
        testThatLines = quote({
            info = 'MV2ModelsomeDiff_logProbTRUE'
            test_that('single', expect_identical(as.numeric(mv[['x0']][[3]]), as.numeric(m[['x0']]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d0']][[3]]), as.numeric(m[['d0']]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['x1']][[3]])[1:2], as.numeric(m[['x1']])[3:4], info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d1']][[3]])[1:2], as.numeric(m[['d1']])[3:4], info = info ))
            test_that('single', expect_identical(as.numeric(mv[['x2']][[3]][1:2, 1:2]), as.numeric(m[['x2']][3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d2']][[3]][1:2, 1:2]), as.numeric(m[['d2']][3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['x3']][[3]][1:2, 1:2, 1:2]), as.numeric(m[['x3']][3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d3']][[3]][1:2, 1:2, 1:2]), as.numeric(m[['d3']][3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['x4']][[3]][1:2, 1:2, 1:2, 1:2]), as.numeric(m[['x4']][3:4, 3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d4']][[3]][1:2, 1:2, 1:2, 1:2]), as.numeric(m[['d4']][3:4, 3:4, 3:4, 3:4]), info = info ))

            test_that('single', expect_identical(as.numeric(mv[['x5']][[3]][1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(m[['x5']][3:4, 3:4, 3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d5']][[3]][1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(m[['d5']][3:4, 3:4, 3:4, 3:4, 3:4]), info = info ))

            test_that('single', expect_identical(as.numeric(mv[['x6']][[3]][1:2, 1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(m[['x6']][3:4, 3:4, 3:4, 3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['d6']][[3]][1:2, 1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(m[['d6']][3:4, 3:4, 3:4, 3:4, 3:4, 3:4]), info = info ))

            test_that('single', expect_identical(as.numeric(mv[['v1']][[3]][1:2]), as.numeric(m[['v1']][3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['w1']][[3]][1:2, 1:2]), as.numeric(m[['w1']][3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_x0']][[3]]), as.numeric(m[['logProb_x0']]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_x1']][[3]])[1:2], as.numeric(m[['logProb_x1']])[3:4], info = info ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_x2']][[3]][1:2, 1:2]), as.numeric(m[['logProb_x2']][3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_x3']][[3]][1:2, 1:2, 1:2]), as.numeric(m[['logProb_x3']][3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_x4']][[3]][1:2, 1:2, 1:2, 1:2]), as.numeric(m[['logProb_x4']][3:4, 3:4, 3:4, 3:4]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_x5']][[3]][1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(m[['logProb_x5']][3:4, 3:4, 3:4, 3:4, 3:4]) ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_x6']][[3]][1:2, 1:2, 1:2, 1:2, 1:2, 1:2]), as.numeric(m[['logProb_x6']][3:4, 3:4, 3:4, 3:4, 3:4, 3:4]) ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_v1']][[3]][1]), as.numeric(m[['logProb_v1']][1]), info = info ))
            test_that('single', expect_identical(as.numeric(mv[['logProb_w1']][[3]][2, 1:2]), as.numeric(m[['logProb_w1']][2, 3:4]), info = info ))
        })
    )
)

## This is code that works for all cases to check that the state of compiled and uncompiled models are identical
compareCompiledToUncompiled_code <- quote({
    for(oneName in m$getVarNames(includeLogProb = TRUE)) {
        info <- paste("compareCompiledToUncompiled",oneName)
        test_that('single', expect_identical(as.numeric(origM[[oneName]]), as.numeric(m[[oneName]]), info = info ))
    }
    for(oneName in mv$getVarNames(includeLogProb = TRUE)) {
        info <- paste("compareCompiledToUncompiled",oneName)
        test_that('single', expect_identical(as.numeric(origMV[[oneName]][[3]]), as.numeric(mv[[oneName]][[3]], info = info )))
    }
})


## Function to iterate through test cases
runCopyTests <- function(testCaseList = copyTestCaseList, testModelCode = copyTestModelCode, testNFcodeList = copyTestNFcodeList, dirName = NULL, verbose = nimbleOptions()$verbose) {
    for(copyTestCase in testCaseList) {
        runOneCopyTest(copyTestCase, testModelCode = testModelCode, testNFcodeList = testNFcodeList, dirName = dirName, verbose = verbose)
    }
}

## Function to handle one test case
runOneCopyTest <- function(copyTestCase, testModelCode, testNFcodeList, dirName = NULL, verbose = nimbleOptions()$verbose) {
    if(verbose) writeLines(paste0('Testing ', copyTestCase$label))
    compileVec <- copyTestCase$compile
    for(compile in compileVec) {
        if(verbose) writeLines(paste0('COMPILE = ', compile))
        m <- nimbleModel(testModelCode, constants = copyTestConstants, data = copyTestData)
        mv <- modelValues(m, 3)
        set.seed(copyTestCase$seed)
        eval(copyTestCase$initCode)
        nf <- eval(testNFcodeList[[ copyTestCase$nfName ]])
        nfM <- eval(copyTestCase$nfMcode)
        if(compile) {
            cm <- compileNimble(m, nfM, dirName = dirName)
            m <- cm$m
            nfM <- cm$nfM
            mv <- mv$CobjectInterface
        }
        runans <- try(nfM$run())
        if(copyTestCase$compareRtoCpp & !compile) {
            origM <- m
            origMV <- mv
        }
        if(inherits(runans, 'try-error')) stop(paste('Error executing a copy test from case ', copyTestCase$label))
        eval(copyTestCase$testThatLine)
        if(copyTestCase$compareRtoCpp & compile) {
            eval(compareCompiledToUncompiled_code)
        }
    }
}

###############
## Master call:
runCopyTests(copyTestCaseList, dirName = getwd(), verbose = nimbleOptions()$verbose)
runCopyTests(copyTestCaseListMVtoModel, dirName = getwd(), verbose = nimbleOptions()$verbose)

###################
### Testing for values() and values()<-
### The general layout that follows is similar to above but is completely separate 

## List of nimbleFunction definitions to use
copyTestNFcodeListValues <- list(
    nfGetValues = quote({
        nimbleFunction(
            setup = function(model, nodes){},
            run = function() {
                P <- values(model, nodes)
                return(P)
                returnType(double(1))
            })
    }),
    nfSetValues = quote({
        nimbleFunction(
            setup = function(model, nodes){},
            run = function(P = double(1)) {
                values(model, nodes) <<- P
            })
    })
)

## List of comparison cases, similar to above but without comparing compiled to uncompiled (not applicable)
copyTestCaseListValues <- list(
    getValues = list(
        label = 'getValues',
        nfName = 'nfGetValues',
        seed = round(runif(1, 1, 10000)),
        initCode = quote({simulate(m); calculate(m);
                          nodes <- c('x0','d0','x1[2:3]','d1[2:3]','x2[2:3,2:3]','d2[2:3, 2:3]',
                                       'x3[2:3,2:3,2:3]','d3[2:3,2:3,2:3]', 'x4[2:3,2:3,2:3,2:3]','d4[2:3,2:3,2:3,2:3]',
                                   'v1[2:3]', 'w1[2:3,2:3]')}),
        compile = c(FALSE, TRUE),
        nfMcode = quote({nf(m, nodes = nodes)}),
        testThatLines = quote({
            P <- nfM$run()
            checkP <- numeric()
            for(i in nodes) checkP <- c(checkP, m[[i]])
            test_that('getValues', expect_identical(P, checkP, info = "getValues copyTestCaseListValues"))
        })
    ),
    setValues = list(
        label = 'setValues',
        nfName = 'nfSetValues',
        seed = round(runif(1, 1, 10000)),
        initCode = quote({simulate(m); calculate(m);
                          nodes <- c('x0','d0','x1[2:3]','d1[2:3]','x2[2:3,2:3]','d2[2:3, 2:3]',
                                     'x3[2:3,2:3,2:3]','d3[2:3,2:3,2:3]', 'x4[2:3,2:3,2:3,2:3]','d4[2:3,2:3,2:3,2:3]',
                                     'v1[2:3]', 'w1[2:3,2:3]')}),
        compile = c(FALSE, TRUE),
        nfMcode = quote({nf(m, nodes = nodes)}),
        testThatLines = quote({
            P <- numeric()
            for(i in nodes) P <- c(P, m[[i]]) ## gets checkP up to the right length
            P[] <- rnorm(length(P))
            nfM$run(P)
            i <- 0;
            for(oneName in nodes) {
                modelP <- m[[oneName]]
                checkP <- P[i + 1:length(modelP)]
                i <- i + length(modelP)
                test_that('setValues', expect_identical(as.numeric(modelP), as.numeric(checkP ), info = "setValues copyTestCaseListValues"))
            }
        })
        )
)

## Iterate through the value tests
runValuesTests <- function(testCaseList = copyTestCaseListValues, testModelCode = copyTestModelCode, testNFcodeList = copyTestNFcodeListValues, testNFconstantsList = copyTestConstants, testNFdataList = copyTestData, dirName = NULL, verbose = nimbleOptions()$verbose) {
    for(copyTestCase in testCaseList) {
        runOneValuesTest(copyTestCase, testModelCode = testModelCode, testNFcodeList = testNFcodeList, testNFconstantsList = testNFconstantsList, testNFdataList = testNFdataList, dirName = dirName, verbose = verbose)
    }
}

## run one value test case
runOneValuesTest <- function(copyTestCase, testModelCode, testNFcodeList, testNFconstantsList, testNFdataList, dirName = NULL, verbose = nimbleOptions()$verbose) {
    if(verbose) writeLines(paste0('Testing ', copyTestCase$label))
    compileVec <- copyTestCase$compile
    for(compile in compileVec) {
        if(verbose) writeLines(paste0('COMPILE = ', compile))
        m <- nimbleModel(testModelCode, constants = testNFconstantsList, data = testNFdataList)
        set.seed(copyTestCase$seed)
        eval(copyTestCase$initCode)
        nf <- eval(testNFcodeList[[ copyTestCase$nfName ]])
        nfM <- eval(copyTestCase$nfMcode)
        if(compile) {
            cm <- compileNimble(m, nfM, dirName = dirName)
            m <- cm$m
            nfM <- cm$nfM
        }
        eval(copyTestCase$testThatLine)
    }
}

## Master call for value() and value()<- tests
runValuesTests(copyTestCaseListValues)

###################
### Testing for values() and values()<-  with indexing of node vector inside values()
### The general layout that follows is similar to above but is completely separate 


## List of comparison cases, similar to above but without comparing compiled to uncompiled (not applicable)
copyTestCaseListValuesIndexed <- list(
    getValues = list(
        label = 'getValues',
        nfName = 'nfGetValues',
        seed = round(runif(1, 1, 10000)),
        initCode = quote({simulate(m); calculate(m);
                          nodes <- c('x0','d0','x1[2:3]','d1[2:3]','x2[2:3,2:3]','d2[2:3, 2:3]',
                                       'x3[2:3,2:3,2:3]','d3[2:3,2:3,2:3]', 'x4[2:3,2:3,2:3,2:3]','d4[2:3,2:3,2:3,2:3]',
                                   'v1[2:3]', 'w1[2:3,2:3]')}),
        compile = c(FALSE, TRUE),
        nfMcode = quote({nf(m, nodes = nodes)}),
        testThatLines = quote({
            for(i in seq_along(nodes)) {
                P <- nfM$run(i)
                test_that('getValues', expect_identical(P, as.numeric(m[[nodes[i]]]), info = "getValues copyTestCaseListValuesIndexed"))
            }
        })
    ),
    setValues = list(
        label = 'setValues',
        nfName = 'nfSetValues',
        seed = round(runif(1, 1, 10000)),
        initCode = quote({simulate(m); calculate(m);
                          nodes <- c('x0','d0','x1[2:3]','d1[2:3]','x2[2:3,2:3]','d2[2:3, 2:3]',
                                     'x3[2:3,2:3,2:3]','d3[2:3,2:3,2:3]', 'x4[2:3,2:3,2:3,2:3]','d4[2:3,2:3,2:3,2:3]',
                                     'v1[2:3]', 'w1[2:3,2:3]')}),
        compile = c(FALSE, TRUE),
        nfMcode = quote({nf(m, nodes = nodes)}),
        testThatLines = quote({
            for(i in seq_along(nodes)) {
                oneName <- nodes[i]
                P <- rnorm(length(m[[oneName]]))
                nfM$run(P, i)
                test_that('setValues', expect_identical(as.numeric(m[[oneName]]), P, info = "setValues copyTestCaseListValuesIndexed"))
            }
        })
    )
)

## List of nimbleFunction definitions to use
copyTestNFcodeListValuesIndexed <- list(
    nfGetValues = quote({
        nimbleFunction(
            setup = function(model, nodes){},
            run = function(index = double(0)) {
                P <- values(model, nodes[index])
                return(P)
                returnType(double(1))
            })
    }),
    nfSetValues = quote({
        nimbleFunction(
            setup = function(model, nodes){},
            run = function(P = double(1), index = double(0)) {
                values(model, nodes[index]) <<- P
            })
    })
)

runValuesTests(copyTestCaseListValuesIndexed, testNFcodeList = copyTestNFcodeListValuesIndexed)


## List of comparison cases, similar to above but without comparing compiled to uncompiled (not applicable)
copyTestCaseListValuesIndexedLoop <- list(
    getValues = list(
        label = 'getValues',
        nfName = 'nfGetValues',
        seed = round(runif(1, 1, 10000)),
        initCode = quote({simulate(m); calculate(m);
                          nodes <- m$expandNodeNames('mu')}),
        compile = c(FALSE, TRUE),
        nfMcode = quote({nf(m, nodes = nodes)}),
        testThatLines = quote({
            P <- nfM$run()
            test_that('getValues', expect_identical(as.numeric(P), as.numeric(m[['mu']]), info = "getValues copyTestCaseListValuesIndexedLoop"))
        })
    ),
    setValues = list(
        label = 'setValues',
        nfName = 'nfSetValues',
        seed = round(runif(1, 1, 10000)),
        initCode = quote({simulate(m); calculate(m);
                          nodes <- m$expandNodeNames('mu')}),
        compile = c(FALSE, TRUE),
        nfMcode = quote({nf(m, nodes = nodes)}),
        testThatLines = quote({
            P <- rnorm(length(m$expandNodeNames('mu')))
            nfM$run(P)
            test_that('setValues', expect_identical(as.numeric(P), as.numeric(m[['mu']]), info = "setValues copyTestCaseListValuesIndexedLoop"))
        })
        )
)

copyTestNFcodeListValuesIndexedLoop <- list(
    nfGetValues = quote({
        nimbleFunction(
            setup = function(model, nodes){
                nn <- length(nodes)
            },
            run = function() {
                P <- numeric(nn)
                for(i in 1:nn)
                    P[i] <- values(model, nodes[i])[1]
                return(P)
                returnType(double(1))
            })
    }),
    nfSetValues = quote({
        nimbleFunction(
            setup = function(model, nodes){
                nn <- length(nodes)
            },
            run = function(P = double(1)) {
                tmp <- numeric(1)
                for(i in 1:nn) {
                    tmp[1] <- P[i] 
                    values(model, nodes[i]) <<- tmp
                }
            })
    })
)


copyTestModelCode <- nimbleCode({
    for(i in 1:n) {
        y[i] ~ dnorm(mu[i], 1)
        mu[i] ~ dnorm(0, 1)
    }
})

n <- 10
copyTestConstants <- list(n = n)
copyTestData <- list(y = rnorm(10))

runValuesTests(copyTestCaseListValuesIndexedLoop, testNFcodeList = copyTestNFcodeListValuesIndexedLoop)

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)

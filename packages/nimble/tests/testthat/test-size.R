## This aims to test whether our size error checking is properly reporting
## problems with model specification, including cases that will fail to compile.
## It's rather complicated because of the different ways that
## incorrect/inconsistent sizes can cause problems at model building,
## compilation, and run-time. For now, test_size is somewhat inconsistent
## with how expect_failure and expect_error are used in our other tests.
## Note that numerous error messages are expected here; check for test failures not denoted with 'KNOWN PROBLEM'

source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)


context("Testing of size/dimension checks in NIMBLE code")

goldFileName <- 'sizeTestLog_Correct.Rout'
tempFileName <- 'sizeTestLog.Rout'
generatingGoldFile <- !is.null(nimbleOptions('generateGoldFileForSizeTesting'))
outputFile <- if(generatingGoldFile) file.path(nimbleOptions('generateGoldFileForSizeTesting'), goldFileName) else tempFileName

## capture warnings in gold file
sink_with_messages(outputFile)

vec2 <- c(1,1)
mat2 <- diag(rep(1, 2))
vec3 <- rep(1, 3)
mat3 <- diag(rep(1, 3))
p3 <- rep(1/3, 3)

fit_model <- function(input, check = TRUE, useConst = FALSE) {
    if(useConst) 
            m <- nimbleModel(code = input$expr, data = input$data, constants = input$inits, check = check)
    else m <- nimbleModel(code = input$expr, data = input$data, inits = input$inits, check = check)
}

   
### basic tests of scalar distribution

# use knownProblem = TRUE to indicate case where we know the test does not pass because of shortcomings in Nimble's size checking

# set expectPass based on whether the syntax is allowed in Nimble or not and therefore whether we hope the test passes or not

testsScalar <- list(
    list(name = 'scalar stochastic, basic', expectPass = TRUE,
         expr = quote({y ~ dnorm(mu, sd = sig)}), 
         inits = list(mu = 0, sig = 1) ),
    list(name = 'scalar stochastic, parameter expression', expectPass = TRUE,
         expr = quote({y ~ dnorm(mu1 + mu2, sd = sig)}), 
         inits = list(mu1 = 0, mu2 = 0, sig = 1) ),
    list(name = 'scalar stochastic, parameter expression non-scalar, no indices',
         expectPass = FALSE, expectPassWithConst = TRUE,
         knownProblem = FALSE, knownProblemWithConst = FALSE, expectWarnWithConst = TRUE,
         expr = quote({y ~ dnorm(mu1%*%mu2, sd = sig)}), 
         inits = list(mu1 = vec2, mu2 = vec2, sig = 1) ),
    
    list(name = 'scalar stochastic, parameter expression non-scalar, RHS index',
         expectPass = FALSE, expectPassWithConst = FALSE,
         knownProblem = FALSE, knownProblemWithConst = FALSE, expectWarnWithConst = TRUE,
         expr = quote({y ~ dnorm((mu1%*%mu2)[1,1], sd = sig)}), 
         inits = list(mu1 = vec2, mu2 = vec2, sig = 1) ),
    
    list(name = 'scalar stochastic, parameter expression non-scalar',
         expectPass = FALSE, expectPassWithConst = TRUE,
         knownProblem = TRUE, knownProblemWithConst = FALSE,
         expr = quote({y ~ dnorm(mu1[1:2]%*%mu2[1:2], sd = sig)}), 
         inits = list(mu1 = vec2, mu2 = vec2, sig = 1) ),
    # passes when shouldn't for RHS init (compileNimble puts one into a browser)
    
    list(name = 'scalar stochastic, parameter expression non-scalar, demotion', expectPass = TRUE,
         expr = quote({y ~ dnorm((mu1[1:2]%*%mu2[1:2])[1,1], sd = sig)}), 
         inits = list(mu1 = vec2, mu2 = vec2, sig = 1) ),
    # gives warning with RHS const

    list(name = 'scalar stochastic, non-scalar parameter', expectPass = TRUE,
         expr = quote({y ~ dnorm((mu1[1:2, 1:2]%*%mu2[1:2])[1,1], sd = sig)}), 
         inits = list(mu1 = mat2, mu2 = vec2, sig = 1) ),
    # gives warning with RHS const
    
    list(name = 'scalar stochastic, non-scalar parameter, demotion',
         expectPass = FALSE, expectPassWithConst = TRUE,
         knownProblem = TRUE, knownProblemWithConst = FALSE,
         expr = quote({y ~ dnorm((mu1[1:2, 1:2]%*%mu2[1:2])[1], sd = sig)}), 
         inits = list(mu1 = mat2, mu2 = vec2, sig = 1) ),
    # passes test with RHS inits when it shouldn't, though compileNimble does give helpful error msg; warning with RHS const

    list(name = 'scalar stochastic, non-scalar value', expectPass = FALSE,
         expr = quote({y[1:2] ~ dnorm(mu, sd = sig)}), 
         inits = list(mu = 1, sig = 1) ),
    list(name = 'scalar stochastic, non-scalar value and parameter', expectPass = FALSE,
         expr = quote({y[1:2] ~ dnorm(mu1[1:2, 1:2]%*%mu2[1:2], sd = sig)}), 
         inits = list(mu1 = mat2, mu2 = vec2, sig = 1) ),

    list(name = 'scalar stochastic, scalar within multivar variable', expectPass = TRUE,
    	 expectPassWithConst = FALSE,
         expr = quote({for(i in 1:1)
                           for(j in 1:1)
                               y[i,j] ~ dnorm(mu[i,j], sd = sig)}), 
         inits = list(mu = 0, sig = 1) ),

    list(name = 'scalar stochastic, scalar within multivar variable 2', expectPass = TRUE,
         expr = quote({for(i in 1:1)
                           for(j in 1:1)
                               y[i,j] ~ dnorm(mu[i,j], sd = sig)}), 
         inits = list(mu = matrix(0, 1,1), sig = 1) ),

    list(name = 'scalar stochastic, scalar within multivar variable 3', expectPass = TRUE,
         expr = quote({for(i in 1:1)
                           for(j in 1:1)
                               y[i,j] ~ dnorm(mu[i,j], sd = sig)}), 
         inits = list(mu = mat2, sig = 1) )
)

testsMultivarParam <- list(
    list(name = 'mv param stochastic, basic', expectPass = TRUE,
         expr = quote({y ~ dcat(p[1:3])}), 
         inits = list(p = p3) ),
    list(name = 'mv param stochastic, no indices', expectPass = FALSE,
         expr = quote({y ~ dcat(p)}), 
         inits = list(p = p3) ),
    # with p as constant ERRORS in compileNimble(); not caught in size check, but replaceConstantsRecurse warning regarding dimensionality is given
    # this now has error that dcat is getting prob = model$p[1] so sum of probs is not 1; this error makes sense given sizes problem.
    
    list(name = 'mv param stochastic, scalar param', expectPass = FALSE,
         expr = quote({y ~ dcat(p[1])}), 
         inits = list(p = 1) ),
    
    list(name = 'mv param stochastic, scalar param, no indices', expectPass = FALSE,
         expr = quote({y ~ dcat(p)}), 
         inits = list(p = 1) ),
    list(name = 'mv param stochastic, parameter expression', expectPass = FALSE,
         expr = quote({y ~ dcat(p1[1:3] + p2[1:3])}), 
         inits = list(p1 = p3/2, p2 = p3/2)),
    list(name = 'mv param stochastic, parameter expression, matrices', expectPass = FALSE,
         expr = quote({y ~ dcat(p1[1:3,1:3] %*% p2[1:3])}), 
         inits = list(p1 = matrix(rep(1, 9)/9, 3), p2 = rep(1,3))),
    # parameter is one-column matrix so fails check (and fails in compileNimble)

    list(name = 'mv param stochastic, parameter expression, matrices subindexing', expectPass = FALSE,
         expr = quote({y ~ dcat((p1[1:3,1:3] %*% p2[1:3])[1:3,1])}), 
         inits = list(p1 = matrix(rep(1, 9)/9, 3), p2 = rep(1,3))),
    # note that if one does not check for multivariate param expressions, this errors out in genVarInfo3() - Error in BUGSdecl$symbolicParentNodes[[iV]] : subscript out of bounds

    list(name = 'mv param stochastic, non-scalar value', expectPass = FALSE,
         expr = quote({y[1:2] ~ dcat(p[1:3])}), 
         inits = list(p = p3)),

    list(name = 'mv param stochastic, scalar value within multivar variable ', expectPass = TRUE,
         expr = quote({for(i in 1:1)
                           for(j in 1:1)
                               y[i,j] ~ dcat(p[1:3])}),
         inits = list(p = p3)), 
    list(name = 'mv param stochastic, vector variable within matrix variable ', expectPass = TRUE,
         expr = quote({for(i in 1:1)
                               y[i] ~ dcat(p[1:3, i])}),
         inits = list(p = matrix(p3, ncol = 1)))
)

testsDeterm <- list(
    list(name = 'deterministic, basic', expectPass = TRUE,
         expr = quote({y <- a + b}),
         inits = list(a = 3, b = 3 )),
    
    list(name = 'deterministic, basic with node in variable, incorrect init size', expectPass = TRUE,
    	 expectPassWithConst = FALSE, 
         expr = quote({for(i in 1:1)
                           for(j in 1:1)
                               y[i,j] <- a[i,j] + b}),
         inits = list(a = 3, b = 3 )),
    
    list(name = 'deterministic, basic with node in variable', expectPass = TRUE,
         expr = quote({for(i in 1:1)
                           for(j in 1:1)
                               y[i,j] <- a[i,j] + b}),
         inits = list(a = matrix(3, 1, 1), b = 3 )),

    list(name = 'deterministic, basic with node in variable 2', expectPass = TRUE,
         expr = quote({for(i in 1:2)
                           for(j in 1:2)
                               y[i,j] <- a[i,j] + b}),
         inits = list(a = mat2, b = 3 )),

    list(name = 'deterministic, basic with node in variable  3', expectPass = TRUE,
         expr = quote({for(i in 1:1)
                           for(j in 1:1)
                               y[i,j] <- a[i,j] + b}),
         inits = list(a = mat2, b = 3 )),
    
    list(name = 'deterministic, non-scalar expression, no indices',
         expectPass = FALSE, expectPassWithConst = TRUE, expectWarnWithConst = TRUE,
         expr = quote({y <- a %*% b}),
         inits = list(a = vec2, b = vec2 )),
    # this compiles fine with RHS const, though warning is given during model building
    
    list(name = 'deterministic, non-scalar expression',
         expectPass = FALSE, expectPassWithConst = TRUE,
         knownProblem = TRUE, knownProblemWithConst = FALSE,
         expr = quote({y <- a[1:2] %*% b[1:2]}),
         inits = list(a = vec2, b = vec2 )),

    ## Different from old devel to newNodeFxns for case of constants = input$inits:
    ## In old system, a[1:2] %*% b[1:2] would be evaluated at model definition time so this would pass
    ## In newNodeFxns, non-scalare constants are not baked in or evaluated, so this does not pass.
    ## KNOWN ISSUE shows correctly as TRUE
    list(name = 'deterministic, non-scalar expression with LHS indexing',
         expectPass = FALSE, expectPassWithConst = TRUE,
         knownProblem = TRUE, knownProblemWithConst = FALSE,
         expr = quote({y[1] <- a[1:2] %*% b[1:2]}),
         inits = list(a = vec2, b = vec2 )),
    
    list(name = 'deterministic, non-scalar expression with RHS indexing', expectPass = TRUE,
         expr = quote({y <- (a[1:2] %*% b[1:2])[1]}),
         inits = list(a = vec2, b = vec2 )),
    
    list(name = 'deterministic, non-scalar expression, dimension mismatch',
         expectPass = FALSE,
         knownProblem = TRUE,
         expectWarn = TRUE,
         expr = quote({y <- a[1:2,1:2] %*% b[1:2]}),
         inits = list(a = mat2, b = vec2 )),

    list(name = 'deterministic, vector value, dimension mismatch',
         expectPass = FALSE,
         knownProblem = TRUE,
         expr = quote({y[1:2] <- a + b}),
         inits = list(a = 3, b = 3)),

    list(name = 'deterministic, vector value, missing indices', expectPass = FALSE,
         expectWarnWithConst = TRUE,
         expr = quote({y[1:2] <- a + b}),
         inits = list(a = vec2, b = vec2)),
    # errors for RHS const but not in model_check()


    list(name = 'deterministic, basic vector', expectPass = TRUE,
         expr = quote({y[1:2] <- a[1:2,1:2] %*% b[1:2]}),
         inits = list(a = mat2, b = vec2 )),

    list(name = 'deterministic, basic vector, missing indices',
         expectPass = FALSE, expectWarnWithConst = TRUE,
         expr = quote({y[1:2] <- a %*% b[1:2]}),
         inits = list(a = mat2, b = vec2 )),
    # errors for RHS const but not in model_check()
    
    list(name = 'deterministic, basic vector, dimension mismatch',
         expectPass = FALSE,
         knownProblem = TRUE,
         expr = quote({y[1:2] <- a %*% b}),
         inits = list(a = 3, b = 2 )),
    list(name = 'deterministic, basic vector, RHS dimension mismatch',
         expectPass = FALSE,
         knownProblem = TRUE,
         expectWarn = TRUE,
         expr = quote({y[1:2] <- a[1:2] + b[1:2, 1:2]}),
         inits = list(a = vec2, b = mat2 )),
    list(name = 'deterministic, basic vector, size mismatch',
         expectPass = FALSE,
         knownProblem = TRUE,
         expectWarn = TRUE,
         expr = quote({y[1:3] <- a[1:2,1:2] %*% b[1:2]}),
         inits = list(a = mat2, b = vec2 )),
    list(name = 'deterministic, basic vector, size mismatch 2',
         expectPass = FALSE,
         knownProblem = TRUE,
         expectWarn = TRUE,
         expr = quote({y[1:2] <- a[1:3,1:3] %*% b[1:3]}),
         inits = list(a = mat3, b = rep(1,3) )),
    list(name = 'deterministic, basic vector dimension mismatch',
         expectPass = FALSE,
         knownProblem = TRUE,
         expectWarn = TRUE,
         expr = quote({y[1:3] <- a[1:3,1:3] %*% b[1:3,1:3]}),
         inits = list(a = mat3, b = mat3 )),
    list(name = 'deterministic, basic matrix', expectPass = TRUE,
         expr = quote({y[1:3, 1:3] <- a[1:3,1:3] %*% b[1:3,1:3]}),
         inits = list(a = mat3, b = mat3 )),

    # nested nodes
    list(name = 'deterministic, nodes within multivar variables 1', expectPass = TRUE,
         expr = quote({
             for(i in 1:1)
                 y[1:2, i] <- a[1:2,1:2] %*% b[1:2, i]}),
         inits = list(a = mat2, b = mat2 )),

    ## difference from old devel to newNodeFxns:
    ## This used to generate errors, so expectPass=FALSE was the correct expectation
    ## On newNodeFxns it generates only warnings, so expect=FALSE is the incorrect expectation
    ## It is still dubious in any case, so rather than change expectPass to TRUE I will add knownProblem = TRUE
    list(name = 'deterministic, nodes within multivar variables, size mismatch', expectPass = FALSE,
         knownProblem = TRUE,
         expectWarn = TRUE,
         expr = quote({
             for(i in 1:1)
                 y[1:3] <- a[1:2,1:2] %*% b[1:2, i]}),
         inits = list(a = mat2, b = mat2 )),
    
    list(name = 'deterministic, nodes within multivar variables, dim mismatch', expectPass = FALSE,
         expr = quote({
             for(i in 1:1)
                 y[1:2, i] <- a[1:2,1:2] %*% b[1:2, 1:2]}),
         inits = list(a = mat2, b = mat2 )),

    list(name = 'deterministic, nodes within multivar variables, input wrong dimension',
         expectPass = FALSE,
         expr = quote({
             for(i in 1:1)
                 y[1:2, i] <- a[1:2,1:2] %*% b[1:2, 1:2]}),
         inits = list(a = mat2, b = vec2 ))
)

testsMultivar <- list(
    list(name = 'multivar, basic', expectPass = TRUE,
         expr = quote({
             y[1:2] ~ dmnorm(mu[1:2], prec[1:2, 1:2])
         }),
         inits = list(mu = vec2, prec = mat2)),
    list(name = 'multivar, expression as param', expectPass = FALSE,
         expr = quote({
             y[1:2] ~ dmnorm(mu1[1:2] + mu2[1:2], prec[1:2, 1:2])
         }),
         inits = list(mu1 = vec2, mu2 = vec2, prec = mat2)),
    list(name = 'multivar, nested', expectPass = TRUE,
         expr = quote({
             for(i in 1:1)
                 y[1:2, i] ~ dmnorm(mu[1:2, i], prec[1:2, 1:2])
         }),
         inits = list(mu = mat2, prec = mat2)),
    list(name = 'multivar, value missing index', expectPass = FALSE,
         expr = quote({
                 y ~ dmnorm(mu[1:2], prec[1:2, 1:2])
         }),
         inits = list(mu = vec2, prec = mat2)),
    list(name = 'multivar, param missing index', expectPass = FALSE,
         expr = quote({
                 y[1:2] ~ dmnorm(mu, prec[1:2, 1:2])
         }),
         inits = list(mu = vec2, prec = mat2)),
    # passes when RHS provided as constants, though warning is given
    
    list(name = 'multivar, param wrong dimension', expectPass = FALSE,
         expr = quote({
                 y[1:2] ~ dmnorm(mu[1:2], prec[1:2])
         }),
         inits = list(mu = vec2, prec = vec2)),
    list(name = 'multivar, value wrong size', expectPass = FALSE,
         expr = quote({
                 y[1:3] ~ dmnorm(mu[1:2], prec[1:2, 1:2])
         }),
         inits = list(mu = vec2, prec = mat2)),
    list(name = 'multivar, param wrong size', expectPass = FALSE,
         expr = quote({
                 y[1:2] ~ dmnorm(mu[1:3], prec[1:2, 1:2])
         }),
         inits = list(mu = vec3, prec = mat2)),
    list(name = 'matrix dist, basic', expectPass = TRUE,
         expr = quote({
                 y[1:2,1:2] ~ dwish(R[1:2,1:2], df)
         }),
         inits = list(df = 4, R = mat2)),
    list(name = 'matrix dist, param wrong size', expectPass = FALSE,
         expr = quote({
                 y[1:2,1:2] ~ dwish(R[1:3,1:3], df)
         }),
         inits = list(df = 4, R = mat3)),
    list(name = 'matrix dist, param wrong size 2', expectPass = FALSE,
         expr = quote({
                 y[1:2,1:2] ~ dwish(R[1:2, 1:3], df)
         }),
         inits = list(df = 4, R = matrix(1, 2, 3))),
    list(name = 'matrix dist, value missing index', expectPass = FALSE,
         expr = quote({
                 y ~ dwish(R[1:2, 1:2], df)
         }),
         inits = list(df = 4, R = mat2))
)

testsTrunc <- list(
    list(name = 'multivar, trunc', expectPass = TRUE,
         expr = quote({
                 y ~ T(dnorm(mu, sd = sd), -1, 1)
         }),
         inits = list(mu = 0, sd = 1), data = list(y = 0)),
        list(name = 'multivar, trunc, extra index', expectPass = FALSE,
         expr = quote({
                 y ~ T(dnorm(mu[1:2], sd = sd), 0, 1)
         }),
         inits = list(mu = vec2, sd = 1), data = list(y = 0))
)


testsLHSRHSmismatch <- list(
    list(name = 'LHS 0 dim, RHS 1 dim, determ',
         correctErrorMsg = 'which does not match its usage',
         expr = quote({
             beta ~ dnorm(0, 1)
             x <- beta[1] * z + beta[2] * z
             z ~ dnorm(0, 1)
         })
         ),
    list(name = 'LHS 0 dim, RHS 1 dim, stoch', 
         correctErrorMsg = 'which does not match its usage',
         expr = quote({
             beta ~ dnorm(0, 1)
             x ~ dnorm(beta[1] * z + beta[2] * z, 1)
             z ~ dnorm(0, 1)
         })
         ),
    list(name = 'LHS 0 dim, RHS 2 dim, determ', 
         correctErrorMsg = 'which does not match its usage',
         expr = quote({
             beta ~ dnorm(0, 1)
             x <- beta[1, 1] * z + beta[2] * z
             z ~ dnorm(0, 1)
         })
         ),
    list(name = 'LHS 0 dim, RHS 2 dim, stoch', 
         correctErrorMsg = 'which does not match its usage',
         expr = quote({
             beta ~ dnorm(0, 1)
             x ~ dnorm(beta[1, 1] * z + beta[2] * z, 1)
             z ~ dnorm(0, 1)
         })
         ),
    list(name = 'LHS 1 dim, RHS 2 dim, determ', 
         correctErrorMsg = 'which does not match its usage',
         expr = quote({
             beta[1:3] ~ dmnorm(mu[1:3], prec[1:3, 1:3])
             x <- beta[1, 1] * z + beta[2] * z
             z ~ dnorm(0, 1)
         }),
         inits = list(mu = vec3, prec = mat3)
         ),
    list(name = 'LHS 1 dim, RHS 2 dim, stoch', 
         correctErrorMsg = 'which does not match its usage',
         expr = quote({
             beta[1:3] ~ dmnorm(mu[1:3], prec[1:3, 1:3])
             x ~ dnorm(beta[1, 1] * z + beta[2] * z, 1)
             z ~ dnorm(0, 1)
         }),
         inits = list(mu = vec3, prec = mat3)
         ),
    list(name = 'LHS 1 dim, RHS 0 dim, determ', 
         correctErrorMsg = 'which does not match its usage',
         expr = quote({
             x <- beta * z 
             beta[1:3] ~ dmnorm(mu[1:3], prec[1:3, 1:3])
             z ~ dnorm(0, 1)
         }),
         inits = list(mu = vec3, prec = mat3)
         ),
    list(name = 'LHS 1 dim, RHS 0 dim, stoch', 
         correctErrorMsg = 'which does not match its usage',
         expr = quote({
             x ~ dnorm(beta * z, 1)
             beta[1:3] ~ dmnorm(mu[1:3], prec[1:3, 1:3])
             z ~ dnorm(0, 1)
         }),
         inits = list(mu = vec3, prec = mat3)
         ),
    list(name = 'LHS 2 dim, RHS 0 dim, determ', 
         correctErrorMsg = 'which does not match its usage',
         expr = quote({
             y[1:2, 1:2] ~ dwish(R[1:2, 1:2], df)
             x <- y + 3
         }),
         inits = list(df = 4, R = mat2)
         ),  
    list(name = 'LHS 2 dim, RHS 0 dim, stoch', 
         correctErrorMsg = 'which does not match its usage',
         expr = quote({
             y[1:2, 1:2] ~ dwish(R[1:2, 1:2], df)
             x ~ dnorm( y + 3, 1)
         }),
         inits = list(df = 4, R = mat2)
         )    ,
    list(name = 'LHS 2 dim, RHS 1 dim, determ', 
         correctErrorMsg = 'which does not match its usage',
         expr = quote({
             y[1:2, 1:2] ~ dwish(R[1:2, 1:2], df)
             x[1:2] <- y[1:2] + 3
         }),
         inits = list(df = 4, R = mat2)
         ),  
    list(name = 'LHS 2 dim, RHS 1 dim, stoch', 
         correctErrorMsg = 'which does not match its usage',
         expr = quote({
             y[1:2, 1:2] ~ dwish(R[1:2, 1:2], df)
             x[1:2] ~ dmnorm(y[1:2], R[1:2, 1:2])
         }),
         inits = list(df = 4, R = mat2)
         )    
)

# for these we simply want errors caught at model-building not compilation
out <- sapply(testsScalar, test_size)
out <- sapply(testsMultivarParam, test_size)
out <- sapply(testsDeterm, test_size)
out <- sapply(testsMultivar, test_size)
out <- sapply(testsTrunc, test_size)

# for these we want errors caught by NIMBLE error-checking not by R errors at model-building
out <- sapply(testsLHSRHSmismatch, test_size_specific_error)

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

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of getParam")

gpScalar <- nimbleFunction(
    setup = function(model, node, param) {},
    run = function() {
        ans1 <- model$getParam(node, param)
        ans2 <- getParam(model, node, param) ## to become model$getParam(node, param)
        if(ans1 != ans2) stop('oops, ans1 != ans2')
        return(ans1)
        returnType(double())
    })

testGetParam <- function(distCall) {
    dist <- nimble:::getDistributionInfo(as.character(distCall[[1]]))
    code <- substitute({x ~ DISTCALL}, list(DISTCALL = distCall))
    m <- nimbleModel( code = code )
    cm <- compileNimble(m)
    gpFuns <- list()
    expectedResults <- list()
    altParams <- dist$altParams
    altParamNames <- names(altParams)
    distCallText <- deparse(distCall)

    reqdArgs <- dist$reqdArgs ## these are canonical
    exprs <- dist$exprs
    alts <- dist$alts
    providedArgs <- names(distCall)
    providedArgs <- providedArgs[providedArgs != ""]
    whichExpr <- NULL
    ## figure out which way arguments were provided in distCall
    for(i in seq_along(exprs)) {
        if(all(providedArgs %in% alts[[i]])) whichExpr <- i
    }
    if(is.null(whichExpr)) {
        if(all(providedArgs %in% reqdArgs)) whichExpr <- 0
    }
    test_that(paste(distCallText, 'args found'), expect_equal(is.null(whichExpr), FALSE))

    
    ## exprs give expressions for calculating reqdArgs from alts

    ## altParams give expressions for calculating individual alts from reqdArgs
    
    ## put reqd in evalEnv, which means using exprs for the alts as needed
    ## if testing on something provided, grab what was provided.
    ## if testing on something not provided, if it is reqd then use it directly
    ## otherwise calculate it from altParams

    evalEnv <- new.env()

    for(i in seq_along(distCall)) {
        if(names(distCall)[i] != "") assign(names(distCall)[i], distCall[[i]], envir = evalEnv)
    }
    if(whichExpr > 0) {  ## what was provided was not canonical
        for(i in seq_along(exprs[[whichExpr]])) {
            assign(names(exprs[[whichExpr]])[i], eval(exprs[[whichExpr]][[i]], envir = evalEnv), envir = evalEnv)
        }
    }

    ## check recovery of alternative param names from what was provided
    for(i in seq_along(altParamNames)) {
        gpFuns[[i]] <- gpScalar(m, 'x', altParamNames[i])
        if(altParamNames[i] %in% providedArgs) ## it was provided so simply eval the name
            expectedResults[[i]] <- eval(as.name(altParamNames[i]), envir = evalEnv)
        else  ## it wasn't provided so eval the expression to calculate it from reqdArgs
            expectedResults[[i]] <- eval(altParams[[i]], envir = evalEnv)
        test_that(paste(distCallText, 'uncompiled', altParamNames[i]), expect_equal(gpFuns[[i]]$run(), expectedResults[[i]]))
        # test use of getParam in R
        test_that(paste(distCallText, 'Rmodel', altParamNames[i]), expect_equal(m$getParam('x', altParamNames[i]), expectedResults[[i]]))
        test_that(paste(distCallText, 'Cmodel', altParamNames[i]), expect_equal(cm$getParam('x', altParamNames[i]), expectedResults[[i]]))
    }

    resultsNames <- altParamNames
    nextI <- length(expectedResults)+1
    for(i in seq_along(reqdArgs)) {
        gpFuns[[nextI]] <- gpScalar(m, 'x', reqdArgs[i])
        expectedResults[[nextI]] <- eval(as.name(reqdArgs[i]), envir = evalEnv) ## it was already calculated into evalEnv above
        test_that(paste(distCallText, 'uncompiled reqd', reqdArgs[i]), expect_equal(gpFuns[[nextI]]$run(), expectedResults[[nextI]]))
        # test use of getParam in R
        test_that(paste(distCallText, 'Rmodel reqd', reqdArgs[i]), expect_equal(m$getParam('x', reqdArgs[i]), expectedResults[[nextI]]))
        test_that(paste(distCallText, 'Cmodel reqd', reqdArgs[i]), expect_equal(cm$getParam('x', reqdArgs[i]), expectedResults[[nextI]]))
        resultsNames[nextI] <- reqdArgs[i]
        nextI <- nextI + 1
    }
    
    compiled <- do.call('compileNimble', c(list(m), gpFuns, list(resetFunctions = TRUE)))
    for(i in seq_along(expectedResults)) {
        test_that(paste(distCallText, 'compiled', resultsNames[i]), expect_equal(compiled[[i+1]]$run(), expectedResults[[i]]))
    }
}

testGetParam(quote(dbern(prob = 0.2)))
testGetParam(quote(dbin(prob = 0.2, size = 3)))
##testGetParam(quote(dbinom(prob = 0.2, size = 3)))
testGetParam(quote(dnegbin(prob = 0.2, size = 3)))
##testGetParam(quote(dnbinom(prob = 0.2, size = 3)))
testGetParam(quote(dpois(lambda = 2.5)))
testGetParam(quote(dbeta(shape1 = 1.5, shape2 = 2.5)))
testGetParam(quote(dbeta(mean = .6, sd = .05)))
testGetParam(quote(dchisq(df = 3)))
testGetParam(quote(dexp(rate = .3)))
testGetParam(quote(dexp(scale = 3)))
testGetParam(quote(dgamma(shape = 2, scale = 1.5)))
testGetParam(quote(dgamma(shape = 2, rate = 3.0)))
testGetParam(quote(dgamma(mean = 2.0, sd = 1.5)))
testGetParam(quote(dlnorm(meanlog = 2.0, taulog = 1.5)))
testGetParam(quote(dlnorm(meanlog = 2.0, sdlog = .8)))
testGetParam(quote(dlnorm(meanlog = 2.0, varlog = .6)))
testGetParam(quote(dlogis(location = 1.5, rate = .2)))
testGetParam(quote(dlogis(location = 1.5, scale = 5.0)))
testGetParam(quote(dnorm(mean = 10.5, sd = 1.5)))
testGetParam(quote(dnorm(mean = 10.5, var = 1.5)))
testGetParam(quote(dnorm(mean = 10.5, tau = 1.5)))
testGetParam(quote(dt(df = 3, mu = 1.5, tau = 0.9)))
testGetParam(quote(dt(df = 3, mu = 1.5, sigma2 = 1.1)))
testGetParam(quote(dt(df = 3, mu = 1.5, sigma = 1.2)))
testGetParam(quote(dunif(min = 1.2, max = 1.3)))
testGetParam(quote(dweib(shape = 1.2, scale = 1.3)))
testGetParam(quote(dweib(shape = 1.2, rate = 1.3)))
testGetParam(quote(dweib(shape = 1.2, lambda = 1.3)))

## We haven't written an extensive version of testing getParam for non-scalar parameters
## However the following covers testing that the size processing and eigenization steps work with getParam.

testCode <- nimbleCode({
    for(i in 1:3) x[i] ~ dnorm(0, 1)
    y[1:3] ~ dmnorm(x[1:3], mycov[1:3, 1:3])
})

y <- rnorm(3)
mycov <- diag(3)
testModel <- nimbleModel(testCode, data = list(y = y, mycov = diag(3)))
x <- rnorm(3)
testModel$x <- x

nf <- nimbleFunction(
    setup = function(model, mvNode){},
    run = function() {
        ans <- model$getParam(mvNode, 'mean')
        return(ans)
        returnType(double(1))
    },
    methods = list(
        test2 = function() {
            ans <- 1.1 + model$getParam(mvNode, 'mean')
            return(ans)
            returnType(double(1))
        },
        test3 = function(z = double(1)) {
            ans <- z + model$getParam(mvNode, 'mean')
            return(ans)
            returnType(double(1))
        })
)

nf1 <- nf(testModel, 'y[1:3]')
test_that('multivar 1', expect_equivalent(nf1$run(), testModel$x))
test_that('multivar 2', expect_equivalent(nf1$test2(), testModel$x + 1.1))
test_that('multivar 3', expect_equivalent(nf1$test3(11:13), testModel$x + 11:13))

Ctest <- compileNimble(testModel, nf1)
test_that('multivar 4', expect_equivalent(Ctest$nf1$run(), Ctest$testModel$x))
test_that('multivar 5', expect_equivalent(Ctest$nf1$test2(), Ctest$testModel$x + 1.1))
test_that('multivar 6', expect_equivalent(Ctest$nf1$test3(11:13), Ctest$testModel$x + 11:13))

# basic non-scalar test
code = nimbleCode({
    a[1:3] ~ dmnorm(mu[1:3],pr[1:3,1:3])
})
pr1 = diag(3)
pr1[1,2]=pr1[2,1]=.3
pr2 <- pr1
pr1[1,2]=pr1[2,1]=.5

m = nimbleModel(code, inits =list(mu=rep(1,3), pr = pr1))
cm = compileNimble(m)

cm$pr <- pr2
cm$calculate(cm$getDependencies('pr'))

test_that('non-scalar 1', expect_equal(pr1, m$getParam('a', 'prec')))
test_that('non-scalar 2', expect_equal(pr2, cm$getParam('a', 'prec')))
test_that('non-scalar 3', expect_equal(solve(pr1), m$getParam('a', 'cov')))
test_that('non-scalar 4', expect_equal(solve(pr2), cm$getParam('a', 'cov')))


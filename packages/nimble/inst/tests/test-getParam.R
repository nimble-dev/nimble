source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of getParam")

gpScalar <- nimbleFunction(
    setup = function(model, node, param) {},
    run = function() {
        ans1 <- getParam(model, node, param)
        ans2 <- getParam(model, node, param) ## to become model$getParam(node, param)
        if(ans1 != ans2) stop('oops, ans1 != ans2')
        return(ans1)
        returnType(double())
    })

testGetParam <- function(distCall) {
    dist <- nimble:::getDistribution(as.character(distCall[[1]]))
    code <- substitute({x ~ DISTCALL}, list(DISTCALL = distCall))
    m <- nimbleModel( code = code )
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
    test_that(paste(distCallText, 'args found'), expect_that(is.null(whichExpr), equals(FALSE)))

    
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
        test_that(paste(distCallText, 'uncompiled', altParamNames[i]), expect_that(gpFuns[[i]]$run(), equals(expectedResults[[i]])))
    }

    resultsNames <- altParamNames
    nextI <- length(expectedResults)+1
    for(i in seq_along(reqdArgs)) {
        gpFuns[[nextI]] <- gpScalar(m, 'x', reqdArgs[i])
        expectedResults[[nextI]] <- eval(as.name(reqdArgs[i]), envir = evalEnv) ## it was already calculated into evalEnv above
        test_that(paste(distCallText, 'uncompiled reqd', reqdArgs[i]), expect_that(gpFuns[[nextI]]$run(), equals(expectedResults[[nextI]])))
        resultsNames[nextI] <- reqdArgs[i]
        nextI <- nextI + 1
    }
    
    compiled <- do.call('compileNimble', c(list(m), gpFuns, list(resetFunctions = TRUE)))
    for(i in seq_along(expectedResults)) {
        test_that(paste(distCallText, 'compiled', resultsNames[i]), expect_that(compiled[[i+1]]$run(), equals(expectedResults[[i]])))
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


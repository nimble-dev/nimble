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
    
    evalEnv <- new.env()
    for(i in seq_along(distCall)) {
        if(names(distCall)[i] != "") assign(names(distCall[i]), distCall[[i]], envir = evalEnv)
    }

    ## check recovery of alternative param names from what was provided
    for(i in seq_along(altParamNames)) {
        gpFuns[[i]] <- gpScalar(m, 'x', altParamNames[i])
        expectedResults[[i]] <- eval(altParams[[i]], envir = evalEnv)
        test_that(paste(distCallText, 'uncompiled', altParamNames[i]), expect_that(gpFuns[[i]]$run(), equals(expectedResults[[i]])))
    }

    reqdArgs <- dist$reqdArgs ## these are canonical
    exprs <- dist$exprs
    alts <- dist$alts
    providedArgs <- ls(evalEnv)
    whichExpr <- NULL
    ## figure out which way arguments were provided in distCall
    for(i in seq_along(exprs)) {
        if(all(providedArgs %in% alts[[i]])) whichExpr <- i
    }
    if(is.null(whichExpr)) {
        if(all(providedArgs %in% reqdArgs)) whichExpr <- 0
    }
    test_that(paste(distCallText, 'args found'), expect_that(is.null(whichExpr), equals(FALSE)))

    resultsNames <- altParamNames
    ## The chosenExpr contains expressions to get canonicals from what was provided
    if(whichExpr != 0) chosenExpr <- exprs[[whichExpr]]
    else chosenExpr <- list()
    nextI <- length(expectedResults)+1
    for(i in seq_along(reqdArgs)) {
        gpFuns[[nextI]] <- gpScalar(m, 'x', reqdArgs[i])
        if(reqdArgs[i] %in% names(chosenExpr)) { ## There is a calculation 
            expectedResults[[nextI]] <- eval(chosenExpr[[i]], envir = evalEnv)
        } else { ## otherwise the canonical is part of what was provided
            expectedResults[[nextI]] <- eval(as.name(reqdArgs[i]), envir = evalEnv)
        }
        test_that(paste(distCallText, 'uncompiled reqd', reqdArgs[i]), expect_that(gpFuns[[nextI]]$run(), equals(expectedResults[[nextI]])))
        resultsNames[nextI] <- reqdArgs[i]
        nextI <- nextI + 1
    }
    
    compiled <- do.call('compileNimble', c(list(m), gpFuns, list(resetFunctions = TRUE)))
    for(i in seq_along(expectedResults)) {
        test_that(paste(distCallText, 'compiled', resultsNames[i]), expect_that(compiled[[i+1]]$run(), equals(expectedResults[[i]])))
    }
}

testGetParam(quote(debern(prob = 0.2)))
debug(testGetParam)
testGetParam(quote(dgamma(shape = 2, scale = 1.5)))
testGetParam(quote(dgamma(shape = 2, rate = 3.0)))

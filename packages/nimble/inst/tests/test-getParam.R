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
    
    for(i in seq_along(altParamNames)) {
        gpFuns[[i]] <- gpScalar(m, 'x', altParamNames[i])
        expectedResults[[i]] <- eval(altParams[[i]], envir = evalEnv)
        test_that(paste(distCallText, 'uncompiled', altParamNames[i]), expect_that(gpFuns[[i]]$run(), equals(expectedResults[[i]])))
    }

    compiled <- do.call('compileNimble', c(list(m), gpFuns, list(resetFunctions = TRUE)))
    for(i in seq_along(altParamNames)) {
        test_that(paste(distCallText, 'compiled', altParamNames[i]), expect_that(compiled[[i+1]]$run(), equals(expectedResults[[i]])))
    }    
}

testGetParam(quote(dgamma(shape = 2, scale = 1.5)))

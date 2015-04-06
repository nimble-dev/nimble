## a "context" means the for-loops in which the BUGS line was nested
BUGSsingleContextClass <- setRefClass('BUGSsingleContextClass',
                                      fields = list(
                                          indexVarExpr = 'ANY',
                                          indexRangeExpr = 'ANY',
                                          forCode = 'ANY'
                                      )
)

## The singleContexts field is a list of BUGSsingleContextClass objects
BUGScontextClass <- setRefClass('BUGScontextClass',
                                
                                fields = list(
                                    replacementsEnv = 'ANY',
                                    ##### all fields are set in setup(), and never change
                                    singleContexts = 'ANY', ## a list of BUGSsingleContextClass objects
                                    indexVarExprs = 'ANY', ## a list of index variable expressions
                                    indexVarNames = 'ANY' ## vector of index variable names (character)
                                    ),
                                
                                methods = list(
                                    setup             = function() {},
                                    genIndexVarValues = function() {}
                                )
)

## master entry function.
## sets all fields, which never change.
BUGScontextClass$methods(setup = function(singleContexts) {
    singleContexts <<- singleContexts
    indexVarExprs <<- lapply(singleContexts, function(x) x$indexVarExpr)
    indexVarNames <<- if(length(indexVarExprs)>0)      unlist(lapply(indexVarExprs, as.character))      else      character(0)
})

BUGScontextClass$methods(genIndexVarValues = function(constantsEnvCopy) {
    genIndexVarValues_recurse(singleContexts, constantsEnvCopy)
})

genIndexVarValues_recurse <- function(singleContexts, constantsEnvCopy) {
    if(length(singleContexts) == 0)   return(list(list()))
    
    indexExpr <- singleContexts[[1]]$indexVarExpr
    indexName <- as.character(indexExpr)
    rangeExpr <- singleContexts[[1]]$indexRangeExpr
    
    ## this changes the behaviour foe expnding looping ranges (L:U),
    ## in particular in the strange cases when L > U
    if(is.call(rangeExpr) && rangeExpr[[1]] == ':') {
        rangeValueL <- eval(rangeExpr[[2]], envir = constantsEnvCopy)
        rangeValueU <- eval(rangeExpr[[3]], envir = constantsEnvCopy)
        if(rangeValueL <= rangeValueU) {
            rangeValues <- rangeValueL:rangeValueU
        } else {
            optionValue <- nimbleOptions$processBackwardsModelIndexRanges
            if(optionValue)  { rangeValues <- rangeValueU:rangeValueL      ## for(i in 9:7) --> for(i in c(9, 8, 7))
            } else           { rangeValues <- numeric(0)               }   ## for(i in 9:7) --> for(i in numeric(0))
        }
    } else {
        warning('loop range expression in BUGS model not of the form (L):(U).  This may be ok; depending on how we\'ve extended NIMBLE')
        rangeValues <- eval(rangeExpr, envir = constantsEnvCopy)
    }
    
    indexVarValues <- list()
    for(value in rangeValues) {
        assign(indexName, value = value, envir = constantsEnvCopy)
        indexVarValuesNew <- genIndexVarValues_recurse(singleContexts[-1], constantsEnvCopy)
        indexVarValuesNew <- lapply(indexVarValuesNew, function(l) { l<-rev(l); l[[indexName]]<-value;  rev(l) })
        indexVarValues <- c(indexVarValues, indexVarValuesNew)
    }
    return(indexVarValues)
}

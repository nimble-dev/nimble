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
                                    genIndexVarValues = function() {},
                                    embedCodeInForLoop = function() {}
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

BUGScontextClass$methods(embedCodeInForLoop = function(innerLoopCode, useContext = NULL, allowNegativeIndexSequences = NULL) {
    ## innerLoopCode is code to be embedded in (possibly nested) for-loops from this context
    ## useContext is an optional logical vector of which contexts to include
    ## allowNegativeIndexSequences: if TRUE, for(i in 2:1) results in iterating over c(2,1), as R would.  If FALSE (Default),
    ##      behavior is like BUGS: for(i in 2:1) results in no iteration.
    if(is.null(allowNegativeIndexSequences))
        allowNegativeIndexSequences <- if(is.null(getNimbleOption('processBackwardsModelIndexRanges'))) TRUE
                                       else getNimbleOption('processBackwardsModelIndexRanges')
 
    if(is.null(useContext)) {
        useContext <- rep(TRUE, length(singleContexts))
    }
    iContext <- length(singleContexts)
    while(iContext >= 1) {
        if(useContext[iContext]) {
            newCode <- singleContexts[[iContext]]$forCode
            if(!allowNegativeIndexSequences) {
                indexRangeCode <- newCode[[3]]
                isColonExpr <- if(is.name(indexRangeCode[[1]])) if(as.character(indexRangeCode[[1]])==':') TRUE else FALSE else FALSE
                if(isColonExpr) newCode[[3]] <- as.call(list(as.name('nm_seq_noDecrease'), indexRangeCode[[2]], indexRangeCode[[3]]))
            }
            newCode[[4]] <- innerLoopCode
            innerLoopCode <- newCode
        }
        iContext <- iContext - 1
    }
    innerLoopCode
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
            optionValue <- getNimbleOption('processBackwardsModelIndexRanges')
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

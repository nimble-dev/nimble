## The BUGSdeclClass contains the pulled-apart content of a BUGS declaration line

## nimbleOrRfunctionNames is used to determine what (in BUGS code) can be evaluated in R if every argument is known OR in C++ (nimble) if arguments are other nodes
nimblePreevaluationFunctionNames <- c('+',
                                      '-',
                                      '/',
                                      '*',
                                      'exp',
                                      'log',
                                      'pow',
                                      'pow_int',
                                      '^',
                                      '%%',
                                      'equals',
                                      'nimEquals',
                                      'sqrt',
                                      'logit',
                                      'expit',
                                      'ilogit',
                                      'probit',
                                      'iprobit',
                                      'phi',
                                      'cloglog',
                                      'icloglog',
                                      'step',
                                      'nimStep',
                                      'sin',
                                      'cos',
                                      'tan',
                                      'asin',
                                      'acos',
                                      'atan',
                                      'cosh',
                                      'sinh',
                                      'tanh',
                                      'asinh',
                                      'acosh',
                                      'atanh',
                                      'cube',
                                      'abs',
                                      'lgamma',
                                      'loggam',
                                      'log1p',
                                      'lfactorial',
                                      'besselK',
                                      'ceiling',
                                      'floor',
                                      'round',
                                      'nimRound',
                                      'trunc',
                                      '>',
                                      '<',
                                      '>=',
                                      '<=',
                                      '==',
                                      '!=',
                                      '[',
                                      '(',
                                      '%*%',
                                      't',
                                      'inprod',
                                      'optim',
                                      'nimOptim',
                                      'optimDefaultControl',
                                      'nimOptimDefaultControl',
                                      'mean',
                                      'sum',
                                      'sd',
                                      'var',
                                      'max',
                                      'min',
                                      'pmin',
                                      'pmax',
                                      'prod',
                                      'asRow',
                                      'asCol',
                                      'logdet',    
                                      'chol',
                                      'inverse',
                                      'forwardsolve',
                                      'backsolve',
                                      'solve',
                                      'nimEigen',
                                      'nimSvd',  
                                      '&',
                                      '|',
                                      '$',
                                      det_distributionFuns,
                                        # these are allowed in DSL as special
                                        # cases even though exp_nimble and
                                        # t_nonstandard are the canonical NIMBLE
                                        # distribution functions
                                      paste0(c('d','q','p'), 't'),
                                      paste0(c('d','q','p'), 'exp'),
                                      'nimC', 'nimRep', 'nimSeq', 'diag',
                                      'nimNumeric','nimMatrix','nimArray',
                                      'length'
                                      )

nimbleOrRfunctionNames <- c(nimblePreevaluationFunctionNames,
                            distribution_rFuns,
                            paste0(c('r'), 't'),
                            paste0(c('r'), 'exp')
                            )

functionsThatShouldNeverBeReplacedInBUGScode <- c(':','nimC','nimRep','nimSeq', 'diag',
                                                  'nimNumeric', 'nimMatrix', 'nimArray')

#' BUGSdeclClass contains the information extracted from one BUGS declaration
BUGSdeclClass <- setRefClass(
    'BUGSdeclClass',
    
    fields = list(###### the following are set in setup(), and never change.
        contextID = 'ANY', sourceLineNumber = 'ANY', code = 'ANY',        ## original BUGS code line
        type = 'ANY', distributionName = 'ANY', targetExpr = 'ANY',   ## LHS of code
        valueExpr = 'ANY',    ## RHS of code
        transExpr = 'ANY', indexExpr = 'ANY', targetVarExpr = 'ANY',
        targetNodeExpr = 'ANY', targetVarName = 'ANY',
        targetNodeName = 'ANY', ## bounds/truncation information
        truncated = 'ANY', boundExprs = 'ANY', ## set in setIndexVariableExprs(), and never changes.
        indexVariableExprs = 'ANY',
        
        ## set in genSymbolicParentNodes(), and never changes.
        symbolicParentNodes = 'ANY',
        
##### the following are set in genReplacementsAndCodeReplaced(), and
##### never change, with one exception
        replacements = 'ANY',
        codeReplaced = 'ANY',  ## MODIFIED in
                               ## genAltParamsModifyCodeReplaced() to
                               ## remove the .params
        replacementNameExprs = 'ANY',
        logProbNodeExpr = 'ANY',
        
##### the following is set in genAltParamsModifyCodeReplaced(), and
##### never changes
        altParamExprs = 'ANY',   ## contains the "smart" expression
                                 ## from each altParam: from original
                                 ## parameterization, whenever
                                 ## possible
        
        ## the following are set in modelDefClass$genNodeInfo(), and
        ## never change:
        indexedNodeInfo = 'ANY',
        indexedNodeNames = 'ANY',
        
        targetExprReplaced = 'ANY',
        valueExprReplaced = 'ANY',
        symbolicParentNodesReplaced = 'ANY',
        rhsVars = 'ANY',
        targetIndexNamePieces = 'ANY',
        parentIndexNamePieces = 'ANY',
        replacementsEnv = 'ANY',
        nodeFunctionNames = 'ANY',
        
        dynamicIndexInfo = 'ANY', ## store info on dynamicIndex
                                  ## expressions for use in
                                  ## restricting to valid index range
                                  ## in nodeFunctions
        
        outputSize = 'ANY', ## should match nrow(unrolledIndicesMatrix)
        origIDs = 'ANY',
        graphIDs = 'ANY',
        unrolledIndicesMatrix = 'ANY',
        numUnrolledNodes = 'ANY', ## differs from outputSize ONLY for a
                                  ## no-context singleton, a ~ dnorm(b,
                                  ## c), so numUnrolledNodes is 1 but
                                  ## outputSize is 0
        envir = 'ANY' ## environment from which nimbleModel called,
                      ## set here rather than in modelDef because
                      ## environment is used in calls to functions from BUGSdeclClass.
    ),   
    
    methods = list(
        setup                          = function() {},
        setIndexVariableExprs          = function() {},
        genSymbolicParentNodes         = function() {},
        genReplacementsAndCodeReplaced = function() {},
        genAltParamsModifyCodeReplaced = function() {},
        genBounds                      = function() {},
        
        genReplacedTargetValueAndParentInfo = function() {},
        
        allParentVarNames = function() {
            unlist(lapply(symbolicParentNodes, function(x) if(length(x) == 1) as.character(x) else if(x[[1]] == '[') as.character(x[[2]]) else stop('Error in allParentVarNames')))
        },
        
        allTargetNodeNames = function() {
            unlist(lapply(indexedNodeInfo, `[[`, 'targetNodeName'))
        },
        
        allParentNodeNames = function() {
            unique(unlist(lapply(indexedNodeInfo, `[[`, 'parentNodeNames')))
        },
        
        allParentNodeExprs = function() {
            unique(unlist(lapply(indexedNodeInfo, `[[`, 'parentNodeExprs')))
        },
        allEdges = function() {
            edgesIn  <- unlist(
                lapply(
                    indexedNodeInfo, 
                    function(x) {
                        expandedNodeIndices <-
                            nl_vectorizedExpandNodeIndexExprs(x$parentNodeExprs)
                        if(length(expandedNodeIndices) > 0) 
                            strsplit(paste(expandedNodeIndices, 
                                           x$targetNodeName, 
                                           sep = '&',
                                           collapse = '&'),
                                     '&')
                        else NULL
                    }
                )
            )
            edgesOut <- unlist(
                lapply(indexedNodeInfo, 
                       function(x)
                           if(is.vectorized(x$targetNodeName)) 
                               strsplit(
                                   paste(x$targetNodeName, 
                                         nl_expandNodeIndexExpr(x$targetNodeExpr), 
                                         sep = '&',
                                         collapse = '&'),
                                   '&')
                           else NULL))
            return(c(edgesIn,
                     edgesOut))
        },
        getDistributionName = function() {
            return(distributionName) 
        },
        isTruncated = function() {
            return(truncated) 
        }
    )
)


BUGSdeclClass$methods(
    setup = function(code,
                     contextID,
                     sourceLineNum,
                     truncated = FALSE,
                     boundExprs = NULL,
                     userEnv = .GlobalEnv) {
        ## This is the master entry function.
        ## Argument 'contextID' is used to set field: contextID.
        ## Argument 'code' is used to set the fields:
        ##  code
        ##  targetExpr, valueExpr
        ##  targetVarExpr, targetNodeExpr
        ##  targetVarName, targetNodeName
        
        contextID <<- contextID
        sourceLineNumber <<- sourceLineNum
        code <<- code
        truncated <<- truncated
        boundExprs <<- boundExprs
        envir <<- userEnv
        
        if(code[[1]] == '~') {
            type <<- 'stoch'
            
            if(!is.call(code[[3]]) ||
               (!any(code[[3]][[1]] == getAllDistributionsInfo('namesVector')) &&
                code[[3]][[1]] != "T" &&
                code[[3]][[1]] != "I"))
                stop(
                    paste0('Improper syntax for stochastic declaration: ',
                           safeDeparse(code))
                )
        } else if(code[[1]] == '<-') {
            type <<- 'determ'
        } else {
            stop(paste0('Improper syntax for declaration: ',
                        safeDeparse(code))
                 )
        }
        
        targetExpr <<- code[[2]]
        valueExpr <<- code[[3]]

        if(type == 'stoch')
            distributionName <<- as.character(valueExpr[[1]])
        else
            distributionName <<- NA

        transExpr <<- NULL
        indexExpr <<- NULL
        
        if(length(targetExpr) > 1) {
            ## There is a tranformation and/or a subscript
            if(targetExpr[[1]] == '[') {
                ## It is a subscript only
                indexExpr <<- as.list(targetExpr[-c(1,2)]) 
                targetVarExpr <<- targetExpr[[2]]
                targetNodeExpr <<- targetExpr
            } else {
                ## There is a transformation, possibly with a subscript
                transExpr <<- targetExpr[[1]]
                targetNodeExpr <<- targetExpr[[2]]
                if(length(targetNodeExpr)>1) {
                    ## There are subscripts inside the transformation
                    if(targetNodeExpr[[1]] != '[') {
                        print(paste("Invalid subscripting for",
                                    safeDeparse(targetExpr)))
                    }
                    indexExpr <<- as.list(targetNodeExpr[-c(1,2)])
                    targetVarExpr <<- targetNodeExpr[[2]]
                } else {
                    targetVarExpr <<- targetNodeExpr
                }
            }
        } else {
            ## no tranformation or subscript
            targetVarExpr <<- targetExpr
            targetNodeExpr <<- targetVarExpr
        }
        
        targetVarName <<- safeDeparse(targetVarExpr, warn = TRUE)
        targetNodeName <<- safeDeparse(targetNodeExpr, warn = TRUE)
    }
)


BUGSdeclClass$methods(
    setIndexVariableExprs = function(exprs) {
        indexVariableExprs <<- exprs
    }
)


BUGSdeclClass$methods(
    genSymbolicParentNodes =
        function(constantsNamesList,
                 context,
                 nimFunNames,
                 unknownIndexDeclInfo = NULL,
                 contextID = NULL,
                 buildDerivs = FALSE) {
    ## sets the field symbolicparentNodes
            symbolicParentNodes <<-
                unique(
                    getSymbolicParentNodes(valueExpr,
                                           constantsNamesList,
                                           context$indexVarExprs,
                                           nimFunNames,
                                           contextID = contextID,
                                           envir = envir,
                                           buildDerivs = buildDerivs)
                ) 
        }
)

stripParentheses <- function(code) {
    if(is.call(code)) {
        if(code[[1]] == "(")
            return(stripParentheses(code[[2]]))
    }
    code
}

## move this to a util file when everything is working.  It is convenient here for now
makeIndexNamePieces <- function(indexCode) {
    indexCode <- stripParentheses(indexCode)
    if(getNimbleOption('allowDynamicIndexing')) {
        if(length(indexCode) == 1)
            return(
                if(is.numeric(indexCode))
                    indexCode
                else
                    as.character(indexCode))
        ## It is easiest to have indexNamePieces be NA when
        ## dynamically indexed rather than retaining the indexing
        ## code.
        if(length(indexCode) == 2 &&
           indexCode[[1]] == ".DYN_INDEXED")
            return(as.numeric(NA))
    } else
        if(length(indexCode) == 1)
            return(
                if(is.numeric(indexCode))
                    indexCode
                else
                    as.character(indexCode))
    ## Diagnostic for messed up indexing here
    if(as.character(indexCode[[1]] != ':'))
        stop(paste0("Error processing model: something is wrong with the index ",
                    safeDeparse(indexCode),
                    ".\nIndexing in model code requires this syntax: '(start expression):(end expression)'."),
             call. = FALSE)
    p1 <- indexCode[[2]]
    p2 <- indexCode[[3]]
    list(
        if(is.numeric(p1))
            p1
        else
            as.character(p1),
        if(is.numeric(p2))
            p2
        else
            as.character(p2))
    ## e.g. makeIndexNamePieces(quote(i))
    ##      makeIndexNamePieces(quote(i:100))
}

BUGSdeclClass$methods(
    genReplacedTargetValueAndParentInfo = function(constantsNamesList,
                                                   context,
                                                   nimFunNames,
                                                   contextID = NULL,
                                                   buildDerivs = FALSE) {
        ## This assumes codeReplaced is there.
        ## Generate hasBracket info:
        targetExprReplaced <<- codeReplaced[[2]]
        ## targetExprReplaced shouldn't have any link functions at this point.
        valueExprReplaced <<- codeReplaced[[3]]
        if(type == 'stoch')
            distributionName <<- as.character(valueExprReplaced[[1]])
        else
            distributionName <<- NA
    
        symbolicParentNodesReplaced <<-
            unique(
                getSymbolicParentNodes(valueExprReplaced,
                                       constantsNamesList,
                                       c(context$indexVarExprs,
                                         replacementNameExprs),
                                       nimFunNames,
                                       contextID = contextID,
                                       envir = envir,
                                       buildDerivs = buildDerivs)
            )
    if(!nimbleOptions()$allowDynamicIndexing) {
        rhsVars <<-
            unlist(
                lapply(
                    symbolicParentNodesReplaced,
                    function(x) 
                        if(length(x) == 1)
                            as.character(x)
                        else
                            as.character(x[[2]])
                )
            )
    } else {
        ## This use of symbolicParentNodes and not
        ## symbolicParentNodesReplaced deals with fact that 'd' is
        ## inserted in front of digits in symbolicParentNodesReplaced
        ## in naming when we have something like k[9-i] but not in
        ## symbolicParentNodes or in varInfo names.
        rhsVars <<- unlist(
            lapply(
                symbolicParentNodes[seq_along(symbolicParentNodesReplaced)],
                function(x) {
                    x <- stripIndexWrapping(x) ## handles dynamic index wrapping
                    if(length(x) == 1) as.character(x) else as.character(x[[2]])
                })
        )
    }

    ## Note that makeIndexNamePieces is designed only for indices that
    ##     are a single name or number, a `:` operator with single
    ##     name or number for each argument, or an NA (for a dynamic
    ##     index).  This relies on the fact that any expression will
    ##     have been lifted by this point and what it has been
    ##     replaced with is simply a name.  This means
    ##     makeIndexNamePieces can include a diagnostic.
        targetIndexNamePieces <<-
            try(
                if(length(targetExprReplaced) > 1)
                    lapply(targetExprReplaced[-c(1,2)],
                           makeIndexNamePieces)
                else
                    NULL
            )
        if(inherits(targetIndexNamePieces, 'try-error'))
            stop(paste('Error occurred defining ',
                       safeDeparse(targetExprReplaced)),
                 call. = FALSE)
    if(!nimbleOptions()$allowDynamicIndexing) {
        parentIndexNamePieces <<-
            lapply(symbolicParentNodesReplaced,
                   function(x)
                       if(length(x) > 1)
                           lapply(x[-c(1,2)],
                                  makeIndexNamePieces)
                       else
                           NULL
                   )
    } else
        parentIndexNamePieces <<-
            lapply(symbolicParentNodesReplaced,
                   function(x) {
                       x <- stripIndexWrapping(x)
                       if(length(x) > 1)
                           lapply(x[-c(1,2)], makeIndexNamePieces)
                       else
                           NULL
                   }
                   )
    NULL
    }
)

BUGSdeclClass$methods(
    genReplacementsAndCodeReplaced = function(constantsNamesList,
                                              context,
                                              nimFunNames) {
        replacementsAndCode <-
            genReplacementsAndCodeRecurse(code,
                                          c(constantsNamesList,
                                            context$indexVarExprs),
                                          nimFunNames,
                                          envir = envir)
        replacements <<- replacementsAndCode$replacements
        codeReplaced <<- replacementsAndCode$codeReplaced
        
        if(type == 'determ')
            logProbNodeExpr <<- NULL
        if(type == 'stoch') {
            logProbNodeExprAndReplacements <-
                genLogProbNodeExprAndReplacements(code,
                                                  codeReplaced,
                                                  context$indexVarExprs)
            logProbNodeExpr <<-
                logProbNodeExprAndReplacements$logProbNodeExpr
            replacements <<-
                c(replacements,
                  logProbNodeExprAndReplacements$replacements)
    }
    
        replacementNameExprs <<-
            lapply(
                as.list(names(replacements)),
                as.name
            )
        names(replacementNameExprs) <<- names(replacements)
})

## genAltParamsModifyCodeReplaced only affects stochastic nodes.  It
## removes any params in codeReplaced which begin with '.'.  It
## generates the altParamExprs list, which contains the expression for
## each alternate parameter, which is taken from the .param expression
## (no longer taken from getDistributionInfo(distName)$altParams,
## ever).
BUGSdeclClass$methods(
    genAltParamsModifyCodeReplaced = function() {
        
        altParamExprs <<- list()
        
        if(type == 'stoch') {
            RHSreplaced <- codeReplaced[[3]]
            if(length(RHSreplaced) > 1) { ## It actually has argument(s)
                paramNamesAll <- names(RHSreplaced)
                paramNamesDotLogicalVector <- grepl('^\\.', paramNamesAll)
                ## remove all parameters whose name begins with '.'
                ## from distribution:
                RHSreplacedWithoutDotParams <-
                    RHSreplaced[!paramNamesDotLogicalVector]
                codeReplaced[[3]] <<- RHSreplacedWithoutDotParams
                
                altParamExprs <<-
                    if(any(paramNamesDotLogicalVector))
                        as.list(RHSreplaced[paramNamesDotLogicalVector])
                    else
                        list()
                ## remove the '.' from each name:
                names(altParamExprs) <<-
                    gsub('^\\.', '', names(altParamExprs))    
        }
    }
})

## genBounds only affects stochastic nodes.  It generates the
## boundExprs list, which contains the expression for both lower and
## upper, which is taken from the lower and upper expression.  If not
## truncated, it removes 'lower' and 'upper' from codeReplaced as not
## needed for standard nodeFunctions (and in nodeFunction creation, it
## checks for presence of lower,upper to indicate truncation since
## only has access to RHS code not to full declInfo)
BUGSdeclClass$methods(
    genBounds = function() {
        boundExprs <<- list()
        if(type == 'stoch') {
            RHSreplaced <- codeReplaced[[3]]
            if(length(RHSreplaced) > 1) { ## It actually has argument(s)
                boundNames <- c('lower_', 'upper_')
                boundExprs <<- as.list(RHSreplaced[boundNames])
                if(truncated) {  # check for user-provided constant bounds inconsistent with distribution range
                    distName <- as.character(RHSreplaced[[1]])
                    distRange <- getDistributionInfo(distName)$range
                    if(is.numeric(boundExprs$lower_) &&
                       is.numeric(distRange$lower) &&
                       is.numeric(boundExprs$upper_) &&
                       is.numeric(distRange$upper) &&
                       boundExprs$lower_ <= distRange$lower &&
                       boundExprs$upper_ >= distRange$upper)  # user specified bounds irrelevant
                        truncated <<- FALSE
                    
                    if(is.numeric(boundExprs$lower_) &&
                       is.numeric(boundExprs$upper_) &&
                       boundExprs$lower_ >= boundExprs$upper_)
                        warning(paste0("Lower bound is greater than or equal to upper bound in ",
                                       safeDeparse(codeReplaced),
                                       "; proceeding anyway, but this is likely to cause numerical issues."))
                    if(is.numeric(boundExprs$lower_) &&
                       is.numeric(distRange$lower) &&
                       boundExprs$lower_ < distRange$lower) {
                        warning(paste0("Lower bound is less than or equal to distribution lower bound in ",
                                       safeDeparse(codeReplaced), "; ignoring user-provided lower bound."))
                        boundExprs$lower_ <<- distRange$lower
                        codeReplaced[[3]]['lower_'] <<- distRange$lower
                    }
                    if(is.numeric(boundExprs$upper_) &&
                       is.numeric(distRange$upper) &&
                       boundExprs$upper_ > distRange$upper) {
                        warning(paste0("Upper bound is greater than or equal to distribution upper bound in ",
                                       safeDeparse(codeReplaced),
                                       "; ignoring user-provided upper bound."))
                        boundExprs$upper_ <<- distRange$upper
                        codeReplaced[[3]]['upper_'] <<- distRange$upper
                    }
                }
                if(!truncated) {
                    boundNamesLogicalVector <-
                        names(RHSreplaced) %in% boundNames
                    RHSreplacedWithoutBounds <-
                        RHSreplaced[!boundNamesLogicalVector]    
                    codeReplaced[[3]] <<- RHSreplacedWithoutBounds
                }
            }
        }
    }
)

getSymbolicParentNodes <- function(code,
                                   constNames = list(),
                                   indexNames = list(),
                                   nimbleFunctionNames = list(),
                                   addDistNames = FALSE,
                                   contextID = NULL,
                                   envir = .GlobalEnv,
                                   buildDerivs = FALSE) {
    if(addDistNames)
        nimbleFunctionNames <- c(nimbleFunctionNames,
                                 getAllDistributionsInfo('namesExprList'))
    ans <- getSymbolicParentNodesRecurse(code,
                                         constNames,
                                         indexNames,
                                         nimbleFunctionNames,
                                         contextID,
                                         envir,
                                         buildDerivs = buildDerivs)
    return(ans$code)
}

getSymbolicParentNodesRecurse <- function(code, constNames = list(), indexNames = list(),
                                          nimbleFunctionNames = list(), contextID = NULL,
                                          envir = .GlobalEnv, buildDerivs = FALSE) {
    ## This takes as input some code and returns the variables in it.
    ## It expects one line of code, not a '{' expression.
    ##
    ## However, indexNames (from for-loop indices) and constNames are
    ## not identified as separate variables.  e.g. x[i] is returned as
    ## 'x' and 'i' or as 'x[i]' if i is an indexName
    ##
    ## indexNames and constNames can be substituted at compile time,
    ## such as a block index variable.
    ##
    ## Every function EXCEPT those in nimbleFunctionNames can be
    ## evaluated at compile time.
    ##
    ## constNames, indexNames and nimbleFunctionNames should be lists
    ## of names.
    ##
    ## details: each recursion returns a list with:
    ## - code: a list of symbolicParentExprs
    ## - replaceable: logical of whether it can be part of a partially
    ##                evaluated expression.
    ##                This includes numbers, constants, indices and
    ##                functions that can be evaluated in R.
    ##                replacements aren't actually done but are used to
    ##                decide handling.
    ##                Something replaceable doesn't need to become a
    ##                symbolicParentNode.
    ##                Something replaceable in an index represents static
    ##                indexing, not dynamic indexing
    ## - hasIndex: is there an index inside
    ## numeric constant
    if(is.numeric(code) || is.logical(code) || 
       (nimbleOptions()$allowDynamicIndexing &&
                       length(code) > 1 &&
                       code[[1]] == ".DYN_INDEXED")
       ) 
        ## Check for .DYN_INDEXED deals with processing of code when
        ## we add unknownIndex declarations.
        return(list(code = NULL,
                    replaceable = TRUE,
                    hasIndex = FALSE))
    cLength <- length(code)
    ## a single name:
    if(cLength == 1) {
        if(is.name(code)) {
            ## Is this for a blank index?, e.g. from first index of
            ## x[, j].  At this point indices have been filled so
            ## there shouldn't be blanks.
            if(code == ''){
                return(list(code = NULL,
                            replaceable = TRUE,
                            hasIndex = FALSE))
            }
            ## an index name
            if(any(code == indexNames)) {
                return(list(code = NULL,
                            replaceable = TRUE,
                            hasIndex = TRUE))
            }
            ## a constant name
            if(any(code == constNames)) {
                return(list(code = NULL,
                            replaceable = TRUE,
                            hasIndex = FALSE))
            }
            ## just something regular: not constant or index
            return(list(code = list(code),
                        replaceable = FALSE,
                        hasIndex = FALSE))
        }
    }

    ## a call:
    if(is.call(code)) {
        indexingBracket <- code[[1]] == '['
        if(indexingBracket) {
            if(is.call(code[[2]])){
                ## a case like foo(x)[i] (when will this occur in BUGS?), 
              indexingBracket <- FALSE
            } 
        }
        if(indexingBracket) {
            ## recurse on the index arguments
            contents <-
                lapply(code[-c(1,2)],
                       function(x)
                           getSymbolicParentNodesRecurse(x,
                                                         constNames,
                                                         indexNames,
                                                         nimbleFunctionNames,
                                                         contextID,
                                                         envir,
                                                         buildDerivs = buildDerivs)
                       )
            ## unpack the codes returned from recursion
            contentsCode <-
                unlist(
                    lapply(contents,
                           function(x)
                               x$code),
                    recursive = FALSE)
            ## unpack whether each index has an index
            contentsHasIndex <-
                unlist(lapply(contents,
                              function(x) x$hasIndex))
            ## unpack whether each index is replaceable
            contentsReplaceable <-
                unlist(lapply(contents,
                              function(x) x$replaceable))
            ## recuse on the variable, e.g. mu in mu[i]
            variable <-
                getSymbolicParentNodesRecurse(code[[2]],
                                              constNames,
                                              indexNames,
                                              nimbleFunctionNames,
                                              contextID,
                                              envir,
                                              buildDerivs = buildDerivs)
            
            ## error if it looks like mu[i][j] where i is a for-loop index
            if(variable$hasIndex)
                stop('Error: Variable',
                     safeDeparse(code[[2]]),
                     'on outside of [ contains a BUGS code index.')
            
            if(variable$replaceable) {
                ## a case like x[ block[i] ], dealing with the
                ## block[i], so block is replaceable
                if(!all(contentsReplaceable)) 
                    ## dynamic index on a constant
                    stop('getSymbolicParentNodesRecurse: dynamic indexing of constants is not allowed in ', safeDeparse(code), '. Try adding the dynamically-indexed constant as data instead (using the data argument of nimbleModel).')
                boolIndexingBlock <-
                    unlist(
                        lapply(code[-c(1,2)],
                               function(x)
                                   if(length(x) > 1)
                                       if(x[[1]] == ':')
                                           TRUE
                                       else
                                           FALSE
                                   else
                                       FALSE)
                    )
                
                if(any(boolIndexingBlock)) {
                    return(list(code = c(contentsCode, code),
                                replaceable = FALSE,
                                hasIndex = any(contentsHasIndex)))
                } else {
                    return(list(code = contentsCode, 
                                replaceable = all(contentsReplaceable),
                                hasIndex = any(contentsHasIndex)))
                }
            } else {
                ## x[i] with x a variable and no dynamic indices
                if(all(contentsReplaceable)) {
                    return(list(code = c(contentsCode, list(code)),
                                replaceable = FALSE,
                                hasIndex = any(contentsHasIndex)))
                } else { ## non-replaceable indices are dynamic indices (or constant vectors, which are not allowed)
                    if(!nimbleOptions()$allowDynamicIndexing) {
                        message("  [Note] It appears you are trying to use dynamic indexing (i.e., the index of a variable is determined by something that is not a constant) in: `",
                                safeDeparse(code),
                                "`. Please set `nimbleOptions(allowDynamicIndexing = TRUE)`.")
                        dynamicIndexParent <- code[[2]]
                    } else {
                        if(isTRUE(nimbleOptions("doADerrorTraps")))
                          if(isTRUE(buildDerivs))
                            message("  [Warning] Derivatives cannot currently be built for models that include dynamic indexing (found in `", safeDeparse(code), "`).  Please set 'nimbleOptions(buildDerivs = FALSE)' to proceed with this model.")
                      
                        if(any(
                            sapply(contentsCode,
                                   detectNonscalarIndex))
                           )
                            stop("getSymbolicParentNodesRecurse: only scalar indices are allowed; vector indexing found in ",
                                 safeDeparse(code))
                        indexedVariable <- safeDeparse(code[[2]], warn = TRUE)
                        dynamicIndexParent <-
                            addUnknownIndexToVarNameInBracketExpr(code, contextID)
                        ## Instead of inserting NA, leave indexing
                        ## code but with indication it is a dynamic
                        ## index so we can detect that later. We need
                        ## the indexing code so we can add it to
                        ## declInfo$dynamicIndexInfo for range
                        ## checking.
                        dynamicIndexParent[-c(1, 2)][ !contentsReplaceable ] <- 
                            lapply(dynamicIndexParent[-c(1, 2)][ !contentsReplaceable ],
                                   addDynamicallyIndexedWrapping)
                        contentsCode = lapply(
                            contentsCode,
                            addIndexWrapping)
                    }
                    return(list(code = c(contentsCode,
                                         list(dynamicIndexParent)),
                                replaceable = FALSE,
                                hasIndex = any(contentsHasIndex)))
                }
            }
        } else {
            ## a regular call like foo(x)
            if(cLength > 1) {
                if(code[[1]] == '$') ## a$x: recurse on a 
                    contents <- lapply(
                        code[2],
                        function(x)
                            getSymbolicParentNodesRecurse(x,
                                                          constNames,
                                                          indexNames,
                                                          nimbleFunctionNames,
                                                          contextID,
                                                          envir,
                                                          buildDerivs = buildDerivs)
                    )
                else ## foo(x): recurse on x
                    contents <- lapply(
                        code[-1],
                        function(x)
                            getSymbolicParentNodesRecurse(x,
                                                          constNames,
                                                          indexNames,
                                                          nimbleFunctionNames,
                                                          contextID,
                                                          envir,
                                                          buildDerivs = buildDerivs)
                    )
                ## unpack results of recursion
                contentsCode <- unlist(
                    lapply(contents,
                           function(x) x$code),
                    recursive = FALSE)
                ## unpack hasIndex entries
                contentsHasIndex <- unlist(
                    lapply(contents,
                           function(x) x$hasIndex)
                )
                ## unpack replaceable entries
                contentsReplaceable <- unlist(
                    lapply(contents,
                           function(x) x$replaceable)
                )
                allContentsReplaceable <- all(contentsReplaceable)
            } else { ## no arguments: foo()
                contentsCode <- NULL
                contentsHasIndex <- FALSE
                allContentsReplaceable <- TRUE
            }
            ## check if the function can be called only in R, not NIMBLE
            isRfunction <- !any(code[[1]] == nimbleFunctionNames)
            funName <- safeDeparse(code[[1]], warn = TRUE)
            isRonly <- isRfunction &
                (!checkNimbleOrRfunctionNames(funName, envir))
            ## if it can be called only in R but not all contents are replaceable, generate error:
            if(isRonly & !allContentsReplaceable) {
                if(!exists(funName, envir))
                    stop("R function '", funName,"' in the code '", safeDeparse(code), "' does not exist.")
                if(funName == ":") ## dynamic indexing in a vector of indices
                    stop("Dynamic indexing found in a vector of indices, ", safeDeparse(code), ". Only scalar indices, such as 'idx' in 'x[idx]', can be dynamic. One can instead use dynamic indexing in a vector of indices inside a nimbleFunction.") 
                unreplaceable <-
                    sapply(contents[!contentsReplaceable],
                           function(x) as.character(x$code)
                           )
                stop("R function '",
                     funName,
                     "' has arguments that cannot be evaluated in the code '", safeDeparse(code), "'. Either the function must be a nimbleFunction or values for the following inputs must be specified as constants in the model: ",
                     paste(unreplaceable, collapse = ","),
                     ".")
            }
            return(list(code = contentsCode,
                        replaceable = allContentsReplaceable & isRfunction,
                        hasIndex = any(contentsHasIndex)))
        }
    }
    stop(paste('Something went wrong in getSymbolicVariablesRecurse with',
               safeDeparse(code)))
}

checkNimbleOrRfunctionNames <- function(functionName, envir) {
    if(any(functionName == nimbleOrRfunctionNames))
        return(TRUE)
    if(exists(functionName, envir) &&
       is.rcf(get(functionName, envir)))
        return(TRUE)  ## Would like to do this by R's scoping rules here and in genCpp_sizeProcessing but that is problematic
    return(FALSE)
}


## The replaceVariableLHS arg avoids replacement of x[i] if x[i] is on the LHS of <- or ~
genReplacementsAndCodeRecurse <- function(code,
                                          constAndIndexNames,
                                          nimbleFunctionNames,
                                          replaceVariableLHS = TRUE,
                                          debug = FALSE,
                                          envir = .GlobalEnv) {
    if(debug) browser()
    if(is.numeric(code) || is.logical(code) ||
       (nimbleOptions()$allowDynamicIndexing &&
                       length(code) > 1 &&
                       code[[1]] == '.DYN_INDEXED')
       )
        ## Check for .DYN_INDEXED deals with processing of code when
        ## we add unknownIndex declarations.
        return(list(codeReplaced = code,
                    replacements = list(),
                    replaceable = TRUE))
    cLength <- length(code)
    if(cLength == 1) {
        if(is.name(code)) {
            if(any(code == constAndIndexNames) &
               replaceVariableLHS)
                return(replaceAllCodeSuccessfully(code))
            else
                return(list(codeReplaced = code,
                            replacements = list(),
                            replaceable = FALSE))
        }
    }
    if(is.call(code)) {
        indexingBracket <- code[[1]] == '['
        if(indexingBracket) {
            if(is.call(code[[2]])) indexingBracket <- FALSE ## treat like any other function
        }
        if(indexingBracket) { 
            contents <- lapply(
                code[-c(1,2)],
                function(x)
                    genReplacementsAndCodeRecurse(x,
                                                  constAndIndexNames,
                                                  nimbleFunctionNames,
                                                  debug = debug,
                                                  envir = envir)
            )
            contentsCodeReplaced <-
                lapply(contents, function(x) x$codeReplaced)
            contentsReplacements <-
                lapply(contents, function(x) x$replacements)
            contentsReplaceable  <-
                unlist(lapply(contents, function(x) x$replaceable))
            if(replaceVariableLHS) {
                variable <-
                    genReplacementsAndCodeRecurse(code[[2]],
                                                  constAndIndexNames,
                                                  nimbleFunctionNames,
                                                  debug = debug,
                                                  envir = envir)
                if(variable$replaceable &&
                   all(contentsReplaceable))
                    return(replaceAllCodeSuccessfully(code))
            }
            return(replaceWhatWeCan(code,
                                    contentsCodeReplaced,
                                    contentsReplacements,
                                    contentsReplaceable,
                                    startingAt=3))
        }
        assignment <- any(code[[1]] == c('<-', '~'))
        if(cLength > 1) {
            if(assignment) {
                ## In an assignment, prevent the outermost variable on
                ## the LHS from being replaced.
                contents <-
                    c(
                        list(
                            genReplacementsAndCodeRecurse(code[[2]],
                                                          constAndIndexNames,
                                                          nimbleFunctionNames,
                                                          replaceVariableLHS = FALSE, debug,
                                                          envir = envir)
                        ),
                        lapply(
                            code[-c(1,2)],
                            function(x)
                                genReplacementsAndCodeRecurse(x,
                                                              constAndIndexNames,
                                                              nimbleFunctionNames,
                                                              debug = debug,
                                                              envir = envir))
                    )
            } else {
                contents <- lapply(
                    code[-1],
                    function(x)
                        genReplacementsAndCodeRecurse(x,
                                                      constAndIndexNames,
                                                      nimbleFunctionNames,
                                                      debug = debug,
                                                      envir = envir))
            }
            contentsCodeReplaced <- lapply(contents, function(x) x$codeReplaced)
            contentsReplacements <- lapply(contents, function(x) x$replacements)
            contentsReplaceable  <-
                unlist(lapply(contents, function(x) x$replaceable))
            allContentsReplaceable <- all(contentsReplaceable)
        } else {
            contentsCodeReplaced <- list()
            contentsReplacements <- list()
            contentsReplaceable  <- list()
            allContentsReplaceable <- TRUE
        }
        ## Do not replace if it is from a special set of functions
        ## or is a nimbleFunction (specifically, an RCfunction)
        if(
        {
            funName <- safeDeparse(code[[1]], warn = TRUE)
            (
                (funName %in% functionsThatShouldNeverBeReplacedInBUGScode) ||
                (exists(funName, envir) && is.rcf(get(funName, envir)))
            )
        }
        )
           return(replaceWhatWeCan(code,
                                    contentsCodeReplaced,
                                    contentsReplacements,
                                    contentsReplaceable,
                                    startingAt=2))
        if(assignment)
            return(replaceWhatWeCan(code,
                                    contentsCodeReplaced,
                                    contentsReplacements,
                                    contentsReplaceable,
                                    startingAt=2))
        isRfunction <- !any(code[[1]] == nimbleFunctionNames)
        isRonly <-
            isRfunction &
            !checkNimbleOrRfunctionNames(safeDeparse(code[[1]], warn = TRUE), envir)
        if(safeDeparse(code[[1]], warn = TRUE) == '$')
            isRonly <- FALSE
        if(isRonly & !allContentsReplaceable)
            stop(paste0('Error, R function \"',
                        safeDeparse(code[[1]]),
                        '\" has non-replaceable node values as arguments.  Must be a nimble function.')
                 )
        if(isRfunction & allContentsReplaceable)
            return(replaceAllCodeSuccessfully(code))
        return(replaceWhatWeCan(code,
                                contentsCodeReplaced,
                                contentsReplacements,
                                contentsReplaceable,
                                startingAt=2))
    }
    stop(paste('Something went wrong in genReplacementsAndCodeRecurse with',
               safeDeparse(code)))
}

replaceAllCodeSuccessfully <- function(code) {
    deparsedCode <- Rname2CppName(code, colonsOK = TRUE)
    replacements <- list()
    replacements[[deparsedCode]] <- code
    return(list(codeReplaced = as.name(deparsedCode),
                replacements = replacements,
                replaceable = TRUE))
}
replaceWhatWeCan <- function(code,
                             contentsCodeReplaced,
                             contentsReplacements,
                             contentsReplaceable,
                             startingAt,
                             replaceable=FALSE) {
    replacements <- list()
    codeReplaced <- code
    if(length(code) >= startingAt)
        for(i in seq_along(contentsReplaceable)) {
            replacements <- c(replacements,
                              contentsReplacements[[i]])
            codeReplaced[[i+startingAt-1]] <- contentsCodeReplaced[[i]]
        }
    replacements <- replacements[unique(names(replacements))]
    list(codeReplaced = codeReplaced,
         replacements = replacements,
         replaceable = replaceable)
}

genLogProbNodeExprAndReplacements <- function(code,
                                              codeReplaced,
                                              indexVarExprs) {
    logProbNodeExpr <- codeReplaced[[2]]   ## Initially, use the replaced version
    replacements <- list()
    
    if(length(logProbNodeExpr) == 1) {
        ## no indexing present
        logProbNodeExpr <- as.name(makeLogProbName(logProbNodeExpr))
        
    } else {
        ## indexing on the LHS node
        if(logProbNodeExpr[[1]] != '[')
            stop('something wrong')
        logProbNodeExpr[[2]] <- as.name(makeLogProbName(logProbNodeExpr[[2]]))
        
        origLHS <- code[[2]]
        for(i in seq_along(origLHS)[-c(1,2)]) {
            origIndex <- origLHS[[i]]
            if(is.vectorized(origIndex)) {
                if(any(indexVarExprs %in% all.vars(origIndex))) {
                    ## Rhe vectorized index includes a loop-indexing
                    ## variable; we will create a replacement, for a
                    ## memberData, for each nodeFunction.
                    replacementExpr <- substitute(min(EXPR),
                                                  list(EXPR=origIndex))
                    replacementName <- Rname2CppName(replacementExpr,
                                                     colonsOK = TRUE)
                    logProbNodeExpr[[i]] <- as.name(replacementName)
                    replacements[[replacementName]] <- replacementExpr
                } else {
                    ## no loop-indexing variables present in the vectorized index.  this index should be constant for all instances of this nodeFunction
                    logProbIndexValue <- as.numeric(min(eval(origIndex)))   ## if this eval() causes an error... then is some case I haven't thought of.
                    logProbNodeExpr[[i]] <- logProbIndexValue
                }
            }
        }
    }
    list(logProbNodeExpr = logProbNodeExpr,
         replacements = replacements)
}

addIndexWrapping <- function(expr) {
    if(length(expr) > 1 &&
       expr[[1]] == '.USED_IN_INDEX') ## nested random indexing
        return(expr)
    substitute(.USED_IN_INDEX(EXPR),
               list(EXPR = expr))
}

addDynamicallyIndexedWrapping <- function(expr) {
    substitute(.DYN_INDEXED(EXPR), list(EXPR = expr))
}

stripDynamicallyIndexedWrapping <- function(expr) {
    if(length(expr) == 1 || !isDynamicIndex(expr))
        return(expr)
    else
        return(expr[[2]])
}

usedInIndex <- function(expr)
    length(expr) > 1 && expr[[1]] == ".USED_IN_INDEX"

isDynamicIndex <- function(expr) {
(length(expr) > 1 && expr[[1]] == ".DYN_INDEXED") ||
    identical(expr, quote(NA_real_))
}

stripIndexWrapping <- function(expr) { 
    if(length(expr) == 1 || !usedInIndex(expr))
        return(expr)
    else
        return(expr[[2]])
}

isVectorIndex <- function(expr) {
    if(isDynamicIndex(expr))
        return(FALSE)
    if(length(expr) > 1 && expr[[1]] == ":")
        return(TRUE)
    return(FALSE)
}

detectNonscalarIndex <- function(expr) {
    if(usedInIndex(expr) || length(expr) == 1)
        return(FALSE)  ## The condition is needed because recursion
                       ## means that we might already have processed
                       ## the dynamic index.
    if(length(expr) == 2) {  ## This can occur if we have mu[k[j[i]]]
        expr <- stripIndexWrapping(expr)
        if(length(expr) <= 2)
            stop("detectNonscalarIndex: unexpected expression ", expr)
    }
    return(
        any(sapply(expr[3:length(expr)],
                   isVectorIndex)
            )
    )
}

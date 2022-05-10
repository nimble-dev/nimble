######	KEYWORD PROCESSING OBJECTS



###		CLASSES

# keywordInfoClass is a class which contains the processor for each keyword
keywordInfoClass <- setRefClass('keywordInfoClass',
                                fields = list(
                                    keyword = 'ANY',
                                    processor = 'ANY'))


# setupCodeTemplateClass is a class that contains the template for generating 
# new setupCode. Objects of this class are used by the function addNecessarySetupCode
setupCodeTemplateClass <- setRefClass('setupCodeTemplateClass',
                                      fields = list(
                                          makeName = 'ANY',
                                          codeTemplate = 'ANY',
                                          makeCodeSubList = 'ANY',
                                          makeOtherNames = 'ANY'),
                                          methods = list(
                                          initialize = function(...){
                                          	makeOtherNames <<- function(name, argList)	return(character(0))
                                          	callSuper(...)
                                          }
                                          ) )


### KEYWORD INFO OBJECTS
		
d_gamma_keywordInfo <- keywordInfoClass(
    keyword = 'dgamma',
    processor = function(code, nfProc, RCfunProc){
        code <- handleScaleAndRateForGamma(code)
	return(code)
    }) 

pq_gamma_keywordInfo <- keywordInfoClass(
    keyword = 'pq_gamma',
    processor = function(code, nfProc, RCfunProc){
        code <- handleScaleAndRateForGamma(code)
	return(code)
    })

rgamma_keywordInfo <- keywordInfoClass(
    keyword = 'rgamma',
    processor = function(code, nfProc, RCfunProc){
        code <- handleScaleAndRateForGamma(code)
        return(code)
    }
)

d_exp_nimble_keywordInfo <- keywordInfoClass(
	keyword = 'dexp_nimble',
	processor = function(code, nfProc, RCfunProc){
		code <- handleScaleAndRateForExpNimble(code)
	return(code)
	}) 

pq_exp_nimble_keywordInfo <- keywordInfoClass(
	keyword = 'pq_exp_nimble',
	processor = function(code, nfProc, RCfunProc){
		code <- handleScaleAndRateForExpNimble(code)
	return(code)
})

rexp_nimble_keywordInfo <- keywordInfoClass(
	keyword = 'rexp_nimble',
	processor = function(code, nfProc, RCfunProc){
		code <- handleScaleAndRateForExpNimble(code)
		return(code)
	}
)

besselK_keywordInfo <- keywordInfoClass(
    keyword = 'besselK',
    processor = function(code, nfProc, RCfunProc) {
        expon.scaledArg <- code$expon.scaled
        if(is.null(expon.scaledArg))
            expon.scaledArg <- FALSE
        if(is.numeric(expon.scaledArg) || is.logical(expon.scaledArg)) {
            code$expon.scaled <- 1 + as.logical(expon.scaledArg)
        } else code$expon.scaled <- substitute(1 + A, list(A = expon.scaledArg))
        return(code)
    }
)


nimSeq_keywordInfo <- keywordInfoClass(
    keyword = 'nimSeq',
    processor = function(code, nfProc, RCfunProc) {
        useBy <- !isCodeArgBlank(code, 'by')
        useLen <- !isCodeArgBlank(code, 'length.out')
        if(useBy && useLen)
            newRunCode <- substitute(nimSeqByLen(FROM, 0, BY, LEN), list(FROM = code$from, BY = code$by, LEN = code$length.out))
        else {
            if(useLen) {
                newRunCode <- substitute(nimSeqLen(FROM, TO, 0, LEN), list(FROM = code$from, TO = code$to, LEN = code$length.out))
            } else {
                byVal <- if(useBy) code$by else 1 ## default by = 1
                newRunCode <- substitute(nimSeqBy(FROM, TO, BY, 0), list(FROM = code$from, TO = code$to, BY = byVal))
            }
        }
        return(newRunCode)
    }
)
    

values_keywordInfo <- keywordInfoClass(
    keyword = 'values',
    processor = function(code, nfProc, RCfunProc){
      if(!isCodeArgBlank(code, 'accessor'))
      	return(code)
      if(isCodeArgBlank(code, 'model'))
      	stop('model argument missing from values call, with no accessor argument supplied')
      
      accessArgList <- list(model = code$model, nodes = code$nodes, logProb = FALSE, logProbOnly = FALSE)

      useAccessorVectorByIndex <- FALSE
      if(hasBracket(accessArgList$nodes)) { 
            useAccessorVectorByIndex <- TRUE
            if(length(accessArgList$nodes) != 3) stop(paste0('Problem with ', deparse(code),'. If you need to index on the nodes argument there should be only one index.'))
            nodesIndexExpr <- accessArgList$nodes[[3]]
            accessArgList$nodes <- accessArgList$nodes[[2]]
            accessArgList$sortUnique <- FALSE   
      }

      accessName <- modelVariableAccessorVector_setupCodeTemplate$makeName(accessArgList)
      addNecessarySetupCode(accessName, accessArgList, modelVariableAccessorVector_setupCodeTemplate, nfProc)
      if(!useAccessorVectorByIndex)
          newRunCode <- substitute(values(accessor = ACCESS_NAME), 
                               list(ACCESS_NAME = as.name(accessName)))
        else
            newRunCode <- substitute(values(accessor = ACCESS_NAME, accessorIndex = ACCESSVECINDEX),
                                 list(ACCESS_NAME = as.name(accessName), ACCESSVECINDEX = nodesIndexExpr))
      return(newRunCode)
    })                                    

getParam_keywordInfo <- keywordInfoClass(
    keyword = 'getParam',
    processor = function(code, nfProc, RCfunProc) {
        if(!isCodeArgBlank(code, 'nodeFunction'))
            return(code)
        errorContext <- deparse(code)
        nodeFunVec_ArgList <- list(model = code$model, nodes = code$node, includeData = TRUE, sortUnique = TRUE, errorContext = errorContext)
        if(!isCodeArgBlank(code, 'nodeFunctionIndex')) { ## new case: calculate(myNodeFunctionVector, nodeFunctionIndex = i), if myNodeFunctionVector was hand-created in setup code
            if(!isCodeArgBlank(code, 'nodes'))
                stop('nodes argument cannot be provided to getParam if nodeFunctionIndex is specified')
            return(code) ## no modification needed!
        }

        if(isCodeArgBlank(code, 'model'))
            stop('model argument missing from getParam, with no accessor argument supplied')
        if(isCodeArgBlank(code, 'node'))
            stop('node argument missing from getParam, with no accessor argument supplied')

        useNodeFunctionVectorByIndex <- FALSE
        if(hasBracket(nodeFunVec_ArgList$nodes)) { ## like calculate(model, nodes[i]), which could have started as model$calculate(nodes[i])
            useNodeFunctionVectorByIndex <- TRUE
            if(length(nodeFunVec_ArgList$nodes) != 3) stop(paste0('Problem with ', deparse(code),'. If you need to index on the nodes argument there should be only one index.'))
            nodesIndexExpr <- nodeFunVec_ArgList$nodes[[3]]
            nodeFunVec_ArgList$nodes <- nodeFunVec_ArgList$nodes[[2]]
            nodeFunVec_ArgList$sortUnique <- FALSE
        }

        nodeFunName <- nodeFunctionVector_SetupTemplate$makeName(nodeFunVec_ArgList)

        if(isCodeArgBlank(code, 'param'))
            stop("'param' argument missing from 'getParam', with no accessor argument supplied")
        paramInfo_ArgList <- list(model = code$model, node = nodeFunVec_ArgList$nodes, param = code$param, hasIndex = useNodeFunctionVectorByIndex) ## use nodeFunVec_ArgList$nodes instead of code$node because nodeFunVec_ArgList$nodes may have been updated if code$nodes has a run-time index.  In that case the paramID will be vector
        paramInfoName <- paramInfo_SetupTemplate$makeName(paramInfo_ArgList)
        paramIDname <- paramInfo_SetupTemplate$makeOtherNames(paramInfoName, paramInfo_ArgList)

        addNecessarySetupCode(nodeFunName, nodeFunVec_ArgList, nodeFunctionVector_SetupTemplate, nfProc)
        addNecessarySetupCode(paramInfoName, paramInfo_ArgList, paramInfo_SetupTemplate, nfProc)
        if(!useNodeFunctionVectorByIndex)
            newRunCode <- substitute(getParam(nodeFunction = NODEFUNVEC_NAME, paramID = PARAMID_NAME, paramInfo = PARAMINFO_NAME),
                                     list(NODEFUNVEC_NAME = as.name(nodeFunName), PARAMID_NAME = as.name(paramIDname), PARAMINFO_NAME = as.name(paramInfoName)))
        else
            newRunCode <- substitute(getParam(nodeFunction = NODEFUNVEC_NAME, paramID = PARAMID_NAME, paramInfo = PARAMINFO_NAME, nodeFunctionIndex = NODEFUNVECINDEX),
                                     list(NODEFUNVEC_NAME = as.name(nodeFunName), PARAMID_NAME = as.name(paramIDname), PARAMINFO_NAME = as.name(paramInfoName), NODEFUNVECINDEX = nodesIndexExpr))
        
        return(newRunCode)
    }
)

getBound_keywordInfo <- keywordInfoClass(
    keyword = 'getBound',
    processor = function(code, nfProc, RCfunProc) {
        if(!isCodeArgBlank(code, 'nodeFunction'))
            return(code)
        errorContext <- deparse(code)
        nodeFunVec_ArgList <- list(model = code$model, nodes = code$node, includeData = TRUE, sortUnique = TRUE, errorContext = errorContext)
        if(!isCodeArgBlank(code, 'nodeFunctionIndex')) { ## new case: calculate(myNodeFunctionVector, nodeFunctionIndex = i), if myNodeFunctionVector was hand-created in setup code
            if(!isCodeArgBlank(code, 'nodes'))
                stop('nodes argument cannot be provided to getParam if nodeFunctionIndex is specified')
            return(code) ## no modification needed!
        }

        if(isCodeArgBlank(code, 'model'))
            stop('model argument missing from getParam, with no accessor argument supplied')
        if(isCodeArgBlank(code, 'node'))
            stop('node argument missing from getParam, with no accessor argument supplied')

        useNodeFunctionVectorByIndex <- FALSE
        if(hasBracket(nodeFunVec_ArgList$nodes)) { ## like calculate(model, nodes[i]), which could have started as model$calculate(nodes[i])
            useNodeFunctionVectorByIndex <- TRUE
            if(length(nodeFunVec_ArgList$nodes) != 3) stop(paste0('Problem with ', deparse(code),'. If you need to index on the nodes argument there should be only one index.'))
            nodesIndexExpr <- nodeFunVec_ArgList$nodes[[3]]
            nodeFunVec_ArgList$nodes <- nodeFunVec_ArgList$nodes[[2]]
            nodeFunVec_ArgList$sortUnique <- FALSE
        }

        nodeFunName <- nodeFunctionVector_SetupTemplate$makeName(nodeFunVec_ArgList)

        if(isCodeArgBlank(code, 'bound'))
            stop("'bound' argument missing from 'getBound', with no accessor argument supplied")
        boundInfo_ArgList <- list(model = code$model, node = nodeFunVec_ArgList$nodes, bound = code$bound) ## use nodeFunVec_ArgList$nodes instead of code$node because nodeFunVec_ArgList$nodes may have been updated if code$nodes has a run-time index.  In that case the boundID will be vector
        boundInfoName <- boundInfo_SetupTemplate$makeName(boundInfo_ArgList)
        boundIDname <- boundInfo_SetupTemplate$makeOtherNames(boundInfoName, boundInfo_ArgList)

        addNecessarySetupCode(nodeFunName, nodeFunVec_ArgList, nodeFunctionVector_SetupTemplate, nfProc)
        addNecessarySetupCode(boundInfoName, boundInfo_ArgList, boundInfo_SetupTemplate, nfProc)
        if(!useNodeFunctionVectorByIndex)
            newRunCode <- substitute(getBound(nodeFunction = NODEFUNVEC_NAME, boundID = BOUNDID_NAME, boundInfo = BOUNDINFO_NAME),
                                     list(NODEFUNVEC_NAME = as.name(nodeFunName), BOUNDID_NAME = as.name(boundIDname), BOUNDINFO_NAME = as.name(boundInfoName)))
        else
            newRunCode <- substitute(getBound(nodeFunction = NODEFUNVEC_NAME, boundID = BOUNDID_NAME, boundInfo = BOUNDINFO_NAME, nodeFunctionIndex = NODEFUNVECINDEX),
                                     list(NODEFUNVEC_NAME = as.name(nodeFunName), BOUNDID_NAME = as.name(boundIDname), BOUNDINFO_NAME = as.name(boundInfoName), NODEFUNVECINDEX = nodesIndexExpr))
                                               
        return(newRunCode)
    }
)


calculate_keywordInfo <- keywordInfoClass(
    keyword = 'calculate',
    processor = function(code, nfProc, RCfunProc){
        if(!isCodeArgBlank(code, 'nodeFxnVector'))
            return(code)
        errorContext <- deparse(code)

        enableDerivs <- FALSE
        withDerivsOutputOnly <- FALSE
        if(isTRUE(nimbleOptions("enableDerivs")) && isTRUE(nimbleOptions("buildDerivs"))) {
          derivControl <- environment(nfProc$nfGenerator)$enableDerivs[[RCfunProc$name]]
          ## There are two cases.
          ## If enableDerivs has an entry, then derivatives are enabled either
          ## explicitly or because of a nimDerivs(model$calculate), i.e. a direct line into derivs for just a model calculation.
          ##
          ## These need to be disambiguated for use below to choose the right setup template case
          enableDerivs <- !is.null(derivControl) | !is.null(code$wrt)
          withDerivsOutputOnly <- enableDerivs & is.null(code$wrt) ## This is the case of *not* nimDerivs(model$calculate)
        }
        
        nodeFunVec_ArgList <- list(model = code$model, nodes = code$nodes, wrtNodes = code$wrt,
                                   includeData = TRUE, sortUnique = TRUE, errorContext = errorContext)

        ##
        if(!isCodeArgBlank(code, 'nodeFunctionIndex')) { ## new case: calculate(myNodeFunctionVector, nodeFunctionIndex = i), if myNodeFunctionVector was hand-created in setup code
            if(!isCodeArgBlank(code, 'nodes'))
                stop('nodes argument cannot be provided to calculate if nodeFunctionIndex is specified')
            return(code) ## no modification needed!
        }
        if(isCodeArgBlank(code, 'model'))
            stop('model argument missing from calculate, with no accessor argument supplied')
        if(isCodeArgBlank(code, 'nodes')){
            LHSnodes_ArgList <- list(model = code$model)
            LHSnodes_name <- allLHSNodes_SetupTemplate$makeName(LHSnodes_ArgList)
            addNecessarySetupCode(LHSnodes_name, LHSnodes_ArgList, allLHSNodes_SetupTemplate, nfProc, allowToCpp = FALSE)
            nodeFunVec_ArgList$nodes = as.name(LHSnodes_name)
        }
        useNodeFunctionVectorByIndex <- FALSE
        if(hasBracket(nodeFunVec_ArgList$nodes)) { ## like calculate(model, nodes[i]), which could have started as model$calculate(nodes[i])
            if(enableDerivs)
                stop(paste0('Derivatives for an indexed node are not allowed.\n',
                            'Use nimDerivs(model$calculate(nodes), ...) instead of \n',
                            'nimDerivs(model$calculate(nodes[i]).\n'))
            useNodeFunctionVectorByIndex <- TRUE
            if(length(nodeFunVec_ArgList$nodes) != 3)
                stop(paste0('Problem with ',
                            deparse(code),
                            '. If you need to index on the nodes argument there should be only one index.'))
            nodesIndexExpr <- nodeFunVec_ArgList$nodes[[3]]
            nodeFunVec_ArgList$nodes <- nodeFunVec_ArgList$nodes[[2]]
            nodeFunVec_ArgList$sortUnique <- FALSE
        }

        if(!withDerivsOutputOnly) { ## This is regular mode, including without derivs at all and with enableDerivs but not nimDerivs(model$calculate...)
          nodeFunName <- nodeFunctionVector_SetupTemplate$makeName(nodeFunVec_ArgList)	
          addNecessarySetupCode(nodeFunName, nodeFunVec_ArgList, nodeFunctionVector_SetupTemplate, nfProc)
        } else {
          nodeFunName <- nodeFunctionVector_WithDerivsOutput_SetupTemplate$makeName(nodeFunVec_ArgList)	
          addNecessarySetupCode(nodeFunName, nodeFunVec_ArgList, nodeFunctionVector_WithDerivsOutput_SetupTemplate, nfProc)
        }
        if(!useNodeFunctionVectorByIndex){
            newRunCode <- substitute(calculate(nodeFxnVector = NODEFUNVEC_NAME),
                                     list(NODEFUNVEC_NAME = as.name(nodeFunName)))
        }
        else{
            newRunCode <- substitute(
                calculate(nodeFxnVector = NODEFUNVEC_NAME,
                          nodeFunctionIndex = NODEFUNVECINDEX),
                list(NODEFUNVEC_NAME = as.name(nodeFunName),
                     NODEFUNVECINDEX = nodesIndexExpr))
        }
        return(newRunCode)
    }
)

calculateDiff_keywordInfo <- keywordInfoClass(
    keyword = 'calculateDiff',
    processor = function(code, nfProc, RCfunProc){
        if(!isCodeArgBlank(code, 'nodeFxnVector'))
            return(code)
        errorContext <- deparse(code)
        nodeFunVec_ArgList <- list(model = code$model, nodes = code$nodes, includeData = TRUE, sortUnique = TRUE, errorContext = errorContext)
        if(!isCodeArgBlank(code, 'nodeFunctionIndex')) { ## new case: calculate(myNodeFunctionVector, nodeFunctionIndex = i), if myNodeFunctionVector was hand-created in setup code
            if(!isCodeArgBlank(code, 'nodes'))
                stop('nodes argument cannot be provided to calculateDiff if nodeFunctionIndex is specified')
            return(code) ## no modification needed!
        }
        
        if(isCodeArgBlank(code, 'model'))
            stop('model argument missing from calculateDiff, with no accessor argument supplied')
        if(isCodeArgBlank(code, 'nodes')){
            LHSnodes_ArgList <- list(model = code$model)
            LHSnodes_name <- allLHSNodes_SetupTemplate$makeName(LHSnodes_ArgList)
            addNecessarySetupCode(LHSnodes_name, LHSnodes_ArgList, allLHSNodes_SetupTemplate, nfProc, allowToCpp = FALSE)
            nodeFunVec_ArgList$nodes = as.name(LHSnodes_name)
        }
        useNodeFunctionVectorByIndex <- FALSE
        if(hasBracket(nodeFunVec_ArgList$nodes)) { ## like calculate(model, nodes[i]), which could have started as model$calculate(nodes[i])
            useNodeFunctionVectorByIndex <- TRUE
            if(length(nodeFunVec_ArgList$nodes) != 3) stop(paste0('Problem with ', deparse(code),'. If you need to index on the nodes argument there should be only one index.'))
            nodesIndexExpr <- nodeFunVec_ArgList$nodes[[3]]
            nodeFunVec_ArgList$nodes <- nodeFunVec_ArgList$nodes[[2]]
            nodeFunVec_ArgList$sortUnique <- FALSE
        }
        
        nodeFunName <- nodeFunctionVector_SetupTemplate$makeName(nodeFunVec_ArgList)	
        addNecessarySetupCode(nodeFunName, nodeFunVec_ArgList, nodeFunctionVector_SetupTemplate, nfProc)
        if(!useNodeFunctionVectorByIndex)
            newRunCode <- substitute(calculateDiff(nodeFxnVector = NODEFUNVEC_NAME),
                                     list(NODEFUNVEC_NAME = as.name(nodeFunName)))
        else
            newRunCode <- substitute(calculateDiff(nodeFxnVector = NODEFUNVEC_NAME, nodeFunctionIndex = NODEFUNVECINDEX),
                                 list(NODEFUNVEC_NAME = as.name(nodeFunName), NODEFUNVECINDEX = nodesIndexExpr))
        return(newRunCode)	
    }
)


simulate_keywordInfo <- keywordInfoClass(
    keyword = 'simulate',
    processor = function(code, nfProc, RCfunProc){
        if(!isCodeArgBlank(code, 'nodeFxnVector')){
            return(substitute(simulate(nodeFxnVector = NODEFXNVECTOR), list(NODEFXNVECTOR = code$nodeFxnVector) ) )
        }
        if(!isCodeArgBlank(code, 'INDEXEDNODEINFO_'))
            return(code)
        errorContext <- deparse(code)
        nodeFunVec_ArgList <- list(model = code$model, nodes = code$nodes, includeData = code$includeData, sortUnique = TRUE, errorContext = errorContext)
        if(!isCodeArgBlank(code, 'nodeFunctionIndex')) { ## new case: calculate(myNodeFunctionVector, nodeFunctionIndex = i), if myNodeFunctionVector was hand-created in setup code
            if(!isCodeArgBlank(code, 'nodes'))
                stop('nodes argument cannot be provided to simulate if nodeFunctionIndex is specified')
            return(code) ## no modification needed!
        }
        if(isCodeArgBlank(code, 'model'))
            stop('model argument missing from simulate, with no accessor argument supplied')
        if(isCodeArgBlank(code, 'nodes')){
            LHSnodes_ArgList <- list(model = code$model)
            LHSnodes_name <- allLHSNodes_SetupTemplate$makeName(LHSnodes_ArgList)
            addNecessarySetupCode(LHSnodes_name, LHSnodes_ArgList, allLHSNodes_SetupTemplate, nfProc, allowToCpp = FALSE)
            nodeFunVec_ArgList$nodes = as.name(LHSnodes_name)
        }
        useNodeFunctionVectorByIndex <- FALSE
        if(hasBracket(nodeFunVec_ArgList$nodes)) { ## like calculate(model, nodes[i]), which could have started as model$calculate(nodes[i])
            useNodeFunctionVectorByIndex <- TRUE
            if(length(nodeFunVec_ArgList$nodes) != 3) stop(paste0('Problem with ', deparse(code),'. If you need to index on the nodes argument there should be only one index.'))
            nodesIndexExpr <- nodeFunVec_ArgList$nodes[[3]]
            nodeFunVec_ArgList$nodes <- nodeFunVec_ArgList$nodes[[2]]
            nodeFunVec_ArgList$sortUnique <- FALSE # If includeData = FALSE, this can trigger error from nodeFunctionVector if nodes does contain data
        }
  
        nodeFunName <- nodeFunctionVector_SetupTemplate$makeName(nodeFunVec_ArgList)	
        addNecessarySetupCode(nodeFunName, nodeFunVec_ArgList, nodeFunctionVector_SetupTemplate, nfProc)
        if(!useNodeFunctionVectorByIndex)
            newRunCode <- substitute(simulate(nodeFxnVector = NODEFUNVEC_NAME),
                                     list(NODEFUNVEC_NAME = as.name(nodeFunName)))
        else
            newRunCode <- substitute(simulate(nodeFxnVector = NODEFUNVEC_NAME, nodeFunctionIndex = NODEFUNVECINDEX),
                                 list(NODEFUNVEC_NAME = as.name(nodeFunName), NODEFUNVECINDEX = nodesIndexExpr))
        
        return(newRunCode)	
    }
)

getLogProb_keywordInfo <- keywordInfoClass(
    keyword = 'getLogProb',
    processor = function(code, nfProc, RCfunProc){
        if(!isCodeArgBlank(code, 'nodeFxnVector'))
            return(code)
        errorContext <- deparse(code)
        nodeFunVec_ArgList <- list(model = code$model, nodes = code$nodes, includeData = TRUE, sortUnique = TRUE, errorContext = errorContext)
        if(!isCodeArgBlank(code, 'nodeFunctionIndex')) { ## new case: calculate(myNodeFunctionVector, nodeFunctionIndex = i), if myNodeFunctionVector was hand-created in setup code
            if(!isCodeArgBlank(code, 'nodes'))
                stop('nodes argument cannot be provided to getLogProb if nodeFunctionIndex is specified')
            return(code) ## no modification needed!
        }
        if(isCodeArgBlank(code, 'model'))
            stop('model argument missing from getLogProb, with no accessor argument supplied')
        if(isCodeArgBlank(code, 'nodes')){
            LHSnodes_ArgList <- list(model = code$model)
            LHSnodes_name <- allLHSNodes_SetupTemplate$makeName(LHSnodes_ArgList)
            addNecessarySetupCode(LHSnodes_name, LHSnodes_ArgList, allLHSNodes_SetupTemplate, nfProc, allowToCpp = FALSE)
            nodeFunVec_ArgList$nodes = as.name(LHSnodes_name)
        }
        useNodeFunctionVectorByIndex <- FALSE
        if(hasBracket(nodeFunVec_ArgList$nodes)) { ## like calculate(model, nodes[i]), which could have started as model$calculate(nodes[i])
            useNodeFunctionVectorByIndex <- TRUE
            if(length(nodeFunVec_ArgList$nodes) != 3) stop(paste0('Problem with ', deparse(code),'. If you need to index on the nodes argument there should be only one index.'))
            nodesIndexExpr <- nodeFunVec_ArgList$nodes[[3]]
            nodeFunVec_ArgList$nodes <- nodeFunVec_ArgList$nodes[[2]]
            nodeFunVec_ArgList$sortUnique <- FALSE
        }

        nodeFunName <- nodeFunctionVector_SetupTemplate$makeName(nodeFunVec_ArgList)	
        addNecessarySetupCode(nodeFunName, nodeFunVec_ArgList, nodeFunctionVector_SetupTemplate, nfProc)
        if(!useNodeFunctionVectorByIndex)
            newRunCode <- substitute(getLogProb(nodeFxnVector = NODEFUNVEC_NAME),
                                     list(NODEFUNVEC_NAME = as.name(nodeFunName)))
        else
            newRunCode <- substitute(getLogProb(nodeFxnVector = NODEFUNVEC_NAME, nodeFunctionIndex = NODEFUNVECINDEX),
                                     list(NODEFUNVEC_NAME = as.name(nodeFunName), NODEFUNVECINDEX = nodesIndexExpr))
        return(newRunCode)	
    }
)    

nimCopy_keywordInfo <- keywordInfoClass(
	keyword = 'nimCopy',
    processor = function(code, nfProc, RCfunProc){
        if(is.null(nfProc)) stop("Can\'t call copy (nimCopy) from a nimbleFunction without setup code")
		possibleObjects <- c('symbolModel', 'symbolModelValues', 'symbolModelVariableAccessorVector', 'symbolModelValuesAccessorVector')
		modelValuesTypes <- c('symbolModelValues', 'symbolModelValuesAccessorVector')
		accessTypes <- c('symbolModelVariableAccessorVector', 'symbolModelValuesAccessorVector')
		from_ArgList <- list(name = code$from, class = symTypeFromSymTab(code$from, nfProc$setupSymTab, options = possibleObjects))
		to_ArgList <- list(name = code$to, class = symTypeFromSymTab(code$to, nfProc$setupSymTab, options = possibleObjects))
		if(is.null(from_ArgList$class)) 
			stop("Error in nimCopy: '", code$from, "' is not a recognized model or modelValues object.") 		
		if(is.null(to_ArgList$class)) 
			stop("Error in nimCopy: '", code$to, "' is not a recognized model or modelValues object.") 		
		if(from_ArgList$class %in% modelValuesTypes){
			if(isCodeArgBlank(code, 'row'))		stop('row argument missing in copy call')
			from_ArgList$row = code$row
		}
		if(to_ArgList$class %in% modelValuesTypes){
			if(isCodeArgBlank(code, 'rowTo')){
				if(isCodeArgBlank(code, 'row'))		stop('row argument missing in copy call')
				else								to_ArgList$row = code$row
			}
			else		to_ArgList$row = code$rowTo
		}
		if(isCodeArgBlank(code, 'nodes')){
			if(from_ArgList$class == 'symbolModel'){
				node_ArgList <- list(model = from_ArgList$name)
				allNodes_name <- allModelNodes_SetupTemplate$makeName( node_ArgList )
				addNecessarySetupCode(allNodes_name, node_ArgList, allModelNodes_SetupTemplate, nfProc, allowToCpp = FALSE)
			}
			else if(from_ArgList$class == 'symbolModelValues'){
				from_ArgList$row = code$row
				mvVar_ArgList <- list(modelValues = from_ArgList$name)
				allNodes_name <- allModelValuesVars_SetupTemplate$makeName(mvVar_ArgList)
				addNecessarySetupCode(allNodes_name, mvVar_ArgList, allModelValuesVars_SetupTemplate, nfProc, allowToCpp = FALSE)
			}
			from_ArgList$nodes <- as.name(allNodes_name)
		}
		else	from_ArgList$nodes <- code$nodes
		
		if(isCodeArgBlank(code, 'nodesTo'))		to_ArgList$nodes <- from_ArgList$nodes
		else									to_ArgList$nodes <- code$nodesTo
				
		if(from_ArgList$class == 'symbolModel'){
                    isMVfrom <- 0 
                    accessFrom_ArgList <- list(model = code$from, nodes = from_ArgList$nodes, logProb = code$logProb, logProbOnly = code$logProbOnly)
                    accessFrom_name <- modelVariableAccessorVector_setupCodeTemplate$makeName(accessFrom_ArgList)
                    addNecessarySetupCode(accessFrom_name, accessFrom_ArgList, modelVariableAccessorVector_setupCodeTemplate, nfProc)
		}
		else if(from_ArgList$class == 'symbolModelValues'){
                    isMVfrom <- 1
                    accessFrom_ArgList <- list(modelValues = code$from, nodes = from_ArgList$nodes, logProb = code$logProb, logProbOnly = code$logProbOnly, row = from_ArgList$row)
                    accessFrom_name <- modelValuesAccessorVector_setupCodeTemplate$makeName(accessFrom_ArgList)
                    addNecessarySetupCode(accessFrom_name, accessFrom_ArgList, modelValuesAccessorVector_setupCodeTemplate, nfProc)
		}
		else if(from_ArgList$class %in% accessTypes) {
                    isMVfrom <- as.integer(from_ArgList$class == 'symbolModelValuesAccessorVector') 
                    accessFrom_name <- as.character(code$from)
                }
        
		if(to_ArgList$class == 'symbolModel'){
                    isMVto <- 0
                    accessTo_ArgList <- list(model = code$to, nodes = to_ArgList$nodes, logProb = code$logProb, logProbOnly = code$logProbOnly)
                    accessTo_name <- modelVariableAccessorVector_setupCodeTemplate$makeName(accessTo_ArgList)
                    addNecessarySetupCode(accessTo_name, accessTo_ArgList, modelVariableAccessorVector_setupCodeTemplate, nfProc)
		}
		else if(to_ArgList$class == 'symbolModelValues'){
                    isMVto <- 1
                    accessTo_ArgList <- list(modelValues = code$to, nodes = to_ArgList$nodes, logProb = code$logProb, logProbOnly = code$logProbOnly, row = to_ArgList$row)
                    accessTo_name <- modelValuesAccessorVector_setupCodeTemplate$makeName(accessTo_ArgList)
                    addNecessarySetupCode(accessTo_name, accessTo_ArgList, modelValuesAccessorVector_setupCodeTemplate, nfProc)
		}
		else if(to_ArgList$class %in% accessTypes) {
                    isMVto <- as.integer(to_ArgList$class == 'symbolModelValuesAccessorVector') 
                    accessTo_name <- as.character(code$to) 
                }
        if(nimbleOptions()$useNewNimCopy) {
            copierVector_ArgList <- list(accessFrom_name = accessFrom_name, accessTo_name = accessTo_name, isMVto = isMVto, isMVfrom = isMVfrom)
            copierVector_name <- copierVector_setupCodeTemplate$makeName(copierVector_ArgList)
            addNecessarySetupCode(copierVector_name, copierVector_ArgList, copierVector_setupCodeTemplate, nfProc) 
        }
        
        if(!nimbleOptions()$useNewNimCopy) {
            ##What happens below is a bit convoluted and really for backwards compatibility 	
            runCode <- substitute(nimCopy(from = FROM_ACCESS, rowFrom = NA, to = TO_ACCESS, rowTo = NA), 
                                  list(FROM_ACCESS = as.name(accessFrom_name), TO_ACCESS = as.name(accessTo_name)))
            if(from_ArgList$class %in% modelValuesTypes)
                runCode$rowFrom = from_ArgList$row
            if(to_ArgList$class %in% modelValuesTypes)
                runCode$rowTo = to_ArgList$row
        } else {
            rowFromArg <- if(from_ArgList$class %in% modelValuesTypes) from_ArgList$row else NA
            rowToArg <- if(to_ArgList$class %in% modelValuesTypes) {
                if(identical(rowFromArg, NA)) {rowFromArg <- 0; unusedArg <- NA} else unusedArg <- 0
                to_ArgList$row
            } else {
                unusedArg <- NA
                NA
            }
            runCode <- substitute(nimCopy(copierVector = COPIER_VECTOR, rowFrom = ROWFROM, rowTo = ROWTO, unused = UNUSED), 
                                  list(COPIER_VECTOR = as.name(copierVector_name),
                                       ROWFROM = rowFromArg, ROWTO = rowToArg, UNUSED  = unusedArg))
        }
        runCode <- runCode[as.character(runCode) != 'NA']
        return(runCode)
    })


doubleBracket_keywordInfo <- keywordInfoClass(
	keyword = '[[', 
    processor = function(code, nfProc, RCfunProc){
        if(is.null(nfProc)) stop("No allowed use of [[ in a nimbleFunction without setup code.")
        possibleObjects <- c('symbolModel', 'symbolNimPtrList', 'symbolNimbleFunctionList', 'symbolNimbleList')
        class = symTypeFromSymTab(code[[2]], nfProc$setupSymTab, options = possibleObjects)
        if(is.null(class)){  ##assume that an element of a run-time provided nimbleList is being accessed
          nl_charName <- as.character(callerCode)
          nl_fieldName <-as.character(code[[3]])
          newRunCode <- substitute(nfVar(NIMBLELIST, VARNAME), list(NIMBLELIST = as.name(nl_charName), VARNAME = nl_fieldName))
          return(newRunCode)
        }
        if(class == 'symbolNimPtrList' || class == 'symbolNimbleFunctionList')
            return(code)
        if(class == 'symbolNimbleList'){
          #	Code is of the form 
          #  myNimbleList[['myVar']]
          nl_charName <- as.character(callerCode)
          nl_fieldName <-as.character(code[[3]])
          newRunCode <- substitute(nlVar(NIMBLELIST, VARNAME), list(NIMBLELIST = as.name(nl_charName), VARNAME = nl_fieldName))
          return(newRunCode)
        }
        if(class == 'symbolModel'){
            singleAccess_ArgList <- list(code = code, model = code[[2]], nodeExpr = code[[3]])
            nodeArg <- code[[3]]
            if(is.character(nodeArg)){
                if(length(nodeArg) > 1)
                    stop(paste0("Problem in ",
                                deparse(code), ". ",
                                deparse(code[[3]]),
                                " is too long.  It can only have one element."),
                         call. = FALSE)
                varAndIndices <- nimbleInternalFunctions$getVarAndIndices(nodeArg)
                
                allNDims <- lapply(nfProc$instances, function(x) {
                    model <- eval(singleAccess_ArgList$model,
                                  envir = x)
                    if(length(nodeArg) != 1)
                        stop(paste0("Length of ",
                                    nodeArg,
                                    " requested from ",
                                    deparse(singleAccess_ArgList$model),
                                    " using '[[' is ",
                                    length(nodeArg),
                                    ". It must be 1.")
                           , call. = FALSE)
                    determineNdimFromOneCase(model, varAndIndices)
                } )

                if(length(unique(allNDims)) > 1)
                    stop(paste0('Error for ', deparse(code),
                                '. Inconsistent numbers of dimensions for different instances.'),
                         call. = FALSE)
                nDim <- allNDims[[1]]
                useMap <- nDim > 0
                
                ## ##
                ## ## If input is of the form model[['a']]
                ## ## and a is non-scalar,
                ## ## we treat it like model$a, which handles either
                if(useMap & length(varAndIndices$indices) == 0) {
                    return(
                        keywordList[['$']]$processor(code, nfProc, RCfunProc)
                    )
                }
                ## ## Following line adds up 0 for each scalar index
                ## ## and 1 for each non-scalar index to determine if the
                ## ## accessed node is scalar
                ## nDim <- sum(1 - unlist(lapply(varAndIndices$indices, is.numeric) ) )
                ## useMap <- nDim > 0
            }
            else{
                allNDims <- determineNdimsFromNfproc(singleAccess_ArgList$model, nodeArg, nfProc)
                if(length(unique(allNDims)) > 1) stop(paste0('Error for ', deparse(code), '. Inconsistent numbers of dimensions for different instances.'), call. = FALSE)
                nDim <- allNDims[[1]]
                useMap <- nDim > 0
            }
            if(useMap){
                accessName <- map_SetupTemplate$makeName(singleAccess_ArgList)
                addNecessarySetupCode(accessName, singleAccess_ArgList, map_SetupTemplate, nfProc)	
                ans <- makeMapAccessExpr(accessName, as.name(accessName), nDim)
            }
            else{
                accessName <- singleModelIndexAccess_SetupTemplate$makeName(singleAccess_ArgList)
                addNecessarySetupCode(accessName, singleAccess_ArgList, singleModelIndexAccess_SetupTemplate, nfProc)
                                        #ans <- substitute(ACCESSNAME[MFLATINDEX], list(ACCESSNAME = as.name(accessName), MFLATINDEX = as.name(paste0(accessName, '_flatIndex'))))
                ans <- makeSingleIndexAccessExpr(accessName, as.name(accessName))
            }
            return(ans)			
        }
        stop("Incorrect use of double brackets in: '", deparse(code), "'.")
    })

modelMemberFun_keywordInfo <- keywordInfoClass(
    keyword = 'multiple',
    processor = function(code, nfProc, RCfunProc) {
        ## if we get here it must be model$member(args)
        ## We will turn it into member(model, args)
        newRunCode <- do.call("call",
                              c(list(as.character(code[[1]][[3]]),
                                     code[[1]][[2]]),
                                as.list(code[-1])),
                              quote = TRUE)
        return(newRunCode)
    })

dollarSign_keywordInfo <- keywordInfoClass(
    keyword = '$',
    processor = function(code, nfProc, RCfunProc){
        callerCode <- code[[2]]
        
        if(is.null(nfProc)) { 

            nl_fieldName <-as.character(code[[3]])
            newRunCode <- substitute(nfVar(NIMBLELIST, VARNAME), list(NIMBLELIST = callerCode, VARNAME = nl_fieldName))
            return(newRunCode)
        }
        
        doubleBracketCase <- FALSE
        if(length(callerCode) > 1) {
            if(deparse(callerCode[[1]] == '[[')) {
                doubleBracketCase <- TRUE
                symObj <- getSymObj_recurse(callerCode[[2]], nfProc$setupSymTab)
            }
        }
        if(!doubleBracketCase)
            symObj <- getSymObj_recurse(callerCode, nfProc$setupSymTab)

        class <- class(symObj)[1] ## symObj is allowed to be NULL       

                                        #	This extracts myNimbleFunction from the expression myNimbleFunction$foo()
        if(length(callerCode) > 1){
            if(callerCode[[1]] == '$'){ ## nested NL or NF case
                callerCode <- processKeyword(callerCode, nfProc, RCfunProc)
            }
        }
                                        #       This extracts myNimbleFunctionList from the expression myNimbleFunctionList[[i]]
                                        #       May be a better way to do this
        

        if(is.null(class) || class == 'NULL'){  ##assume that an element of a run-time provided nimbleList is being accessed
            nl_fieldName <-as.character(code[[3]])
            newRunCode <- substitute(nfVar(NIMBLELIST, VARNAME), list(NIMBLELIST = callerCode, VARNAME = nl_fieldName))
            return(newRunCode)				
        }
        if(class == 'symbolNimbleFunctionList'){
             nf_fieldName <-as.character(code[[3]])
            newRunCode <- substitute(nfMethod(NIMBLEFXN, METHODNAME), list(NIMBLEFXN = callerCode, METHODNAME = nf_fieldName))
            return(newRunCode)
        }
        if(class == 'symbolModel'){
            singleAccess_ArgList <- list(code = code, model = callerCode, var = as.character(code[[3]]) )
            accessName <- singleVarAccess_SetupTemplate$makeName(singleAccess_ArgList)
            addNecessarySetupCode(accessName, singleAccess_ArgList, singleVarAccess_SetupTemplate, nfProc)
            return(as.name(accessName))
        }
        if(class == 'symbolNimbleFunction'){
            
                                        #	Code is of the form myNimbleFunction$myMethod
                                        #   or myNimbleFunction$myVar
            
            
                                        #	Note that we have cut off '()' in the case of myMethod, so we must inspect the
                                        #   nested symbol for myMethod to determine whether it is a method or variable
            
            nf_fieldName <-as.character(code[[3]])
            objectSymbol = symObj$nfProc$setupSymTab$getSymbolObject(nf_fieldName)
            if(class(objectSymbol)[[1]] == 'symbolMemberFunction'){
                newRunCode <- substitute(nfMethod(NIMBLEFXN, METHODNAME), list(NIMBLEFXN = callerCode, METHODNAME = nf_fieldName))				
                return(newRunCode)
            }
            else {
                                        # We *assume* that if its not a member function, it should be treated with 
                                        # nfVar
                newRunCode <- substitute(nfVar(NIMBLEFXN, VARNAME), list(NIMBLEFXN = callerCode, VARNAME = nf_fieldName))
                return(newRunCode)
            }
        }
        if(class == 'symbolNimbleList'){
                                        #	Code is of the form 
                                        #  myNimbleList$myVar
            nl_fieldName <-as.character(code[[3]])
            newRunCode <- substitute(nfVar(NIMBLELIST, VARNAME), list(NIMBLELIST = callerCode, VARNAME = nl_fieldName))
            return(newRunCode)
        }
        if(class == 'symbolNimbleFunctionList'){
            
                                        #	Code is of the form myNimbleFunctionList[[i]]$foo	(foo should be a method)
                                        #	At this point, we cannot access variables of a nimble function list, ie
                                        #	myNimbleFunctionList[[i]]$myVariable is not allowed
                                        #	If we add this functionality, we will need to look up what foo as we do
                                        #	for the nimbleFunction case above
            
            nf_name <-code[[2]]
            nf_fieldName <- as.character(code[[3]])
            newRunCode <- substitute(nfMethod(NIMBLEFXN, METHODNAME), list(NIMBLEFXN = nf_name, METHODNAME = nf_fieldName))
            return(newRunCode)				
        }
    }
)
    
singleBracket_keywordInfo <- keywordInfoClass(
	keyword = '[',
    processor = function(code, nfProc, RCfunProc){
        if(is.null(nfProc)) return (code)
        class <- symTypeFromSymTab(code[[2]], nfProc$setupSymTab)
        if(class == 'symbolModelValues'){
            if(length(code) < 4) {
                stop(paste0("incorrect syntax for accessing element of modelValues: ", deparse(code)))
            }
            singleMVAccess_ArgList <- list(code = code, modelValues = code[[2]], var = code[[3]], row = code[[4]])
            accessName <- singleModelValuesAccessor_SetupTemplate$makeName(singleMVAccess_ArgList)
            addNecessarySetupCode(accessName, singleMVAccess_ArgList, singleModelValuesAccessor_SetupTemplate, nfProc)
            if(length(code) == 4)
                indexExpr = code[[4]]
            else
                indexExpr = substitute(1)
            
            return(substitute(ACCESS[INDEX], list(ACCESS = as.name(accessName), INDEX = indexExpr) ) )
        }
        return(code)
    }
)    

length_char_keywordInfo <- keywordInfoClass(
    keyword = 'length',
    processor = function(code, nfProc, RCfunProc) {
        if(is.null(nfProc)) return(code)
        if(length(code) < 2) stop('length() used without an argument')
        class <- symTypeFromSymTab(code[[2]], nfProc$setupSymTab)
        if(class == "symbolString") {
            length_char_ArgList <- list(code = code, string = code[[2]])
            accessName <- length_char_SetupTemplate$makeName(length_char_ArgList)
            addNecessarySetupCode(accessName, length_char_ArgList, length_char_SetupTemplate, nfProc)
            return(substitute(LENGTHNAME, list(LENGTHNAME = as.name(accessName))))
        }
        return(code)
    })

nimOptim_model_keywordInfo <- keywordInfoClass(
    keyword = "nimOptim_model",
    processor = function(code, nfProc, RCfunProc) {
        nimOptim_model_keywordInfo_impl(code, nfProc, RCfunProc)
    }
)

nimOptim_model_keywordInfo_impl <- function(code, nfProc, RCfunProc) {
    wrt_arg <- code[['wrt']]
    nodes_arg <- code[['nodes']]
    model_arg <- code[['model']]

    ## This will treat the line of code as if it is calculate(model, nodes, wrt),
    ## which in turn sets up the nodeFxnVector with derivs info.
    newCode <- calculate_keywordInfo$processor(code, nfProc, RCfunProc)
    ## new code has only calculate(nodeFxnVector), so we need to add other argument back in
    newCode[[1]] <- as.name('nimOptim_model')
    newCode$use.gr <- code[['use.gr']]
    newCode$method <-code[['method']]
    newCode$lower <-code[['lower']]
    newCode$upper <-code[['upper']]
    newCode$control <-code[['control']]
    newCode$hessian <-code[['hessian']]
    return(newCode)
}

nimDerivs_keywordInfo <- keywordInfoClass(
  keyword = 'nimDerivs',
  processor = function(code, nfProc, RCfunProc) {
    wrtArgs <- code[['wrt']]
    ## First check to see if nimFxn argument is a method.
    fxnCall <- code[[2]][[1]]
    order <- code[['order']]

    calculateCase <- FALSE
    if(deparse(fxnCall) == 'calculate'){
        calculateCase <- TRUE
    } 
    else if(length(fxnCall) == 3 &&
            (deparse(fxnCall[[1]]) == '$' &&
             deparse(fxnCall[[3]]) == 'calculate')){
        ## re-arrange nimDerivs(model$calculate(...), ...) to nimDerivs(calculate(model, ...), ...)
        code[[2]] <- modelMemberFun_keywordInfo$processor(code[[2]], nfProc)
        code[[2]] <- matchKeywordCode(code[[2]], nfProc)
        calculateCase <- TRUE
    }
    if(calculateCase) {
        innerCode <- code[[2]]
        if(is.null(code$wrt))
            stop("Derivatives of a call to 'calculate()' must have 'wrt' argument specified.")
        innerCode$wrt <- code$wrt
        newCode <- calculate_keywordInfo$processor(innerCode, nfProc, RCfunProc)
        newCode[[1]] <- as.name('nimDerivs_calculate')
        newCode$orderVector <- code$order
        newCode$reset <- if(is.null(code[['reset']])) FALSE else code[['reset']]
        return(newCode)
    }
    ## Not a calculate case:
    if(!is.null(nfProc$origMethods[[deparse(fxnCall)]])) {
      derivMethod <- nfProc$origMethods[[deparse(fxnCall)]]
      derivMethodArgs <- derivMethod$getArgInfo()

      ## Only make the wrt substitution if the names are baked in as character
      wrtArg <- code$wrt
      doPreprocess <- FALSE
      if(is.numeric(wrtArg) | is.logical(wrtArg)) {
        if(any(is.na(wrtArg[1]))) { ## wrt = NULL (default), set to NA, which will become -1 by convertWrtArgToIndices and then c(-1, -1) in the setup code
          doPreprocess <- TRUE
        }
      } else if(is.name(wrtArg)) { ## wrt = some_index_var_maybe_from_setup
        ## do nothing
      } else if(deparse(wrtArg[[1]]) == 'nimC') { ## wrt = c('x', ...)
        if(is.character(wrtArg[[2]]))
          doPreprocess <- TRUE
      } else if(is.character(wrtArg)) { ## wrt = 'x'
        doPreprocess <- TRUE
      } ## In any other setting, it must be an index vector and we leave it alone here.
      ##
      if(doPreprocess) {
        wrt_argList <- list(fxnCall = fxnCall,
                            wrt = code$wrt,
                            #                   vector = wrtArgIndices
                            derivMethodArgs = derivMethodArgs)
        accessName <- wrtVector_setupCodeTemplate$makeName(wrt_argList)
        addNecessarySetupCode(accessName, wrt_argList, wrtVector_setupCodeTemplate, nfProc)
        code[['wrt']] <- substitute(VECNAME,
                                    list(VECNAME = as.name(accessName)))
      }
      ## Check if model and updateNodes are provided.
      ## If so, create a nodeFxnVector for them.
      modelProvided <- !is.null(code[['model']])
      updateNodesProvided <- !is.null(code[['updateNodes']]) ## This can be a list of both updateNodes and constantNodes or a vector of node names
      constantNodesProvided <- !is.null(code[['constantNodes']]) ## This can be a vector of node names or provided in updateNodes
      if(xor(modelProvided, updateNodesProvided || constantNodesProvided))
        stop('nimDerivs arguments model and at least one of updateNodes or constantNodes must be provided together or not at all.')
      if(modelProvided) { ## At this point that means both were provided
        DMUNargList <- list(model = code$model, updateNodes = code[["updateNodes"]], constantNodes = code[["constantNodes"]])
        accessName <- nodeFunctionVector_DerivsModelUpdateNodes_SetupTemplate$makeName(DMUNargList)
        addNecessarySetupCode(accessName, DMUNargList, nodeFunctionVector_DerivsModelUpdateNodes_SetupTemplate, nfProc)
        code[['model']] <- NULL
        if(updateNodesProvided) code[['updateNodes']] <- NULL
        if(constantNodesProvided) code[['constantNodes']] <- NULL
        code[['updateNodesName']] <- accessName
      }
    }
    return(code)
  }
)

derivInfo_keywordInfo <- keywordInfoClass(
  keyword = 'derivInfo',
  processor = function(code, nfProc, RCfunProc) {
    fxnCall <- code[[2]][[1]]
    calculateCase <- FALSE
    if(deparse(fxnCall) == 'calculate'){
      calculateCase <- TRUE
    } 
    else if(length(fxnCall) == 3 &&
              (deparse(fxnCall[[1]]) == '$' &&
                 deparse(fxnCall[[3]]) == 'calculate')){
      ## re-arrange nimDerivs(model$calculate(...), ...) to nimDerivs(calculate(model, ...), ...)
      code[[2]] <- modelMemberFun_keywordInfo$processor(code[[2]], nfProc)
      code[[2]] <- matchKeywordCode(code[[2]], nfProc)
      calculateCase <- TRUE
    }
    if(!calculateCase) {
      stop('The first argument to derivInfo must be a call to calculate(model, ...) or model$calculate(...)')
    }
    innerCode <- code[[2]]
    if(is.null(code$wrt))
      stop("derivInfo must provide 'wrt' argument.")
    innerCode$wrt <- code$wrt
    newCode <- calculate_keywordInfo$processor(innerCode, nfProc, RCfunProc)
    return(newCode)
  }
)

#	KeywordList
keywordList <- new.env()
keywordList[['nimSeq']] <- nimSeq_keywordInfo
keywordList[['getParam']] <- getParam_keywordInfo
keywordList[['getBound']] <- getBound_keywordInfo
keywordList[['values']] <- values_keywordInfo
keywordList[['calculate']] <- calculate_keywordInfo
keywordList[['calculateDiff']] <- calculateDiff_keywordInfo
keywordList[['simulate']] <- simulate_keywordInfo
keywordList[['getLogProb']] <- getLogProb_keywordInfo
keywordList[['nimCopy']] <- nimCopy_keywordInfo
keywordList[['[[']] <- doubleBracket_keywordInfo
keywordList[['$']] <- dollarSign_keywordInfo
keywordList[['[']] <- singleBracket_keywordInfo
keywordList[['besselK']] <- besselK_keywordInfo
keywordList[['dgamma']] <- d_gamma_keywordInfo
keywordList[['pgamma']] <- pq_gamma_keywordInfo
keywordList[['qgamma']] <- pq_gamma_keywordInfo
keywordList[['rgamma']] <- rgamma_keywordInfo
keywordList[['dinvgamma']] <- d_gamma_keywordInfo
keywordList[['pinvgamma']] <- pq_gamma_keywordInfo
keywordList[['qinvgamma']] <- pq_gamma_keywordInfo
keywordList[['rinvgamma']] <- rgamma_keywordInfo
keywordList[['dsqrtinvgamma']] <- d_gamma_keywordInfo
keywordList[['psqrtinvgamma']] <- pq_gamma_keywordInfo
keywordList[['qsqrtinvgamma']] <- pq_gamma_keywordInfo
keywordList[['rsqrtinvgamma']] <- rgamma_keywordInfo
keywordList[['ddexp']] <- d_gamma_keywordInfo
keywordList[['pdexp']] <- pq_gamma_keywordInfo
keywordList[['qdexp']] <- pq_gamma_keywordInfo
keywordList[['rdexp']] <- rgamma_keywordInfo
# can be handled the same as gamma, so include although we have dexp_nimble too
keywordList[['dexp']] <- d_gamma_keywordInfo
keywordList[['pexp']] <- pq_gamma_keywordInfo
keywordList[['qexp']] <- pq_gamma_keywordInfo
keywordList[['rexp']] <- rgamma_keywordInfo

keywordList[['dexp_nimble']] <- d_exp_nimble_keywordInfo
keywordList[['pexp_nimble']] <- pq_exp_nimble_keywordInfo
keywordList[['qexp_nimble']] <- pq_exp_nimble_keywordInfo
keywordList[['rexp_nimble']] <- rexp_nimble_keywordInfo

keywordList[['length']] <- length_char_keywordInfo ## active only if argument has type character

keywordList[['nimDerivs']] <- nimDerivs_keywordInfo
keywordList[['nimOptim_model']] <- nimOptim_model_keywordInfo 
keywordList[['derivInfo']] <- derivInfo_keywordInfo

keywordListModelMemberFuns <- new.env()
keywordListModelMemberFuns[['calculate']] <- modelMemberFun_keywordInfo
keywordListModelMemberFuns[['simulate']] <- modelMemberFun_keywordInfo
keywordListModelMemberFuns[['calculateDiff']] <- modelMemberFun_keywordInfo
keywordListModelMemberFuns[['getLogProb']] <- modelMemberFun_keywordInfo
keywordListModelMemberFuns[['getParam']] <- modelMemberFun_keywordInfo
keywordListModelMemberFuns[['getBound']] <- modelMemberFun_keywordInfo


matchFunctions <- new.env()
matchFunctions[['setSize']] <- function(var, ..., copy = TRUE, fillZeros = TRUE){} 
matchFunctions[['nimC']] <- nimC
matchFunctions[['nimRep']] <- function(x, times = 1, length.out, each = 1) {}
matchFunctions[['nimSeq']] <- nimSeq
matchFunctions[['nimNumeric']] <- nimNumeric 
matchFunctions[['nimInteger']] <- nimInteger 
matchFunctions[['nimLogical']] <- nimLogical 
matchFunctions[['nimMatrix']] <- nimMatrix 
matchFunctions[['nimArray']] <- nimArray 
matchFunctions[['values']] <- function(model, nodes, accessor){}
matchFunctions[['getParam']] <- getParam
matchFunctions[['getBound']] <- getBound
matchFunctions[['calculate']] <- calculate
matchFunctions[['calculateDiff']] <- calculateDiff
matchFunctions[['simulate']] <- simulate
matchFunctions[['getLogProb']] <- getLogProb
matchFunctions[['nimCopy']] <- function(from, to, nodes, nodesTo, row, rowTo, logProb = FALSE, logProbOnly = FALSE){}
matchFunctions[['double']] <- function(nDim, dim, default, ...){}
matchFunctions[['int']] <- function(nDim, dim, default, ...){}
matchFunctions[['nimOptim']] <- nimOptim
matchFunctions[['nimOptimDefaultControl']] <- nimOptimDefaultControl
matchFunctions[['nimEigen']] <- function(squareMat, symmetric = FALSE, only.values = FALSE){}
matchFunctions[['nimSvd']] <- function(mat, vectors = 'full'){}
matchFunctions[['nimOptim_model']] <- function(model, wrt, nodes, use.gr = TRUE, method = "BFGS", lower = -Inf, upper = Inf, control = nimOptimDefaultControl(), hessian = FALSE) {} ## Any changes here need to be reflected in the keyword processor, which has to re-insert arguments to a modified call.
matchFunctions[['nimDerivs']] <- function(nimFxn = NA, order = nimC(0,1,2), dropArgs = NA, wrt = NA, calcNodes = NA, static = FALSE,
                                          model, updateNodes, constantNodes, reset = FALSE) {} ## Avoid NULLs in R default args.
matchFunctions[['derivInfo']] <- derivInfo
matchFunctions[['besselK']] <- function(x, nu, expon.scaled = FALSE){}
matchFunctions[['dgamma']] <- function(x, shape, rate = 1, scale, log = FALSE){}
matchFunctions[['rgamma']] <- function(n, shape, rate = 1, scale){}
matchFunctions[['qgamma']] <- function(p, shape, rate = 1, scale, lower.tail = TRUE, log.p = FALSE){}
matchFunctions[['pgamma']] <- function(q, shape, rate = 1, scale, lower.tail = TRUE, log.p = FALSE){}
matchFunctions[['dinvgamma']] <- matchFunctions[['dsqrtinvgamma']] <- function(x, shape, scale = 1, rate, log = FALSE){}
matchFunctions[['rinvgamma']] <- matchFunctions[['rsqrtinvgamma']] <- function(n, shape, scale = 1, rate){}
matchFunctions[['qinvgamma']] <- matchFunctions[['qsqrtinvgamma']] <- function(p, shape, scale = 1, rate, lower.tail = TRUE, log.p = FALSE){}
matchFunctions[['pinvgamma']] <- matchFunctions[['psqrtinvgamma']] <- function(q, shape, scale = 1, rate, lower.tail = TRUE, log.p = FALSE){}
matchFunctions[['dexp']] <- function(x, rate = 1, log = FALSE){}
matchFunctions[['rexp']] <- function(n, rate = 1){}
matchFunctions[['qexp']] <- function(p, rate = 1, lower.tail = TRUE, log.p = FALSE){}
matchFunctions[['pexp']] <- function(q, rate = 1, lower.tail = TRUE, log.p = FALSE){}
matchFunctions[['dexp_nimble']] <- function(x, rate, scale = 1, log = FALSE){}
matchFunctions[['rexp_nimble']] <- function(n, rate, scale = 1){}
matchFunctions[['qexp_nimble']] <- function(p, rate, scale = 1, lower.tail = TRUE, log.p = FALSE){}
matchFunctions[['pexp_nimble']] <- function(q, rate, scale = 1, lower.tail = TRUE, log.p = FALSE){}

matchFunctions[['ddexp']] <- function(x, location, scale = 1, rate, log = FALSE){}
matchFunctions[['rdexp']] <- function(n, location, scale = 1, rate){}
matchFunctions[['qdexp']] <- function(p, location, scale = 1, rate, lower.tail = TRUE, log.p = FALSE){}
matchFunctions[['pdexp']] <- function(q, location, scale = 1, rate, lower.tail = TRUE, log.p = FALSE){}

matchModelMemberFunctions <- new.env()
matchModelMemberFunctions[['calculate']] <- function(nodes) {}
matchModelMemberFunctions[['calculateDiff']] <- function(nodes) {}
matchModelMemberFunctions[['getLogProb']] <- function(nodes) {}
matchModelMemberFunctions[['simulate']] <- function(nodes, includeData = FALSE) {}
matchModelMemberFunctions[['getParam']] <- function(node, param) {}
matchModelMemberFunctions[['getBound']] <- function(node, bound) {}

# remove ncp from signatures
stripArgs <- function(fname, argNames) {
    if(exists(fname)) {
        args <- formals(eval(as.name(fname)))
        args <- args[-which(names(args) %in% argNames)]
        template <- function() {}
        formals(template) <- args
        return(template)
    } else return(NULL)
}

for(distfun in paste0(c('d','p','q','r'), 'beta'))
    matchFunctions[[distfun]] <- stripArgs(distfun, 'ncp')
for(distfun in paste0(c('d','p','q','r'), 'chisq'))
    matchFunctions[[distfun]] <- stripArgs(distfun, 'ncp')
for(distfun in paste0(c('d','p','q','r'), 't'))
    matchFunctions[[distfun]] <- stripArgs(distfun, 'ncp')
for(distfun in paste0(c('d','p','q','r'), 'nbinom'))
    matchFunctions[[distfun]] <- stripArgs(distfun, 'mu')


# the following are standard in terms of both matchFunctions and keywordList
matchDistList <- list('binom', 'cat', 'dirch', 'interval', 'lnorm', 'logis', 'multi', 'mnorm_chol', 'mvt_chol', 'norm', 'pois', 't_nonstandard', 'unif', 'weibull', 'wish_chol', 'invwish_chol', 'car_normal', 'car_proper')

addDistList2matchFunctions <- function(distList, matchFunEnv){
	for(thisDist in distList){
		pFun <- paste0('p', thisDist)
		qFun <- paste0('q', thisDist)
		rFun <- paste0('r', thisDist)
		dFun <- paste0('d', thisDist)
		
                eval(substitute(matchFunctions[[dFun]] <- DFUN, list(DFUN = as.name(dFun))))
                eval(substitute(matchFunctions[[rFun]] <- RFUN, list(RFUN = as.name(rFun))))
                if(exists(qFun))
                    eval(substitute(matchFunctions[[qFun]] <- QFUN, list(QFUN = as.name(qFun))))
                if(exists(pFun))
                    eval(substitute(matchFunctions[[pFun]] <- PFUN, list(PFUN = as.name(pFun))))
	}
}
          

addDistList2matchFunctions(matchDistList, matchFunctions)

#	processKeyword function to be called by nfProc
processKeyword <- function(code, nfProc, RCfunProc){
  thisKeywordInfo <- keywordList[[ as.character(code[[1]]) ]]
  if(!is.null(thisKeywordInfo))
    return(thisKeywordInfo$processor(code, nfProc, RCfunProc))
  return(code)
}

processKeywordCodeMemberFun <- function(code, nfProc, RCfunProc) { ## handle cases like a$b(c) as one unit
    ## this includes nf$method()
    ## nfList[[i]]$method
    ## model$calculate(nodes)
    dollarSignPart <- code[[1]]
    objectPart <- dollarSignPart[[2]]
    
    isModel <- FALSE
    if(length(objectPart) != 1) isModel <- FALSE ## a case like a[[i]]$b(), which can only be a nimbleFunction list
    else {
        symObj <- nfProc$setupSymTab$getSymbolObject(as.character(objectPart))
        if(is.null(symObj)) stop(paste0("In processKeywordCodeMemberFun: not sure what to do with ", deparse(code)))
        if(inherits(symObj, 'symbolModel'))
            isModel <- TRUE
    }
    if(isModel) {
        thisKeywordInfo <- keywordListModelMemberFuns[[ as.character(dollarSignPart[[3]]) ]]
        if(is.null(thisKeywordInfo)) stop(paste0("In processKeywordCodeMemberFun, don't know what do with: ", deparse(code)))
        rearrangedCode <- thisKeywordInfo$processor(code, nfProc)
        rearrangedCode <- matchKeywordCode(rearrangedCode, nfProc)
        return(processKeywords_recurse(rearrangedCode, nfProc, RCfunProc))
    } else {
        ## same as processKeywords_recurse
        ## first line here creates something like nfMethod(model, method)(args)
        ## which is handled as a chainedCall in later processing
        code[[1]] <- processKeywords_recurse(code[[1]], nfProc, RCfunProc)
        cl <- length(code)
        if(cl >= 2) {
            for(i in 2:cl) {
                code[[i]] <- processKeywords_recurse(code[[i]], nfProc, RCfunProc)
            }
        }
        return(code)
    }
}

processKeywords_recurse <- function(code, nfProc = NULL, RCfunProc) {
    cl = length(code)
    if(cl == 1) {
        if(is.call(code)) {
            if(length(code[[1]]) > 1)
                if(deparse(code[[1]][[1]] == '$'))
                    code <- processKeywordCodeMemberFun(code, nfProc, RCfunProc)
                else
                    code[[1]] <- processKeywords_recurse(code[[1]], nfProc, RCfunProc)
        }
        return(code)
    }
    
    if(length(code[[1]]) == 1) {
        code <- processKeyword(code, nfProc, RCfunProc)
    }
    
    cl = length(code)
    
    if(is.call(code)) {
        if(length(code[[1]]) > 1) {
            if(deparse(code[[1]][[1]] == '$')) {
                code <- processKeywordCodeMemberFun(code, nfProc, RCfunProc) ## case like model$calculate(nodes)
                return(code) ## don't recurse on arguments of anything in this category
            }
            code[[1]] <- processKeywords_recurse(code[[1]], nfProc, RCfunProc)
        }
        if(cl >= 2) {
            for(i in 2:cl) {
                code[[i]] <- processKeywords_recurse(code[[i]], nfProc, RCfunProc)
            }
        }
    }
    return(code)
}

#####	SETUPCODE TEMPLATES

wrtVector_setupCodeTemplate <- setupCodeTemplateClass(
    makeName = function(argList){Rname2CppName(paste0('wrtVec_',
                                                      deparse(argList$fxn),
                                                      paste(deparse(argList$wrt), sep='_'),
                                                      '_'))},
  codeTemplate = quote(
  {
      WRTVEC <- nimble:::convertWrtArgToIndices(WRT, #code$wrt
                                                 DERIVMETHODARGS,
                                                 FXNNAME)
      if(length(WRTVEC) == 1)
          WRTVEC <- c(WRTVEC, -1)
  }
  ),
  makeCodeSubList = function(resultName, argList){
      list(WRTVEC = as.name(resultName),
           WRT = argList$wrt
          ,DERIVMETHODARGS = argList$derivMethodArgs
          ,FXNNAME = deparse(argList$fxnCall)
           )##argList$vector)
  }
)

length_char_SetupTemplate <- setupCodeTemplateClass(
    makeName = function(argList) {Rname2CppName(paste0(paste("length", deparse(argList$string), sep='_'), '_KNOWN_'))},
    codeTemplate = quote(LENGTHNAME <- CODE),
    makeCodeSubList = function(resultName, argList) {
        list(LENGTHNAME = as.name(resultName),
             CODE = argList$code)
    })

modelVariableAccessorVector_setupCodeTemplate <- setupCodeTemplateClass(
	#Note to programmer: required fields of argList are model, nodes and logProb
    makeName = function(argList) {Rname2CppName(paste(deparse(argList$model), deparse(argList$nodes), 'access_logProb', deparse(argList$logProb),'LPO', deparse(argList$logProbOnly), sep = '_'))},
    codeTemplate = quote( ACCESSNAME <- nimble:::modelVariableAccessorVector(MODEL, NODES, logProb = LOGPROB, logProbOnly = LOGPROBONLY) ),
    makeCodeSubList = function(resultName, argList) {
        list(ACCESSNAME = as.name(resultName),
             MODEL = argList$model,
             NODES = argList$nodes,
             LOGPROB = argList$logProb,
             LOGPROBONLY = argList$logProbOnly)
    })

copierVector_setupCodeTemplate <- setupCodeTemplateClass(
    makeName = function(argList) {Rname2CppName(paste0(argList$accessFrom_name, '_', argList$accessTo_name))},
    codeTemplate = quote( COPIERNAME <- nimble:::copierVector(ACCESS_FROM, ACCESS_TO, ISMVFROM, ISMVTO) ),
    makeCodeSubList = function(resultName, argList) {
        list(COPIERNAME = as.name(resultName),
             ACCESS_FROM = as.name(argList$accessFrom_name),
             ACCESS_TO   = as.name(argList$accessTo_name),
             ISMVFROM    = as.integer(argList$isMVfrom),
             ISMVTO      = as.integer(argList$isMVto)) 
    })
    

modelValuesAccessorVector_setupCodeTemplate <- setupCodeTemplateClass(
	#Note to programmer: required fields of argList are model, nodes and logProb

    makeName = function(argList) {Rname2CppName(paste(deparse(argList$model), deparse(argList$nodes), 'access_logProb', deparse(argList$logProb), 'LPO', deparse(argList$logProbOnly), deparse(argList$row), sep = '_'))},
    codeTemplate = quote( ACCESSNAME <- nimble:::modelValuesAccessorVector(MODEL, NODES, logProb = LOGPROB, logProbOnly = LOGPROBONLY) ),
	makeCodeSubList = function(resultName, argList) {
        list(ACCESSNAME = as.name(resultName),
             MODEL = argList$model,
             NODES = argList$nodes,
             LOGPROB = argList$logProb,
             LOGPROBONLY = argList$logProbOnly)
    })


nodeFunctionVector_WithDerivsOutput_SetupTemplate <- setupCodeTemplateClass(
  makeName = function(argList){
    Rname2CppName(paste(deparse(argList$model),
                        deparse(argList$nodes),
                        'nodeFxnVector_WithDerivsOutput', ## The "_derivs_" tag is used later to determine the right class - klugey
                        deparse(argList$includeData),
                        if(argList$sortUnique) "SU" else "notSU",
                        '_derivs_',
                        sep = '_')
                  )
  },
  codeTemplate = quote(
    NODEFXNVECNAME <- nimble:::nodeFunctionVector_WithDerivsOutputNodes(
      model = MODEL,
      calcNodes = CALCNODES,
      excludeData = EXCLUDEDATA,
      sortUnique = SORTUNIQUE)
  ), 
  makeCodeSubList = function(resultName, argList){
    list(NODEFXNVECNAME = as.name(resultName),
         MODEL = argList$model,
         CALCNODES = argList$nodes,
         EXCLUDEDATA = !argList$includeData,
         SORTUNIQUE = argList$sortUnique
         )
  })


nodeFunctionVector_DerivsModelUpdateNodes_SetupTemplate <- setupCodeTemplateClass(
  makeName = function(argList){
    Rname2CppName(paste(deparse(argList$model),
                        deparse(argList[['updateNodes']]),
                        deparse(argList[['constantNodes']]),
                        'nodeFxnVector_DerivsModelUpdateNodes_derivs_', ## The "_derivs_" tag is used later to determine the right class - klugey
                        sep = '_')
                  )
  },
  codeTemplate = quote(
    NODEFXNVECNAME <- nimble:::nodeFunctionVector_DerivsModelUpdateNodes(
      model = MODEL,
      updateNodes = UPDATENODES,
      constantNodes = CONSTANTNODES)
  ), 
  makeCodeSubList = function(resultName, argList){
    list(NODEFXNVECNAME = as.name(resultName),
         MODEL = argList$model,
         UPDATENODES = argList[["updateNodes"]],
         CONSTANTNODES = argList[["constantNodes"]]
         )
  })

nodeFunctionVector_SetupTemplate <- setupCodeTemplateClass(
                                        #Note to programmer: required fields of argList are model, nodes and includeData
    
    makeName = function(argList){
        Rname2CppName(paste(deparse(argList$model),
                            deparse(argList$nodes),
                            'nodeFxnVector_includeData',
                            deparse(argList$includeData),
                            if(argList$sortUnique) "SU" else "notSU",
                            if(is.null(argList$wrtNodes)) '' else '_derivs_', sep = '_')
                      )
    },
    codeTemplate = quote(
        NODEFXNVECNAME <- nimble:::nodeFunctionVector(model = MODEL,
                                                      nodeNames = NODES,
                                                      wrtNodes = WRTNODES,
                                                      excludeData = EXCLUDEDATA,
                                                      sortUnique = SORTUNIQUE,
                                                      errorContext = ERRORCONTEXT)
    ), 
    makeCodeSubList = function(resultName, argList){
        list(NODEFXNVECNAME = as.name(resultName),
             MODEL = argList$model,
             NODES = argList$nodes,
             WRTNODES = argList$wrtNodes,
             EXCLUDEDATA = !argList$includeData,
             SORTUNIQUE = argList$sortUnique,
             ERRORCONTEXT = argList$errorContext
             )
    })


paramInfo_SetupTemplate <- setupCodeTemplateClass(
    #Note to programmer: required fields of argList are model, node and param
    makeName = function(argList){Rname2CppName(paste(deparse(argList$model), deparse(argList$node), deparse(argList$param), 'paramInfo', sep='_'))},
    makeOtherNames = function(name,argList) {Rname2CppName(paste0(name,'_ID'))},
    codeTemplate = quote({
        PARAMINFONAME <- nimble:::makeParamInfo(model = MODEL, nodes = NODE, param = PARAM, vector = HASINDEX )
        PARAMIDNAME <- PARAMINFONAME$paramID
        PARAMINFONAME$paramID <- NULL
       }),
    makeCodeSubList = function(resultName, argList){
        list(PARAMINFONAME = as.name(resultName),
             PARAMIDNAME = as.name(paste0(resultName,'_ID')),
             MODEL = argList$model,
             NODE = argList$node,
             PARAM = argList$param,
             HASINDEX = argList$hasIndex)
    })

boundInfo_SetupTemplate <- setupCodeTemplateClass(
    #Note to programmer: required fields of argList are model, node and param
    makeName = function(argList){Rname2CppName(paste(deparse(argList$model), deparse(argList$node), deparse(argList$bound), 'boundInfo', sep='_'))},
    makeOtherNames = function(name,argList) {Rname2CppName(paste0(name,'_ID'))},
    codeTemplate = quote({
        BOUNDINFONAME <- nimble:::makeBoundInfo(MODEL, NODE, BOUND)
        BOUNDIDNAME <- BOUNDINFONAME$boundID
        BOUNDINFONAME$boundID <- NULL
       }),
    makeCodeSubList = function(resultName, argList){
        list(BOUNDINFONAME = as.name(resultName),
             BOUNDIDNAME = as.name(paste0(resultName,'_ID')),
             MODEL = argList$model,
             NODE = argList$node,
             BOUND = argList$bound)
    })

allLHSNodes_SetupTemplate <- setupCodeTemplateClass(
	#Note to programmer: required fields of argList are model

	makeName = function(argList){
		Rname2CppName(paste('allLHSnodes', deparse(argList$model), sep = '_'))
	},
	codeTemplate = quote(NODENAMES <- MODEL$getMaps('nodeNamesLHSall')),
	makeCodeSubList = function(resultName, argList){
		list(NODENAMES = as.name(resultName),
			MODEL = argList$model)
	})
	
allModelNodes_SetupTemplate <- setupCodeTemplateClass(
	#Note to programmer: required fields of argList are model

	makeName = function(argList){
		Rname2CppName(paste('allModelNodes', deparse(argList$model), sep = '_'))
	},
	codeTemplate = quote(NODENAMES <- MODEL$getVarNames()),
	makeCodeSubList = function(resultName, argList){
		list(NODENAMES = as.name(resultName),
			MODEL = argList$model)
	})	
	
allModelValuesVars_SetupTemplate <- setupCodeTemplateClass(
	#Note to programmer: required fields of argList are modelValues

	makeName = function(argList){
		Rname2CppName(paste('allMVVars', deparse(argList$modelValues), sep = '_'))
	},
	codeTemplate = quote(NODENAMES <- MODELVALUES$getVarNames(includeLogProb = FALSE)),	
		
	makeCodeSubList = function(resultName, argList){
		list(NODENAMES = as.name(resultName),
			MODELVALUES = argList$modelValues)
	})	
	
code2Name_fromArgList <- function(argList)
	Rname2CppName(deparse(argList$code))	
	
singleVarAccess_SetupTemplate <- setupCodeTemplateClass(
	#Note to progammer: required fields of argList are 'code' (raw code to be processed), model and var

	makeName = code2Name_fromArgList,

	codeTemplate = quote(SINGLEACCESSOR <- nimble:::singleVarAccess(MODEL, VAR)),

	makeCodeSubList = function(resultName, argList){
		list(SINGLEACCESSOR = as.name(resultName),
			MODEL = argList$model,
			VAR = argList$var)	
	})	
	
singleModelIndexAccess_SetupTemplate <- setupCodeTemplateClass(
	#Note to progammer: required fields of argList are code, varAndIndices, node (character) and model(expression)
	makeOtherNames = function(name, argList){ paste0(name, '_flatIndex')},
	makeName = code2Name_fromArgList,
	
	codeTemplate = quote({
		VARANDINDICES <- nimble:::nimbleInternalFunctions$getVarAndIndices(NODEVARNAME)
		NEWVARNAME <- as.character(VARANDINDICES$varName)
		MFLATINDEX <- nimble:::varAndIndices2flatIndex(VARANDINDICES, MODELVAREXPR$getVarInfo(NEWVARNAME))
		VARACCESSOR <- nimble:::singleVarAccess(MODELVAREXPR, NEWVARNAME, useSingleIndex = TRUE)
	}),
	makeCodeSubList = function(resultName, argList){
		list(VARACCESSOR = as.name(resultName),
			VARANDINDICES = as.name(paste0(resultName, '_varAndIndices') ),
			NEWVARNAME = as.name(paste0(resultName, '_newVarName')),
			NODEVARNAME = argList$nodeExpr,
			MFLATINDEX = as.name(paste0(resultName, '_flatIndex')),
			MODELVAREXPR = argList$model)
	})
	
map_SetupTemplate <- setupCodeTemplateClass(
	#Note to programmer: required fields of argList are code, model
	makeName  = code2Name_fromArgList,
	makeOtherNames = function(name, argList){
		output <- character()
		output[1] = paste0(name, '_strides')
		output[2] = paste0(name, '_offset')
		output[3] = paste0(name, '_sizes')
		return(output)
	},
	codeTemplate = quote({
		VARANDINDICES <- nimble:::nimbleInternalFunctions$getVarAndIndices(NODEVARNAME)
		NEWVARNAME <- as.character(VARANDINDICES$varName)
                map_SetupTemplate_vi <- MODEL$getVarInfo(NEWVARNAME)
		map_SetupTemplate_mapParts <- nimble:::varAndIndices2mapParts(VARANDINDICES, map_SetupTemplate_vi$maxs, map_SetupTemplate_vi$nDim)
		MSTRIDES <- map_SetupTemplate_mapParts$strides
		MOFFSET <- map_SetupTemplate_mapParts$offset
		MSIZES <- map_SetupTemplate_mapParts$sizes
		VARACCESSOR <- nimble:::singleVarAccess(model, NEWVARNAME)
	}),
	makeCodeSubList = function(resultName, argList){
		list(VARACCESSOR = as.name(resultName),
		NODEVARNAME =	argList$nodeExpr,
		NEWVARNAME = as.name(paste0(resultName, '_newVarName')),
		VARANDINDICES = as.name(paste0(resultName, '_varAndIndices')),
		MODEL = argList$model,
		MSTRIDES = as.name(paste0(resultName, '_strides')),
		MOFFSET = as.name(paste0(resultName, '_offset')),

                     MSIZES = as.name(paste0(resultName, '_sizes')))
	})
	
singleModelValuesAccessor_SetupTemplate <- setupCodeTemplateClass(
	#Note to programmer: required fields of argList are modelValues, var, row, code
	makeName = code2Name_fromArgList,
	codeTemplate = quote({
		MVACCESS <- nimble:::singleModelValuesAccess(MODELVALUES, VAR)
	}),
	makeCodeSubList = function(resultName, argList){
		list(MVACCESS = as.name(resultName),
		MODELVALUES = argList$modelValues,
		VAR = argList$var)
	}
)





#### KEYWORD PROCESSING UTILITIES


isCodeArgBlank <- function(code, arg){
	#return(nchar(code[[arg]]) == 0)
	return(is.null(code[[arg]]))
}

# Utility functions to make things a little neater
isSetupCodeGenerated <- function(name, nfProc)
    name %in% nfProc$newSetupOutputNames
addSetupCodeNames <- function(name, otherNames, nfProc)
    nfProc$newSetupOutputNames <- c(name, otherNames, nfProc$newSetupOutputNames)
addBlockFromCppName <- function(name, nfProc)
    nfProc$blockFromCppNames <- c(name, nfProc$blockFromCppNames)
addNewCode <- function(name, subList, template, nfProc)
    nfProc$newSetupCode[[name]] <- eval(substitute(substitute(TEMPLATE, subList), list(TEMPLATE = template$codeTemplate) ) )

addNecessarySetupCode <- function(name, argList, template, nfProc, allowToCpp = TRUE){
    if(is.null(nfProc)) stop("Trying to add setup code for a nimbleFunction with no setup code.")
    if(!isSetupCodeGenerated(name, nfProc)){
        addSetupCodeNames(name, template$makeOtherNames(name, argList), nfProc)
        subList <- template$makeCodeSubList(name, argList)
        addNewCode(name, subList, template, nfProc)
        if(!allowToCpp) addBlockFromCppName(name, nfProc) ## ignores makeOtherNames for now
    }
}

getSymObj_recurse <- function(code, symTab, recurse = FALSE) { #code will be like a$b$c with expectation all are NF or NL
    if(length(code) > 1) {
        if(deparse(code[[1]]) != '$') return(NULL) ## This is valid if we have makeNewNimbleListObject(...)$member or foo()$member
        firstArg <- code[[2]]
        memberName <- deparse(code[[3]])
        symTab <- getSymObj_recurse(firstArg, symTab, recurse = TRUE) ## when recursing, return the symTab
    } else {
        memberName <- deparse(code)
    }
    symObject <- if(is.null(symTab)) NULL else symTab$getSymbolObject(memberName)
    if(recurse) {
        if(inherits(symObject, 'symbolNimbleFunction')) return(symObject$nfProc$setupSymTab)
        if(inherits(symObject, 'symbolNimbleList')) return(symObject$nlProc$symTab) ## can only be known if it was created in setup code
        if(is.null(symObject)) return(NULL) ## will assume to be nimbleList, only kind of nested data structure that could be of unknown type right now
    } else {
        return(symObject)
    }
    stop(paste0('Problem (ii) working through ', deparse(code)), .call = FALSE)
}

symTypeFromSymTab <- function(codeName, symTab, options = character(0) ){
    if(is.language(codeName))
        codeName <- as.character(codeName)
    if(length(codeName) > 1)
        return('NULL')
    class <- class(symTab$symbols[[codeName]])[1]
    if(length(options) == 0)
        return(class)
    if(!(class %in% options))
        return(NULL)  ## nimbleList that was not constructed in setup code 
    return(class)
}

isSymbolType <- function(symTab, varName, symType)
	inherits(symTab$symbols[[varName]], symType)

matchAndFill.call <- function(def, call){
    ##   Note re: matchAndFill.call(function(a = 1, ..., b = 2){}, quote(foo(b = 1, 2, 3))): behavior on this is due to match.call
    theseFormals <- formals(def)
    formalNames <- names(theseFormals) # formalArgs are the arguments that are defined, i.e. does NOT include anything that is from the args "..."
    theseFormals <- theseFormals[nchar(theseFormals) > 0]
    matchedCall <- match.call(def, call) # problem with match.call for our needs is it omits formals that were not provided
    missingArgs <- which(!(names(theseFormals) %in% names(matchedCall)))
    for(ind in missingArgs){ ## this puts back in anything omitted, but order may become wrong
        name <- names(theseFormals)[ind]
        matchedCall[[name]] <- (if(is.null(theseFormals[[name]])) list(NULL) else theseFormals[[name]])
    }

    newCall <- matchedCall[1]

    if(is.null(names(matchedCall))) names(matchedCall) <- c("CALL_", rep("", length(matchedCall) - 1)) ## strangely assigning all "" values results in NULL
    indexAdditionalArgs <- which(!(names(matchedCall)[-1] %in% formalNames))
    
    for(thisArgName in formalNames){					# This is to get the order of the arguments correctly and to insert unmatched arguemnts to ... location if appropriate
        if(thisArgName == '...') {
            for(thisIndex in indexAdditionalArgs) {
                thisName <- names(matchedCall)[thisIndex+1]
                if(thisName=="")
                    newCall[[length(newCall) + 1]] <- matchedCall[[thisIndex + 1]]
                else {
                    newCall[[thisName]] <- matchedCall[[thisName]]
                }
            }
        } else {        
            thisArg <- matchedCall[[thisArgName]]
            if(!is.null(thisArg))
                newCall[[thisArgName]] <- thisArg
        }
    }
    
    return(newCall)
}


determineNdimsFromNfproc <- function(modelExpr, varOrNodeExpr, nfProc) {
    allNDims <- lapply(nfProc$instances, function(x) {
        model <- eval(modelExpr, envir = x)
        if(!exists(as.character(varOrNodeExpr), x, inherits = FALSE) ) {
            stop(paste0('Problem accessing node or variable ', deparse(varOrNodeExpr), '.'), call. = FALSE)
        }
        lab <- eval(varOrNodeExpr, envir = x)
        if(length(lab) != 1)
            stop(paste0("Length of ",
                        deparse(varOrNodeExpr),
                        " requested from ",
                        deparse(modelExpr),
                        " using '[[' is ",
                        length(lab),
                        ". It must be 1." )
               , call. = FALSE)
        varAndIndices <- nimbleInternalFunctions$getVarAndIndices(lab)
        determineNdimFromOneCase(model, varAndIndices)
    } )
    return(allNDims)
}

## from a$b(), goal is to get symbolObject for a
## from a$b$c(), goal is to get symbolObject for b
matchKeywordCodeMemberFun <- function(code, nfProc) {  ## handles cases like a$b(c) as one unit so the member function template for b can be looked up
    dollarSignPart <- code[[1]] ## we already checked that code[[1]][[1]] is '$' before calling this function
    leftSide <- dollarSignPart[[2]]  ## could have further recursing needed. at the moment, nimbleFunction list can only appear at the end of nesting (member data are not part of inheritance)
    rightSide <- dollarSignPart[[3]]
    memFunName <- deparse(rightSide)
    nestedLeftSide <- FALSE

    ## Step 1: Find the symbol object for whatever is on left-hand side of $
    ## This may involve recursion
    ## Only case where the symbol object can be missing (NULL) is nimbleListDef$new() if nimbleListDef if from global env
    if(!is.null(nfProc)) { ## If nfProc is null, we are in an RCfunction and the only valid case is nimbleListDef$new()
        if(length(leftSide) != 1) {
            if(deparse(leftSide[[1]])=='$') {
                symObj <- getSymObj_recurse(leftSide, nfProc$setupSymTab)
                nestedLeftSide <- TRUE
            } else {
                nfNestedCall <- leftSide[[1]]
                if(length(nfNestedCall) != 1) stop(paste0("Cannot handle this expression: ", deparse(code)))
                if(deparse(nfNestedCall) != '[[') stop(paste0("Cannot handle this expression: ", deparse(code)))
                leftLeftSide <- leftSide[[2]]
                if(length(leftLeftSide) != 1) {
                    if(deparse(leftLeftSide[[1]]) == '$') { ##a$b[[i]]$foo()
                        symObj <- getSymObj_recurse(leftLeftSide, nfProc$setupSymTab)
                        nestedLeftSide <- TRUE
                    } else stop('Problem processing something like a$b[[i]]$foo()')
                } else {
                    nfListName <- deparse(leftSide[[2]])
                    if(nfProc$setupSymTab$symbolExists(nfListName)) { ## look in symbol table
                        symObj <- nfProc$setupSymTab$getSymbolObject(nfListName)
                    }
                }
                
            }
        } else {
            symTab <- nfProc$setupSymTab
            symObj <- symTab$getSymbolObject(deparse(leftSide))
        }
    } else {
        symTab <- symObj <- NULL
    }
    
    if(memFunName=='new') { ## this is unique because in non-nested mode, this can be looking for a nlDef in global environment (or possibly elsewhere, but not dealt with)
        ## symObj can be null here
        if(is.null(symObj)) {
            if(nestedLeftSide) stop('Cannot find nested nimbleList definition')
            nlGenName <- deparse(leftSide)
            if(exists(nlGenName, where = globalenv())) {
                possibleNLgen <- get(nlGenName, envir = globalenv())

                if(is.nlGenerator(possibleNLgen)) {
                    thisFunctionMatch <- makeNimbleListTemplateWithBlankFirstArg(nl.getListDef(possibleNLgen))
                } else {
                    stop(paste0('problem with ', deparse(code)))
                }
            } else {
                stop(paste0('problem with ', deparse(code)))
            }
        } else {
            thisFunctionMatch <- makeNimbleListTemplateWithBlankFirstArg(nl.getListDef(symObj$nlProc$nlGenerator))
        }
        for(i in seq_along(code[-1])) code[[i+1]] <- matchKeywords_recurse(code[[i+1]], nfProc)
        code[[ length(code) + 1]] <- leftSide      ## add '.LEFTSIDE = leftSide' arg to code
        names(code)[length(code)] <- '.LEFTSIDE'
        code[[1]] <- as.name('makeNewNimbleListObject') ## modify first arg of code to be desired name of call
        return(matchAndFill.call(thisFunctionMatch, code ) ) ## should create makeNewNimbleListObject( nimbleListCallMaybeNested, var1, var2, etc.)
    }

    if(is.null(symObj)) stop('Problem looking up object')
    
    if(symObj$type == 'nimbleFunction') {
        thisRCfunProc <- symObj$nfProc$RCfunProcs[[memFunName]] 
        if(is.null(thisRCfunProc)) stop(paste0("Cannot handle this expression (member function may not exist): ", deparse(code)), call. = FALSE)
        thisFunctionMatch <- thisRCfunProc$RCfun$template
        return(matchAndFill.call(thisFunctionMatch, code ) )
    } 
    if(inherits(symObj, 'symbolModel')) {
        if(nestedLeftSide) stop('Access to a model cannot be nested.')
        thisFunctionMatch <- matchModelMemberFunctions[[ memFunName ]]
        if(is.null(thisFunctionMatch)) stop(paste0("Cannot handle this expression (looks like a model with an invalid member function call?): ", deparse(code)))
        return(matchAndFill.call(thisFunctionMatch, code) )
    }
    if(inherits(symObj, 'symbolNimbleFunctionList')) {
        thisBaseClass <- symObj$baseClass
        thisFunctionMatch <- environment(symObj$baseClass)$methodList[[memFunName]]$template
        return(matchAndFill.call(thisFunctionMatch, code ) )
    }
    stop(paste0("Cannot handle this expression: ", deparse(code))) 
}


matchKeywordCode <- function(code, nfProc){
    callName <- as.character(code[[1]])
    thisFunctionMatch <- matchFunctions[[ callName ]]
    ## see if this is a member function of an nf object
    if(!is.null(nfProc)) {
        modCallName <- callName
        if(nfProc$setupSymTab$symbolExists(modCallName)) {
            symObj <- nfProc$setupSymTab$getSymbolObject(modCallName)
            if(inherits(symObj, "symbolMemberFunction")) {
                thisRCfunProc <- nfProc$RCfunProcs[[modCallName]]
                if(is.null(thisRCfunProc)) stop(paste0("Cannot handle this expression (looks like a member function but something is wrong): ", deparse(code)), call. = FALSE)
                thisFunctionMatch <- thisRCfunProc$RCfun$template
                return(matchAndFill.call(thisFunctionMatch, code ) )
            }
        }
    }
    
    ## see if this is a call to an RCfunction
    if(is.null(thisFunctionMatch)) {
        if(exists(callName)) {
            callObj <- get(callName)
            if(is.rcf(callObj)) {
                thisFunctionMatch <- callObj
            }
        }
    }
    
    if(!is.null(thisFunctionMatch))
        return(matchAndFill.call(thisFunctionMatch, code ) )
    return(code)
}

matchKeywords_recurse <- function(code, nfProc = NULL) {
    cl = length(code)
    if(cl == 1){ ## There are no arguments
        if(is.call(code)){  
            if(length(code[[1]]) > 1) {
                if(deparse(code[[1]][[1]]) == '$') code <- matchKeywordCodeMemberFun(code, nfProc)
                else
                    code[[1]] <- matchKeywords_recurse(code[[1]], nfProc) ## recurse on the "a$b" part of a$b() (or the "a(b)" part of a(b)()), etc
            } else
                code <- matchKeywordCode(code, nfProc)
        }
        return(code)
    }
    if(length(code[[1]]) == 1) ## a simple call like a(b,c), not a$b(c)
        code <- matchKeywordCode(code, nfProc)
    
    if(is.call(code)) {
        if(length(code[[1]]) > 1) {
            if(deparse(code[[1]][[1]]) == '$') code <- matchKeywordCodeMemberFun(code, nfProc) ## handle a$b(c) as one unit
            else code[[1]] <- matchKeywords_recurse(code[[1]], nfProc) ## handle "a(b)" part of a(b)(c), which is probably *never* triggered
        }
        if(cl >= 2) { ## recurse through arguments
            for(i in 2:cl) {
                code[[i]] <- matchKeywords_recurse(code[[i]], nfProc)
            }
        }
    }
    return(code)
}


makeSingleIndexAccessExpr <- function(newName, newNameExpr) {
    codeNames <- makeSingleIndexAccessCodeNames(newName)
    subList <- lapply(codeNames, as.name)
    ans <- substitute( name[MflatIndex], c(list(name = newNameExpr), subList))
    ans
}

## want map(name, nDim, offset, sizelist, stridelist)
## this is a unique case, where sizelist and stridelist are just lists
## stuck in there
makeMapAccessExpr <- function(newName, newNameExpr, nDim) { ## newNameExpr not used any more!
    codeNames <- makeMapSetupCodeNames(newName)
    subList <- lapply(codeNames, as.name)
    if(nDim == 0) { ## not sure this can happen
        sizeExprList <- strideExprList <- list()
    }
    if(nDim == 1) {
        sizeExprList <- list(substitute(Msizes, subList))
        strideExprList <- list(substitute(Mstrides, subList))
    }
    if(nDim >= 2) {
        sizeExprList <- rep(list( substitute(Msizes[1], subList)), nDim)
        for(i in 1:nDim) sizeExprList[[i]][[3]] <- i
        strideExprList <- rep(list( substitute(Mstrides[1], subList)), nDim)
        for(i in 1:nDim) strideExprList[[i]][[3]] <- i
    }
    ans <- substitute(map( name, nDim, Moffset, sizes, strides),
                      c(subList, list(nDim = nDim, name = newName, sizes = sizeExprList, strides = strideExprList)))
    ans
}

determineNdimFromOneCase <- function(model, varAndIndices) {
    varInfo <- try(model$getVarInfo(as.character(varAndIndices$varName)))
    if(inherits(varInfo, 'try-error')) stop(paste0('In determineNdimFromOneCase: error in extracting varInfo for ', varAndIndices$varName), call. = FALSE)
    varNdim <- varInfo$nDim
    if(length(varAndIndices$indices) == 0) return(varNdim)
    if(length(varAndIndices$indices) != varNdim) {
        stop(paste0('Error, wrong number of dimensions in a node label for ', varAndIndices$varName, '.  Expected ',varNdim,' indices but got ', length(varAndIndices$indices),'.'))
    }
    dropNdim <- sum(unlist(lapply(varAndIndices$indices, is.numeric)))
    return(varNdim - dropNdim)
}



## steps here are similar to makeMapExprFromBrackets, but that uses exprClasses

varAndIndices2mapParts <- function(varAndIndices, sizes, nDim) {
    indices <- varAndIndices$indices
    ## put together offsetExpr, sizeExprs, strideExprs
    ## need sizes to get strides
    if(length(sizes) == 0) sizes <- 1
    if(nDim > 0 & length(indices)==0) {
        blockBool <- rep(TRUE, nDim)
        firstIndexRexprs <- rep(list(1), nDim)
        targetSizes <- sizes
    } else {
        nDim <- length(indices) ## may be redundant/moot
        firstIndexRexprs <- vector('list', nDim)
        targetSizes <- integer(nDim)
        blockBool <- rep(FALSE, nDim)
        for(i in seq_along(indices)) {
            if(is.blank(indices[[i]])) {
                blockBool[i] <- TRUE
                firstIndexRexprs[[i]] <- 1
                targetSizes[i] <- sizes[i]
            }
            else if(is.numeric(indices[[i]])) {
                firstIndexRexprs[[i]] <- indices[[i]]
            } else {
                ## better be :
                if(indices[[i]][[1]] != ":") stop("error, expecting : here")
                blockBool[i] <- TRUE
                firstIndexRexprs[[i]] <- indices[[i]][[2]]
                targetSizes[i] <- indices[[i]][[3]] - indices[[i]][[2]] + 1
            }
        }
    }
    strides <- c(1, cumprod(sizes[-length(sizes)]))
    sourceStrideRexprs <- as.list(strides)
    targetOffsetRexpr <- makeOffsetRexpr(firstIndexRexprs, sourceStrideRexprs)
    targetStrides <- strides[blockBool]
    targetSizes <- targetSizes[blockBool]
    list(offset = eval(targetOffsetRexpr),
         sizes = targetSizes,
         strides = targetStrides)
}


getVarAndIndices <- function(code) {
    if(is.character(code)) code <- parse(text = code, keep.source = FALSE)[[1]]
    if(length(code) > 1) {
        if(code[[1]] == '[') {
            varName <- code[[2]]
            indices <- as.list(code[-c(1,2)])
        } else {
            stop(paste('Error:', deparse(code), 'is a malformed node label.'))
        }
    } else {
        varName <- code
        indices <- list()
    }
    list(varName = varName, indices = indices)
}

## This takes the indices field returned by getVarAndIndices and turns it into a matrix
## e.g. from getVarAndIndices('x[1:3, 2:4]'), we have varName = 'x' and indices = list(quote(1:3), quote(2:4))
## indexExprs2matrix takes the indices and makes [1 3; 2 4]

varAndIndices2flatIndex <- function(varAndIndices, varInfo) {
    if(length(varInfo$maxs) == 0) return(1) ## A -1 is done automatically, later, so here we should stay in R's 1-based indexing
    sizes <- varInfo$maxs
    strides <- c(1, cumprod(sizes[-length(sizes)]))
    flatIndex <- 1 + sum((unlist(varAndIndices$indices)-1) * strides)
    flatIndex
}


makeMapSetupCodeNames <- function(baseName) {
    list(Mstrides = paste0(baseName, '_strides'),
         Msizes = paste0(baseName, '_sizes'),
         Moffset = paste0(baseName, '_offset'))
}


makeSingleIndexAccessCodeNames <- function(baseName) {
    list(MflatIndex = paste0(baseName, '_flatIndex'))
}

handleScaleAndRateForGamma <- function(code){
  ## also handles core R dexp, as well as ddexp, and dinvgamma
  scaleArg <- code$scale
  rateArg <- code$rate
  if(is.null(scaleArg) && is.null(rateArg))	stop('neither scale nor rate defined in dgamma, invgamma, dexp, or ddexp')
  codeName <- deparse(code[[1]])
  dist <- substring(codeName, 2, nchar(codeName))
  if(dist == 'invgamma' || dist == 'sqrtinvgamma') {  ## For [drpq]invgamma
    if(is.null(rateArg)) {
      rateArg <- substitute(1.0/(A), list(A = code$scale)) 
      code$scale <- rateArg
      names(code)[which(names(code) == 'scale')] <- 'rate'  # to preserve correct order
    }
    code$scale <- NULL
  } else if(dist == 'dexp') { ## This is for [drpq]dexp
    ## The logic here is trickier.  scale has a default value and is canonical (what is needed for C).
    ## If rate was provided, then by the time we are here, matchFunctions has been used (matchAndFill.call),
    ## so there will also be a scale, since it has a default
    setScaleArg <- FALSE
    if(is.null(scaleArg)) setScaleArg <- TRUE ## scale is not expected to be null ever (see previous comment), but this is defensive.
    if(!is.null(scaleArg) & !is.null(rateArg)) setScaleArg <- TRUE ## Both are there, so set scale from the provided rate
    if(setScaleArg) {
      code$scale <- NULL
      scaleArg <- substitute(1.0/(A), list(A = code$rate)) 
      code$rate <- scaleArg
      names(code)[which(names(code) == 'rate')] <- 'scale'  # to preserve correct order
    }
    code$rate <- NULL
  } else { # dgamma
    if(is.null(scaleArg)) {
      scaleArg <- substitute(1.0/(A), list(A = code$rate)) 
      code$rate <- scaleArg
      names(code)[which(names(code) == 'rate')] <- 'scale'  # to preserve correct order
    }
    code$rate <- NULL
  }
  return(code)
}

handleScaleAndRateForExpNimble <- function(code){
    scaleArg <- code$scale
    rateArg <- code$rate
    if(is.null(scaleArg) && is.null(rateArg))	stop('neither scale nor rate defined in dexp_nimble')
    if(is.null(rateArg)) {
        rateArg <- substitute(1.0/(A), list(A = code$scale)) 
        code$scale <- rateArg
        names(code)[which(names(code) == 'scale')] <- 'rate'  # to preserve correct order
    }
    code$scale <- NULL
    return(code)
}

## wrtNodes are nodes with respect to which derivatives will be taken
##     (i.e. denominator of dy/dx).  wrtNodes may or may not be part of
##     calcNodes
## extraInputNodes are nodes needed in calculation of derivatives but
##     are not wrtNodes.  They include two categories: (i) parents of
##     any calcNodes that are not themselves calcNodes or wrtNodes;
##     (ii) any stochastic calcNodes, because the node values are needed.
##     However, stochastic calcNodes that are in constantNodes are not
##     included in extraInputNodes.
## constantNodes are nodes that are assumed to be constant for all
##     derivative calls throughout the life of the nimbleFunction.
##     They will be baked in to the CppAD tape and cannot then be changed.
##     They will typically be model data.
## modelOutputNodes are nodes whose values are calculated as part of the
##     tape "forward-zero" stage (value calculation) and need to be stored
##     in the model.  It appears more efficient to copy them in to the model
##     than to use the regular model$calculate() itself.
##     modelOutputNodes will include all deterministic nodes in calcNodes
##     as well as the logProb_ node for any stochastic nodes in calcNodes.
## calcNodes is the same as calcNodes for model$calculate(calcNodes).  It is
##     the ordered sequence of nodes to be calculated
nimDerivsInfoClass_init_impl <- function(.self
                                       , wrtNodes
                                       , calcNodes
                                       , model) {
    .self$model <- model

    ## wrt nodes
    wrtNodes <- model$expandNodeNames(wrtNodes, returnScalarComponents = TRUE)
    wrtNodesAccessor <- modelVariableAccessorVector(model,
                                               wrtNodes,
                                               logProb = FALSE)
    .self$wrtMapInfo <- makeMapInfoFromAccessorVectorFaster(wrtNodesAccessor)

    # See comment in makeDerivsInfo for explanation of why both of the next
    # two lines are necessary.
    calcNodes <- model$expandNodeNames(calcNodes)
    calcNodes <- model$expandNodeNames(calcNodes, returnScalarComponents = TRUE)
    derivsInfo <- makeDerivsInfo_impl(model,
                                      wrtNodes,
                                      calcNodes,
                                      dataAsConstantNodes = TRUE)
    constantNodes <- derivsInfo$constantNodes
    updateNodes <- derivsInfo$updateNodes
    
    extraInputNodesAccessor <- modelVariableAccessorVector(model,
                                                           updateNodes,
                                                           logProb = FALSE)
    .self$extraInputMapInfo <-
        makeMapInfoFromAccessorVectorFaster(extraInputNodesAccessor)

    constantNodesAccessor <- modelVariableAccessorVector(model,
                                                         constantNodes,
                                                         logProb = FALSE)
    .self$constantMapInfo <- makeMapInfoFromAccessorVectorFaster(constantNodesAccessor)

    
    ## output nodes: deterministic nodes in calcNodes plus logProb nodes
    ##   but not the actual data nodes.

    modelOutputNodes <- makeOutputNodes(model, calcNodes)
    
    modelOutputNodesAccessor <- modelVariableAccessorVector(model,
                                                            modelOutputNodes,
                                                            logProb = FALSE)
    .self$modelOutputMapInfo <-
        makeMapInfoFromAccessorVectorFaster(modelOutputNodesAccessor)
    NULL
}

makeOutputNodes <- function(model,
                            calcNodes) {
  calcNodeNames <- model$expandNodeNames(calcNodes, returnScalarComponents = TRUE)
  logProbCalcNodeNames <- model$modelDef$nodeName2LogProbName(calcNodeNames)
  isDetermCalcNodes <- model$isDeterm(calcNodeNames)
  modelOutputNodes <- c(calcNodeNames[isDetermCalcNodes],
                        logProbCalcNodeNames)
  modelOutputNodes
}

nimDerivsInfoClass_output_init_impl <- function(.self,
                                                calcNodes,
                                                model) {
  .self$model <- model
  wrtNodesAccessor <- modelVariableAccessorVector(model,
                                                  character(),
                                                  logProb = FALSE)
  .self$wrtMapInfo <- makeMapInfoFromAccessorVectorFaster(wrtNodesAccessor)
  
  constantNodesAccessor <- modelVariableAccessorVector(model,
                                                       character(),
                                                       logProb = FALSE)
  .self$constantMapInfo <- makeMapInfoFromAccessorVectorFaster(constantNodesAccessor)

  extraInputNodesAccessor <- modelVariableAccessorVector(model,
                                                         character(),
                                                         logProb = FALSE)
  .self$extraInputMapInfo <-
    makeMapInfoFromAccessorVectorFaster(extraInputNodesAccessor)

  modelOutputNodes <- makeOutputNodes(model, calcNodes)
  modelOutputNodesAccessor <- modelVariableAccessorVector(model,
                                                          modelOutputNodes,
                                                          logProb = FALSE)
  .self$modelOutputMapInfo <-
    makeMapInfoFromAccessorVectorFaster(modelOutputNodesAccessor)

  NULL
}

makeDerivsInfo <- function(model,
                           wrtNodes,
                           calcNodes,
                           dataAsConstantNodes = TRUE) {
  wrtNodes <- model$expandNodeNames(wrtNodes, returnScalarComponents = TRUE)
  # This ensures that elements of a non-scalar node become the entire non-scalare node
  calcNodes <- model$expandNodeNames(calcNodes)
  # And then this splits into scalar components.
  # E.g. if calcNodes is 'mu[1]' but 'mu[1:3]' is a vector node,
  # the above call gets `mu[1:3]` and then the below call splits it.
  calcNodes <- model$expandNodeNames(calcNodes, returnScalarComponents = TRUE) 
  makeDerivsInfo_impl(model,
                      wrtNodes,
                      calcNodes,
                      dataAsConstantNodes)
}

getImmediateParentNodes <- function(nodes, model) { 
  ## adapted from BUGS_modelDef creation of edgesFrom2To
  maps <- model$modelDef$maps
  maxNodeID <- length(maps$vertexID_2_nodeID) ## should be same as length(maps$nodeNames)
  
  edgesLevels <- if(maxNodeID > 0) 1:maxNodeID else numeric(0)
  fedgesTo <- factor(maps$edgesTo, levels = edgesLevels) ## setting levels ensures blanks inserted into the splits correctly
  edgesTo2From <<- split(maps$edgesFrom, fedgesTo)
  nodeIDs <- model$expandNodeNames(nodes, returnType = "ids")
  fromIDs <- sort(unique(unlist(edgesTo2From[nodeIDs])))
  fromNodes <- maps$graphID_2_nodeName[fromIDs]
  fromNodes
}

makeDerivsInfo_impl <- function(model,
                                wrtNodes,
                                 calcNodes,
                                dataAsConstantNodes = TRUE) {
  nonWrtCalcNodes <- setdiff(calcNodes, wrtNodes)
  nonWrtCalcNodeNames <- model$expandNodeNames(nonWrtCalcNodes, returnScalarComponents = TRUE)
  nonWrtStochCalcNodeNames <- nonWrtCalcNodeNames[ model$isStoch(nonWrtCalcNodeNames) ]  
  ## Some duplication of work in expandNodeNames
  parentNodes <- getImmediateParentNodes(calcNodes, model)
  parentNodes <- model$expandNodeNames(parentNodes, returnScalarComponents = TRUE)
  neededParentNodes <- setdiff(parentNodes, c(wrtNodes, nonWrtCalcNodeNames))
  extraInputNodes <- model$expandNodeNames(c(neededParentNodes,
                                             nonWrtStochCalcNodeNames),
                                           returnScalarComponents = TRUE,
                                           sort = TRUE)
  constantNodes <- character()
  if(dataAsConstantNodes) {
    boolData <- model$isData(extraInputNodes)
    constantNodes <- extraInputNodes[boolData]
    extraInputNodes <- extraInputNodes[!boolData]
  }
  list(updateNodes = extraInputNodes,
       constantNodes = constantNodes)
}

nimDerivsInfoClass_update_init_impl <- function(.self,
                                                updateNodes = NULL,
                                                constantNodes = NULL,
                                                model) {  
  if(is.null(updateNodes)) updateNodes <- character()
  if(is.null(constantNodes)) constantNodes <- character()

  .self$model <- model
  wrtNodesAccessor <- modelVariableAccessorVector(model,
                                                  character(),
                                                  logProb = FALSE)
  .self$wrtMapInfo <- makeMapInfoFromAccessorVectorFaster(wrtNodesAccessor)
  
  constantNodesAccessor <- modelVariableAccessorVector(model,
                                                       constantNodes,
                                                       logProb = FALSE)
  .self$constantMapInfo <- makeMapInfoFromAccessorVectorFaster(constantNodesAccessor)

  extraInputNodesAccessor <- modelVariableAccessorVector(model,
                                                         updateNodes,
                                                         logProb = FALSE)
  .self$extraInputMapInfo <-
    makeMapInfoFromAccessorVectorFaster(extraInputNodesAccessor)
  
  modelOutputNodesAccessor <- modelVariableAccessorVector(model,
                                                          character(),
                                                          logProb = FALSE)
  .self$modelOutputMapInfo <-
    makeMapInfoFromAccessorVectorFaster(modelOutputNodesAccessor)

  NULL
}

nimDerivsInfoClass <- setRefClass(
    'nimDerivsInfoClass',
    fields = list(
        wrtMapInfo = 'ANY'
      , extraInputMapInfo = 'ANY'
      , modelOutputMapInfo = 'ANY'
      , constantMapInfo = 'ANY'
      , model = 'ANY'
    ),
    methods = list(
        initialize = function(wrtNodes = NULL,
                              calcNodes = NULL,
                              updateNodes = NULL,
                              constantNodes = NULL,
                              thisModel = NULL,
                              case = "nimDerivsCalculate",
                              ...) {
          switch(case,
                 nimDerivsCalculate = nimDerivsInfoClass_init_impl(.self = .self,
                                                                   wrtNodes = wrtNodes,
                                                                   calcNodes = calcNodes,
                                                                   model = thisModel),
                 
                 updateOnly = nimDerivsInfoClass_update_init_impl(.self = .self,
                                                                  updateNodes = updateNodes,
                                                                  constantNodes = constantNodes,
                                                                  model = thisModel),
                 
                 outputOnly = nimDerivsInfoClass_output_init_impl(.self = .self,
                                                                  calcNodes = calcNodes,
                                                                  model = thisModel)
                 )
        }
    )
)

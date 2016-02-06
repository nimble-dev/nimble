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

#		Current objects:
#		d_gamma_keywordInfo
#		pq_gamma_keywordInfo
#		rgamma_keywordInfo
#		d_dist_keywordInfo
#		qp_dist_keywordInfo
#		nimOptim_keywordInfo
#		values_keywordInfo
#		calculate_keywordInfo
#		simulate_keywordInfo
#		getLogProb_keywordInfo
#		nimCopy_keywordInfo
#		doubleBracket_keywordInfo
#		dollarSign_keywordInfo
#		singleBracket_keywordInfo
		
		
		
d_gamma_keywordInfo <- keywordInfoClass(
	keyword = 'dgamma',
	processor = function(code, nfProc){
		logArg <- code$log
		if(logArg == TRUE)	code$log <- 1
			else code$log <- 0
		code <- handleScaleAndRateForGamma(code)
	return(code)
	}) 

pq_gamma_keywordInfo <- keywordInfoClass(
	keyword = 'pq_gamma',
	processor = function(code, nfProc){
		lower.tailArg <- code$lower.tail
		if(lower.tailArg == TRUE) code$lower.tail <- 1
			else code$lower.tail <- 0
			
		logArg <- code$log.p
		if(logArg == TRUE)	code$log.p <- 1
			else code$log.p <- 0
		code <- handleScaleAndRateForGamma(code)
	return(code)
})

rgamma_keywordInfo <- keywordInfoClass(
	keyword = 'rgamma',
	processor = function(code, nfProc){
		code <- handleScaleAndRateForGamma(code)
		return(code)
	}
)

d_exp_nimble_keywordInfo <- keywordInfoClass(
	keyword = 'dexp_nimble',
	processor = function(code, nfProc){
		logArg <- code$log
		if(logArg == TRUE)	code$log <- 1
			else code$log <- 0
		code <- handleScaleAndRateForExpNimble(code)
	return(code)
	}) 

pq_exp_nimble_keywordInfo <- keywordInfoClass(
	keyword = 'pq_exp_nimble',
	processor = function(code, nfProc){
		lower.tailArg <- code$lower.tail
		if(lower.tailArg == TRUE) code$lower.tail <- 1
			else code$lower.tail <- 0
			
		logArg <- code$log.p
		if(logArg == TRUE)	code$log.p <- 1
			else code$log.p <- 0
		code <- handleScaleAndRateForExpNimble(code)
	return(code)
})

rexp_nimble_keywordInfo <- keywordInfoClass(
	keyword = 'rexp_nimble',
	processor = function(code, nfProc){
		code <- handleScaleAndRateForExpNimble(code)
		return(code)
	}
)


d_dist_keywordInfo <- keywordInfoClass(
	keyword = 'd',
	processor = function(code, nfProc){
		logArg <- code$log
		if(logArg == TRUE)	code$log <- 1
			else code$log <- 0
			
		return(code)
	}
)

qp_dist_keywordInfo <- keywordInfoClass(	##q and p functions treated the same
	keyword = 'p',
	processor = function(code, nfProc){
		lower.tailArg	<- code$lower.tail
		if(lower.tailArg == TRUE) code$lower.tail <- 1
			else code$lower.tail <- 0
			
		logArg <- code$log.p
		if(logArg == TRUE)	code$log.p <- 1
			else code$log.p <- 0
		return(code)		
	}
	)

nimOptim_keywordInfo <- keywordInfoClass(
	keyword = 'nimOptim',
	processor = function(code, nfProc){
		
		
		rawFunName <- as.character(code$optFun)		#grabs the run code name of the function to be optimized
		optimFunSym <- nfProc$setupSymTab$symbols[[rawFunName]]		#gets the symbol associated with this function

		optimFunName <- funName2OptimFunName(optimFunSym$nfProc$name)		#makes the name of the C++ function
		argInfo <- getArgInfoFromNFSym(optimFunSym)							#extracts the information about the function to be optimized from it's symbol
				
		voidPointerNFName <- makeVoidPointerName_fromObjName(rawFunName)	#Makes the name for the void* that will ultimately be passed to nimble_optim
		voidPointerForNFLine <- paste0(voidPointerNFName, ' <- voidPtr(', rawFunName, ', nimbleFunction)')		#line of run code for making the void pointer
		
		argPointInfo <- makeDSLCallforVoidPtr_fromArgInfAndCall(code, argInfo)			#Makes the lines of code that create the void*'s for the arguments to be passed to optimized function
		argPointerRunCode <- argPointInfo$newRunCode									#run code for void pointers to arguments
		argVPtrNames <- as.character(argPointInfo$pointerNames)							#names of pointers
								
		allNewBuildCode <- c(voidPointerForNFLine, argPointerRunCode)					#combining all the code for wrapping the void pointers
		code$optFun <- parse(text = optimFunName)[[1]]									#replacing the optFun argument with the name of the new wrapper
		
		
		cNimCall_to_use <- parse(text = 'bareBonesOptim')[[1]]	#As we update the flexbility of optim, we are going to be
																#changing the cNimCall_to_use
																#this is currently in place to bypass issues such as 
																#getting the optimAns and optimControl into the DSL
																
																#When we get those working, the C++ function we will use (i.e. no longer bareBonesOptim) 
																#will need different arguments, and so the next few lines of code will probably need expanding

		code[[1]] <-cNimCall_to_use	
		code[[4]] <- parse(text = voidPointerNFName)[[1]]
		names(code)[4] <- 'voidNimFunPtr'
		code[[5]] <- parse(text = length(argPointInfo$pointerNames))[[1]]
		names(code)[5] <- 'numAdditionalArgs'
		
		for(argName in names(argPointInfo$pointerNames) )
			code[[argName]] <- parse(text = argPointInfo$pointerNames[[argName]])[[1]]

		newRunCode <- quote({})
		for(i in seq_along(allNewBuildCode))
			newRunCode[[i+1]] <- parse(text = allNewBuildCode[[i]])[[1]]
		newRunCode[[length(newRunCode) + 1]] <- code
				
		
		setupArgList <- list(name = optimFunName, nimbleFunctionName = rawFunName)			
		addNecessarySetupCode(optimFunName, setupArgList, optimReadyFun_setupCodeTemplate, nfProc)
		return(newRunCode)
	}
)
                            
                                    
values_keywordInfo <- keywordInfoClass(
    keyword = 'values',
    processor = function(code, nfProc){
      if(!isCodeArgBlank(code, 'accessor'))
      	return(code)
      if(isCodeArgBlank(code, 'model'))
      	stop('model argument missing from values call, with no accessor argument supplied')
      
      accessArgList <- list(model = code$model, nodes = code$nodes, logProb = FALSE)
      accessName <- modelVariableAccessorVector_setupCodeTemplate$makeName(accessArgList)
      addNecessarySetupCode(accessName, accessArgList, modelVariableAccessorVector_setupCodeTemplate, nfProc)
      
      ## accessLengthArgList <- list(accessName = accessName)
      ## accessLengthName <- accessorVectorLength_setupCodeTemplate$makeName(accessLengthArgList)
      ## addNecessarySetupCode(accessLengthName, accessLengthArgList, accessorVectorLength_setupCodeTemplate, nfProc)

      newRunCode <- substitute(values(accessor = ACCESS_NAME), 
                               list(ACCESS_NAME = as.name(accessName)))
      return(newRunCode)
    })                                    

get_param_keywordInfo <- keywordInfoClass(
    keyword = 'get_param',
    processor = function(code, nfProc) {

        if(!isCodeArgBlank(code, 'nodeFunction'))
            return(code)
        if(isCodeArgBlank(code, 'model'))
            stop('model argument missing from get_param, with no accessor argument supplied')
        if(isCodeArgBlank(code, 'node'))
            stop('node argument missing from get_param, with no accessor argument supplied')
        nodeFunVec_ArgList <- list(model = code$model, nodes = code$node, includeData = FALSE)
        nodeFunName <- nodeFunctionVector_SetupTemplate$makeName(nodeFunVec_ArgList)

        if(isCodeArgBlank(code, 'param'))
            stop('param argument missing from get_param, with no accessor argument supplied')
        paramInfo_ArgList <- list(model = code$model, node = code$node, param = code$param)
        paramInfoName <- paramInfo_SetupTemplate$makeName(paramInfo_ArgList)
        paramIDname <- paramInfo_SetupTemplate$makeOtherNames(paramInfoName, paramInfo_ArgList)

        addNecessarySetupCode(nodeFunName, nodeFunVec_ArgList, nodeFunctionVector_SetupTemplate, nfProc)
        addNecessarySetupCode(paramInfoName, paramInfo_ArgList, paramInfo_SetupTemplate, nfProc)
        
        newRunCode <- substitute(get_param(nodeFunction = NODEFUNVEC_NAME, paramID = PARAMID_NAME, paramInfo = PARAMINFO_NAME),
                                 list(NODEFUNVEC_NAME = as.name(nodeFunName), PARAMID_NAME = as.name(paramIDname), PARAMINFO_NAME = as.name(paramInfoName)))
        return(newRunCode)
    }
)

calculate_keywordInfo <- keywordInfoClass(
	keyword = 'calculate',
	processor = function(code, nfProc){
		if(!isCodeArgBlank(code, 'nodeFxnVector'))
			return(code)
		nodeFunVec_ArgList <- list(model = code$model, nodes = code$nodes, includeData = TRUE)
		if(isCodeArgBlank(code, 'model'))
			stop('model argument missing from calculate, with no accessor argument supplied')
		if(isCodeArgBlank(code, 'nodes')){
			LHSnodes_ArgList <- list(model = code$model)
			LHSnodes_name <- allLHSNodes_SetupTemplate$makeName(LHSnodes_ArgList)
			addNecessarySetupCode(LHSnodes_name, LHSnodes_ArgList, allLHSNodes_SetupTemplate, nfProc, allowToCpp = FALSE)
			nodeFunVec_ArgList$nodes = as.name(LHSnodes_name)
			}
		nodeFunName <- nodeFunctionVector_SetupTemplate$makeName(nodeFunVec_ArgList)	
		addNecessarySetupCode(nodeFunName, nodeFunVec_ArgList, nodeFunctionVector_SetupTemplate, nfProc)
		newRunCode <- substitute(calculate(nodeFunctionVector = NODEFUNVEC_NAME),
											list(NODEFUNVEC_NAME = as.name(nodeFunName)))
		return(newRunCode)	
		}
	)

calculateDiff_keywordInfo <- keywordInfoClass(
	keyword = 'calculateDiff',
	processor = function(code, nfProc){
		if(!isCodeArgBlank(code, 'nodeFxnVector'))
			return(code)
		nodeFunVec_ArgList <- list(model = code$model, nodes = code$nodes, includeData = TRUE)
		if(isCodeArgBlank(code, 'model'))
			stop('model argument missing from calculateDiff, with no accessor argument supplied')
		if(isCodeArgBlank(code, 'nodes')){
			LHSnodes_ArgList <- list(model = code$model)
			LHSnodes_name <- allLHSNodes_SetupTemplate$makeName(LHSnodes_ArgList)
			addNecessarySetupCode(LHSnodes_name, LHSnodes_ArgList, allLHSNodes_SetupTemplate, nfProc, allowToCpp = FALSE)
			nodeFunVec_ArgList$nodes = as.name(LHSnodes_name)
			}
		nodeFunName <- nodeFunctionVector_SetupTemplate$makeName(nodeFunVec_ArgList)	
		addNecessarySetupCode(nodeFunName, nodeFunVec_ArgList, nodeFunctionVector_SetupTemplate, nfProc)
		newRunCode <- substitute(calculateDiff(nodeFunctionVector = NODEFUNVEC_NAME),
											list(NODEFUNVEC_NAME = as.name(nodeFunName)))
		return(newRunCode)	
		}
	)


simulate_keywordInfo <- keywordInfoClass(
	keyword = 'simulate',
	processor = function(code, nfProc){
		if(!isCodeArgBlank(code, 'nodeFxnVector')){
			return(substitute(simulate(nodeFxnVector = NODEFXNVECTOR), list(NODEFXNVECTOR = code$nodeFxnVector) ) )
			}
		nodeFunVec_ArgList <- list(model = code$model, nodes = code$nodes, includeData = code$includeData)
		if(isCodeArgBlank(code, 'model'))
			stop('model argument missing from simulate, with no accessor argument supplied')
		if(isCodeArgBlank(code, 'nodes')){
			LHSnodes_ArgList <- list(model = code$model)
			LHSnodes_name <- allLHSNodes_SetupTemplate$makeName(LHSnodes_ArgList)
			addNecessarySetupCode(LHSnodes_name, LHSnodes_ArgList, allLHSNodes_SetupTemplate, nfProc, allowToCpp = FALSE)
			nodeFunVec_ArgList$nodes = as.name(LHSnodes_name)
			}
		nodeFunName <- nodeFunctionVector_SetupTemplate$makeName(nodeFunVec_ArgList)	
		addNecessarySetupCode(nodeFunName, nodeFunVec_ArgList, nodeFunctionVector_SetupTemplate, nfProc)
		newRunCode <- substitute(simulate(nodeFunctionVector = NODEFUNVEC_NAME),
											list(NODEFUNVEC_NAME = as.name(nodeFunName)))
											
		return(newRunCode)	
		}
	)

getLogProb_keywordInfo <- keywordInfoClass(
	keyword = 'getLogProb',
	processor = function(code, nfProc){
		if(!isCodeArgBlank(code, 'nodeFxnVector'))
			return(code)
		nodeFunVec_ArgList <- list(model = code$model, nodes = code$nodes, includeData = TRUE)
		if(isCodeArgBlank(code, 'model'))
			stop('model argument missing from getLogProb, with no accessor argument supplied')
		if(isCodeArgBlank(code, 'nodes')){
			LHSnodes_ArgList <- list(model = code$model)
			LHSnodes_name <- allLHSNodes_SetupTemplate$makeName(LHSnodes_ArgList)
			addNecessarySetupCode(LHSnodes_name, LHSnodes_ArgList, allLHSNodes_SetupTemplate, nfProc, allowToCpp = FALSE)
			nodeFunVec_ArgList$nodes = as.name(LHSnodes_name)
			}
		nodeFunName <- nodeFunctionVector_SetupTemplate$makeName(nodeFunVec_ArgList)	
		addNecessarySetupCode(nodeFunName, nodeFunVec_ArgList, nodeFunctionVector_SetupTemplate, nfProc)
		newRunCode <- substitute(getLogProb(nodeFunctionVector = NODEFUNVEC_NAME),
											list(NODEFUNVEC_NAME = as.name(nodeFunName)))
		return(newRunCode)	
		}
	)    

nimCopy_keywordInfo <- keywordInfoClass(
	keyword = 'nimCopy',
    processor = function(code, nfProc){
        if(is.null(nfProc)) stop("Can\'t call copy (nimCopy) from a nimbleFunction without setup code")
		possibleObjects <- c('symbolModel', 'symbolModelValues', 'symbolModelVariableAccessorVector', 'symbolModelValuesAccessorVector')
		modelValuesTypes <- c('symbolModelValues', 'symbolModelValuesAccessorVector')
		accessTypes <- c('symbolModelVariableAccessorVector', 'symbolModelValuesAccessorVector')
		from_ArgList <- list(name = code$from, class = symTypeFromSymTab(code$from, nfProc$setupSymTab, options = possibleObjects))
		to_ArgList <- list(name = code$to, class = symTypeFromSymTab(code$to, nfProc$setupSymTab, options = possibleObjects))
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
                    isMVfrom <- 0 ## for newNimCopy
                        accessFrom_ArgList <- list(model = code$from, nodes = from_ArgList$nodes, logProb = code$logProb)
			accessFrom_name <- modelVariableAccessorVector_setupCodeTemplate$makeName(accessFrom_ArgList)
			addNecessarySetupCode(accessFrom_name, accessFrom_ArgList, modelVariableAccessorVector_setupCodeTemplate, nfProc)
		}
		else if(from_ArgList$class == 'symbolModelValues'){
                    isMVfrom <- 1 ## for newNimCopy
                        accessFrom_ArgList <- list(modelValues = code$from, nodes = from_ArgList$nodes, logProb = code$logProb, row = from_ArgList$row)
			accessFrom_name <- modelValuesAccessorVector_setupCodeTemplate$makeName(accessFrom_ArgList)
			addNecessarySetupCode(accessFrom_name, accessFrom_ArgList, modelValuesAccessorVector_setupCodeTemplate, nfProc)
		}
		else if(from_ArgList$class %in% accessTypes) {
                    isMVfrom <- as.integer(from_ArgList$class == 'symbolModelValuesAccessorVector') 
                    accessFrom_name <- as.character(code$from)
                }
        
		if(to_ArgList$class == 'symbolModel'){
                    isMVto <- 0 ## for newNimCopy
                        accessTo_ArgList <- list(model = code$to, nodes = to_ArgList$nodes, logProb = code$logProb)
			accessTo_name <- modelVariableAccessorVector_setupCodeTemplate$makeName(accessTo_ArgList)
			addNecessarySetupCode(accessTo_name, accessTo_ArgList, modelVariableAccessorVector_setupCodeTemplate, nfProc)
		}
		else if(to_ArgList$class == 'symbolModelValues'){
                    isMVto <- 1 ## for newNimCopy
                        accessTo_ArgList <- list(modelValues = code$to, nodes = to_ArgList$nodes, logProb = code$logProb, row = to_ArgList$row)
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

#	Need to get setupCodeTemplates working first...
doubleBracket_keywordInfo <- keywordInfoClass(
	keyword = '[[', 
    processor = function(code, nfProc){
        if(is.null(nfProc)) stop("No allowed use of [[ in a nimbleFunction without setup code.")
        possibleObjects <- c('symbolModel', 'symbolNimPtrList', 'symbolNimbleFunctionList')
        class = symTypeFromSymTab(code[[2]], nfProc$setupSymTab, options = possibleObjects)
        if(class == 'symbolNimPtrList' || class == 'symbolNimbleFunctionList')
            return(code)
        if(class == 'symbolModel'){
            singleAccess_ArgList <- list(code = code, model = code[[2]], nodeExpr = code[[3]])
            nodeArg <- code[[3]]
            if(is.character(nodeArg)){
                varAndIndices <- getVarAndIndices(nodeArg)
                nDim <- sum(1 - unlist(lapply(varAndIndices$indices, is.numeric) ) )
                useMap <- nDim > 0
            }
            else{
                allNDims <- determineNdimsFromNfproc(singleAccess_ArgList$model, nodeArg, nfProc)
                if(length(unique(allNDims)) > 1) stop(paste0('Error for ', deparse(code), '. Inconsistent numbers of dimensions for different instances.'))
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
        stop(paste('in keywordProcessing of "[[", type not recognized. Code = ', code) )
    })

dollarSign_keywordInfo <- keywordInfoClass(
	keyword = '$',
	processor = function(code, nfProc){
            if(is.null(nfProc)) stop("No legal use of dollar sign in nimbleFunction with no setup code")
			#	First thing we need to do is remove 'run', for backward compatibility, i.e.
			#   replace myNimbleFunction$run() -> myNimbleFunction() or
			#	myNimbleFunList[[i]]$run() -> myNimbleFunList[[i]]()
			#   Probably a better way to handle this
			
					
		if(code[[3]] == 'run'){
			newRunCode <- code[[2]]
			return(newRunCode)
		}
				
		possibleObjects <- c('symbolModel', 'symbolNimPtrList', 'symbolNimbleFunction', 'symbolNimbleFunctionList')

		callerCode <- code[[2]]
		#	This extracts myNimbleFunction from the expression myNimbleFunction$foo()
		
		
		if(length(callerCode) > 1)
			callerCode <- callerCode[[2]]
		#	This extracts myNimbleFunctionList from the expression myNimbleFunctionList[[i]]
		#	May be a better way to do this
		
		class <- symTypeFromSymTab(callerCode, nfProc$setupSymTab, options = possibleObjects)
				
		if(class == 'symbolNimPtrList'){
			return(code)
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
						
			nf_charName <- as.character(callerCode)
			nf_fieldName <-as.character(code[[3]])
			objectSymbol = nfProc$setupSymTab$symbols[[nf_charName]]$nfProc$setupSymTab$symbols[[nf_fieldName]]
			if(class(objectSymbol)[[1]] == 'symbolMemberFunction'){
				newRunCode <- substitute(nfMethod(NIMBLEFXN, METHODNAME), list(NIMBLEFXN = as.name(nf_charName), METHODNAME = nf_fieldName))				
				return(newRunCode)
			}
			else {
				# I *assume* that if its not a member function, it should be treated with 
				# nfVar
				newRunCode <- substitute(nfVar(NIMBLEFXN, VARNAME), list(NIMBLEFXN = as.name(nf_charName), VARNAME = nf_fieldName))
				return(newRunCode)
			}
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
			stop(paste('in keywordProcessing of "$", type not recognized. Code = ', code) )
	}
)
    
singleBracket_keywordInfo <- keywordInfoClass(
	keyword = '[',
    processor = function(code, nfProc){
        if(is.null(nfProc)) return (code)
        class <- symTypeFromSymTab(code[[2]], nfProc$setupSymTab)
        if(class == 'symbolModelValues'){
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



#	KeywordList
keywordList <- new.env()
keywordList[['get_param']] <- get_param_keywordInfo
keywordList[['values']] <- values_keywordInfo
keywordList[['calculate']] <- calculate_keywordInfo
keywordList[['calculateDiff']] <- calculateDiff_keywordInfo
keywordList[['simulate']] <- simulate_keywordInfo
keywordList[['getLogProb']] <- getLogProb_keywordInfo
keywordList[['nimCopy']] <- nimCopy_keywordInfo
keywordList[['[[']] <- doubleBracket_keywordInfo
keywordList[['$']] <- dollarSign_keywordInfo
keywordList[['[']] <- singleBracket_keywordInfo
keywordList[['nimOptim']] <- nimOptim_keywordInfo
keywordList[['dgamma']] <- d_gamma_keywordInfo
keywordList[['pgamma']] <- pq_gamma_keywordInfo
keywordList[['qgamma']] <- pq_gamma_keywordInfo
keywordList[['rgamma']] <- rgamma_keywordInfo
# can be handled the same as gamma, so include although we have dexp_nimble too
keywordList[['dexp']] <- d_gamma_keywordInfo
keywordList[['pexp']] <- pq_gamma_keywordInfo
keywordList[['qexp']] <- pq_gamma_keywordInfo
keywordList[['rexp']] <- rgamma_keywordInfo

keywordList[['dexp_nimble']] <- d_exp_nimble_keywordInfo
keywordList[['pexp_nimble']] <- pq_exp_nimble_keywordInfo
keywordList[['qexp_nimble']] <- pq_exp_nimble_keywordInfo
keywordList[['rexp_nimble']] <- rexp_nimble_keywordInfo

# necessary keywords:
#	calculate 	(done)
#	simulate	(done)
#	getLogProb	(done)
#	values		(done)
#	getValues	(removed from DSL)
#	setValues	(removed from DSL)
#	nimCopy		(done)
#	[[			(done)
#	$			(done)
#	[			(done)
#	resize		(special processing only for numericLists, not really used)
#	setSize		(special processing only for numericLists, not really used)
#	getSize		(special processing only for numericLists, not really used)
#	Also see replaceAccessorsOneFunction



matchFunctions <- new.env()
matchFunctions[['values']] <- function(model, nodes, accessor){}
matchFunctions[['get_param']] <- get_param
matchFunctions[['calculate']] <- calculate		#function(model, nodes, nodeFunctionVector){}
matchFunctions[['calculateDiff']] <- calculateDiff		#function(model, nodes, nodeFunctionVector){}
matchFunctions[['simulate']] <- simulate		#function(model, nodes, includeData = FALSE, nodeFunctionVector){}
matchFunctions[['getLogProb']] <- getLogProb	#function(model, nodes, nodeFunctionVector){}
matchFunctions[['nimCopy']] <- function(from, to, nodes, nodesTo, row, rowTo, logProb = FALSE){}
matchFunctions[['double']] <- function(dim, default){}
matchFunctions[['int']] <- function(dim, default){}
matchFunctions[['nimOptim']] <- function(initPar, optFun, ...){} 
matchFunctions[['dgamma']] <- function(x, shape, rate = 1, scale, log = FALSE){}
matchFunctions[['rgamma']] <- function(n, shape, rate = 1, scale){}
matchFunctions[['qgamma']] <- function(p, shape, rate = 1, scale, lower.tail = TRUE, log.p = FALSE){}
matchFunctions[['pgamma']] <- function(q, shape, rate = 1, scale, lower.tail = TRUE, log.p = FALSE){}
matchFunctions[['dexp']] <- function(x, rate = 1, log = FALSE){}
matchFunctions[['rexp']] <- function(n, rate = 1){}
matchFunctions[['qexp']] <- function(p, rate = 1, lower.tail = TRUE, log.p = FALSE){}
matchFunctions[['pexp']] <- function(q, rate = 1, lower.tail = TRUE, log.p = FALSE){}
matchFunctions[['dexp_nimble']] <- function(x, rate, scale = 1, log = FALSE){}
matchFunctions[['rexp_nimble']] <- function(n, rate, scale = 1){}
matchFunctions[['qexp_nimble']] <- function(p, rate, scale = 1, lower.tail = TRUE, log.p = FALSE){}
matchFunctions[['pexp_nimble']] <- function(q, rate, scale = 1, lower.tail = TRUE, log.p = FALSE){}

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
matchDistList <- list('binom', 'cat', 'dirch', 'interval', 'lnorm', 'logis', 'multi', 'mnorm_chol', 'norm', 'pois', 't_nonstandard', 'unif', 'weibull', 'wish_chol')
# these are standard for keywordList and handled specially above for matchFunctions
keywordOnlyMatchDistList <- list('t', 'beta', 'chisq', 'nbinom')

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

addDistKeywordProcessors <- function(distList, keywordEnv){
		for(thisDist in distList) {
                    pFun <- paste0('p', thisDist)
                    qFun <- paste0('q', thisDist)
                    dFun <- paste0('d', thisDist)
                    
                    keywordEnv[[dFun]] <- d_dist_keywordInfo
                    keywordEnv[[pFun]] <- qp_dist_keywordInfo
                    keywordEnv[[qFun]] <- qp_dist_keywordInfo
		}
            }
          

addDistList2matchFunctions(matchDistList, matchFunctions)
addDistKeywordProcessors(c(matchDistList, keywordOnlyMatchDistList), keywordList)

#	processKeyword function to be called by nfProc
processKeyword <- function(code, nfProc){
  thisKeywordInfo <- keywordList[[ as.character(code[[1]]) ]]
  if(!is.null(thisKeywordInfo))
    return(thisKeywordInfo$processor(code, nfProc))
  return(code)
}







#####	SETUPCODE TEMPLATES

#		Current Available Templates:
#		optimReadyFun_setupCodeTemplate
#		modelVariableAccessorVector_setupCodeTemplate
#		modelValuesAccessorVector_setupCodeTemplate
#		accessorVectorLength_setupCodeTemplate
#		nodeFunctionVector_SetupTemplate
#		allLHSNodes_SetupTemplate
#		allModelNodes_SetupTemplate
#		allModelValuesVars_SetupTemplate
#		singleVarAccess_SetupTemplate
#		singleModelIndexAccess_SetupTemplate
#		map_SetupTemplate
#       singleModelValuesAccessor_SetupTemplate
        
        
        
                                         

optimReadyFun_setupCodeTemplate <- setupCodeTemplateClass(
	makeName = function(argList){Rname2CppName(deparse(argList$name))},
	codeTemplate = quote(OPTIM_FUN <- OptimReadyFunction(name = OPTIM_FUN_INQUOTES, nimbleFunction = NFNAME, localNimbleFunctionName = LOCALORGNAME)),
	makeCodeSubList = function(resultName, argList){
		list(OPTIM_FUN = as.name(argList$name),
			OPTIM_FUN_INQUOTES = argList$name, 
			NFNAME = as.name(argList$nimbleFunctionName),
			LOCALORGNAME = argList$nimbleFunctionName)
	})

                                          
modelVariableAccessorVector_setupCodeTemplate <- setupCodeTemplateClass(
	#Note to programmer: required fields of argList are model, nodes and logProb
    makeName = function(argList) {Rname2CppName(paste(deparse(argList$model), deparse(argList$nodes), 'access_logProb', deparse(argList$logProb), sep = '_'))},
    codeTemplate = quote( ACCESSNAME <- modelVariableAccessorVector(MODEL, NODES, logProb = LOGPROB) ),
    makeCodeSubList = function(resultName, argList) {
        list(ACCESSNAME = as.name(resultName),
             MODEL = argList$model,
             NODES = argList$nodes,
             LOGPROB = argList$logProb)
    })

copierVector_setupCodeTemplate <- setupCodeTemplateClass(
    makeName = function(argList) {Rname2CppName(paste0(argList$accessFrom_name, '_', argList$accessTo_name))},
    codeTemplate = quote( COPIERNAME <- copierVector(ACCESS_FROM, ACCESS_TO, ISMVFROM, ISMVTO) ),
    makeCodeSubList = function(resultName, argList) {
        list(COPIERNAME = as.name(resultName),
             ACCESS_FROM = as.name(argList$accessFrom_name),
             ACCESS_TO   = as.name(argList$accessTo_name),
             ISMVFROM    = as.integer(argList$isMVfrom),
             ISMVTO      = as.integer(argList$isMVto)) 
    })
    

modelValuesAccessorVector_setupCodeTemplate <- setupCodeTemplateClass(
	#Note to programmer: required fields of argList are model, nodes and logProb

    makeName = function(argList) {Rname2CppName(paste(deparse(argList$model), deparse(argList$nodes), 'access_logProb', deparse(argList$logProb), deparse(argList$row), sep = '_'))},
    codeTemplate = quote( ACCESSNAME <- modelValuesAccessorVector(MODEL, NODES, logProb = LOGPROB) ),
	makeCodeSubList = function(resultName, argList) {
        list(ACCESSNAME = as.name(resultName),
             MODEL = argList$model,
             NODES = argList$nodes,
             LOGPROB = argList$logProb)
    })


    
accessorVectorLength_setupCodeTemplate <- setupCodeTemplateClass( ## This is not very nice: modify the accessor to have element 5 as the length name, so that when makeMapInfo...() is called it will set the length variable in the calling environment.  kind of convoluted but doing it for now.
  #Note to programmer: required fields of argList are accessName
 
  makeName = function(argList){ Rname2CppName(paste(deparse(argList$accessName), 'length', sep = '_')) },
    codeTemplate = quote({ACCESSLENGTH <- 0;
        ACCESSNAME[[5]] <- ACCESSLENGTHNAME}),
  makeCodeSubList = function(resultName, argList){
      list(ACCESSNAME = as.name(argList$accessName),
           ACCESSLENGTH = as.name(resultName),
           ACCESSLENGTHNAME = resultName)
  })

## accessorVectorLength_setupCodeTemplate <- setupCodeTemplateClass( ## NEW ACCESSORS
##   #Note to programmer: required fields of argList are accessName 
##   makeName = function(argList){ Rname2CppName(paste(deparse(argList$accessName), 'length', sep = '_')) },
##   codeTemplate = quote(ACCESSLENGTH <- ACCESSNAME$getLength()),
##   makeCodeSubList = function(resultName, argList){
##   	list(ACCESSNAME = as.name(argList$accessName),
##   		ACCESSLENGTH = as.name(resultName) )
##   	})


nodeFunctionVector_SetupTemplate <- setupCodeTemplateClass(
	#Note to programmer: required fields of argList are model, nodes and includeData
	
	makeName = function(argList){Rname2CppName(paste(deparse(argList$model), deparse(argList$nodes), 'nodeFxnVector_includeData', deparse(argList$includeData), sep = '_'))},
	codeTemplate = quote(NODEFXNVECNAME <- nodeFunctionVector(model = MODEL, nodeNames = NODES, excludeData = EXCLUDEDATA)), 
	makeCodeSubList = function(resultName, argList){
		list(NODEFXNVECNAME = as.name(resultName),
			MODEL = argList$model,
			NODES = argList$nodes,
			EXCLUDEDATA = !argList$includeData)
	})

paramInfo_SetupTemplate <- setupCodeTemplateClass(
    #Note to programmer: required fields of argList are model, node and param
    makeName = function(argList){Rname2CppName(paste(deparse(argList$model), deparse(argList$node), deparse(argList$param), 'paramInfo', sep='_'))},
    makeOtherNames = function(name,argList) {Rname2CppName(paste0(name,'_ID'))},
    codeTemplate = quote({
        PARAMINFONAME <- makeParamInfo(MODEL, NODE, PARAM)
        PARAMIDNAME <- PARAMINFONAME$paramID
       }),
    makeCodeSubList = function(resultName, argList){
        list(PARAMINFONAME = as.name(resultName),
             PARAMIDNAME = as.name(paste0(resultName,'_ID')),
             MODEL = argList$model,
             NODE = argList$node,
             PARAM = argList$param)
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

	codeTemplate = quote(SINGLEACCESSOR <- singleVarAccess(MODEL, VAR)),

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
		VARANDINDICES <- getVarAndIndices(NODEVARNAME)
		NEWVARNAME <- as.character(VARANDINDICES$varName)
		MFLATINDEX <- varAndIndices2flatIndex(VARANDINDICES, MODELVAREXPR$getVarInfo(NEWVARNAME))
		VARACCESSOR <- singleVarAccess(MODELVAREXPR, NEWVARNAME, useSingleIndex = TRUE)
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
		VARANDINDICES <- getVarAndIndices(NODEVARNAME)
		NEWVARNAME <- as.character(VARANDINDICES$varName)
                map_SetupTemplate_vi <- MODEL$getVarInfo(NEWVARNAME)
		map_SetupTemplate_mapParts <- varAndIndices2mapParts(VARANDINDICES, map_SetupTemplate_vi$maxs, map_SetupTemplate_vi$nDim)
		MSTRIDES <- map_SetupTemplate_mapParts$strides
		MOFFSET <- map_SetupTemplate_mapParts$offset
		MSIZES <- map_SetupTemplate_mapParts$sizes
		VARACCESSOR <- singleVarAccess(model, NEWVARNAME)
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
		MVACCESS <- singleModelValuesAccess(MODELVALUES, VAR)
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
                                        #	name <- template$makeName(argList)
    if(!isSetupCodeGenerated(name, nfProc)){
        addSetupCodeNames(name, template$makeOtherNames(name, argList), nfProc)
        subList <- template$makeCodeSubList(name, argList)
        addNewCode(name, subList, template, nfProc)
        if(!allowToCpp) addBlockFromCppName(name, nfProc) ## ignores makeOtherNames for now
    }
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
		stop(paste('invalid class for object ', codeName, 'class provided = ', class) )
	return(class)
}

isSymbolType <- function(symTab, varName, symType)
	inherits(symTab$symbols[[varName]], symType)

matchAndFill.call <- function(def, call){
  theseFormals <- formals(def)
  formalNames <- names(theseFormals) # formalArgs are the arguments that are defined, i.e. does NOT include anything that is from the args "..."
  theseFormals <- theseFormals[nchar(theseFormals) > 0]
  matchedCall <- match.call(def, call)
  missingArgs <- which(!(names(theseFormals) %in% names(matchedCall)))
  for(ind in missingArgs){
    name <- names(theseFormals)[ind]
    matchedCall[[name]] <- theseFormals[[name]]    
  }
    
  newCall <- matchedCall[1]

  for(thisArgName in formalNames){					# This is to get the order of the arguments correctly
  	thisArg <- matchedCall[[thisArgName]]
	if(!is.null(thisArg))
	  	newCall[[thisArgName]] <- thisArg
  }
  
  informalArgNames <- names(matchedCall)[!(names(matchedCall) %in% formalNames)]
 		# i.e. are there any "..." args? if so, adds them on in the end
  		# Note: this will preserve arguments EVEN if no '...' is declared, i.e.
  		# dnorm(jnk = 3, x= 10) will turn into dnorm(x = 10, mean = 0, sd = 1, log = FALSE, jnk = 3)
  informalArgNames <- informalArgNames[-1]	#removing "", which is the function call, not an argument 

  for(thisArg in informalArgNames)
  	newCall[[thisArg]] <- matchedCall[[thisArg]]
  return(newCall)
}

## pasteExpr <- function(expr1, expr2)
## 	parse(text=paste0(as.character(expr1), as.character(expr2) ) )[[1]]


determineNdimsFromNfproc <- function(modelExpr, varOrNodeExpr, nfProc) {
    allNDims <- lapply(nfProc$instances, function(x) {
        model <- eval(modelExpr, envir = x)
        if(!exists(as.character(varOrNodeExpr), x, inherits = FALSE) ) {
            stop(paste0('Error, ', as.character(varOrNodeExpr), ' does not exist in an instance of this nimbleFunction.'))
        }
        lab <- eval(varOrNodeExpr, envir = x)
        varAndIndices <- getVarAndIndices(lab)
        determineNdimFromOneCase(model, varAndIndices)
    } )
    return(allNDims)
}

matchKeywordCodeMemberFun <- function(code, nfProc) {  ## handles cases like a$b(c) as one unit so the member function template for b can be looked up
    dollarSignPart <- code[[1]] ## we already checked that code[[1]][[1]] is '$' before calling this function
    nfPart <- dollarSignPart[[2]]
    if(length(nfPart) != 1) { ## It could be a nimbleFunctionList with nfl[[i]]$member(a)
        nfNestedPart <- nfPart[[1]]
        if(length(nfNestedPart) != 1) stop(paste0("Cannot handle this expression: ", deparse(code)))
        if(deparse(nfNestedPart) != '[[') stop(paste0("Cannot handle this expression: ", deparse(code)))
        nfListName <- deparse(nfPart[[2]])
        memFunName <- deparse(dollarSignPart[[3]])
        if(is.null(nfProc)) stop(paste0("Cannot handle what looks like a nimbleFunctionList usage unless it was created in setup code: ", deparse(code)), call. = FALSE)
        
        if(nfProc$setupSymTab$symbolExists(nfListName)) { ## look in symbol table
            symObj <- nfProc$setupSymTab$getSymbolObject(nfListName)
            if(symObj$type == 'nimbleFunctionList') {
                thisBaseClass <- symObj$baseClass
                thisFunctionMatch <- environment(symObj$baseClass)$methodList[[memFunName]]$template
                return(matchAndFill.call(thisFunctionMatch, code ) )
            } else stop(paste0("Syntax looks like a nimbleFunctionList member function, but the object isn\'t the right type: ", deparse(code)), call. = FALSE)
        } else stop(paste0("Syntax looks like a nimbleFunctionList member function, but we can\'t find it in setup: ", deparse(code)), call. = FALSE)
    }
    ## must be nfObj$member(args)
    nfName <- as.character(nfPart) ## nfObj
    memFunName <- deparse(dollarSignPart[[3]])
    if(nfProc$setupSymTab$symbolExists(nfName)) { ## first look in symbolTable
        symObj <- nfProc$setupSymTab$getSymbolObject(nfName)
        if(symObj$type == 'nimbleFunction') {
            thisRCfunProc <- if(memFunName == 'run') symObj$nfProc$RCfunProcs[["operator()"]] else symObj$nfProc$RCfunProcs[[memFunName]] 
            if(is.null(thisRCfunProc)) stop(paste0("Cannot handle this expression (member function may not exist): ", deparse(code)), call. = FALSE)
            thisFunctionMatch <- thisRCfunProc$RCfun$template
            return(matchAndFill.call(thisFunctionMatch, code ) )
        } else stop(paste0("Cannot handle this expression (maybe it's not a nimbleFunction?): ", deparse(code))) 
    }
    ## then look in R
    if(exists(nfName)) {
       stop(paste0("Cannot use a specialized nimbleFunction that is not in setup code (it can be an argument to setup or created in setup): ", deparse(code)), call. = FALSE)
    }
}

matchKeywordCode <- function(code, nfProc){
    callName <- as.character(code[[1]])
    thisFunctionMatch <- matchFunctions[[ callName ]]

    ## see if this is a member function of an nf object
    if(!is.null(nfProc)) {
        modCallName <- if(callName == "run") "operator()" else callName
        if(nfProc$setupSymTab$symbolExists(modCallName)) {
            symObj <- nfProc$setupSymTab$getSymbolObject(modCallName)
            if(class(symObj) == "symbolMemberFunction") {
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
    if(inherits(varInfo, 'try-error')) browser()
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
    ##varName <- varAndIndices$name
    indices <- varAndIndices$indices
    ## put together offsetExpr, sizeExprs, strideExprs
    ## need sizes to get strides
    if(length(sizes) == 0) sizes <- 1
##    sizes <- if(length(varInfo$maxs) > 0) varInfo$maxs else 1 ## would be wierd to be mapping into something with length 1 anyway
##    if(varInfo$nDim > 0 & length(indices)==0) { ## A case like model[[node]] where node == 'x', and we should treat like 'x[,]', e.g.
##        nDim <- varInfo$nDim
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
    ## also handles core R dexp
    scaleArg <- code$scale
    rateArg <- code$rate
    if(is.null(scaleArg) && is.null(rateArg))	stop('neither scale nor rate defined in dgamma or dexp')
    if(is.null(scaleArg)) {
        scaleArg <- substitute(1.0/(A), list(A = code$rate)) ## parse(text = paste0('1/', code$rate))[[1]]
        code$rate <- scaleArg
        names(code)[which(names(code) == 'rate')] <- 'scale'  # to preserve correct order
    }
    code$rate <- NULL
    return(code)
}

handleScaleAndRateForExpNimble <- function(code){
    scaleArg <- code$scale
    rateArg <- code$rate
    if(is.null(scaleArg) && is.null(rateArg))	stop('neither scale nor rate defined in dexp_nimble')
    if(is.null(rateArg)) {
        rateArg <- substitute(1.0/(A), list(A = code$scale)) ##parse(text = paste0('1/', code$scale))[[1]]
        code$scale <- rateArg
        names(code)[which(names(code) == 'scale')] <- 'rate'  # to preserve correct order
    }
    code$scale <- NULL
    return(code)
}

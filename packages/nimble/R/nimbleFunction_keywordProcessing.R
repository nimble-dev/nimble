######	KEYWORD PROCESSING OBJECTS

#Keyword Class and objects






keywordInfoClass <- setRefClass('keywordInfoClass',
                                fields = list(
                                    keyword = 'ANY',
                                    processor = 'ANY'))
                                    
                                    
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
      	
      accessLengthArgList <- list(accessName = accessName)
      accessLengthName <- accessorVectorLength_setupCodeTemplate$makeName(accessLengthArgList)
      addNecessarySetupCode(accessLengthName, accessLengthArgList, accessorVectorLength_setupCodeTemplate, nfProc)

      newRunCode <- substitute(values(accessor = ACCESS_NAME), 
                               list(ACCESS_NAME = as.name(accessName)))
      return(newRunCode)
    })                                    
 
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
			addNecessarySetupCode(LHSnodes_name, LHSnodes_ArgList, allLHSNodes_SetupTemplate, nfProc)
			nodeFunVec_ArgList$nodes = as.name(LHSnodes_name)
			}
		nodeFunName <- nodeFunctionVector_SetupTemplate$makeName(nodeFunVec_ArgList)	
		addNecessarySetupCode(nodeFunName, nodeFunVec_ArgList, nodeFunctionVector_SetupTemplate, nfProc)
		newRunCode <- substitute(calculate(nodeFunctionVector = NODEFUNVEC_NAME),
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
			addNecessarySetupCode(LHSnodes_name, LHSnodes_ArgList, allLHSNodes_SetupTemplate, nfProc)
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
			addNecessarySetupCode(LHSnodes_name, LHSnodes_ArgList, allLHSNodes_SetupTemplate, nfProc)
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
				addNecessarySetupCode(allNodes_name, node_ArgList, allModelNodes_SetupTemplate, nfProc)
			}
			else if(from_ArgList$class == 'symbolModelValues'){
				from_ArgList$row = code$row
				mvVar_ArgList <- list(modelValues = from_ArgList$name)
				allNodes_name <- allModelValuesVars_SetupTemplate$makeName(mvVar_ArgList)
				addNecessarySetupCode(allNodes_name, mvVar_ArgList, allModelValuesVars_SetupTemplate, nfProc)
			}
			from_ArgList$nodes <- as.name(allNodes_name)
		}
		else	from_ArgList$nodes <- code$nodes
		
		if(isCodeArgBlank(code, 'nodesTo'))		to_ArgList$nodes <- from_ArgList$nodes
		else									to_ArgList$nodes <- code$nodesTo
				
		if(from_ArgList$class == 'symbolModel'){
			accessFrom_ArgList <- list(model = code$from, nodes = from_ArgList$nodes, logProb = code$logProb)
			accessFrom_name <- modelVariableAccessorVector_setupCodeTemplate$makeName(accessFrom_ArgList)
			addNecessarySetupCode(accessFrom_name, accessFrom_ArgList, modelVariableAccessorVector_setupCodeTemplate, nfProc)
		}
		else if(from_ArgList$class == 'symbolModelValues'){
			accessFrom_ArgList <- list(modelValues = code$from, nodes = from_ArgList$nodes, logProb = code$logProb)
			accessFrom_name <- modelValuesAccessorVector_setupCodeTemplate$makeName(accessFrom_ArgList)
			addNecessarySetupCode(accessFrom_name, accessFrom_ArgList, modelValuesAccessorVector_setupCodeTemplate, nfProc)
		}
		else if(from_ArgList$class %in% accessTypes)
			accessFrom_name <- as.character(code$from)
		
		if(to_ArgList$class == 'symbolModel'){
			accessTo_ArgList <- list(model = code$to, nodes = to_ArgList$nodes, logProb = code$logProb)
			accessTo_name <- modelVariableAccessorVector_setupCodeTemplate$makeName(accessTo_ArgList)
			addNecessarySetupCode(accessTo_name, accessTo_ArgList, modelVariableAccessorVector_setupCodeTemplate, nfProc)
		}
		else if(to_ArgList$class == 'symbolModelValues'){
			accessTo_ArgList <- list(modelValues = code$to, nodes = to_ArgList$nodes, logProb = code$logProb)
			accessTo_name <- modelValuesAccessorVector_setupCodeTemplate$makeName(accessTo_ArgList)
			addNecessarySetupCode(accessTo_name, accessTo_ArgList, modelValuesAccessorVector_setupCodeTemplate, nfProc)
		}
		else if(to_ArgList$class %in% accessTypes)
			accessTo_name <- as.character(code$to)
			
		#What happens below is a bit convoluted and really for backwards compatibility 	
		runCode <- substitute(nimCopy(from = FROM_ACCESS, rowFrom = NA, to = TO_ACCESS, rowTo = NA), 
							  list(FROM_ACCESS = as.name(accessFrom_name), TO_ACCESS = as.name(accessTo_name)))
		if(from_ArgList$class %in% modelValuesTypes)
			runCode$rowFrom = from_ArgList$row
		if(to_ArgList$class %in% modelValuesTypes)
			runCode$rowTo = to_ArgList$row
		runCode <- runCode[as.character(runCode) != 'NA']
		
		return(runCode)
	})

#	Need to get setupCodeTemplates working first...
doubleBracket_keywordInfo <- keywordInfoClass(
	keyword = '[[', 
	processor = function(code, nfProc){
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
keywordList[['values']] <- values_keywordInfo
keywordList[['calculate']] <- calculate_keywordInfo
keywordList[['simulate']] <- simulate_keywordInfo
keywordList[['getLogProb']] <- getLogProb_keywordInfo
keywordList[['nimCopy']] <- nimCopy_keywordInfo
keywordList[['[[']] <- doubleBracket_keywordInfo
keywordList[['$']] <- dollarSign_keywordInfo
keywordList[['[']] <- singleBracket_keywordInfo
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


#	processKeyword function to be called by nfProc
processKeyword <- function(code, nfProc){
  thisKeywordInfo <- keywordList[[ as.character(code[[1]]) ]]
  if(!is.null(thisKeywordInfo))
    return(thisKeywordInfo$processor(code, nfProc))
  return(code)
}







#####	SETUPCODE TEMPLATES

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
                                          
                                          
modelVariableAccessorVector_setupCodeTemplate <- setupCodeTemplateClass(
	#Note to programmer: required fields of argList are model, nodes and logProb

    makeName = function(argList) {Rname2CppName(paste(argList$model, argList$nodes, 'access_logProb', argList$logProb, sep = '_'))},
    codeTemplate = quote( ACCESSNAME <- modelVariableAccessorVector(MODEL, NODES, logProb = LOGPROB) ),
	makeCodeSubList = function(resultName, argList) {
        list(ACCESSNAME = as.name(resultName),
             MODEL = argList$model,
             NODES = argList$nodes,
             LOGPROB = argList$logProb)
    })

modelValuesAccessorVector_setupCodeTemplate <- setupCodeTemplateClass(
	#Note to programmer: required fields of argList are model, nodes and logProb

    makeName = function(argList) {Rname2CppName(paste(argList$model, argList$nodes, 'access_logProb', argList$logProb, sep = '_'))},
    codeTemplate = quote( ACCESSNAME <- modelValuesAccessorVector(MODEL, NODES, logProb = LOGPROB) ),
	makeCodeSubList = function(resultName, argList) {
        list(ACCESSNAME = as.name(resultName),
             MODEL = argList$model,
             NODES = argList$nodes,
             LOGPROB = argList$logProb)
    })


    
accessorVectorLength_setupCodeTemplate <- setupCodeTemplateClass(
  #Note to programmer: required fields of argList are accessName
 
  makeName = function(argList){ Rname2CppName(paste(argList$accessName, 'length', sep = '_')) },
  codeTemplate = quote(ACCESSLENGTH <- ACCESSNAME$getLength()),
  makeCodeSubList = function(resultName, argList){
  	list(ACCESSNAME = as.name(argList$accessName),
  		ACCESSLENGTH = as.name(resultName) )
  	})


nodeFunctionVector_SetupTemplate <- setupCodeTemplateClass(
	#Note to programmer: required fields of argList are model, nodes and includeData
	
	makeName = function(argList){Rname2CppName(paste(argList$model, argList$nodes, 'nodeFxnVector_includeData', argList$includeData, sep = '_'))},
	codeTemplate = quote(NODEFXNVECNAME <- nodeFunctionVector(model = MODEL, nodeNames = NODES, excludeData = EXCLUDEDATA)), 
	makeCodeSubList = function(resultName, argList){
		list(NODEFXNVECNAME = as.name(resultName),
			MODEL = argList$model,
			NODES = argList$nodes,
			EXCLUDEDATA = !argList$includeData)
	})

allLHSNodes_SetupTemplate <- setupCodeTemplateClass(
	#Note to programmer: required fields of argList are model

	makeName = function(argList){
		Rname2CppName(paste('allLHSnodes', argList$model, sep = '_'))
	},
	codeTemplate = quote(NODENAMES <- MODEL$getMaps('nodeNamesLHSall')),
	makeCodeSubList = function(resultName, argList){
		list(NODENAMES = as.name(resultName),
			MODEL = argList$model)
	})
	
allModelNodes_SetupTemplate <- setupCodeTemplateClass(
	#Note to programmer: required fields of argList are model

	makeName = function(argList){
		Rname2CppName(paste('allModelNodes', argList$model, sep = '_'))
	},
	codeTemplate = quote(NODENAMES <- MODEL$getNodeNames()),
	makeCodeSubList = function(resultName, argList){
		list(NODENAMES = as.name(resultName),
			MODEL = argList$model)
	})	
	
allModelValuesVars_SetupTemplate <- setupCodeTemplateClass(
	#Note to programmer: required fields of argList are modelValues

	makeName = function(argList){
		Rname2CppName(paste('allMVVars', argList$modelValues, sep = '_'))
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
addNewCode <- function(name, subList, template, nfProc)
	nfProc$newSetupCode[[name]] <- eval(substitute(substitute(TEMPLATE, subList), list(TEMPLATE = template$codeTemplate) ) )

addNecessarySetupCode <- function(name, argList, template, nfProc){
#	name <- template$makeName(argList)
	if(!isSetupCodeGenerated(name, nfProc)){
		addSetupCodeNames(name, template$makeOtherNames(name, argList), nfProc)
		subList <- template$makeCodeSubList(name, argList)
		addNewCode(name, subList, template, nfProc)
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
  theseFormals <- theseFormals[nchar(theseFormals) > 0]
  matchedCall <- match.call(def, call)
  missingArgs <- which(!(names(theseFormals) %in% names(matchedCall)))
  for(ind in missingArgs){
    name <- names(theseFormals)[ind]
    matchedCall[[name]] <- theseFormals[[name]]    
  }
  return(matchedCall)
}

pasteExpr <- function(expr1, expr2)
	parse(text=paste0(as.character(expr1), as.character(expr2) ) )[[1]]


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

matchFunctions <- new.env()
matchFunctions[['values']] <- function(model, nodes, accessor){}
matchFunctions[['calculate']] <- calculate		#function(model, nodes, nodeFunctionVector){}
matchFunctions[['simulate']] <- simulate		#function(model, nodes, includeData = FALSE, nodeFunctionVector){}
matchFunctions[['getLogProb']] <- getLogProb	#function(model, nodes, nodeFunctionVector){}
matchFunctions[['nimCopy']] <- function(from, to, nodes, nodesTo, row, rowTo, logProb = FALSE){}

matchKeywordCode <- function(code){
	thisFunctionMatch <- matchFunctions[[ as.character( code[[1]] ) ]]
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

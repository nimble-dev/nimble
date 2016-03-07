##		For the code that generates the run function for nimOptim, see nimOptim_keywordInfo$processor
##		When we update the interface to actually work (i.e. pass the results back to an object that can be manipulated in the DSL)
##		we will need to change $processor so that it changes the run code 

##		This is the class that represents the generated C++ wrapper to a nimbleFunction that fits the optim
##		framework, i.e. three arguments are int n, double* par and void* ex
OptimReadyFunction <- setRefClass("OptimReadyFunction", 
							fields = list(name = 'ANY',
										  nimbleFunction = 'ANY',
										  localNimbleFunctionName = 'ANY',
										  generatorName = 'ANY'),
							methods = list(
							getName = function() return(name),
							getNimbleFunction = function() return(nimbleFunction),
							getOriginalName = function() return(localNimbleFunctionName),
							getGeneratorName = function(){
								if(inherits(generatorName, 'uninitializedField'))	generatorName <<- paste0('OPTIMREADY_',environment(nf_getGeneratorFunction(nimbleFunction) )$name)
								return(generatorName)
							}
							)
)


##		This is where the actual C++ wrapper is actually built
cppOptimObject <- function(name, nfSym){
		
	castingExLine <- cppLiteral('vector<void*>* vP = static_cast<vector<void*>*>(ex);')
	thisNimbleCppTypeNamePtr <- paste0( nfSym$nfProc$name, '*')
	thisFunName <- nfSym$name
	castingNFChar <- paste(thisNimbleCppTypeNamePtr, thisFunName, ' = static_cast<', thisNimbleCppTypeNamePtr, '>( (*vP)[0]);')
	castingNFLine <- cppLiteral(castingNFChar)	
	castingNewParLine <- cppLiteral('NimArr<1, double>* newpars = new NimArr<1, double>(n);')
	copyingNewParLine <- cppLiteral(makeNewParsCharCopyLine())

	argInfo <- getAdditonalArgsListForOptim(nfSym)

	castingArgsChunk <- makeStaticCastingChunk(argInfo)
	callLine <- makeOptimFunctionCallLine(thisFunName, argInfo)	
	fnScaleLine <-cppLiteral(makefnScaleLine()) 
	returnLine <- cppLiteral('return(ans);')
	deleteParsLine <- cppLiteral('delete newpars;')
	newCodeLinesList <- c(list(castingExLine, castingNFLine, castingNewParLine, copyingNewParLine), castingArgsChunk, list(callLine, fnScaleLine, deleteParsLine,  returnLine))
	codeBlock <- putCodeLinesInBrackets(codeLines = newCodeLinesList)
	ans <- cppFunctionDef(name = name,
											args = OptimFun_argDefs,
											code = cppCodeBlock(code = codeBlock, skipBrackets = TRUE),
											returnType = cppDouble(),
											externC = FALSE
											)
	return(ans)
}




OptimFun_argDefs <- list(cppInt('n'),
				cppDouble('par', ptr = 1),
				cppVoid('ex', ptr = 1)
				)
				
funName2OptimFunName <- function(nameChar){
	paste0('OPTIMREADY_', nameChar)
}



argInfo2PointerStaticCast <- function(argInf){
	
	argInf <- matchKeywordCode(argInf)
	charInf <- as.character(argInf)
	if(length(charInf) == 1 || charInf[2] == '0')
		return(paste0(charInf[1], '*'))
	
#	if(length(charInf) > 2)		#actually, args could have defaults...
#		stop('too many arguments in argInfo2PointerStaticCast')
	
	if(charInf != 'int' && charInf != 'double' )
		stop('value type not recognized in argInfo2PointerStaticCast')
	
	ans <-paste0('NimArr<', charInf[2], ',' , charInf[1], '>*')
	return(ans)
}

getArgInfoFromNFSym <- function(nfSym){
	nfSym$nfProc$setupSymTab$getSymbolObject('operator()')$nfMethodRCobj$argInfo
}

nimOptim <- matchFunctions[['nimOptim']]

getArgsFromOptimCallForTargetFun <- function(call){
	optSpecificNames <- names(formals(nimOptim))
	allCallNames <- names(call)[-1]
	argIsForOptim <- allCallNames %in% optSpecificNames
	return(allCallNames[!argIsForOptim])
}

makeVoidPointerName_fromObjName <- function(objName)
	Rname2CppName(paste0('vPtrFor_', objName) )

makeDSLCallforVoidPtr_fromArgInfAndCall <- function(call, argInfo){
	argNames <- names(argInfo)
	namesFromCall <- getArgsFromOptimCallForTargetFun(call)

	unnecessaryCallNames <- setdiff(namesFromCall, argNames)
	if(length(unnecessaryCallNames) > 0){
			stop(paste('argument(s) ', paste(unnecessaryCallNames, collapse = ", "), ' provided to nimOptim, but is not part of function to be optimized'))
	}


	if(argNames[1] %in% namesFromCall){
		errorMessage <- paste0('in call to nimOptim, first argument (', argNames[1], ') included as fixed parameter')
		stop(errorMessage)
	}
	newRunCode <- list()
	pointerNames <- list()
	uniquePointerNames <- character()
	if(length(argNames) > 1){
		for(thisName in argNames[-1]){
			if(thisName %in% namesFromCall)
				objName <- as.character(call[[thisName]])
			else{
				if(!('default' %in% names(argInfo[[thisName]]))){
					errorMsg <- paste0('in call to nimOptim, argument for optimized function ',
					 as.character(call$optFun), ', the argument ',
					 thisName, ' missing with no default')
					stop(errorMsg)
					}
				objName <- as.character(argInfo[[thisName]][['default']])
			}
		vptrName <- makeVoidPointerName_fromObjName(objName)
		pointerNames[[thisName]] <- vptrName
		if(!vptrName %in% uniquePointerNames ){
			uniquePointerNames <- c(vptrName, uniquePointerNames)
			runLine <- paste0(vptrName, ' <- voidPtr(',objName, ', ',deparse(argInfo[[thisName]]), ')' )
			newRunCode[[thisName]] <- runLine
			}
		}
	}
	ans <- list(newRunCode = newRunCode, pointerNames = pointerNames)
	return(ans)
}

callToPackVPtrs <- function(listName, voidPointerNames){
	args <- paste(voidPointerNames, collapse = ', ')
	paste0(listName, " = voidPointerList(", args,  ")")
}
getAdditonalArgsListForOptim <- function(nfSym){
	argInfo <- getArgInfoFromNFSym(nfSym)
	argNames <- names(argInfo)
	if(argInfo[[1]] != substitute(double(1) ) ){
		functionName <- nfSym$name
		warningMessage <- paste('trying to optimize function', functionName, 'but first argument is not a numeric vector')
		stop(warningMessage)
		}
		
	staticCastList <- list()
#	staticCastList[[argNames[1] ]] <- argInfo2PointerStaticCast(argInfo[[1]])		#This is unnecessary
	if(length(argNames) >= 2)
	for(i in 2:length(argNames))
		staticCastList[[ argNames[i] ]] <- argInfo2PointerStaticCast(argInfo[[i]] )
	return(staticCastList)
}

makeNewParsCharCopyLine <- function(){
	#This looks like a pretty dumb function. However, it will be shortly expanded to allow us to provide scaling parameters
	return('for(int i = 0; i < n; i++) {(*newpars)[i] = par[i];}')
}

makefnScaleLine <- function(){
	#Again a dumb function. Right now the default is going to be set fnscale = -1
	return('ans = -ans;')
}

makeStaticCastingChunk <- function(argInfo){
	argNames <- names(argInfo)
	ans <- list()
	if(length(argNames) > 0)
		ans[[1]] <- cppLiteral('vector<void*>* argInfoPtr = static_cast<vector<void*>* > ( (*vP)[1] );')
	for(i in 1:length(argInfo)){
		thisType <- argInfo[[i]]
		thisName <- argNames[i]
		thisLine <- paste(thisType, thisName, ' = static_cast<', thisType, '>( (*argInfoPtr)[', i - 1, ']);')
		ans[[i+1]] <- cppLiteral(thisLine)
	}
	return(ans)
}

makeFixedParameterCallCodeForOptim <- function(argInfo){
	ans <- character()
	argNames <- names(argInfo)
	if(length(argInfo) > 0){
		for(i in 1:length(argInfo))
			ans <- paste0(ans, paste0(', (*', argNames[i], ')'))
	}
	return(ans)
}

makeOptimFunctionCallLine <- function(funName, argInfo){
	callChar <- paste0('double ans = (*', funName,') ( (*newpars) ', makeFixedParameterCallCodeForOptim(argInfo), ');')
	return(cppLiteral(callChar))	
}





#OptimReadyFunction <- function(name, nimbleFunction) new('OptimReadyFunction', name = name, nimbleFunction = nimbleFunction)



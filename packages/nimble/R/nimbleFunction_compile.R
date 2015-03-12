virtualNFprocessing <- setRefClass('virtualNFprocessing',
                                   fields = list(
                                       name = 'ANY', ## character
                                       nfGenerator =  'ANY',		#'function',
                                       compileInfos =  'ANY',		#'list', ## A list of RCfunctionCompileClass objects
                                       origMethods = 'ANY',
                                       RCfunProcs =  'ANY',		#'list', ## A list of RCfunProcessing  or RCvirtualFunProcessing objects
                                      ## RCfuns = 'list', ## A list of RCfun objects
                                       cppDef = 'ANY'
                                       ),
                                   methods = list(
                                       ## setClassName = function(className) {
                                       ##     print('what to do with setClassName'); browser()
                                       ##     name <<- className
                                       ##     assign('CclassName', name, envir = environment(nfGenerator)) ## to avoid warnings on <<-
                                       ## },
                                       show = function() {
                                           writeLines(paste0('virtualNFprocessing object ', name))
                                       },
                                       initialize = function(f = NULL, className, virtual = TRUE) {
                                       		compileInfos <<- list()
                                       		RCfunProcs <<- list()
                                       		
                                           if(!is.null(f)) { ## This allows successful default instantiation by R when defining nfProcessing below -- crazy
                                               ## nfGenerator allowed if it is a nimbleFunctionVirtual
                                               if(is.nf(f) | is.nfGenerator(f)) nfGenerator <<- nf_getGeneratorFunction(f)
                                               else if(inherits(f, 'list')) {
                                                   if(length(unique(lapply(f, nfGetDefVar, 'name'))) != 1)
                                                       stop('Error with list of instances not having same nfGenerator')
                                                   nfGenerator <<- nf_getGeneratorFunction(f[[1]])
                                               }
                                               if(missing(className)) {
                                                   sf <- environment(nfGenerator)$name
                                                   ##setClassName(Rname2CppName(deparse(sf)))
                                                   name <<- Rname2CppName(sf)
                                               } else {
                                                   ##setClassName(className)
                                                   name <<- className
                                               }
                                               origMethods <<- nf_getMethodList(nfGenerator)
                                               RCfunProcs <<- list()
                                               for(i in seq_along(origMethods)) {
                                                   RCname <- names(origMethods)[i]
                                                   if(RCname == 'run') RCname <- 'operator()'
                                                   RCfunProcs[[RCname]] <<- if(virtual) RCvirtualFunProcessing$new(origMethods[[i]], RCname) else RCfunProcessing$new(origMethods[[i]], RCname)
                                               }
                                               compileInfos <<- lapply(RCfunProcs,
                                                                       function(x) x$compileInfo)
                                           }
                                        },
                                       setupLocalSymbolTables = function() {
                                           for(i in seq_along(RCfunProcs)) {
                                               RCfunProcs[[i]]$setupSymbolTables()
                                           }
                                       },
                                       doRCfunProcess = function(control = list(debug = FALSE, debugCpp = FALSE)) {
                                           for(i in seq_along(RCfunProcs)) {
                                               RCfunProcs[[i]]$process(debug = control$debug, debugCpp = control$debugCpp, debugCppLabel = name)
                                           }
                                       },
                                       process = function(control = list(debug = FALSE, debugCpp = FALSE)) {
                                           setupLocalSymbolTables()
                                           doRCfunProcess(control)
                                       }
                                       )
                                   )


nfProcessing <- setRefClass('nfProcessing',
                            contains = 'virtualNFprocessing',
                            fields = list(
                                instances = 'ANY',
                                setupSymTab = 'ANY',
                                neededTypes =  'ANY',		#'list', ## A list of symbolTable entries of non-native types, such as derived models or modelValues, that will be needed
                              neededObjectNames =  'ANY',		#'character', ## a character vector of the names of objects such as models or modelValues that need to exist external to the nimbleFunction object so their contents can be pointed to 
                                newSetupOutputNames =  'ANY',		#'character',
                                newSetupCode =  'ANY',		#'list',
                                newSetupCodeOneExpr = 'ANY',
                                nimbleProject = 'ANY',
                                inModel =  'ANY'		#'logical'
                              ),
                          methods = list(
                              show = function() {
                                  writeLines(paste0('nfProcessing object ', name))
                              },
                              initialize = function(f = NULL, className, fromModel = FALSE, project) {
                              	neededTypes <<- list()
                              	neededObjectNames <<- character()
                              	newSetupCode <<- list()
                                  if(!is.null(f)) {
                                      ## in new system, f must be a specialized nf, or a list of them
                                      nimbleProject <<- project
                                      inModel <<- fromModel
                                      if(missing(className)) {
                                          sf <- if(is.list(f)) nfGetDefVar(f[[1]], 'name') else nfGetDefVar(f, 'name')
                                          name <<- Rname2CppName(sf)
                                      } else {
                                          name <<- className
                                      }
                                      callSuper(f, name, virtual = FALSE)
                                      instances <<- if(inherits(f, 'list')) lapply(f, nf_getRefClassObject) else list(nf_getRefClassObject(f))
                                     
                                      newSetupOutputNames <<- character()
                                      newSetupCode <<- list()
                                  }
                              },


								##NEW PROCESSING TOOLS.   
								processKeywords_all = function(){},
								processKeywords_one = function(){},
								matchKeywords_all = function(){},
								matchKeywords_one = function(){},

  
  
                              doSetupTypeInference_processNF = function() {},
                              makeTypeObject = function() {},
                              replaceCall = function() {},
                              evalNewSetupLines = function(){},
                              makeNewSetupLinesOneExpr = function() {},
                              evalNewSetupLinesOneInstance = function(instances, check = FALSE){},
                              setupTypesForUsingFunction = function(){},
                              doSetupTypeInference = function(){},
                              addMemberFunctionsToSymbolTable = function(){},
                              setupLocalSymbolTables = function() {
                                  for(i in seq_along(RCfunProcs)) {
                                      RCfunProcs[[i]]$setupSymbolTables(parentST = setupSymTab)
                                  }
                              },
                              ## buildRCfun = function(RCname) {
                              ##     ans <- RCfunctionDef$new()
                              ##     ans$buildFunction(RCfunProcs[[RCname]])
                              ##     ans$buildSEXPinterfaceFun(className = name)
                              ##     ans
                              ## },
                              ## buildRCfuns = function() {
                              ##     RCfuns <<- list()
                              ##     for(i in seq_along(RCfunProcs)) {
                              ##         RCname <- names(RCfunProcs)[i]
                              ##         RCfuns[[RCname]] <<- buildRCfun(RCname)
                              ##     }
                              ## },
                              collectRCfunNeededTypes = function() {
                                  for(i in seq_along(RCfunProcs)) {
                                      for(j in names(RCfunProcs[[i]]$neededRCfuns)) {
                                          if(is.null(neededTypes[[j]])) {
                                              neededTypes[[j]] <<- RCfunProcs[[i]]$neededRCfuns[[j]]
                                          }
                                      }
                                      ## could clear RCfunProc[[i]]$neededRCtypes, but instead will prevent them from being used at compilation
                                  }
                              },
                              addBaseClassTypes = function() {
                                  ## If this class has a virtual base class, we add it to the needed types here
                                  contains <- environment(nfGenerator)$contains
                                  if(!is.null(contains)) {
                                      ## generatorFun <- nf_getGeneratorFunction(contains)
                                      className <- environment(contains)$className
                                      nfp <- nimbleProject$setupVirtualNimbleFunction(contains, fromModel = inModel)
                                      newSym <- symbolNimbleFunction(name = name, type = 'nimbleFunctionVirtual', nfProc = nfp) 
                                      if(!(className %in% names(neededTypes))) neededTypes[[className]] <<- newSym
                                  }
                              },
                              process = function(control = list(debug = FALSE, debugCpp = FALSE)) {
                                  ## Modifications to R code
                                  debug <- control$debug
                                  debugCpp <- control$debugCpp
                                  if(!is.null(nimbleOptions$debugNFProcessing)) {
                                      if(nimbleOptions$debugNFProcessing) {
                                          debug <- TRUE
                                          control$debug <- TRUE
                                          writeLines('Debugging nfProcessing (nimbleOptions$debugRCfunProcessing is set to TRUE)') 
                                      }
                                  }
                                  
                                  if(debug) {
                                      print('setupSymTab')
                                      print(setupSymTab)
                                      
                                      writeLines('***** READY FOR replaceModelSingleValues *****')
                                      browser()
                                  }

                                  if(inherits(setupSymTab, 'uninitializedField')) {
                                      ## This step could have already been done if the types were needed by another nimbleFunction
                                      setupTypesForUsingFunction()
                                  }

                                  if(debug) browser()
                                  addBaseClassTypes()
								
									##NEW PROCESSING TOOLS.   
									matchKeywords_all()
									processKeywords_all()



                                  if(debug) browser()

                                  makeNewSetupLinesOneExpr()
                                  
                                  evalNewSetupLines()
								
                                  if(debug) {
                                      print('setupSymTab')
                                      print(setupSymTab) 
                                      print('newSetupOutputNames')
                                      print(newSetupOutputNames)
                                      print('newSetupCode')
                                      print(newSetupCode)
                                      writeLines('***** READY FOR doSetupTypeInference *****')
                                      browser()
                                  }
                                  doSetupTypeInference(setupOrig = FALSE, setupNew = TRUE)  
                                                                    
                                  if(debug) {
                                      writeLines('***** READY FOR setupLocalSymbolTables *****')
                                      browser()
                                  }
                                  setupLocalSymbolTables()

                                  if(debug) {
                                      print('lapply(compileInfos, function(x) print(x$newLocalSymTab))')
                                      lapply(compileInfos, function(x) print(x$newLocalSymTab))
                                      writeLines('**** READY FOR RFfunProcessing *****')
                                      browser()
                                  }

                                  doRCfunProcess(control)

                                  collectRCfunNeededTypes()

##                                  buildRCfuns() ## This should go somewhere else, I think in cppNimbleFunction
                                  
                                  if(debug) {
                                      print('done with RCfunProcessing')
                                      browser()
                                  }
                              }
                              )
                          )

nfProcessing$methods(evalNewSetupLines = function() {
    if(length(instances) == 0)      { warning('No specialized instances of nimble function');   return() }
    for(i in seq_along(instances)) {
        evalNewSetupLinesOneInstance(instances[[i]])
    }
})

nfProcessing$methods(makeNewSetupLinesOneExpr = function() {
    newSetupCodeOneExpr <<- as.call(c(list(as.name('{')), newSetupCode))
})
                     
nfProcessing$methods(evalNewSetupLinesOneInstance = function(instance, check = FALSE) {
    if(is.nf(instance)) instance <- nf_getRefClassObject(instance)
    if(check) {
        go <- if(inherits(instance$.newSetupLinesProcessed, 'uninitializedField')) {
            TRUE
        } else {
            if(length(instance$.newSetupLinesProcessed) == 0) TRUE else !instance$.newSetupLinesProcessed
        }
        if(!go) return(invisible(NULL))
    }
    eval(newSetupCodeOneExpr, envir = instance)
    ## for(j in seq_along(newSetupCode)) {
    ##     ## Warning: this relies on the fact that althought refClass environments are closed, we can
    ##     ## eval in them and create new variables in them that way.
    ##     eval(newSetupCode[[j]], envir = instance)
    ## }
    instance$.newSetupLinesProcessed <- TRUE
})

nfProcessing$methods(setupTypesForUsingFunction = function() {
    doSetupTypeInference(TRUE, FALSE)
    addMemberFunctionsToSymbolTable()
})

nfProcessing$methods(addMemberFunctionsToSymbolTable = function() {
    for(i in seq_along(origMethods)) {
        thisName <- names(origMethods)[i]
        if(thisName == 'run') thisName <- 'operator()'
        newSym <- symbolMemberFunction(name = thisName, nfMethodRCobj = origMethods[[i]])
        setupSymTab$addSymbol(newSym)
    }
})

nfProcessing$methods(doSetupTypeInference = function(setupOrig, setupNew) {
    if(!setupOrig & !setupNew) {
        warning('Weird, doSetupTypeInference was called with both setupOrig and setupNew FALSE.  Nothing to do.', call. = FALSE)
        return(NULL)
    }
    if(length(instances) == 0)      {
        warning('No specialized instances of nimble function', call. = FALSE)
        return(NULL)
    } 
    outputNames <- character()
    if(setupOrig) {
    	setupSymTab <<- symbolTable(parentST = NULL)
    	outputNames <- c(outputNames, nf_getSetupOutputNames(nfGenerator))
    }
    if(setupNew) {
        outputNames <- c(outputNames, newSetupOutputNames)
    }
    doSetupTypeInference_processNF(setupSymTab, outputNames, instances, add = TRUE)   # add info about each setupOutput to symTab
})

nfProcessing$methods(doSetupTypeInference_processNF = function(symTab, setupOutputNames, instances, add = FALSE, firstOnly = FALSE) {
    if(length(instances) == 0) {
        warning('Can not infer setup output types with no instances.')
        return(invisible(NULL))
    }
    for(name in setupOutputNames) {
        symbolRCobject <- makeTypeObject(name, instances, firstOnly)
        if(is.logical(symbolRCobject)) {
            warning(paste0('There is an error involving the type of ', name,'.  Going into browser(). Press Q to exit.'), call. = FALSE)
            browser()
        }
        if(add) symTab$addSymbol(symbolRCobject)
    }
})


nfProcessing$methods(getModelVarDim = function(modelVarName, labelVarName, firstOnly = FALSE) {
    firstNDim <- instances[[1]][[modelVarName]]$modelDef$varInfo[[labelVarName]]$nDim
    if(!firstOnly) {
        if(!all(unlist(lapply(instances, function(x) x[[modelVarName]]$modelDef$varInfo[[labelVarName]]$nDim == firstNdim)))) {
            warning(paste0('Problem: not all instances of label ',labelVarName,' in model ', modelVarName, ' have the same number of dimensions.'))
            return(invisible(NULL))
        }
    }
    return(firstNDim)
})

## firstOnly is supposed to indicate whether we look at only the first instance, or use all of them
## but actually, right now, we use it inconsistently.
## this is a function that could use a lot of polishing, but it's ok for now.
nfProcessing$methods(makeTypeObject = function(name, instances, firstOnly = FALSE) { 
    if(inherits(instances[[1]][[name]], 'nimbleFunctionList')) {
        
        neededObjectNames <<- c(neededObjectNames, name)
        baseClass <- instances[[1]][[name]]$baseClass ## an nfGenerator created by virtualNimbleFunction()
        baseClassName <- environment(baseClass)$className
 
        newSym <- symbolNimbleFunctionList(name = name, type = 'nimbleFunctionList', baseClass = baseClass)
        if(!(baseClassName %in% names(neededTypes))) {
            nfp <- nimbleProject$setupVirtualNimbleFunction(baseClass, fromModel = inModel)
            neededTypeSim <- symbolNimbleFunction(name = baseClassName, type = 'virtualNimbleFunction', nfProc = nfp)
            neededTypes[[baseClassName]] <<- newSym
        }

        allInstances <- unlist(lapply(instances, function(x) x[[name]]$contentsList), recursive = FALSE)
        newNFprocs <- nimbleProject$compileNimbleFunctionMulti(allInstances, initialTypeInference = TRUE)
        ## only types are needed here, not initialTypeInference, because nfVar's from a nimbleFunctionList are not available (could be in future)
        for(nfp in newNFprocs) {
            newTypeName <- environment(nfp$nfGenerator)$name
            neededTypes[[ newTypeName ]] <<- symbolNimbleFunction(name = newTypeName, type = 'nimbleFunction',
                                                                  nfProc = nfp)
        }
        return(newSym)
    }
    if(inherits(instances[[1]][[name]], 'OptimReadyFunction')){
    	thisObj <- instances[[1]][[name]]
    	nfName <- thisObj$getOriginalName()
        nfp <- getnfProcFromProject_fromnfClassName(nfName, nimbleProject) ## will return existing nfProc if it exists
    	newSym <- symbolOptimReadyFunction(name = name, type = 'OptimReadyFunction', nfName = nfName, nfProc = nfp)
    	return(newSym)
    }
    if(is.character(instances[[1]][[name]])) {
        return(symbolBase(name = name, type = 'Ronly'))
    }
    if(is.nf(instances[[1]][[name]])) { ## nimbleFunction
        funList <- lapply(instances, `[[`, name)
        nfp <- nimbleProject$compileNimbleFunction(funList, initialTypeInferenceOnly = TRUE) ## will return existing nfProc if it exists

        className <- class(nf_getRefClassObject(funList[[1]]))
        neededObjectNames <<- c(neededObjectNames, name)
        newSym <- symbolNimbleFunction(name = name, type = 'nimbleFunction', nfProc = nfp)
        if(!(className %in% names(neededTypes))) neededTypes[[className]] <<- newSym
        return(newSym)
    }
    if(inherits(instances[[1]][[name]], 'modelValuesBaseClass')) { ## In some cases these could be different derived classes.  If locally defined they must be the same
        if(!firstOnly) {
            if(!all(unlist(lapply(instances, function(x) inherits(x[[name]], 'modelValuesBaseClass'))))) {
                warning(paste0('Problem: some but not all instances have ', name,' as a modelValues.  Types must be consistent.'))
                return(invisible(NULL))
            }
        }
        ## Generate one set of symbolModelValues objects for the neededTypes, and each of these can have its own mvSpec
        ## Generate another symbolModelValues to return and have in the symTab for this compilation
        ## I don't think that mvSpec gets used, since they all get Values *        
        for(i in seq_along(instances)) {
            className <- class(instances[[i]][[name]])
            if(!(className %in% names(neededTypes))) {
                ## these are used only to build neededTypes
                ntSym <- symbolModelValues(name = name, type = 'Values', mvSpec = instances[[i]][[name]]$mvSpec)
                neededTypes[[className]] <<- ntSym
            }
        }
        ## this is used in the symbol table
        neededObjectNames <<- c(neededObjectNames, name)
        newSym <- symbolModelValues(name = name, type = 'Values', mvSpec = NULL)
        return(newSym)
    }
    if(inherits(instances[[1]][[name]], 'modelBaseClass')) {
        if(!firstOnly) {
            if(!all(unlist(lapply(instances, function(x) inherits(x[[name]], 'modelBaseClass'))))) {
                warning(paste0('Problem: some but not all instances have ', name,' as a model.  Types must be consistent.'))
                return(invisible(NULL))
            }
            if(!all(unlist(lapply(instances, function(x) inherits(x[[name]], 'RModelBaseClass'))))) {
                warning(paste0('Problem: models should be provided as R model objects, not C model objects'))
                return(invisible(NULL))
            }
        }
        return(symbolModel(name = name, type = 'Ronly', className = class(instances[[1]][[name]]))) 
    }

    if(inherits(instances[[1]][[name]], 'NumericListBase')) {
    	    	
    	varinfo <- instances[[1]][[name]]
        if(!firstOnly) {
            if(!all(unlist(lapply(instances, function(x) inherits(x[[name]], 'NumericListBase'))))) {
                warning(paste0('Problem: some but not all instances have ', name,' as a NumericList.  Types must be consistent.'))
                return(invisible(NULL))
            }
        }
        
        return(symbolNumericList(name = name, type = varinfo$listType, nDim = max(varinfo$nDim, 1),  className = class(instances[[1]][[name]]))) 
    }

    if(inherits(instances[[1]][[name]], 'singleVarAccessClass')) {
        ## Keeping this simple: only doing first instance for now
        varInfo <- instances[[1]][[name]]$model$getVarInfo( instances[[1]][[name]]$var )
        ## Maybe we should intercept this case in the model, but for now here:
        if(instances[[1]][[name]]$useSingleIndex) {
            nDim <- 1
            size <- prod(varInfo$maxs)
        } else {
            nDim <- varInfo$nDim
            size <- varInfo$maxs
            if(length(nDim) == 0) browser()
            if(is.na(nDim)) browser()
            if(nDim == 0) {nDim <- 1; size <- 1;} ## There is no such thing as a scalar in a model
        }
        return(symbolNimArrDoublePtr(name = name, type = 'double', nDim = nDim, size = size))
    }

    if(inherits(instances[[1]][[name]], 'singleModelValuesAccessClass')) {
        
        varOrgName <- instances[[1]][[name]]$var
        varSym <- instances[[1]][[name]]$modelValues$symTab$getSymbolObject(varOrgName)
        nDim <- max( c(varSym$nDim, 1) )
        type = instances[[1]][[name]]$modelValues$symTab$symbols[[varOrgName]]$type
        return(symbolVecNimArrPtr(name = name, type = type, nDim = nDim, size = varSym$size))
    }


    if(inherits(instances[[1]][[name]], 'nodeFunctionVector')) { 
        return(symbolNodeFunctionVector(name = name))
    }
    if(inherits(instances[[1]][[name]], 'modelVariableAccessorVector')){
    	return(symbolModelVariableAccessorVector(name = name, lengthName = paste0(name, '_length')) )
    }
    if(inherits(instances[[1]][[name]], 'modelValuesAccessorVector')){
    	return(symbolModelValuesAccessorVector(name = name) )     	
    }    
    if(is.numeric(instances[[1]][[name]]) | is.logical(instances[[1]][[name]])) {
        if(firstOnly) {
            type <- storage.mode(instances[[1]][[name]])
            nDim <- if(is.null(dim(instances[[1]][[name]]))) 1L else length(dim(instances[[1]][[name]]))
            size <- if(length(instances[[1]][[name]])==1) rep(1L, nDim) else rep(as.numeric(NA), nDim)
            return(symbolBasic(name = name, type = type, nDim = nDim, size = size))
        } else {
            instanceObjs <- lapply(instances, `[[`, name)
            types <- unlist(lapply(instanceObjs, storage.mode))
            dims <- lapply(instanceObjs, dim)
            dimsNULL <- unlist(lapply(dims, is.null))
            if(any(dimsNULL)) { ## dimsNULL TRUE means it is a vector
                if(!all(dimsNULL)) {
                    warning(paste0('Problem, dimensions do no all match for ', name))
                    return(NA)
                }
                nDim <- 1
                lengths <- unlist(lapply(instanceObjs, length))
                size <- if(!all(lengths == 1)) as.numeric(NA) else 1L
                if(nimbleOptions$convertSingleVectorsToScalarsInSetupArgs) {
                    if(nDim == 1 & identical(as.integer(size), 1L)) nDim <- 0
                }
            } else {
                ## no dims are null, so everything is matrix or array
                dimsLengths <- unlist(lapply(dims, length))
                if(length(unique(dimsLengths)) > 1) {
                    warning(paste0('Problem, dimensions do no all match for ', name))
                    return(NA)
                    }
                nDim <- dimsLengths[[1]]
                size <- rep(as.numeric(NA), nDim)
            }
            
            if(any(types == 'double')) {
                if(!all(types %in% c('double','integer'))) {
                    warning('Problem: some but not all instances have ', name, ' as double or integer.  Types must be consistent.')
                    return(NA)
                }
                return(symbolBasic(name = name, type = 'double', nDim = nDim, size = size))
            }
            if(any(types == 'integer')) {
                if(!all(types == 'integer')) {
                    warning('Problem: some but not all instances have ', name, ' as integer.  Types must be consistent.')
                    return(NA)
                }
                return(symbolBasic(name = name, type = 'integer', nDim = nDim, size = size))
            }
            if(any(types == 'logical')) {
                if(!all(types == 'logical')) {
                    warning('Problem: some but not all instances have ', name, ' as logical.  Types must be consistent.')
                    return(NA)
                }   
                return(symbolBasic(name = name, type = 'logical', nDim = nDim, size = size))
            }
        }
    }
    return(NA)
})




nfProcessing$methods(determineNdimsFromInstances = function(modelExpr, varOrNodeExpr) {
    allNDims <- lapply(instances, function(x) {
        model <- eval(modelExpr, envir = x)
        if(!exists(as.character(varOrNodeExpr), x, inherits = FALSE) ) {
            stop(paste0('Error, ', as.character(varOrNodeExpr), ' does not exist in an instance of this nimbleFunction.'))
        }
        lab <- eval(varOrNodeExpr, envir = x)
        varAndIndices <- getVarAndIndices(lab)
        determineNdimFromOneCase(model, varAndIndices)
    } )
    return(allNDims)
})




nfProcessing$methods(processKeywords_all = function(){	
	for(i in seq_along(compileInfos)){
		compileInfos[[i]]$newRcode <<- processKeywords_one(compileInfos[[i]]$origRcode)
		}
	})

nfProcessing$methods(processKeywords_one = function(code){
	cl = length(code)
	if(cl == 1){
		if(is.call(code)){
			if(length(code[[1]]) > 1)	code[[1]] <- processKeywords_one(code[[1]])
		}
		return(code)
	}

	if(length(code[[1]]) == 1)
		{
		code <- processKeyword(code, .self)
		}

	cl = length(code)

  if(is.call(code)) {
        if(length(code[[1]]) > 1) code[[1]] <- processKeywords_one(code[[1]])
        if(cl >= 2) {
            for(i in 2:cl) {
                code[[i]] <- processKeywords_one(code[[i]])
            }
        }
    }
	return(code)
})

nfProcessing$methods(matchKeywords_all = function(){
	for(i in seq_along(compileInfos))
	compileInfos[[i]]$origRcode <<- matchKeywords_one(compileInfos[[i]]$origRcode)
})

nfProcessing$methods(matchKeywords_one = function(code){
	cl = length(code)
	if(cl == 1){
		if(is.call(code)){
			if(length(code[[1]]) > 1)	code[[1]] <- matchKeywords_one(code[[1]])
		}
		return(code)
	}
	if(length(code[[1]]) == 1)
		code <- matchKeywordCode(code) 
	if(is.call(code)) {
        if(length(code[[1]]) > 1) code[[1]] <- matchKeywords_one(code[[1]])
        if(cl >= 2) {
            for(i in 2:cl) {
                code[[i]] <- matchKeywords_one(code[[i]])
            }
        }
    }
    return(code)
})






singleVarAccessClass <- setRefClass('singleVarAccessClass',
                                    fields = list(model = 'ANY', var = 'ANY', useSingleIndex = 'ANY'),
                                    methods = list(
                                        show = function() {
                                            writeLines(paste('singleVarAccess for model',model$name,'to var',var))
                                        })
                                    )

singleVarAccess <- function(model, var, useSingleIndex = FALSE) {
    singleVarAccessClass$new(model = model, var = var, useSingleIndex = useSingleIndex)
}


singleModelValuesAccessClass <- setRefClass('singleModelValuesAccessClass',
                                    fields = list(modelValues = 'ANY', var = 'ANY'),
                                    methods = list(
                                        show = function() {
                                            writeLines(paste('singleModelValuesAccess for model to var',var))
                                        })
                                    )

singleModelValuesAccess <- function(modelValues, var) {
    singleModelValuesAccessClass$new(modelValues = modelValues, var = var)
}



getnfProcFromProject_fromnfClassName <- function(nfClassName, nimbleProject){
	compInfoNames <- names(nimbleProject$nfCompInfos)
	i = 0
	foundClass = FALSE
	while(i < length(compInfoNames) & !foundClass){
		i = i+1
		foundClass <- grepl(pattern = compInfoNames[i], x = nfClassName)
	}
	if(!foundClass){
		errorMessage <-paste0("class not found in available names. nfClassName = ", nfClassName, " available names = ", compInfoNames)
		stop(errorMessage)
	}
	nfProc <- nimbleProject$nfCompInfos[[compInfoNames[i] ]]$nfProc
	return(nfProc)
}

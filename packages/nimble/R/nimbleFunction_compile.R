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
                              replaceAccessors = function(){},


								##NEW PROCESSING TOOLS.   
								processKeywords_all = function(){},
								processKeywords_one = function(){},
								matchKeywords_all = function(){},
								matchKeywords_one = function(){},

  
  
                              doSetupTypeInference_processNF = function() {},
                              makeTypeObject = function() {},
                              replaceCall = function() {},
                              replaceAccessorsOneFunction = function() {},
                              evalNewSetupLines = function(){},
                              makeNewSetupLinesOneExpr = function() {},
                              evalNewSetupLinesOneInstance = function(instances, check = FALSE){},
                              setupTypesForUsingFunction = function(){},
                              doSetupTypeInference = function(){},
                              replaceOneCalcGLP = function(){},
                              replaceOneSimulate = function(){},
                              replaceOneValues = function() {},
                              replaceOneGetOrSetValues = function(){},
                              replaceOneModelAccessor = function(){},
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
                                  
                            #      replaceAccessors()
								
								
									##NEW PROCESSING TOOLS.   
									matchKeywords_all()
									processKeywords_all()

							#		print(newSetupCode)


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


makeMapSetupCodeNames <- function(baseName) {
    list(Mstrides = paste0(baseName, '_strides'),
         Msizes = paste0(baseName, '_sizes'),
         Moffset = paste0(baseName, '_offset'))
}


makeSingleIndexAccessCodeNames <- function(baseName) {
    list(MflatIndex = paste0(baseName, '_flatIndex'))
}


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

nfProcessing$methods(replaceOneModelAccessor = function(code) {
    useSingleIndexAccessor <- FALSE
    if(code[[1]] == '[[') {
        if(!is.character(code[[3]])) { ## It is an expression, like model[[var]] or model[[node]], so inspect instances
            allNDims <- determineNdimsFromInstances(code[[2]], code[[3]])
            if(length(unique(allNDims)) > 1) stop(paste0('Error for ', deparse(code), '.  There are inconsistent numbers of dimensions in different instances.'))
            nDim <- allNDims[[1]]
            useMap <- nDim > 0
            useSingleIndexAccessor <- !useMap
            labelGiven <- code[[3]] ## this is not used if useMap == TRUE, so it is only used for model[[var]] case, not model[[node]]
        } else { ## it is something like model[['x[1:3]']], with the node hard-wired
            labelAndIndices <- getVarAndIndices(code[[3]])
            labelGiven <- as.character(labelAndIndices$varName)
            nDim <- sum(1-unlist(lapply(labelAndIndices$indices, is.numeric))) ##length(labelAndIndices$indices)
            useMap <- nDim > 0
        }
    } else {
        labelGiven <- as.character(code[[3]])
        useMap <- FALSE
        ## nDim should be unnecessary
    }
    newName <- Rname2CppName(deparse(code))
    newNameExpr <- as.name(newName)
    if(!(newName %in% newSetupOutputNames)) {
        if(useSingleIndexAccessor) {
            newSetupCodeInfo <- makeSingleIndexAccessorSetupCode(newName, code[[2]], code[[3]])
            newSetupCode[[newName]] <<- newSetupCodeInfo$newSetupCode
            newSetupOutputNames <<- c(newSetupOutputNames, newSetupCodeInfo$newOutputNames)
        } else {
            if(useMap) {
                newSetupCodeInfo <- makeMapSetupCodeExprs(newName, code[[2]], code[[3]])
                newSetupCode[[newName]] <<- newSetupCodeInfo$newSetupCode
                newSetupOutputNames <<- c(newSetupOutputNames, newSetupCodeInfo$newOutputNames)
            } else {
                newSetupOutputNames <<- c(newSetupOutputNames, newName) ## These will have to be evaluated before doing type inference on them
                nsc <- substitute(nn <- singleVarAccess(oc2, oc3), list(nn = newNameExpr, oc2 = code[[2]], oc3 = labelGiven))
                newSetupCode[[ newName ]] <<- nsc
            }
        }
    }
    if(useSingleIndexAccessor) {
        ans <- makeSingleIndexAccessExpr(newName, newNameExpr)
    } else {
        if(useMap) {
            ans <- makeMapAccessExpr(newName, newNameExpr, nDim) 
        } else {
            ans <- newNameExpr
        }
    }
    return(ans)
})

nfProcessing$methods(replaceOneModelValuesAccessorSingleBrackets = function(code) {
    ## upon entry we have already checked that length(code[[2]])==1
    charCode2 <- as.character(code[[2]])
    if( charCode2 %in% setupSymTab$getSymbolNames() ) { 

        sym = setupSymTab$getSymbolObject(charCode2)
        if( inherits(sym, 'symbolModelValues') ) {
            varGiven = code[[3]] ## We decided the syntax is mv['x', i]
            row = code[[4]]
            newName = Rname2CppName(deparse(code) )
            newNameExpr = as.name(newName)
            if(!(newName %in% newSetupOutputNames) ) {
                newSetupOutputNames <<- c(newSetupOutputNames, newName)
                nsc <- substitute(nn <- singleModelValuesAccess(MODELVALUES, VARIABLE), list(nn = newNameExpr, MODELVALUES = code[[2]], VARIABLE = varGiven))
                newSetupCode[[ newName ]] <<- nsc
            }
            ans <- substitute( NAME[ROW], list(NAME = newNameExpr, ROW = row))
            return(ans)
        }
      ##  writeLines(paste("Warning, not sure what to do with", deparse(code), "because the LHS of [ is not a modelValues"))
    } ##else {
      ## This warning is unfounded because it could just be a local variable 
      ##  writeLines(paste("Warning, not sure what to do with", deparse(code), "because the LHS of [ is not in the symbol table"))
##    }
    return(code)
})

nfProcessing$methods(replaceOneValues = function(code) { ## new syntax with values(model, nodes) <- P, or P <- values(model, node)
    if(length(code) != 3) stop("Error in processing code, values() should have 2 arguments", call. = FALSE)
    newRunCode <- code[1:2]
    newName <- paste(code[[2]], code[[3]], "access", sep = "_") 
    newNameExpr <- as.name(newName)
    newRunCode[[2]] = newNameExpr
    if(!(newName %in% newSetupOutputNames)) {
        newSetupOutputNames <<- c(newSetupOutputNames, newName) ## These will have to be evaluated before doing type inference on them
        newSetupLine <- substitute(ACCESSNAME <- modelVariableAccessorVector(MODEL, NODES, logProb = FALSE), list(ACCESSNAME = newNameExpr, MODEL = code[[2]], NODES = code[[3]]))
        newSetupCode[[ newName ]] <<- newSetupLine
    }
    lengthName <- paste0(newName,'_length')
    if(!(lengthName %in% newSetupOutputNames)) {
        newSetupOutputNames <<- c(newSetupOutputNames, lengthName)
        newSetupLine <- substitute({LENGTHNAME <- ACCESSNAME$length; if(length(LENGTHNAME) != 1) stop("Error making a length for a values accessor", call. = FALSE)}, list(LENGTHNAME = as.name(lengthName), ACCESSNAME = newNameExpr))
        newSetupCode[[ lengthName ]] <<- newSetupLine
    }
    return(newRunCode) ## used to be one line up inside the if clause
} )

nfProcessing$methods(replaceOneGetOrSetValues = function(code) {
    if(length(code) == 3) return(code)

    if(length(code) != 4) stop("Error in processing code: unrecognized number of arguments for setValues or getValues")
    newRunCode <- code[1:3]
    newName <- paste(code[[3]], code[[4]], "access", sep = "_") 
    newNameExpr <- as.name(newName)
    newRunCode[[3]] = newNameExpr
    if(!(newName %in% newSetupOutputNames)) {
        newSetupOutputNames <<- c(newSetupOutputNames, newName) ## These will have to be evaluated before doing type inference on them
        newSetupLine = substitute(ACCESSNAME <- modelVariableAccessorVector(MODEL, NODES, logProb = FALSE), list(ACCESSNAME = newNameExpr, MODEL = code[[3]], NODES = code[[4]]))
        newSetupCode[[ newName ]] <<- newSetupLine
    }
    return(newRunCode) ## used to be one line up inside the if clause
} )

nfProcessing$methods(replaceOneSimulate = function(code) {
    newRunCode <- code
    matchCode <- match.call(simulate, code)
    varName = as.character(matchCode[['model']])
    if(!inherits(setupSymTab$symbols[[varName]], "symbolModel")) {
        warning(paste("This call to simulate does not simulate does not seem to provide a valid model argument: ", deparse(code)), call. = FALSE )
        return(newRunCode)
    }
    simData = FALSE
    if('includeData' %in% names(matchCode) ) 
        simData <- matchCode[['includeData']]
    nodeNames <- paste('all', varName, 'nodes', sep = '_')
    if('nodes' %in% names(matchCode) ) {
        nodeNames = as.character(matchCode[['nodes']])
    } else {
        if(!(nodeNames %in% newSetupOutputNames) ){
            newSetupOutputNames <<- c(newSetupOutputNames, nodeNames)
            ## It is usually redundant to call getDependencies here, but it allows the includeData argument
            newSetupLine = substitute(NODESNAME <- MODEL$getDependencies(MODEL$getMaps('nodeNamesLHSall'), includeData = SIMDATA), list(NODESNAME = as.name(nodeNames), MODEL = as.name(varName), SIMDATA = simData ) )
            newSetupCode[[nodeNames]] <<- newSetupLine
        }
    }
    
    if(simData == TRUE)
        newName = paste(varName, nodeNames, 'nodeFxnVector', 'simulateNoData', sep = '_')
    else
        newName = paste(varName, nodeNames, 'nodeFxnVector', 'simulate', sep = '_')
    newNameExpr <- as.name(newName)
    newRunCode = substitute(simulate(NEWNFV), list(NEWNFV = newNameExpr) )
    if(!(newName %in% newSetupOutputNames) ){
        newSetupOutputNames <<- c(newSetupOutputNames, newName) ## These will have to be evaluated before doing type inference on them
        newSetupLine = substitute(NODEFXNVEC <- nodeFunctionVector(MODEL, nl_nodeVectorReadyNodes(MODEL, NODES, SIMDATA) ),  list(NODEFXNVEC = newNameExpr, MODEL = as.name(varName), SIMDATA = simData, NODES = as.name(nodeNames) ) )
        newSetupCode[[ newName ]] <<- newSetupLine
    }
    return(newRunCode)
})

nfProcessing$methods(replaceOneCalcGLP = function(code) {
    if(length(code) < 2) stop(paste("Error, don't know what to do with", deparse(code)))
    
    if(length(code) == 2) { ## Default case like calculate(model), where default is to use all nodes
        newRunCode <- code
        matchCode <- match.call(calculate, code)
        varName = as.character(matchCode[['model']])
        if(!inherits(setupSymTab$getSymbolObject(varName), 'symbolModel')) return(newRunCode)
        nodeNames <- 'all'
        if('nodes' %in% names(matchCode) ) 
            nodeNames <- as.character(matchCode[['nodes']])

        newName <- paste(varName, nodeNames, 'nodeFxnVector', sep = '_')
        nodeArg <- substitute(MODEL$getMaps('nodeNamesLHSall'), list(MODEL = code[[2]]) )    
    } else {## end length(code)==2
        if(length(code) != 3) stop(paste("Error in processing code: unrecognized number of arguments for calculate, simulate or getLogProbs :", deparse(code)))
        newRunCode <- code[1:2]
        newName <- paste(code[[2]], code[[3]], "nodeFxnVector", sep = "_")
        nodeArg <- code[[3]]
    }
    newNameExpr <- as.name(newName)
    newRunCode[[2]] = newNameExpr
    if(!(newName %in% newSetupOutputNames)) {
        newSetupOutputNames <<- c(newSetupOutputNames, newName) ## These will have to be evaluated before doing type inference on them
        newSetupLine = substitute(ACCESSNAME <- nodeFunctionVector(MODEL, NODES), list(ACCESSNAME = newNameExpr, MODEL = code[[2]], NODES = nodeArg))
        newSetupCode[[ newName ]] <<- newSetupLine
    }
    return(newRunCode)
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



nfProcessing$methods(replaceAccessors = function() {
    for(i in seq_along(compileInfos)) {
        compileInfos[[i]]$newRcode <<- replaceAccessorsOneFunction(compileInfos[[i]]$origRcode)
		print(compileInfos[[i]]$newRcode)
    }
})

## turn m$x into access(m_x), with setupCode m_x <- singleVarAccess(m, 'x')
## turn m[[v]] into access(m_v) with setupCode m_v <- singleVarAccess(m, v)
## These are for variables
## For nodes we will use m['x[1]'] or m[n]
nfProcessing$methods(replaceAccessorsOneFunction = function(code) {
    cl <- length(code)	
    if(cl==1) {
        if(is.call(code)) { ## Catches a case like A[[i]]() which has length 1 but A[[i]] should be processed
            if(length(code[[1]]) > 1) code[[1]] <- replaceAccessorsOneFunction(code[[1]])
        }
        return(code)
    }
    if(length(code[[1]]) == 1) {    	
        charCode1 <- as.character(code[[1]])
        if(charCode1 %in% c('$', '[[')) {
            varExpr <- code[[2]]
            varName <- as.character(varExpr)

            if(inherits(setupSymTab$symbols[[varName]], 'symbolNimPtrList')) {
             ##   code[[1]] <- as.name('nimFunListAccess') ## will be handled by `[[` name
                return(code)
            }
            if(inherits(setupSymTab$symbols[[varName]], "symbolModel"))
                varType <- "model"
            else if(inherits(setupSymTab$symbols[[varName]], "symbolNumericList"))
                varType <- "numericList"
            else stop(paste('Error: We do not know what to do with', deparse(code)))
            
            if(varType == "model") {
                return(replaceOneModelAccessor(code))
            } 

            if(varType == "numericList"){
                ans <- substitute( (numListAccess(NAME)[ROW]), list(NAME = as.name(varName), ROW = code[[3]]) )
                return(ans)
            } 
            
            stop("Calling '[[]]' on a data type which does not support this\n")
        } ## end $ or `[[`
        if(charCode1 == '[' & length(code[[2]]) == 1) {
            return(replaceOneModelValuesAccessorSingleBrackets(code))
        } ## end `[`

        if(charCode1 == 'values') {
            return(replaceOneValues(code))
        }
        if(charCode1 == 'setValues' || charCode1 == 'getValues') {
            return(replaceOneGetOrSetValues(code))
        } ## end setValues and getValues
        if(charCode1 %in% c( 'calculate', 'getLogProb') ) {
            return(replaceOneCalcGLP(code))
        } ## end calculate,  getLogProb
        if(charCode1 == 'simulate') {
            return(replaceOneSimulate(code))
        }
        
        if(charCode1 == 'nimCopy'){
        	
            matchCode = match.call(nimCopy, code)
            logProb = FALSE
            if('logProb' %in% names(matchCode) ) 
                logProb = matchCode[['logProb']]
            
            logProbChar = as.character(logProb)				#CHANGE1 (added)
     
            toName = as.character(matchCode[['to']])
            fromName = as.character(matchCode[["from"]])
            fromType = NA
            if(inherits(setupSymTab$symbols[[fromName]], "symbolModel"))
                fromType = "model"
            else if(inherits(setupSymTab$symbols[[fromName]], "symbolModelValues") )
                fromType = "modelValues"        	
            if('nodes' %in% names(matchCode) ){
                nodesFromCode = matchCode[["nodes"]]
                
#               nodesFromName = Rname2CppName(as.character(nodesFromCode) )
                nodesFromName = Rname2CppName(  paste0(as.character(nodesFromCode), "_LP_", logProbChar ) )		#CHANGE2

                }
            else{
            	
#                nodesFromName = paste('All', fromName, 'Nodes', sep = '_')
                nodesFromName = paste('All', paste0(fromName, '_LP_', logProbChar), 'Nodes', sep = '_')			#CHANGE3			

                nodesFromCode = as.name(nodesFromName)
                if(!(nodesFromName) %in% newSetupOutputNames){
                    newSetupOutputNames <<- c(newSetupOutputNames, nodesFromName)
                    newNodesFromCode = substitute(NFN <- allNodeNames(OBJECT), list(NFN = as.name(nodesFromName), OBJECT = as.name(fromName) ) )
                    newSetupCode[[nodesFromName]] <<- newNodesFromCode
                }
            }
            if("nodesTo" %in% names(matchCode) ){ 
                nodesToCode = matchCode[["nodesTo"]]
                
 #           	nodesToName = Rname2CppName(as.character(nodesToCode) )
             	nodesToName = Rname2CppName(	paste0(as.character(nodesToCode), '_LP_', logProbChar ) )			#CHANGE4

            }
            else{
                nodesToName = nodesFromName
                nodesToCode = nodesFromCode
                if(nodesFromName == paste('All', fromName, 'Nodes', sep = '_') ){
                    nodesToName = paste('All', toName, 'Nodes', sep = '_')
                    
                    
#                    if(!(nodesToName) %in% newSetupOutputNames){
#                        newSetupOutputNames <<- c(newSetupOutputNames, nodesToName)
 
                 if(nodesToName == paste('All', paste0(fromName, '_LP_', logProbChar), 'Nodes', sep = '_')){		#CHANGE5
                    nodesToName = paste('All', paste0(toName, '_LP_', logProbChar), 'Nodes', sep = '_')				#CHANGE6

 
                        newNodesToCode  = substitute(NTN <- allNodeNames(OBJECT), list(NTN = as.name(nodesToName), OBJECT = as.name(toName) ) )
                        newSetupCode[[nodesToName]] <<- newNodesToCode
                    }
                }
            }
            toType = NA
            if(inherits(setupSymTab$symbols[[toName]], "symbolModel"))
                toType = "model"
            else if(inherits(setupSymTab$symbols[[toName]], "symbolModelValues"))
                toType = "modelValues"
            toAccessName = paste(toName, nodesToName, "access", sep = "_") 
            if(!(toAccessName) %in% newSetupOutputNames)  {
                newSetupOutputNames <<- c(newSetupOutputNames, toAccessName)
                if(toType == "model")
                    newSetupLine = substitute(ACCESSNAME <- modelVariableAccessorVector(MODEL, NODES, logProb = LOGPROB), list(ACCESSNAME = as.name(toAccessName), MODEL = as.name(toName), NODES = nodesToCode, LOGPROB = logProb ) ) 
                if(toType == "modelValues")
                    newSetupLine = substitute(ACCESSNAME <- modelValuesAccessorVector(MODEL, NODES, logProb = LOGPROB), list(ACCESSNAME = as.name(toAccessName), MODEL = as.name(toName), NODES = nodesToCode, LOGPROB = logProb ) ) 
                newSetupCode[[toAccessName]] <<- newSetupLine	
            }
            fromAccessName = paste(fromName, nodesFromName, "access", sep = "_") 
            if(!(fromAccessName) %in% newSetupOutputNames){
                newSetupOutputNames <<- c(newSetupOutputNames, fromAccessName)
                if(fromType == "model"){
                    newSetupLine = substitute(ACCESSNAME <- modelVariableAccessorVector(MODEL, NODES, logProb = LOGPROB), list(ACCESSNAME = as.name(fromAccessName), MODEL = as.name(fromName), NODES = nodesFromCode, LOGPROB = logProb ) ) 
                    }
                if(fromType == "modelValues")
                    newSetupLine = substitute(ACCESSNAME <- modelValuesAccessorVector(MODEL, NODES, logProb = LOGPROB), list(ACCESSNAME = as.name(fromAccessName), MODEL = as.name(fromName), NODES = nodesFromCode, LOGPROB = logProb ) )
			    newSetupCode[[fromAccessName]] <<- newSetupLine
            }  
            if(toType == "model" & fromType == "model")
                newRunCode <- substitute(nimCopy(ACCESSFROM, ACCESSTO), list(ACCESSFROM = as.name(fromAccessName), ACCESSTO = as.name(toAccessName) ) )
            else if(toType == "model" & fromType == "modelValues")
                {
                    rowFrom = matchCode[["row"]]	
                    newRunCode <- substitute(nimCopy(ACCESSFROM, ROWFROM, ACCESSTO), list(ACCESSFROM = as.name(fromAccessName), ROWFROM = rowFrom, ACCESSTO = as.name(toAccessName) ) ) 		
                }
            else if(toType == "modelValues" & fromType == "model"){
                if( "row" %in% names(matchCode) ) 
                    rowTo = matchCode[["row"]]
                if("rowTo" %in% names(matchCode) ) 
                    rowTo = matchCode[["rowTo"]]
                newRunCode <- substitute(nimCopy(ACCESSFROM, ACCESSTO, ROWTO), list(ACCESSFROM = as.name(fromAccessName), ACCESSTO = as.name(toAccessName), ROWTO = rowTo ) )
            }	
            else if(toType == "modelValues" & fromType == "modelValues"){
                rowFrom = matchCode[["row"]]
                if( "row" %in% names(matchCode) ) 
                    rowTo = matchCode[["row"]]
                if("rowTo" %in% names(matchCode) ) 
                    rowTo = matchCode[["rowTo"]]
                newRunCode <- substitute(nimCopy(ACCESSFROM,ROWFROM, ACCESSTO, ROWTO), list(ACCESSFROM = as.name(fromAccessName), ACCESSTO = as.name(toAccessName), ROWTO = rowTo, ROWFROM = rowFrom ) )
            }	
            return(newRunCode)
    	}

        ## The last two pieces of processing here should happen in nfCodeProc because they do not generate newSetupCode
        else if(as.character(code[[1]]) %in% c('resize') )
            {
                varName = as.character(code[[2]])
                if(inherits(setupSymTab$symbols[[varName]], 'symbolNumericList')){
                    newCode = code
                    newCode[[1]] = substitute(resizeNoPtr)
                    return(newCode)	
                }
                return(code)
            } ## end resize
        else if(as.character(code[[1]]) %in% c('getsize') )
            {
                varName = as.character(code[[2]])
                if(inherits(setupSymTab$symbols[[varName]], 'symbolNumericList')){
                    newCode = code
                    newCode[[1]] = substitute(size)
                    return(newCode)	
                }
                return(code)
            } ## end getsize
        else if(as.character(code[[1]]) %in% c('setSize') )
            {
                varName = as.character(code[[2]])
                if(inherits(setupSymTab$symbols[[varName]], 'symbolNumericList')){
                    newCode = code
                    newCode[[3]] = substitute(ORG_INDEX - 1, list(ORG_INDEX = code[[3]] ) )
                    return(newCode)
                }
                return(code)
            } ## end setSize
    }
    if(is.call(code)) {
        if(length(code[[1]]) > 1) code[[1]] <- replaceAccessorsOneFunction(code[[1]])
        if(cl >= 2) {
            for(i in 2:cl) {
                code[[i]] <- replaceAccessorsOneFunction(code[[i]])
            }
        }
        return(code)
    }
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

## steps here are similar to makeMapExprFromBrackets, but that uses exprClasses

varAndIndices2mapParts <- function(varAndIndices, varInfo) {
    varName <- varAndIndices$name
    indices <- varAndIndices$indices
    ## put together offsetExpr, sizeExprs, strideExprs
    ## need sizes to get strides
    sizes <- if(length(varInfo$maxs) > 0) varInfo$maxs else 1 ## would be wierd to be mapping into something with length 1 anyway
    if(varInfo$nDim > 0 & length(indices)==0) { ## A case like model[[node]] where node == 'x', and we should treat like 'x[,]', e.g.
        nDim <- varInfo$nDim
        blockBool <- rep(TRUE, nDim)
        firstIndexRexprs <- rep(list(1), nDim)
        targetSizes <- sizes
    } else {
        nDim <- length(indices)
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

makeSingleIndexAccessorSetupCode <- function(baseName, modelVarName, nodeVarName) {
    codeNames <- makeSingleIndexAccessCodeNames(baseName)
    newOutputNames <- unlist(codeNames)
    namesToSub <- c('varAndIndices', 'newVarName') 
    namesSubList <- lapply(paste0(baseName, '_', namesToSub), as.name)
    names(namesSubList) <- namesToSub
    namesSubList <- c(namesSubList, lapply(codeNames, as.name))
    namesSubList[['modelVarExpr']] <- as.name(modelVarName)
    namesSubList[['nodeVarName']] <- nodeVarName
    namesSubList[['varAccessor']] <- as.name(baseName)
    newSetupCode <- substitute(
        {
            varAndIndices <- getVarAndIndices(nodeVarName)
            newVarName <- as.character(varAndIndices$varName)
            MflatIndex <- varAndIndices2flatIndex(varAndIndices, modelVarExpr$getVarInfo(newVarName))
            varAccessor <- singleVarAccess(modelVarExpr, newVarName, useSingleIndex = TRUE)
        }, namesSubList )
    list(newSetupCode = newSetupCode, newOutputNames = c(newOutputNames, baseName))
}

makeMapSetupCodeExprs <- function(baseName, modelVarName, nodeVarName) {
    codeNames <- makeMapSetupCodeNames(baseName)

    newOutputNames <- unlist(codeNames)
    namesToSub <- c('varAndIndices', 'newVarName') 
    namesSubList <- lapply(paste0(baseName, '_', namesToSub), as.name)
    names(namesSubList) <- namesToSub
    namesSubList <- c(namesSubList, lapply(codeNames, as.name))
    namesSubList[['modelVarExpr']] <- as.name(modelVarName)
    namesSubList[['nodeVarName']] <- nodeVarName
    namesSubList[['varAccessor']] <- as.name(baseName)

    newSetupCode <- substitute(
        {
            varAndIndices <- getVarAndIndices(nodeVarName)
            newVarName <- as.character(varAndIndices$varName)
            mapParts <- varAndIndices2mapParts(varAndIndices, modelVarExpr$getVarInfo(newVarName))
            Mstrides <- mapParts$strides
            Moffset <- mapParts$offset
            Msizes <- mapParts$sizes
            varAccessor <- singleVarAccess(modelVarExpr, newVarName)
        },
        namesSubList )
    list(newSetupCode = newSetupCode, newOutputNames = c(newOutputNames, baseName))
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

makeNameFromMatchedCall <- function(arg, matchCode, default){
	output = default
	if( as.character(arg) %in% names(matchCode) ) 
		output = as.character(matchCode[[as.character(arg) ]]) 
	return(output)
}

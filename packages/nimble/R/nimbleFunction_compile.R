virtualNFprocessing <- setRefClass('virtualNFprocessing',
    fields = list(
        name = 'ANY',                    ## character
        setupSymTab = 'ANY',             ## symbolTable
        nfGenerator =  'ANY',		## 'function',
        compileInfos =  'ANY',		## list of RCfunctionCompileClass objects
        origMethods = 'ANY',             ## list of original methods
        RCfunProcs =  'ANY',		## list of RCfunProcessing  or RCvirtualFunProcessing objects
        nimbleProject = 'ANY',           ## nimbleProjectclass object
        cppDef = 'ANY',                  ## cppNimbleFunctionClass or cppVirtualNimbleFunctionClass object
        isNode = 'ANY'                  ## logical, is it a nodeFunction?
    ),
    methods = list(
        show = function() {
            writeLines(paste0('virtualNFprocessing object ', name))
        },
        initialize = function(f = NULL, className, virtual = TRUE, isNode = FALSE, project = NULL) {
            nimbleProject <<- project
            compileInfos <<- list()
            RCfunProcs <<- list()
            
            isNode <<- isNode
    
            if(!is.null(f)) { ## This allows successful default instantiation by R when defining nfProcessing below -- crazy.
                ## nfGenerator is allowed if it is a nimbleFunctionVirtual.
                if(is.nf(f) | is.nfGenerator(f)) nfGenerator <<- nf_getGeneratorFunction(f)
                else if(inherits(f, 'list')) {
                    if(length(unique(lapply(f, nfGetDefVar, 'name'))) != 1)
                        stop('Error with list of instances not having same nfGenerator')
                    nfGenerator <<- nf_getGeneratorFunction(f[[1]])
                }
                if(missing(className)) {
                    sf <- environment(nfGenerator)$name
                    name <<- Rname2CppName(sf)
                } else {
                    name <<- className
                }
                origMethods <<- nf_getMethodList(nfGenerator)
                RCfunProcs <<- list()
                for(i in seq_along(origMethods)) {
                    RCname <- names(origMethods)[i]
                    if(isNode && strsplit(RCname, '_', fixed = TRUE)[[1]][1] == getCalcADFunName()) constFlag <- FALSE
                    else constFlag <- isNode
                    RCfunProcs[[RCname]] <<- if(virtual) RCvirtualFunProcessing$new(origMethods[[i]], RCname, const = constFlag) else RCfunProcessing$new(origMethods[[i]], RCname, const = constFlag)
                }
                compileInfos <<- lapply(RCfunProcs,
                                        function(x) x$compileInfo)
            }
         },
        setupLocalSymbolTables = function() {
            for(i in seq_along(RCfunProcs)) {
                RCfunProcs[[i]]$setupSymbolTables(parentST = setupSymTab, neededTypes = list(), nimbleProject = nimbleProject)
            }
        },
        doRCfunProcess = function(control = list(debug = FALSE, debugCpp = FALSE)) {
            for(i in seq_along(RCfunProcs)) {
                RCfunProcs[[i]]$process(debug = control$debug, debugCpp = control$debugCpp, debugCppLabel = name, doKeywords = FALSE, nimbleProject = nimbleProject)
            }
        },
        addMemberFunctionsToSymbolTable = function() {
            for(i in seq_along(origMethods)) {
                thisName <- names(origMethods)[i]
                newSym <- symbolMemberFunction(name = thisName, nfMethodRCobj = origMethods[[i]], RCfunProc = RCfunProcs[[i]])
                setupSymTab$addSymbol(newSym)
            }
        },
        process = function(control = list(debug = FALSE, debugCpp = FALSE)) {
            setupSymTab <<- symbolTable(parentST = NULL)
            addMemberFunctionsToSymbolTable()
            setupLocalSymbolTables()
            doRCfunProcess(control)
        }
    )
)


nfProcessing <- setRefClass('nfProcessing',
    contains = 'virtualNFprocessing',
    fields = list(
        instances = 'ANY',           ## list of instances of the nimbleFunction to used for setup types and receive newSetupCode
        neededTypes =  'ANY',	     ## list of symbols for non-trivial types that will be needed for compilation, such as derived models or modelValues
        neededObjectNames =  'ANY',     ## character vector of the names of objects such as models or modelValues that need to exist during C++ instantiation and population so their contents can be pointed to 
        newSetupOutputNames =  'ANY', ## character vector of names of objects created by newSetupCode from "keyword processing"
        blockFromCppNames = 'ANY',    ## character vector of names of setup outputs that should not be propagated to C++
        newSetupCode =  'ANY',	      ## list of lines of setup code populated by keyword processing
        newSetupCodeOneExpr = 'ANY',  ## all lines of new setup code put into one expression for evaluation
        inModel =  'ANY'	      ## logical: whether this nfProcessing object is for a nodeFunction in a model
    ),
    methods = list(
        show = function() {
            writeLines(paste0('nfProcessing object ', name))
        },
        initialize = function(f = NULL, className, fromModel = FALSE, project, isNode = FALSE) {
            neededTypes <<- list()
            neededObjectNames <<- character()
            newSetupCode <<- list()
            if(!is.null(f)) {
                ## f must be a specialized nf, or a list of them
                inModel <<- fromModel
                if(missing(className)) {
                    sf <- if(is.list(f)) nfGetDefVar(f[[1]], 'name') else nfGetDefVar(f, 'name')
                    name <<- Rname2CppName(sf)
                } else {
                    name <<- className
                }
                callSuper(f, name, virtual = FALSE, isNode = isNode, project = project)
                instances <<- if(inherits(f, 'list')) lapply(f, nf_getRefClassObject) else list(nf_getRefClassObject(f))
               
                newSetupOutputNames <<- character()
                blockFromCppNames <<- character()
                newSetupCode <<- list()
            }
        },

        getSymbolTable = function() setupSymTab,
        getMethodInterfaces = function() origMethods,
        
        processKeywords_all = function(){},
        matchKeywords_all = function(){},
        
        doSetupTypeInference_processNF = function() {},
        makeTypeObject = function() {},
        replaceCall = function() {},
        evalNewSetupLines = function(){},
        makeNewSetupLinesOneExpr = function() {},
        evalNewSetupLinesOneInstance = function(instances, check = FALSE){},
        setupTypesForUsingFunction = function(){},
        doSetupTypeInference = function(){},
        clearSetupOutputs = function() {},
        
        setupLocalSymbolTables = function() {
            for(i in seq_along(RCfunProcs)) {
                RCfunProcs[[i]]$setupSymbolTables(parentST = setupSymTab, neededTypes = neededTypes, nimbleProject = nimbleProject)
            }
        },
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
        collect_nimDerivsCalculate_info = function() {
            newEnableDerivs <- list()
            for(i in seq_along(RCfunProcs)) {
                numNimDerivsCalculate <- RCfunProcs[[i]]$compileInfo$typeEnv[['numNimDerivsCalculate']]
                if(!is.null(numNimDerivsCalculate)) {
                    methodName <- RCfunProcs[[i]]$name
                    if(is.character(methodName)) ## not sure when it wouldn't be; this is defensive
                        newEnableDerivs[[ methodName ]] <- list(numNimDerivsCalculate = numNimDerivsCalculate)
                }
            }
            if(length(newEnableDerivs)) {
                environment(nfGenerator)$enableDerivs <<- c(newEnableDerivs,
                                                            environment(nfGenerator)$enableDerivs)
            }
        },
        addBaseClassTypes = function() {
            ## If this class has a virtual base class, we add it to the needed types here
            contains <- environment(nfGenerator)$contains
            if(!is.null(contains)) {
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
            if(!is.null(nimbleOptions()$debugNFProcessing)) {
                if(nimbleOptions()$debugNFProcessing) {
                    debug <- TRUE
                    control$debug <- TRUE
                    writeLines('Debugging nfProcessing (nimbleOptions()$debugRCfunProcessing is set to TRUE)') 
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
                print('lapply(compileInfos, function(x) print(x$newLocalSymTab))')
                lapply(compileInfos, function(x) print(x$newLocalSymTab))
                writeLines('**** READY FOR RFfunProcessing *****')
                browser()
            }

            doRCfunProcess(control)

            collectRCfunNeededTypes()

            if(isTRUE(nimbleOptions('experimentalEnableDerivs'))) {
                collect_nimDerivsCalculate_info()
            }
            
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

nfProcessing$methods(clearSetupOutputs = function(inst) {
    for(i in nf_getSetupOutputNames(nfGenerator)) {
        inst[[i]] <- NULL
    }
    for(i in newSetupOutputNames) {
        inst[[i]] <- NULL
    }
    NULL
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
    ## Warning: this relies on the fact that although refClass environments are closed, we can
    ## eval in them and create new variables in them that way.
    eval(newSetupCodeOneExpr, envir = instance)
    instance$.newSetupLinesProcessed <- TRUE
})

nfProcessing$methods(setupTypesForUsingFunction = function() {
    if(inherits(setupSymTab, 'uninitializedField')) {
        doSetupTypeInference(TRUE, FALSE)
        addMemberFunctionsToSymbolTable()
        addBaseClassTypes()						
        matchKeywords_all()
        processKeywords_all()
        setupLocalSymbolTables()
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
        setupSymTab$addSymbol(symbolNimbleFunctionSelf(name = ".self",
                                                       nfProc = .self) )
    	outputNames <- c(outputNames, nf_getSetupOutputNames(nfGenerator))
    	if(length(outputNames)>0) outputNames <- unique(outputNames)
    }
    if(setupNew) {
        ## Kluge that results from adding string handling to the compiler:
        ## Previously any character objects were assigned a symbol object with
        ## type 'Ronly'.  In later processing all 'Ronly' types are filtered out of
        ## propagation to C++.
        ## Now that we have added string handling, character objects are assigned
        ## a symbolString symbol with type "character" and not automatically filtered.
        ## Unfortunately this means that vectors of node names that are only used
        ## in lines like calculate(model, nodeNames), which undergoes keyword processing
        ## would be propogated to C++ wastefully.
        ## As a kluge, we will step in here, during second round of setup type inference
        ## to re-assign type 'Ronly' to any symbols that, as a result of
        ## keyword processing, we can now see are not needed
        ## We also need the section added below to filter out newSetupOutputs
        ## that are really created as intermediates for others that are really needed
        ## during the keyword processing, the newSetupOutputNames is used for
        ## bookkeeping, so it would not be trivial to remove them at an earlier stage.
        
        origSetupOutputs <- nf_getSetupOutputNames(nfGenerator)
        declaredSetupOutputs <- getFunctionEnvVar(nfGenerator, 'declaredSetupOutputNames')
        origSetupOutputs <- setdiff(origSetupOutputs, declaredSetupOutputs)
        newRcodeList <- lapply(compileInfos, `[[`, 'newRcode')
        allNamesInCodeAfterKeywordProcessing <- unique(unlist(lapply(newRcodeList, all.names)))
        origSetupOutputNamesToKeep <- intersect(allNamesInCodeAfterKeywordProcessing, origSetupOutputs) ## this loses mv!
        origSetupOutputNamesNotNeeded <- setdiff(origSetupOutputs,origSetupOutputNamesToKeep) ## order matters
        for(nameNotNeeded in origSetupOutputNamesNotNeeded) {
            thisSym <- setupSymTab$getSymbolObject(nameNotNeeded)
            if(!is.null(thisSym))  if(!thisSym$type == 'Values') thisSym$type <- 'Ronly' ## must keep modelValues, nimbleFunctions, possibly others
        }

        outputNames <- c(outputNames, newSetupOutputNames)
    }
    doSetupTypeInference_processNF(setupSymTab, outputNames, instances, add = TRUE)   # add info about each setupOutput to symTab

    if(setupNew) {
        ## This is the second part of the kluge.
        ## Probably it would be ok to never add these to the symbol table in the first place
        ## but right now I am doing it this way to minimize unforeseen consequences by more closely mimicing what would have been created prior to adding string support
        ## This is trickier because keyword processing can create objects for propogation to C++ that never appear in method code (e.g. manyVariableAccessors used to construct copierVectors)
        ## So I added a blockFromCppNames that is populated during keyword processing
        
        for(nameNotNeeded in blockFromCppNames) {
            thisSym <- setupSymTab$getSymbolObject(nameNotNeeded)
            if(!is.null(thisSym)) thisSym$type <- 'Ronly'
        }
    }
})

nfProcessing$methods(doSetupTypeInference_processNF = function(symTab, setupOutputNames, instances, add = FALSE, firstOnly = FALSE) {
    if(length(instances) == 0) {
        warning('Can not infer setup output types with no instances.')
        return(invisible(NULL))
    }
    for(name in setupOutputNames) {
        symbolRCobject <- makeTypeObject(name, instances, firstOnly)
        if(is.null(symbolRCobject)) next
        if(is.logical(symbolRCobject)) {
            stop(paste0('There is an error involving the type of ', name,'.'), call. = FALSE)
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
    makeTypeObj_impl(.self, name, instances, firstOnly)
})

makeTypeObj_impl <- function(.self, name, instances, firstOnly) {
  isNLG <- FALSE
  if(is.nlGenerator(instances[[1]][[name]])){
    nlGen <- instances[[1]][[name]] 
    isNLG <- TRUE
  } else if(exists(name, envir = globalenv())) {
      foundObject <- get(name, envir = globalenv())
      if(is.nlGenerator(foundObject)) {
          nlGen <- foundObject 
          isNLG <- TRUE
      }
  }
  if(isNLG){
    nlp <- .self$nimbleProject$compileNimbleList(nlGen, initialTypeInferenceOnly = TRUE)
    className <- nl.getListDef(nlGen)$className
    newSym <- symbolNimbleList(name = name, nlProc = nlp)
    .self$neededTypes[[className]] <- newSym  ## if returnType is a NLG, this will ensure that it can be found in argType2symbol()
    returnSym <- symbolNimbleListGenerator(name = name, nlProc = nlp)
    return(returnSym)
  }
  if(is.nl(instances[[1]][[name]])) {
    ## This case mimics the nimbleFunction case below (see is.nf)

    ## We need all instances created in setup code from all instances
    nlList <- lapply(instances, `[[`, name)
    ## trigger initial procesing to set up an nlProc object
    ## that will have a symbol table.
    ## Issue: We may also need to trigger this step from run code
    nlp <- .self$nimbleProject$compileNimbleList(nlList, initialTypeInferenceOnly = TRUE)
    ## get the unique name that we use to generate a unique C++ definition
    className <- nlList[[1]]$nimbleListDef$className
    ## add the setupOutput name to objects that we need to instantiate and point to
    .self$neededObjectNames <- c(.self$neededObjectNames, name)
    
    ## create a symbol table object
    newSym <- symbolNimbleList(name = name, nlProc = nlp)
    
    ## If this is the first time this type is encountered,
    ## add it to the list of types whose C++ definitions will need to be generated
    if(!(className %in% names(.self$neededTypes))) .self$neededTypes[[className]] <- newSym
    return(newSym)
  }
  if(inherits(instances[[1]][[name]], 'indexedNodeInfoTableClass')) {
      return(symbolIndexedNodeInfoTable(name = name, type = 'symbolIndexedNodeInfoTable')) ## the class type will get it copied but the Ronly will make it skip a type declaration, which is good since it is in the nodeFun base class.
  }
  if(inherits(instances[[1]][[name]], 'nimbleFunctionList')) {
      
      .self$neededObjectNames <- c(.self$neededObjectNames, name)
      baseClass <- instances[[1]][[name]]$baseClass ## an nfGenerator created by virtualNimbleFunction()
      baseClassName <- environment(baseClass)$className
      
      if(!(baseClassName %in% names(.self$neededTypes))) {
          nfp <- .self$nimbleProject$setupVirtualNimbleFunction(baseClass, fromModel = .self$inModel)
          newSym <- symbolNimbleFunctionList(name = name, type = 'nimbleFunctionList', baseClass = baseClass, nfProc = nfp)
          neededTypeSim <- symbolNimbleFunction(name = baseClassName, type = 'virtualNimbleFunction', nfProc = nfp)
          .self$neededTypes[[baseClassName]] <- newSym
      } else {
          newSym <- .self$neededTypes[[baseClassName]]
      }
    
      allInstances <- unlist(lapply(instances, function(x) x[[name]]$contentsList), recursive = FALSE)
      newNFprocs <- .self$nimbleProject$compileNimbleFunctionMulti(allInstances, initialTypeInference = TRUE)
      ## only types are needed here, not initialTypeInference, because nfVar's from a nimbleFunctionList are not available (could be in future)
      for(nfp in newNFprocs) {
          newTypeName <- environment(nfp$nfGenerator)$name
          .self$neededTypes[[ newTypeName ]] <- symbolNimbleFunction(name = newTypeName, type = 'nimbleFunction',
                                                                nfProc = nfp)
      }
    return(newSym)
  }
  if(is.nf(instances[[1]][[name]])) { ## nimbleFunction
    funList <- lapply(instances, `[[`, name)
    nfp <- .self$nimbleProject$compileNimbleFunction(funList, initialTypeInferenceOnly = TRUE) ## will return existing nfProc if it exists
    className <- class(nf_getRefClassObject(funList[[1]]))
    .self$neededObjectNames <- c(.self$neededObjectNames, name)
    newSym <- symbolNimbleFunction(name = name, type = 'nimbleFunction', nfProc = nfp)
    if(!(className %in% names(.self$neededTypes))) .self$neededTypes[[className]] <- newSym
    
    return(newSym)
  }
  if(inherits(instances[[1]][[name]], 'modelValuesBaseClass')) { ## In some cases these could be different derived classes.  If locally defined they must be the same
    if(!firstOnly) {
      if(!all(unlist(lapply(instances, function(x) inherits(x[[name]], 'modelValuesBaseClass'))))) {
        warning(paste0('Problem: some but not all instances have ', name,' as a modelValues.  Types must be consistent.'))
        return(invisible(NULL))
      }
    }
    ## Generate one set of symbolModelValues objects for the neededTypes, and each of these can have its own mvConf
    ## Generate another symbolModelValues to return and have in the symTab for this compilation
    ## I don't think that mvConf gets used, since they all get Values *        
    for(i in seq_along(instances)) {
      className <- class(instances[[i]][[name]])
      if(!(className %in% names(.self$neededTypes))) {
        ## these are used only to build neededTypes
        ntSym <- symbolModelValues(name = name, type = 'Values', mvConf = instances[[i]][[name]]$mvConf)
        .self$neededTypes[[className]] <- ntSym
      }
    }
    ## this is used in the symbol table
    .self$neededObjectNames <- c(.self$neededObjectNames, name)
    newSym <- symbolModelValues(name = name, type = 'Values', mvConf = NULL)
    return(newSym)
  }
  if(inherits(instances[[1]][[name]], 'modelBaseClass')) {
    if(!firstOnly) {
      if(!all(unlist(lapply(instances, function(x) inherits(x[[name]], 'modelBaseClass'))))) {
        warning(paste0('Problem: some but not all instances have ', name,' as a model.  Types must be consistent.'))
        return(invisible(NULL))
      }
      if(!all(unlist(lapply(instances, function(x) inherits(x[[name]], 'RmodelBaseClass'))))) {
        warning(paste0('Problem: models should be provided as R model objects, not C model objects'))
        return(invisible(NULL))
      }
    }
    return(symbolModel(name = name, type = 'Ronly', className = class(instances[[1]][[name]]))) 
  }
  if(isTRUE(getNimbleOption('experimentalEnableDerivs'))) {
      if(inherits(instances[[1]][[name]], 'ADproxyModelClass')) {
          return(symbolModel(name = name, type = 'Ronly', className = class(instances[[1]][[name]]$model))) 
      }
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
  
  if(inherits(instances[[1]][[name]], 'copierVectorClass')) {
    newSym <- symbolCopierVector(name = name, type = 'symbolCopierVector')
    return(newSym)
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
  if(inherits(instances[[1]][[name]], 'nodeFunctionVector_nimDerivs')) { 
      return(symbolNodeFunctionVector_nimDerivs(name = name))
  }
  if(inherits(instances[[1]][[name]], 'modelVariableAccessorVector')){
    return(symbolModelVariableAccessorVector(name = name, lengthName = paste0(name, '_length')) )
  }
  if(inherits(instances[[1]][[name]], 'modelValuesAccessorVector')){
    return(symbolModelValuesAccessorVector(name = name) )     	
  }
  if(inherits(instances[[1]][[name]], 'getParam_info')) { ## the paramInfo in an instance is allowed to be NULL (see GitHub Issue #327). Hence we search for the first valid case and default to double()
      iInst <- 1
      paramInfo <- instances[[iInst]][[name]]
      while(is.na(paramInfo$type) & iInst < length(instances)) {
          iInst <- iInst + 1
          paramInfo <- instances[[iInst]][[name]]
      }
      if(is.na(paramInfo$type)) paramInfo <- defaultParamInfo()
      return(symbolGetParamInfo(name = name, paramInfo = paramInfo))
  }
  if(inherits(instances[[1]][[name]], 'getBound_info')) {
    return(symbolGetBoundInfo(name = name, boundInfo = instances[[1]][[name]]))
  }
  if(is.character(instances[[1]][[name]])) {
    if(firstOnly) {
      nDim <- if(is.null(dim(instances[[1]][[name]]))) 1L else length(dim(instances[[1]][[name]]))
      if(nDim > 1) {
        warning('character object with nDim > 1 being handled as a vector')
        nDim <- 1
      }
      size <- if(length(instances[[1]][[name]])==1) 1L else as.numeric(NA)
      if(nimbleOptions()$convertSingleVectorsToScalarsInSetupArgs) {
        if(nDim == 1 & identical(as.integer(size), 1L)) nDim <- 0
      }
      return(symbolString(name = name, type = 'character', nDim = nDim, size = size))
    } else {
      instanceObjs <- lapply(instances, `[[`, name)
      types <- unlist(lapply(instanceObjs, storage.mode))
      if(!all(types == 'character')) stop(paste('Inconsistent types for setup variable', name))
      dims <- lapply(instanceObjs, dim)
      dimsNULL <- unlist(lapply(dims, is.null))
      if(any(dimsNULL)) { ## dimsNULL TRUE means it is a vector
        if(!all(dimsNULL)) {
          warning(paste0('Dimensions do no all match for ', name, 'but they will be treated as all vectors anyway.'))
        }
      }
      nDim <- 1
      lengths <- unlist(lapply(instanceObjs, length))
      size <- if(!all(lengths == 1)) as.numeric(NA) else 1L
      if(nimbleOptions()$convertSingleVectorsToScalarsInSetupArgs) {
        if(nDim == 1 & identical(as.integer(size), 1L)) nDim <- 0
      }
      return(symbolString(name = name, type = 'character', nDim = nDim, size = size))
    }
  }
  if(is.numeric(instances[[1]][[name]]) | is.logical(instances[[1]][[name]])) {
    if(firstOnly) {
      type <- storage.mode(instances[[1]][[name]])
      nDim <- if(is.null(dim(instances[[1]][[name]]))) 1L else length(dim(instances[[1]][[name]]))
      size <- if(length(instances[[1]][[name]])==1) rep(1L, nDim) else rep(as.numeric(NA), nDim)
      if(nimbleOptions()$convertSingleVectorsToScalarsInSetupArgs) {
        if(nDim == 1 & identical(as.integer(size), 1L)) nDim <- 0
      }
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
        if(nimbleOptions()$convertSingleVectorsToScalarsInSetupArgs) {
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




nfProcessing$methods(processKeywords_all = function(){	
    for(i in seq_along(compileInfos)){
        RCfunProcs[[i]]$processKeywords(.self)
    }
})

nfProcessing$methods(matchKeywords_all = function(){
	for(i in seq_along(compileInfos))
            RCfunProcs[[i]]$matchKeywords(.self)
})


#' Class \code{singleVarAccessClass}
#' @aliases singleVarAccessClass
#' @export
#' @description
#' Classes used internally in NIMBLE and not expected to be called directly by users.
singleVarAccessClass <- setRefClass('singleVarAccessClass',
                                       methods = list(
                                           initialize = function() cat('Oops: building a singleVarAccessClass refClass -- should be defunct\n')
                                       ))

singleVarAccess <- function(model, var, useSingleIndex = FALSE) {
    ans <- list(model = model, var = var, useSingleIndex = useSingleIndex)
    class(ans) <- 'singleVarAccessClass'
    ans
}


# singleModelValuesAccessClass and singleModelValuesAccess are exported (and 'documented' in nimble-internals.Rd) based on prep_pkg; doing it here causes R CMD check issues with argument names
singleModelValuesAccessClass <- setRefClass('singleModelValuesAccessClass',
                                     methods = list(
                                           initialize = function() cat('Oops: building a singleModelValuesAccessClass refClass -- should be defunct\n')
                                     ))

singleModelValuesAccess <- function(modelValues, var) {
    ans <- list(modelValues = modelValues, var = var)
    class(ans) <- 'singleModelValuesAccessClass'
    ans
}




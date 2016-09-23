
###		buildModelInterface(refName, symbolTable, basePtrCall)
###					Builds a newreference class around a function which returns a
###					pointer to a Model.
###					refname will be the name of the new class
###
###					symbolTable will be a used to determine element types
###					basePtrCall is the C++ call which builds the C++ object
###					i.e. in R .Call(basePtrCall) should return the pointer
###					to our NimbleFunction Class
###					IMPORTANT NOTE: symbolTable MUST match form object pointed to by
###					basePtrCall! In particular, if the symbolTable says an element is a
###					scalar (i.e. nDims = 0) and it's really a NimArr (or vice versa), 
###					attempting to access this element will cause R to crash!

###					Once the new class is built, a new object is create via refname$new()
###					Access to the elements is provided through nimbleModel$Varname
###					Access to the modelValue elements is provided through nimbleModel$modelValues$Varname
###					In this interface, it is assumed that the pointers associated with nimbleModel$Varname
###					are double pointers, i.e. pointers to pointers to NimArr's
###					IMPORTANT NOTE: the function which builds the C++ object (i.e. basePtrCall)
###					MUST initialize the pointers correctly, or accessing will cause R to crash!
###					Default initialization is assumed that elements are pointing to row 1 of
###					the modelValues, but this is not necessary (as long as it is a double pointer to something

getMVptr <- function(rPtr, dll)
  .Call( nimbleUserNamespace$sessionSpecificDll$getModelValuesPtrFromModel, rPtr)
getMVName <- function(modelValuePtr, dll)
  .Call( nimbleUserNamespace$sessionSpecificDll$getMVBuildName, modelValuePtr)

# for now export this as R<3.1.2 give warnings if don't

#' Class \code{CmodelBaseClass}
#' @aliases CmodelBaseClass 
#' @export
#' @description
#' Classes used internally in NIMBLE and not expected to be called directly by users.
CmodelBaseClass <- setRefClass('CmodelBaseClass',
                               contains = 'modelBaseClass',
                               fields = list(
                                   .basePtr = 'ANY',
                                   dll = 'ANY',
                                   Rmodel = 'ANY',
                                   cppNames = 'ANY',
                                   cppCopyTypes = 'ANY', ## At the given moment these will all be 'numeric', but the system allows more flexibility
                                   ##CnodeFunClasses = 'list',
                                   compiledModel = 'ANY',
                                   ##.nodeFxnPointers_byGID = 'ANY',
                                   .nodeFxnPointers_byDeclID = 'ANY', ## Added for newNodeFxns
                                   ##.nodeValPointers_byGID = 'ANY',
                                   ##.nodeLogProbPointers_byGID = 'ANY',
                                   nodeFunctions = 'ANY' ## Added for newNodeFxns, so we can access nodeFunctions by declID. Could be migrated up to modelBaseClass.
                                   ),
                               methods = list(
                                   show = function() {
                                       cat('CmodelBaseClass object\n')
                                   },
                                   finalize = function() {
                                     ##  Rmodel$CobjectInterface <<- NULL
                                       for(vn in cppNames) {
                                           vPtrName <- paste(".", vn, "_Ptr", sep = "")
                                           assign(vPtrName, NULL, inherits = TRUE)
                                       }
                                       browser()
                                       for(i in ls(Rmodel$nodes)) {
                                           if(is.list(nodes[[thisNodeFunName]]))
                                               nodes[[thisNodeFunName]][[1]]$finalize(nodes[[thisNodeFunName]][[2]])
                                           else
                                               nodes[[thisNodeFunName]]$finalize()
                                       }
                                       .basePtr <<- NULL
                                       .nodeFxnPointers_byDeclID$finalize()
                                       nimbleProject <<- NULL
                                   },
                                   setModel = function(model) {
                                       ## This is creating a circular reference, so be careful with show(), and with Rstudio
                                       Rmodel <<- model 
                                       model$CobjectInterface <- .self
                                   },
                                   copyFromModel = function(model) { ## Could potentially be generated customized for each model
                                       if(missing(model)) model <- Rmodel
                                       for(v in cppNames) {
                                           if(cppCopyTypes[[v]] == 'numeric') {
                                               .self[[v]] <<- model[[v]]
                                               next
                                           }
                                       }
                                   },
                                   setupNodes = function(where = classEnvironment, dll = NULL) {
                                       ## assumes we have the model set
                                       ## there is a field nodes in the modelBaseClass
                                       ##
                                       ## 1. generate CnodeFunClasses
                                       ##     - by iterating through the nodeGenerators in the Rmodel
                                       nodesEnv <- new.env()
                                       asTopLevel <- getNimbleOption('buildInterfacesForCompiledNestedNimbleFunctions')
                                       nodeFunctions <<- vector('list', length(Rmodel$nodeFunctions))
                                       for(i in seq_along(Rmodel$nodeFunctions)) {
                                           thisNodeFunName <- names(Rmodel$nodeFunctions)[i]
                                           nodesEnv[[thisNodeFunName]] <- nimbleProject$instantiateNimbleFunction(Rmodel$nodes[[thisNodeFunName]], dll = dll, asTopLevel = asTopLevel)
                                           nodeFunctions[[i]] <<- nodesEnv[[thisNodeFunName]]
                                       }
                                       nodes <<- nodesEnv
                                       names(nodeFunctions) <<- names(Rmodel$nodeFunctions)
                                       
                                       ##.nodeFxnPointers_byGID <<- new('numberedObjects')
                                       .nodeFxnPointers_byDeclID <<- new('numberedObjects', dll = dll) 
                                       ##maxID = length(modelDef$maps$graphIDs)
                                       maxID = length(modelDef$declInfo)
                                       ##.nodeFxnPointers_byGID$resize(maxID)
                                       .nodeFxnPointers_byDeclID$resize(maxID)
                                       ## for(nodeName in ls(nodes)) {
                                       ##     gID <- modelDef$nodeName2GraphIDs(nodeName)
                                       ##     basePtr <- if(is.list(nodes[[nodeName]])) nodes[[nodeName]][[1]]$basePtrList[[ nodes[[nodeName]][[2]] ]] else nodes[[nodeName]]$.basePtr
                                       ##     .self$.nodeFxnPointers_byGID[gID] <- basePtr ## nodes[[nodeName]]$.basePtr
                                       ## }

                                       for(declID in seq_along(nodes)) {
                                           thisNodeFunctionName <- names(Rmodel$nodeFunctions)[declID]
                                           basePtr <- if(is.list(nodes[[thisNodeFunctionName]])) ## it's a multiInterface
                                                          nodes[[thisNodeFunctionName]][[1]]$basePtrList[[ nodes[[thisNodeFunctionName]][[2]] ]]
                                                      else ## it's a direct interface
                                                          nodes[[thisNodeFunctionName]]$.basePtr
                                           .self$.nodeFxnPointers_byDeclID[declID] <- basePtr ## nodes[[nodeName]]$.basePtr
                                       }

                                       ## maxGraphID <- length(modelDef$maps$graphIDs)
                                       ## .nodeValPointers_byGID <<- new('numberedModelVariableAccessors')
                                       ## .nodeValPointers_byGID$resize(maxGraphID)
                                       ## .nodeLogProbPointers_byGID <<- new('numberedModelVariableAccessors')
                                       ## .nodeLogProbPointers_byGID$resize(maxGraphID)
                                
                                       ## for(vName in Rmodel$getVarNames()){
                                       ## 		flatIndices = 1
                                       ## 		if(length(vars[vName]) > 0)
	                               ##         		flatIndices = 1:prod(unlist(vars[vName]))
                                       		
                                       ## 		gIDs_withNAs = unlist(sapply(vName, parseEvalNumeric, env = Rmodel$modelDef$maps$vars2GraphID_values, USE.NAMES = FALSE))
                                       ## 		validIndices = which(!is.na(gIDs_withNAs))
                                       ## 		gIDs = gIDs_withNAs[validIndices] 
											
                                       ## 		.Call('populateNumberedObject_withSingleModelVariablesAccessors', .basePtr , vName, as.integer(gIDs), as.integer(validIndices), .nodeValPointers_byGID$.ptr)
                                       ## 		logVNames <- modelDef$nodeName2LogProbName(vName)
                                   ##    		if(length(logVNames) > 0){
                                  ##     			logVName <- nl_getVarNameFromNodeName(logVNames[1])
                                  ##     			LP_gIDs_withNAs =  unlist(sapply(vName, parseEvalNumeric, env = Rmodel$modelDef$maps$vars2LogProbID, USE.NAMES = FALSE))
	                          ##             		validIndices = which(!is.na(LP_gIDs_withNAs) ) 
	                          ##             		l_gIDs = Rmodel$modelDef$nodeName2LogProbID(vName)
	                           ##            		.Call('populateNumberedObject_withSingleModelVariablesAccessors', .basePtr, logVName, as.integer(l_gIDs), as.integer(validIndices), .nodeLogProbPointers_byGID$.ptr)
                                                    ##	                                       		}
                                   ##}
                                   }
                                   )
                               )

makeModelCppCopyTypes <- function(symTab) {
    ans <- list()
    for(s in symTab$symbols) {
        ans[[s$name]] <- 'numeric' ## only one option right now, in a model
    }
    ans
}



makeModelBindingFields <- function(symTab) {
  fieldList = list(.modelValues_Ptr = "ANY", .DUMMY = "ANY")  
  vNames = names(symTab$symbols)
  for(vn in vNames){
    ptrName = paste(".", vn, "_Ptr", sep = "")
    fieldList[[ptrName]] <- "ANY"
    eval(substitute( fieldList$VARNAME <- function(x){
      if(missing(x) ) 
        nimbleInternalFunctions$getNimValues(VPTR, 2, dll = dll)
      else 
        nimbleInternalFunctions$setNimValues(VPTR, x, 2, allowResize = FALSE, dll = dll)
    }, list(VPTR = as.name(ptrName), VARNAME = vn) ) )
  }
  return(fieldList)
}

buildModelInterface <- function(refName, compiledModel, basePtrCall, project = NULL, dll = NULL, where = globalenv()){
    defaults <- list()
    if(inherits(compiledModel, 'symbolTable')) {
        symTab <- compiledModel
        defaults$cm <- NULL
        warning('No compiled model provided, so interface will be incomplete')
    } else {
        symTab <- compiledModel$model$modelDef$symTab
        defaults$cm <- compiledModel
    }
    defaults$cppCT <- makeModelCppCopyTypes(symTab)
    defaults$project <- project
    
    # resolve basePtrCall rather than looking it up later.
  if(!is.null(dll) && is.character(basePtrCall))
     basePtrCall = getNativeSymbolInfo(basePtrCall, dll)
  else
     warning("creating an initialization method that calls a C routine without any DLL information")

    # add basePtrCall as arg for initialize and substitute on parsed text string to avoid CRAN issues with .Call registration
  eval(substitute(      newClass <-  setRefClass(refName,
                                            fields = FIELDS,
                                            contains = "CmodelBaseClass",
                                            methods = list(initialize = function(model, defaults, basePtrCall, ..., dll = NULL) { ## model is an optional R model
                                                nodes <<- list()
                                                isDataEnv <<- new.env()
                                                classEnvironment <<- new.env()
                                                
                                                callSuper(dll = dll, ...)

                                        # avoid R CMD check problem with registration
                                                ## notice that the following line appears a few lines up:basePtrCall = getNativeSymbolInfo(basePtrCall, dll)
                                                .basePtr <<- eval(parse(text = ".Call(basePtrCall)"))
                                                .Call(nimbleUserNamespace$sessionSpecificDll$register_namedObjects_Finalizer, .basePtr, dll)
                                                # .basePtr <<- .Call(BPTRCALL)
                                                .modelValues_Ptr <<- nimbleInternalFunctions$getMVptr(.basePtr, dll = dll)
                                                defaultModelValues <<- nimbleInternalFunctions$CmodelValues$new(existingPtr = .modelValues_Ptr, buildCall = nimbleInternalFunctions$getMVName(.modelValues_Ptr, dll), initialized = TRUE, dll = dll )
                                                modelDef <<- model$modelDef
                                                graph <<- model$graph
                                                vars <<- model$vars
                                                isDataVars <<- model$isDataVars
                                                nimbleProject <<- defaults$project
                                                for(v in ls(model$isDataEnv)) isDataEnv[[v]] <<- model$isDataEnv[[v]]
                                                setData(modelDef$constantsList, warnAboutMissingNames = FALSE)
                                                cppNames <<- .Call( nimbleUserNamespace$sessionSpecificDll$getAvailableNames, .basePtr) ## or could get this from R objects
                                                cppCopyTypes <<- defaults$cppCT
                                                compiledModel <<- defaults$cm
                                                for(vn in cppNames)
                                                    {
                                                        vPtrName <- paste(".", vn, "_Ptr", sep = "")
                                                     	.self[[vPtrName]] <<- nimbleInternalFunctions$newObjElementPtr(.basePtr, vn, dll = dll)
                                                    }      
                                                if(!missing(model)) {
                                                    setModel(model)
                                                    copyFromModel()
                                                    setupNodes(dll = dll)
                                                }
                                            },
                                                show = function() {
                                                    writeLines(paste0("Derived CmodelBaseClass created by buildModelInterface for model ", modelDef$name))
                                                }),
                                            where = where
  ), list(FIELDS = makeModelBindingFields(symTab), where = where ) ) )

    ans <- function(model, where = globalenv(), dll = NULL, ...) {
        newClass$new(model, defaults, basePtrCall, classEnvironment = where, dll = dll, ...) ## this means defaults will be in the closure of this function and hence found
    } 
    formals(ans)$where = where
    formals(ans)$dll = dll
    ans
}

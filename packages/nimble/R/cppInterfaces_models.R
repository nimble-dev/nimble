
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
  eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getModelValuesPtrFromModel, rPtr))
getMVName <- function(modelValuePtr, dll)
  eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getMVBuildName, modelValuePtr))

# for now export this as R<3.1.2 give warnings if don't

#' Class \code{CmodelBaseClass}
#' @aliases CmodelBaseClass 
#' @export
#' @description
#' Classes used internally in NIMBLE and not expected to be called directly by users.
CmodelBaseClass <- setRefClass('CmodelBaseClass',
                               contains = 'modelBaseClass',
                               fields = list(
                                   .basePtr = 'ANY', ## pointer to derived model C++ class (backwards terminology due to history, unfortunately)
                                   .namedObjectsPtr = 'ANY', ## pointer to base NamedObjects C++ base blass
                                   .ModelBasePtr = 'ANY',
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
                                   finalizeInternal = function() {
                                       for(vn in cppNames) {
                                           vPtrName <- paste(".", vn, "_Ptr", sep = "")
                                           assign(vPtrName, NULL, inherits = TRUE)
                                       }
                                       finalize()
                                       .basePtr <<- NULL
                                       .namedObjectsPtr <<- NULL
                                       .ModelBasePtr <<- NULL
                                       .nodeFxnPointers_byDeclID <<- NULL
                                       nimbleProject <<- NULL
                                   },
                                   finalize = function() {
                                       for(i in ls(Rmodel$nodes)) {
                                           if(is.null(nodes[[i]])) next
                                           if(is.list(nodes[[i]]))
                                               nodes[[i]][[1]]$finalizeInstance(nodes[[i]][[2]])
                                           else
                                               nodes[[i]]$finalize()
                                           nodes[[i]] <<- NULL
                                       }
                                       if(!is.null(.nodeFxnPointers_byDeclID))
                                           .nodeFxnPointers_byDeclID$finalize()
                                       if(!is.null(.namedObjectsPtr)) ## .basePtr
                                           nimbleInternalFunctions$nimbleFinalize(.namedObjectsPtr) ##.basePtr
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
                                       
                                       .nodeFxnPointers_byDeclID <<- new('numberedObjects', dll = dll) 
                                       maxID = length(modelDef$declInfo)
                                       .nodeFxnPointers_byDeclID$resize(maxID)

                                       for(declID in seq_along(nodes)) {
                                           thisNodeFunctionName <- names(Rmodel$nodeFunctions)[declID]
                                           basePtr <- if(is.list(nodes[[thisNodeFunctionName]])) ## it's a multiInterface
                                                          nodes[[thisNodeFunctionName]][[1]]$basePtrList[[ nodes[[thisNodeFunctionName]][[2]] ]]
                                                      else ## it's a direct interface
                                                          nodes[[thisNodeFunctionName]]$.basePtr
                                           .self$.nodeFxnPointers_byDeclID[declID] <- basePtr 
                                       }
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
    defaults$extPtrTypeIndex <- compiledModel$getExtPtrTypeIndex()
    
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
                                                newPtrPair <- eval(parse(text = ".Call(basePtrCall)"))
                                                .basePtr <<- newPtrPair[[1]]
                                                .ModelBasePtr <<- newPtrPair[[ defaults$extPtrTypeIndex['ModelBase'] ]]
                                                .namedObjectsPtr <<- newPtrPair[[ defaults$extPtrTypeIndex['NamedObjects'] ]]
                                                eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$register_namedObjects_Finalizer,
                                                          .namedObjectsPtr, ##.basePtr,
                                                          dll[['handle']], model$name))

                                                .modelValues_Ptr <<- nimbleInternalFunctions$getMVptr(.ModelBasePtr, dll = dll) ## this is a Values*
                                                defaultModelValues <<- nimbleInternalFunctions$CmodelValues$new(existingPtr = .modelValues_Ptr,
                                                                                                                buildCall = nimbleInternalFunctions$getMVName(.modelValues_Ptr, dll),
                                                                                                                initialized = TRUE, dll = dll )
                                                modelDef <<- model$modelDef
                                                graph <<- model$graph
                                                vars <<- model$vars
                                                isDataVars <<- model$isDataVars
                                                nimbleProject <<- defaults$project
                                                for(v in ls(model$isDataEnv)) isDataEnv[[v]] <<- model$isDataEnv[[v]]
                                                setData(modelDef$constantsList, warnAboutMissingNames = FALSE)
                                                cppNames <<- eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getAvailableNames, .namedObjectsPtr)) ## or could get this from R objects
                                                cppCopyTypes <<- defaults$cppCT
                                                compiledModel <<- defaults$cm
                                                for(vn in cppNames)
                                                    {
                                                        vPtrName <- paste(".", vn, "_Ptr", sep = "")
                                                     	.self[[vPtrName]] <<- nimbleInternalFunctions$newObjElementPtr(.namedObjectsPtr, vn, dll = dll)
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

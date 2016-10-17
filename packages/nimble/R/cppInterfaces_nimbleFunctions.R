
###		buildNimbleFxnInterface(refname, symbolTable, basePtrCall)
###						This is a function which builds a new reference class around 
###						a function which returns a pointer to a NimbleFxn. 
###						refname will be the name of the new class
###						symbolTable will be a used to determine element types
###						basePtrCall is the C++ call which builds the C++ object
###						i.e. in R .Call(basePtrCall) should return the pointer
###						to our NimbleFunction Class
###						IMPORTANT NOTE: symbolTable MUST match form object pointed to by
###						basePtrCall! In particular, if the symbolTable says an element is a
###						scalar (i.e. nDims = 0) and it's really a NimArr (or vice versa), 
###						attempting to access this element will cause R to crash!

###						Once the new class is built, a new object is create via refname$new()
###						Access to the elements is provided through nimbleFxnObject$Varname

makeNFBindingFields <- function(symTab, cppNames) {
    fieldList = list(.DUMMY = "ANY")  # We use this .DUMMY field to trick R into not mistakenly printing
                                        # an error. See initialization function inside buildNimbleFxnInterface
    vNames <- if(missing(cppNames)) names(symTab$symbols) else cppNames
    for(vn in vNames) {
        thisSymbol <- symTab$getSymbolObject(vn)
        if(is.null(thisSymbol)) next
        if(thisSymbol$type == 'model' || thisSymbol$type == 'symbolNodeFunctionVector' || thisSymbol$type == 'symbolModelVariableAccessorVector' ||thisSymbol$type == 'symbolModelValuesAccessorVector' || thisSymbol$type == 'symbolCopierVector' || thisSymbol$type == 'symbolIndexedNodeInfoTable') next ## skip models and NodeFunctionVectors and modelVariableAccessors      
        ptrName = paste0(".", vn, "_Ptr")
        fieldList[[ptrName]] <- "ANY" ## "ANY" 
        ## Model variables:
        if(inherits(thisSymbol, 'symbolOptimReadyFunction'))	next
        if(inherits(thisSymbol,'symbolNimArrDoublePtr')) { ## copy type 'modelVar'
            fieldList[[vn]] <- eval(substitute(
                function(x) {
                    if(missing(x))
                        VPTR
                    else {
                        if(!inherits(x, 'externalptr')) stop(paste('Nimble compilation error initializing ptr ', VPTRname, '.'), call. = FALSE)
                        ##message('setting a modelVar')
                        nimbleInternalFunctions$setDoublePtrFromSinglePtr(VPTR, x, dll = dll)   
                    }
                }, list(VPTR = as.name(ptrName), VPTRname = ptrName ) ) )
            next
        }
        if(inherits(thisSymbol, 'symbolVecNimArrPtr')){ ## copy type 'modelValuesPtr'
            fieldList[[vn]] <- eval(substitute(
                function(x) {
                    if(missing(x))
                        VPTR
                    else {
                        if(!inherits(x, 'externalptr')) stop(paste('Nimble compilation error initializing ptr ', VPTRname, '.'), call. = FALSE)
                        ##message('setting a modelValuesPtr')
                        nimbleInternalFunctions$setDoublePtrFromSinglePtr(VPTR, x, dll = dll)
                    }
                }, list(VPTR = as.name(ptrName), VPTRname = ptrName ) ) )
            next      	
        }
        
        if(inherits(thisSymbol, 'symbolNumericList') ) { ## copy type 'numericList' -- has fallen out of support
            fieldList[[vn]] <- 'ANY'
            next
        }
        
        if(inherits(thisSymbol, 'symbolNimbleFunction')) { ## copy type 'nimbleFunction'
            nfName <- paste0(".",vn,"_CnimbleFunction")
            fieldList[[nfName]] <- "ANY" ## This will have the ref class object that interfaces to the C++ nimbleFunction
            fieldList[[vn]] <- eval(substitute(
                function(x) {
                    if(missing(x))
                        NFNAME
                    else {
                        if(is.list(x)) { ## can be a list with first element a CmultiNimbleFunction object and second element an index
                            if(!inherits(x[[1]], 'CmultiNimbleFunctionClass')) stop(paste('Nimble compilation error initializing pointer for nimbleFunction from a CmultiNimbleFunction object ', NFNAMECHAR, '.'), call. = FALSE)
                            basePtr <- x[[1]]$basePtrList[[ x[[2]] ]]
                            ##message('setting a nimbleFunction for a multiInterface')
                            nimbleInternalFunctions$setDoublePtrFromSinglePtr(VPTR, basePtr, dll = dll) ## check field name
                        } else {
                            if(!inherits(x, 'CnimbleFunctionBase')) stop(paste('Nimble compilation error initializing nimbleFunction ', NFNAMECHAR, '.'), call. = FALSE)
                            if(!inherits(x$.basePtr, 'externalptr')) stop(paste('Nimble compilation error initializing pointer for nimbleFunction ', NFNAMECHAR, '.'), call. = FALSE)
                            ##message('setting a nimbleFunction for a regular interface')
                            nimbleInternalFunctions$setDoublePtrFromSinglePtr(VPTR, x$.basePtr, dll = dll) ## check field name
                        }
                        assign(NFNAMECHAR, x, inherits = TRUE) ## avoids <<- warnings
                       }
                }, list(VPTR = as.name(ptrName), NFNAME = as.name(nfName), NFNAMECHAR = nfName) ) )
            next
        }
        if(inherits(thisSymbol, 'symbolModelValues')) { ## copy type 'modelValues' ## similar behavior to symbolNimArrDoublePtr
            mvName <- paste0(".", vn, "_CmodelValues")
            fieldList[[mvName]] <- "ANY" ## This will have the CmodelValues object
            fieldList[[vn]] <- eval(substitute(
                function(x) {
                    if(missing(x))
                        MVNAME
                    else {
                        if(!inherits(x, 'CmodelValues')) stop(paste('Nimble compilation error initializing modelVaues ', MVNAMECHAR, '.'), call. = FALSE)
                        if(!inherits(x$extptr, 'externalptr')) stop(paste('Nimble compilation error initializing pointer for modelValues ', MVNAMECHAR, '.'), call. = FALSE)
                        ##message('setting a modelValues')
                        nimbleInternalFunctions$setDoublePtrFromSinglePtr(VPTR, x$extptr, dll = dll)
                        assign(MVNAMECHAR, x, inherits = TRUE) ## THIS WAY TO AVOID refClass "<<-" warnings 
                    }
                }, list(VPTR = as.name(ptrName), MVNAME = as.name(mvName), MVNAMECHAR = mvName ) ) )
            next
        }
        if(inherits(thisSymbol, 'symbolNimPtrList')) { ## copy type 'nimPtrList' ## for nimbleFunctionList, set up with some partial generality but not all the way
            nflName <- paste0(".",vn,"_CnimbleFunctionList")
            accessorPtrName <- paste0(".", vn, "_setter_Ptr") ## "_setter" part must match nimbleDSL_class_symbolTable symbolNimPtrList 
            fieldList[[accessorPtrName]] <- "ANY"
            fieldList[[nflName]] <- "ANY"
            fieldList[[vn]] <- eval(substitute(
                function(x) {
                    if(missing(x))
                        NFLNAME
                    else {
                        if(!inherits(x, 'nimPointerList')) stop(paste('Nimble compilation error initializing nimPointerList (nimbleFunctionList) ', NFLNAMECHAR, '.'), call. = FALSE)
                        if(!nimbleInternalFunctions$checkNimbleFunctionListCpp(x)) stop(paste('Nimble compilation error initializing nimbleFunctionList ', NFLNAMECHAR, '.  Something is not valid in this list.  It may be the contains (base class) value of one or more functions in the list.'), call. = FALSE)
                        nimbleInternalFunctions$setPtrVectorOfPtrs(ACCESSPTR, CONTENTSPTR, length(x$contentsList), dll = dll)
                        for(i in seq_along(x$contentsList)) {
                            if(is.list(x[[i]])) { ## case of list(CmultiNimbleFunction, index)
                                basePtr <- x[[i]][[1]]$basePtrList[[ x[[i]][[2]] ]]
                            } else {              ## case of CnimbleFunction
                                basePtr <- x[[i]]$.basePtr
                            }
                            if(!inherits(basePtr, 'externalptr')) stop(paste('Nimble compilation error initializing pointer ', i, ' of nimPointerList (nimbleFunctionList) ', NFLNAMECHAR, '.'), call. = FALSE)
                            ##message('setting a nimbleFunctionList')
                            nimbleInternalFunctions$setOnePtrVectorOfPtrs(ACCESSPTR, i, basePtr, dll = dll)
                        }
                        assign(NFLNAMECHAR, x, inherits = TRUE)
                    }
                }, list(NFLNAME = as.name(nflName), NFLNAMECHAR = nflName, CONTENTSPTR = as.name(ptrName), ACCESSPTR = as.name(accessorPtrName) ) ) )
            ## getter: return the nimPointerList set up with the list of CinterfaceObjects
            ## setter: call setPtrVectorOfPtrs(accessorExtPtr, contentsExtrPtr.  Then iterate and call setOnePtrVectorOfPtrs(accessorPtr, i, nimPtrList[[i]]$.basePtr)
            next
        }
        if(thisSymbol$type == "character") { ## cpp copy type 'character'  : 2 sub-cases (vector and scalar)
            if(thisSymbol$nDim > 0) {   ## character vector (nDim can only be 0 or 1)
                eval(substitute( fieldList$VARNAME <- function(x){
                    if(missing(x) ) 
                        nimbleInternalFunctions$getCharacterVectorValue(VPTR, dll = dll)
                    else
                        nimbleInternalFunctions$setCharacterVectorValue(VPTR, x, dll = dll)
                }, list(VPTR = as.name(ptrName), VARNAME = vn) ) )
            next
            } else {                    ## character scalar
                eval(substitute( fieldList$VARNAME <- function(x){
                    if(missing(x) ) 
                        nimbleInternalFunctions$getCharacterValue(VPTR, dll = dll)
                    else
                        nimbleInternalFunctions$setCharacterValue(VPTR, x, dll = dll)
                }, list(VPTR = as.name(ptrName), VARNAME = vn) ) )
                next
            }
        }
        if(inherits(thisSymbol, 'symbolBase')) { ## All numeric and logical cases  ## cpp copy type 'numeric': 4 sub-cases
            if(thisSymbol$nDim > 0) {            ## Anything vector
                eval(substitute( fieldList$VARNAME <- function(x){
                    
                    if(missing(x) ) 
                        nimbleInternalFunctions$getNimValues(VPTR, dll = dll)
                    else {
                        ##message('setting a numeric via setNimValues')
                        nimbleInternalFunctions$setNimValues(VPTR, x, dll = dll)
                    }
                }, list(VPTR = as.name(ptrName), VARNAME = vn) ) )
                next
            }
            
            if(thisSymbol$type == "double"){     ## Scalar double
                eval(substitute( fieldList$VARNAME <- function(x){
                    if(missing(x) ) 
                        nimbleInternalFunctions$getDoubleValue(VPTR, dll = dll)
                    else
                       nimbleInternalFunctions$setDoubleValue(VPTR, x, dll = dll)
                        
                }, list(VPTR = as.name(ptrName), VARNAME = vn) ) )
                next
            }
            if(thisSymbol$type == "integer"){    ## Scalar int
                eval(substitute( fieldList$VARNAME <- function(x){
                    if(missing(x) ) 
                        nimbleInternalFunctions$getIntValue(VPTR, dll = dll)
                    else
                        nimbleInternalFunctions$setIntValue(VPTR, x, dll = dll)                        
                }, list(VPTR = as.name(ptrName), VARNAME = vn) ) )
                next
            }
            if(thisSymbol$type == "logical"){    ## Scalar logical
                eval(substitute( fieldList$VARNAME <- function(x){
                    if(missing(x) ) 
                        nimbleInternalFunctions$getBoolValue(VPTR, dll = dll)
                    else
                        nimbleInternalFunctions$setBoolValue(VPTR, x, dll = dll)
                }, list(VPTR = as.name(ptrName), VARNAME = vn) ) )
                next
            }
            warning("Warning: scalar datatype not current supported for variable ", vn, "\n", call. = FALSE)
            browser()
            return(NULL)
        }
        warning('Warning in trying to build interface: symbol ', vn, ' not understood.', call. = FALSE)
    }
    return(fieldList)
}

#### functions that are similar to what is created in makeNFBindingFields but are standalone and look up pointers each time
getSetModelVarPtr <- function(name, value, basePtr, dll) { ## This only deals with a pointer member data.  It doesn't return or set the model's actual values.
    vptr <- nimbleInternalFunctions$newObjElementPtr(basePtr, name, dll = dll)
    if(missing(value)) {
        return(vptr)
    } else {
        ##message('setting from getSetModelVarPtr')
        nimbleInternalFunctions$setDoublePtrFromSinglePtr(vptr, value, dll = dll)
    }
}

getSetModelValuesPtr <- function(name, value, basePtr, dll) {
    ##message('setting from getSetModelValuesPtr')
    getSetModelVarPtr(name, value, basePtr, dll = dll)
}

getSetNimbleFunction <- function(name, value, basePtr, dll) {
    if(missing(value)) {
        warning('getSetNimbleFunction does not work for getting but was called without value.', call. = FALSE)
        return(NULL)
    } else {
        vptr <- nimbleInternalFunctions$newObjElementPtr(basePtr, name, dll = dll)
        ##message('setting from getSetNimbleFunction')
        nimbleInternalFunctions$setDoublePtrFromSinglePtr(vptr, value, dll = dll) ## previously value$.basePtr
    }
}

getSetModelValues <- function(name, value, basePtr, dll) {
      if(missing(value)) {
        warning('getSetModelValues does not work for getting but was called without value.', call. = FALSE)
        return(NULL)
    } else {
        vptr <- nimbleInternalFunctions$newObjElementPtr(basePtr, name, dll = dll)
        ##message('setting from getSetModelValues')
        nimbleInternalFunctions$setDoublePtrFromSinglePtr(vptr, value$extptr, dll = dll)
    }  
}

getSetNimPtrList <- function(name, value, basePtr, dll) {
    if(missing(value)) {
        warning('getSetNimPtrList does not work for getting but was called without value.', call. = FALSE)
        return(NULL)
    } else {
        if(!inherits(value, 'nimPointerList')) stop(paste('Nimble compilation error initializing nimPointerList (nimbleFunctionList) ', name, '.'), call. = FALSE)
        if(!nimbleInternalFunctions$checkNimbleFunctionListCpp(value)) stop(paste('Nimble compilation error initializing nimbleFunctionList ', name, '.  Something is not valid in this list.  It may be the contains (base class) value of one or more functions in the list.'), call. = FALSE)
        vptr <- nimbleInternalFunctions$newObjElementPtr(basePtr, name, dll = dll)
        accessptr <- nimbleInternalFunctions$newObjElementPtr(basePtr, paste0(name, '_setter'), dll = dll)
        nimbleInternalFunctions$setPtrVectorOfPtrs(accessptr, vptr, length(value$contentsList), dll = dll)
        for(i in seq_along(value$contentsList)) {
            if(is.list(value[[i]])) { ## case of list(CmultiNimbleFunction, index)
                basePtr <- value[[i]][[1]]$basePtrList[[ value[[i]][[2]] ]]
            } else {              ## case of CnimbleFunction
                basePtr <- value[[i]]$.basePtr
            }
            if(!inherits(basePtr, 'externalptr')) stop(paste('Nimble compilation error initializing pointer ', i, ' of nimPointerList (nimbleFunctionList) ', name, '.'), call. = FALSE)
            ##message('setting from getSetNimPtrList')
            nimbleInternalFunctions$setOnePtrVectorOfPtrs(accessptr, i, basePtr, dll = dll)
        }
    }
}

getSetCharacterVector <- function(name, value, basePtr, dll) {
    vptr <- nimbleInternalFunctions$newObjElementPtr(basePtr, name, dll = dll)
     if(missing(value)) 
         nimbleInternalFunctions$getCharacterVectorValue(vptr, dll = dll)
     else
         nimbleInternalFunctions$setCharacterVectorValue(vptr, value, dll = dll)
}

getSetCharacter <- function(name, value, basePtr, dll) {
    vptr <- nimbleInternalFunctions$newObjElementPtr(basePtr, name, dll = dll)
     if(missing(value)) 
         nimbleInternalFunctions$getCharacterValue(vptr, dll = dll)
     else
         nimbleInternalFunctions$setCharacterValue(vptr, value, dll = dll)
}

getSetNumericVector <- function(name, value, basePtr, dll) {
      vptr <- nimbleInternalFunctions$newObjElementPtr(basePtr, name, dll = dll)
     if(missing(value)) 
         nimbleInternalFunctions$getNimValues(vptr, dll = dll)
     else {
         ##message('setting from getSetNumericVector')
         nimbleInternalFunctions$setNimValues(vptr, value, dll = dll)
     }
}

getSetDoubleScalar <- function(name, value, basePtr, dll) {
      vptr <- nimbleInternalFunctions$newObjElementPtr(basePtr, name, dll = dll)
     if(missing(value)) 
         nimbleInternalFunctions$getDoubleValue(vptr, dll = dll)
     else
         nimbleInternalFunctions$setDoubleValue(vptr, value, dll = dll)
}

getSetIntegerScalar <- function(name, value, basePtr, dll) {
      vptr <- nimbleInternalFunctions$newObjElementPtr(basePtr, name, dll = dll)
     if(missing(value)) 
         nimbleInternalFunctions$getIntValue(vptr, dll = dll)
     else
         nimbleInternalFunctions$setIntValue(vptr, value, dll = dll)
}

getSetLogicalScalar <- function(name, value, basePtr, dll) {
      vptr <- nimbleInternalFunctions$newObjElementPtr(basePtr, name, dll = dll)
     if(missing(value)) 
         nimbleInternalFunctions$getBoolValue(vptr, dll = dll)
     else
         nimbleInternalFunctions$setBoolValue(vptr, value, dll = dll)
}


#' Class \code{CnimbleFunctionBase}
#' @aliases CnimbleFunctionBase
#' @export
#' @description
#' Classes used internally in NIMBLE and not expected to be called directly by users.
CnimbleFunctionBase <- setRefClass('CnimbleFunctionBase',
                                   fields = list(
                                       dll = "ANY",
                                       .basePtr = 'ANY',
                                       compiledNodeFun = 'ANY',
                                       Robject = 'ANY', ## this should be the refClassObject, not the function
                                       cppNames = 'ANY',
                                       cppCopyTypes = 'ANY',
                                       neededObjects = 'ANY', ## A list of things like modelValues objects, if they don't already exist
                                       nimbleProject = 'ANY'
                                       ),
                                   methods = list(
                                       finalizeInternal = function() {
                                           neededObjects <<- nimbleInternalFunctions$clearNeededObjects(Robject, compiledNodeFun, neededObjects)
                                           vPtrNames <- 	paste0(".", cppNames, "_Ptr")	
                                           for(vn in seq_along(cppNames) ){
                                               .self[[vPtrNames[vn]]] <- NULL
                                           }
                                           finalize()
                                           .basePtr <<- NULL
                                           nimbleProject <<- NULL
                                       },
                                       finalize = function() {
                                           nimbleInternalFunctions$nimbleFinalize(.basePtr)
                                       },
                                       initialize = function(dll = NULL, project = NULL, test = TRUE, ...) {
                                       		neededObjects <<- list()
                                           if(!test) {
                                               dll <<- dll
                                               if(is.null(project)) {
                                                   warning('Missing project argument in CnimbleFunctionBase, from a nimbleFunction C++ interface. Could crash if there are member objects with project information needed.', call.=FALSE)
                                               }
                                               nimbleProject <<- project
                                           }
                                                callSuper(...)
                                            },
                                       getDefinition = function()
                                           nimble:::getDefinition(.self),
                                       setRobject = function(Robj) {
                                           if(is.nf(Robj)) Robject <<- nimble:::nf_getRefClassObject(Robj)
                                           else Robject <<- Robj
                                           Robject$.CobjectInterface <<- .self 
                                       },
                                       lookupSymbol = function(symname) {
                                           if(is.null(dll))
                                               stop("No DLL for this object")
                                           
                                           getNativeSymbolInfo(symname, dll)
                                       }
                                       ))

makeNimbleFxnCppCopyTypes <- function(symTab, cppNames) {
    ans <- list()
    vNames <- if(missing(cppNames)) names(symTab$symbols) else cppNames
    for(vn in vNames) {
        thisSymbol <- symTab$getSymbolObject(vn)
        if(is.null(thisSymbol)) next
        if(thisSymbol$type == 'Ronly') next ## skip models
        if(inherits(thisSymbol, 'symbolIndexedNodeInfoTable')) {ans[[thisSymbol$name]] <- 'indexedNodeInfoTable'; next}
        if(inherits(thisSymbol, 'symbolNimArrDoublePtr')) {ans[[thisSymbol$name]] <- 'modelVar'; next}
        if(inherits(thisSymbol, 'symbolNodeFunctionVector'))  { ans[[thisSymbol$name]] <- 'nodeFxnVec'; next}
        if(inherits(thisSymbol, 'symbolModelVariableAccessorVector')) {ans[[thisSymbol$name]] <- 'modelVarAccess';next}
        if(inherits(thisSymbol, 'symbolModelValuesAccessorVector')) {ans[[thisSymbol$name]] <- 'modelValuesAccess';next}
        if(inherits(thisSymbol, 'symbolModelValues')) {ans[[thisSymbol$name]] <- 'modelValues'; next}
        if(inherits(thisSymbol, 'symbolNimbleFunction')) {ans[[thisSymbol$name]] <- 'nimbleFunction'; next}
        if(inherits(thisSymbol, 'symbolVecNimArrPtr')) {ans[[thisSymbol$name]] <- 'modelValuesPtr'; next} ## from a singleModelValuesAccessClass, from e.g. mv[i, 'x']
        if(inherits(thisSymbol, 'symbolNumericList')) {ans[[thisSymbol$name]] <- 'numericList'; next}
        if(inherits(thisSymbol, 'symbolNimPtrList')) {ans[[thisSymbol$name]] <- 'nimPtrList'; next}
        if(inherits(thisSymbol, 'symbolCopierVector')) {ans[[thisSymbol$name]] <- 'copierVector'; next}
        if(inherits(thisSymbol, 'symbolString')) {
            if(thisSymbol$nDim > 0)
                ans[[thisSymbol$name]] <- 'characterVector'
            else
                ans[[thisSymbol$name]] <- 'characterScalar'
            next
        }
        if(thisSymbol$nDim > 0) {ans[[thisSymbol$name]] <- 'numericVector'; next}
        if(thisSymbol$type == 'double') {ans[[thisSymbol$name]] <- 'doubleScalar'; next}
        if(thisSymbol$type == 'integer') {ans[[thisSymbol$name]] <- 'integerScalar'; next}
        if(thisSymbol$type == 'logical') {ans[[thisSymbol$name]] <- 'logicalScalar'; next}
        stop(paste0('Confused in assigning a cpp copy type for ',thisSymbol$name), call. = FALSE) 
    }
    ans
}

makeNimbleFxnInterfaceCallMethodCode <- function(compiledNodeFun, includeDotSelfAsArg = FALSE, embedInBrackets = FALSE) {
    numFuns <- length(compiledNodeFun$RCfunDefs)
    ans <- if(!embedInBrackets) quote(list()) else quote({}) 
    if(numFuns == 0) return(ans)
    funNames <- names(compiledNodeFun$RCfunDefs)
    funNames[funNames == 'operator()'] <- 'run'
    for(i in seq_along(compiledNodeFun$RCfunDefs)) {
        ## note that the className is really used as a boolean: any non-NULL value triggers treatment as a class, but name isn't needed
        ans[[i+1]] <- compiledNodeFun$RCfunDefs[[i]]$buildRwrapperFunCode(className = compiledNodeFun$nfProc$name, includeLHS = FALSE, returnArgsAsList = FALSE, includeDotSelf = '.basePtr', includeDotSelfAsArg = includeDotSelfAsArg)
        if(embedInBrackets) ans[[i+1]] <- substitute(THISNAME <- FUNDEF, list(THISNAME = as.name(funNames[i]), FUNDEF = ans[[i+1]]))
    }
    if(!embedInBrackets) names(ans) <- c('', funNames)
    ans
}

clearNeededObjects <- function(Robj, compiledNodeFun, neededObjects) {
    for(iName in compiledNodeFun$nfProc$neededObjectNames) {
        thisObj <- Robj[[iName]]

        if(inherits(thisObj, 'modelValuesBaseClass')) {
            ## if(!(inherits(thisObj$CobjectInterface, 'uninitializedField') || is.null(thisObj$CobjectInterface))) {
            ##     neededObjects[[iName]] <- nimbleProject$clearCmodelValues(thisObj) ## should be redundant since CmodelValues are added to the nimbleProject earlier and so will be cleared anyway via nimbleProject$clearCompiled()
            ## }

            ## We skip over CmodelValues for now because they should be cleared by nimbleProject$clearCompiled
            ## otherwise we don't have a good way to look them up from the project except by index, which isn't known
            next
        }
        if(is.nf(thisObj)) {
            RCO <- nf_getRefClassObject(thisObj)
            if(!(inherits(RCO$.CobjectInterface, 'uninitializedField') || is.null(RCO$.CobjectInterface))) {
                if(is.list(RCO$.CobjectInterface))
                    RCO$.CobjectInterface[[1]]$finalizeInstance(RCO$.CobjectInterface[[2]])
                else
                    RCO$.CobjectInterface$finalize()
                neededObjects[[iName]] <- NULL
            }
            next
        }
        if(inherits(thisObj, 'nimbleFunctionList')) {
            for(i in seq_along(thisObj$contentsList)) {
                RCO <- nf_getRefClassObject(thisObj[[i]])
                if(!(inherits(RCO$.CobjectInterface, 'uninitializedField') || is.null(RCO$.CobjectInterface))) {
                    if(is.list(RCO$.CobjectInterface))
                        RCO$.CobjectInterface[[1]]$finalizeInstance(RCO$.CobjectInterface[[2]])
                    else
                        RCO$.CobjectInterface$finalize()
                    neededObjects[[iName]][[i]] <- NULL
                } 
            }
            thisObj$CobjectInterface <- NULL
            neededObjects[[iName]] <- NULL
            next
        }
    }
    neededObjects
}

buildNeededObjects <- function(Robj, compiledNodeFun, neededObjects, dll, nimbleProject) {
    for(iName in compiledNodeFun$nfProc$neededObjectNames) {
        thisObj <- Robj[[iName]]
        if(inherits(thisObj, 'modelValuesBaseClass')) {
            if(inherits(thisObj$CobjectInterface, 'uninitializedField') || is.null(thisObj$CobjectInterface)) {
                neededObjects[[iName]] <- nimbleProject$instantiateCmodelValues(thisObj, dll)
            }
            next
        }
        if(is.nf(thisObj)) {
            RCO <- nf_getRefClassObject(thisObj)
            if(inherits(RCO$.CobjectInterface, 'uninitializedField') || is.null(RCO$.CobjectInterface)) {
                neededObjects[[iName]] <- nimbleProject$instantiateNimbleFunction(thisObj, dll, asTopLevel = getNimbleOption('buildInterfacesForCompiledNestedNimbleFunctions'))
            }
            next
        }
        if(inherits(thisObj, 'nimbleFunctionList')) {
            neededObjects[[iName]] <- nimPointerList(thisObj$baseClass, length(thisObj$contentsList))
            for(i in seq_along(thisObj$contentsList)) {
                RCO <- nf_getRefClassObject(thisObj[[i]])
                if(inherits(RCO$.CobjectInterface, 'uninitializedField') || is.null(RCO$.CobjectInterface)) {
                    neededObjects[[iName]][[i]] <- nimbleProject$instantiateNimbleFunction(thisObj[[i]], dll, asTopLevel = getNimbleOption('buildInterfacesForCompiledNestedNimbleFunctions'))
                } else {
                    neededObjects[[iName]][[i]] <- RCO$.CobjectInterface ## either CnimbleFunction or list(CmultiNimbleFunction, index)
                }
            }
            names(neededObjects[[iName]]$contentsList) <- names(thisObj$contentsList)
            thisObj$CobjectInterface <- neededObjects[[iName]]
            next
        }
        warning('Warning: object ',iName,' not of a type that can be built.', call. = FALSE)
    }
    neededObjects
}

copyFromRobjectViaActiveBindings = function(Robj, cppNames, cppCopyTypes, .self, dll) {
    for(v in cppNames) {
        if(is.null(cppCopyTypes[[v]])) next
        if(is.null(Robj[[v]])) {
            warning("Problem in copyFromRobject.  There is an object to be copied that is NULL.  Going to browser.", call. = FALSE)
            browser()
        }
        if(cppCopyTypes[[v]] == 'modelVar') {
            modelVar <- Robj[[v]] ## this is a singleVarAccessClass created by replaceModelSingles
            Cmodel <- modelVar$model$CobjectInterface
            varName <- modelVar$var
            .self[[v]] <- eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getModelObjectPtr, Cmodel$.basePtr, varName))
            next
        }
        else if(cppCopyTypes[[v]] == 'nimbleFunction') {
            modelVar <- Robj[[v]]
            Cnf <- nf_getRefClassObject(modelVar)$.CobjectInterface ##environment(modelVar)$.CobjectInterface
            ## Cnf coule be old format (CnimbleFunction) or a list(CmultiNimbleFunction, index)
            .self[[v]] <- Cnf
            next
        }
        else if(cppCopyTypes[[v]] == 'nimPtrList') {
            if(is.null(Robj[[v]]$contentsList)) {
                warning('Problem in copying a nimPtrList to C++ object. The contentsList is NULL. Going to browser', call. = FALSE)
                browser()
            }
            if(any(unlist(lapply(Robj[[v]]$contentsList, is.null)))) {
                warning('Problem in copying a nimPtrList to C++ object. The contentsList is NULL')
                browser()
            }
            modelVar <- Robj[[v]] ## This is a nimPtrList 
            Cmv <- modelVar$CobjectInterface ## This was created above in build neededObjects
            .self[[v]] <- Cmv
        }
        else if(cppCopyTypes[[v]] == 'modelValues') { ## somewhat similar to modelVar
            rModelValues <- Robj[[v]]
            Cmv <- rModelValues$CobjectInterface
            if(!Cmv$initialized) {
                ##cat('We are copying modelValues during copyFromRobjectViaActiveBindings\n')
                k = getsize(rModelValues)
                resize(Cmv, k)
                vNames = rModelValues[['varNames']]
                for(vN in vNames)
                    Cmv[vN,] <- rModelValues[vN,]
                Cmv$symTab <- rModelValues$symTab
                Cmv$initialized <- TRUE
            }
            .self[[v]] <- Cmv
            next
        }
        else if(cppCopyTypes[[v]] == 'nodeFxnVec') {
            ##populateNodeFxnVec(fxnPtr = .self$.basePtr, Robject = Robj, fxnVecName = v)
            populateNodeFxnVecNew(fxnPtr = .self$.basePtr, Robject = Robj, fxnVecName = v, dll = dll) 
            next
        }
        else if(cppCopyTypes[[v]] == 'modelVarAccess'){
            populateManyModelVarMapAccess(fxnPtr = .self$.basePtr, Robject = Robj, manyAccessName = v, dll = dll)
            next
        }
        else if(cppCopyTypes[[v]] == 'modelValuesAccess'){
            populateManyModelValuesMapAccess(fxnPtr = .self$.basePtr, Robject = Robj, manyAccessName = v, dll = dll)
            next
        }
        else if(cppCopyTypes[[v]] == "modelValuesPtr"){
            curObj <- Robj[[v]]
            mvPtr = curObj$modelValues$CobjectInterface$componentExtptrs[[curObj$var]]
            .self[[v]] <- mvPtr
            next
        }
        else if(cppCopyTypes[[v]] == 'numericList'){
            stop('numericList is not working\n')
            ## rawPtr = .Call('getModelObjectPtr', .self$.basePtr, v)
            ## .self[[v]] <- numericList(buildType = 'C', extPtr = rawPtr)
            ## nRows = Robj[[v]]$nRow
            ## resize(.self[[v]], nRows)
            ## for(i in 1:nRows){
            ##     copyDims = max(c(1, dimOrLength[[v]][[i]]) )
            ##     d1 = copyDims[1]
            ##     d2 = copyDims[2]
            ##     d3 = copyDims[3]
            ##     setSize(.self[[v]], row = i, d1, d2, d3)
            ##     .self[[v]][[i]] <- Robj[[v]][[i]]
            ## }
            next
        }
        else if(cppCopyTypes[[v]] == 'indexedNodeInfoTable') {
            populateIndexedNodeInfoTable(fxnPtr = .self$.basePtr, Robject = Robj, indexedNodeInfoTableName = v, dll = dll)
        }
        else if(cppCopyTypes[[v]] == 'characterVector' || cppCopyTypes[[v]] == 'characterScalar') {
            .self[[v]] <- Robj[[v]]
        }
        else if(cppCopyTypes[[v]] %in% c('numericVector','doubleScalar','integerScalar','logicalScalar')) {
            .self[[v]] <- Robj[[v]]
        }
        else if(!(cppCopyTypes[[v]] %in% c('copierVector'))) {
            warning(paste0("Note: cppCopyTypes not recognized. Type = ", cppCopyTypes[[v]], "\n"), call. = FALSE)
        }
    }
    ## second pass is for initializations that require everything from first pass be done
    for(v in cppNames) {
        if(is.null(cppCopyTypes[[v]])) next
        if(cppCopyTypes[[v]] == 'copierVector') {
            populateCopierVector(fxnPtr = .self$.basePtr, Robject = Robj, vecName = v, dll = dll)
        }
    }
}

copyFromRobject <- function(Robj, cppNames, cppCopyTypes, basePtr, dll) {
    for(v in cppNames) {
        if(is.null(cppCopyTypes[[v]])) next
        if(is.null(Robj[[v]])) {
            warning("Problem in copyFromRobject.  There is an object to be copied that is NULL.  Going to browser.", call. = FALSE)
            browser()
        }
        if(cppCopyTypes[[v]] == 'modelVar') {
            modelVar <- Robj[[v]] ## this is a singleVarAccessClass created by replaceModelSingles
            Cmodel <- modelVar$model$CobjectInterface
            varName <- modelVar$var
            ##message('copying modelVar (copyFromRobject)')
            getSetModelVarPtr(v, eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getModelObjectPtr, Cmodel$.basePtr, varName)), basePtr, dll = dll)
            next
        }
        else if(cppCopyTypes[[v]] == 'nimbleFunction') {
            modelVar <- Robj[[v]]
            Cnf <- nf_getRefClassObject(modelVar)$.CobjectInterface ##environment(modelVar)$.CobjectInterface
            if(is.list(Cnf)) {
                valueBasePtr <- Cnf[[1]]$basePtrList[[ Cnf[[2]] ]]
            } else {
                valueBasePtr <- Cnf$.basePtr
            }
            ##message('copying nimbleFunction (copyFromRobject)')
            getSetNimbleFunction(v, valueBasePtr, basePtr, dll = dll)
            ## .self[[v]] <- Cnf
            next
        }
        else if(cppCopyTypes[[v]] == 'nimPtrList') {
            if(is.null(Robj[[v]]$contentsList)) {
                warning('Problem in copying a nimPtrList to C++ object. The contentsList is NULL. Going to browser', call. = FALSE)
                browser()
            }
            if(any(unlist(lapply(Robj[[v]]$contentsList, is.null)))) {
                warning('Problem in copying a nimPtrList to C++ object. The contentsList is NULL')
                browser()
            }
            modelVar <- Robj[[v]] ## This is a nimPtrList 
            Cmv <- modelVar$CobjectInterface ## This was created above in build neededObjects
            ##message('copying nimPtrList (copyFromRobject)')
            getSetNimPtrList(v, Cmv, basePtr, dll = dll)
            ##.self[[v]] <- Cmv
        }
        else if(cppCopyTypes[[v]] == 'modelValues') { ## somewhat similar to modelVar
            rModelValues <- Robj[[v]]
            Cmv <- rModelValues$CobjectInterface
            if(!Cmv$initialized) {
                ##cat('We are copying modelValues during copyFromRobject\n')
                k = getsize(rModelValues)
                resize(Cmv, k)
                vNames = rModelValues[['varNames']]
                for(vN in vNames)
                    Cmv[vN,] <- rModelValues[vN,]
                Cmv$symTab <- rModelValues$symTab
                Cmv$initialized <- TRUE
            }
            ##message('copying modelValues (copyFromRobject)')
            getSetModelValues(v, Cmv, basePtr, dll = dll)
            ##            .self[[v]] <- Cmv
            next
        }
        else if(cppCopyTypes[[v]] == 'nodeFxnVec') {
            populateNodeFxnVecNew(fxnPtr = basePtr, Robject = Robj, fxnVecName = v, dll = dll)
            ## populateNodeFxnVec(fxnPtr = .self$.basePtr, Robject = Robj, fxnVecName = v) 
            next
        }
        else if(cppCopyTypes[[v]] == 'modelVarAccess'){
            populateManyModelVarMapAccess(fxnPtr = basePtr, Robject = Robj, manyAccessName = v, dll = dll)
            ## populateManyModelVarMapAccess(fxnPtr = .self$.basePtr, Robject = Robj, manyAccessName = v)
            next
        }
        else if(cppCopyTypes[[v]] == 'modelValuesAccess'){
            populateManyModelValuesMapAccess(fxnPtr = basePtr, Robject = Robj, manyAccessName = v, dll = dll)
            ## populateManyModelValuesMapAccess(fxnPtr = .self$.basePtr, Robject = Robj, manyAccessName = v)
            next
        }
        else if(cppCopyTypes[[v]] == "modelValuesPtr"){
            curObj <- Robj[[v]]
            mvPtr = curObj$modelValues$CobjectInterface$componentExtptrs[[curObj$var]]
            getSetModelValuesPtr(v, mvPtr, basePtr, dll = dll)
            ## .self[[v]] <<- mvPtr
            next
        }
        else if(cppCopyTypes[[v]] == 'numericList'){
            stop('numericList is not working\n')
            ## rawPtr = .Call('getModelObjectPtr', .self$.basePtr, v)
            ## .self[[v]] <<- numericList(buildType = 'C', extPtr = rawPtr)
            ## nRows = Robj[[v]]$nRow
            ## resize(.self[[v]], nRows)
            ## for(i in 1:nRows){
            ##     copyDims = max(c(1, dimOrLength[[v]][[i]]) )
            ##     d1 = copyDims[1]
            ##     d2 = copyDims[2]
            ##     d3 = copyDims[3]
            ##     setSize(.self[[v]], row = i, d1, d2, d3)
            ##     .self[[v]][[i]] <<- Robj[[v]][[i]]
            ## }
            next
        }
        else if(cppCopyTypes[[v]] == 'indexedNodeInfoTable') {
            populateIndexedNodeInfoTable(fxnPtr = basePtr, Robject = Robj, indexedNodeInfoTableName = v, dll = dll)
        }
        else if(cppCopyTypes[[v]] == 'characterVector') {
            getSetCharacterVector(v, Robj[[v]], basePtr, dll = dll)
            ##.self[[v]] <<- Robj[[v]]
            next
        }
        else if(cppCopyTypes[[v]] == 'characterScalar') {
            getSetCharacter(v, Robj[[v]], basePtr, dll = dll)
            next
        }
        else if(cppCopyTypes[[v]] == 'numericVector') {
            ##message('copying numericVector (copyFromRobject)')
            getSetNumericVector(v, Robj[[v]], basePtr, dll = dll)
            ##.self[[v]] <<- Robj[[v]]
            next
        }
        else if(cppCopyTypes[[v]] == 'doubleScalar') {
            getSetDoubleScalar(v, Robj[[v]], basePtr, dll = dll)
            next
        }
        else if(cppCopyTypes[[v]] == 'integerScalar') {
            getSetIntegerScalar(v, Robj[[v]], basePtr, dll = dll)
            next
        }
        else if(cppCopyTypes[[v]] == 'logicalScalar') {
            getSetLogicalScalar(v, Robj[[v]], basePtr, dll = dll)
            next
        }
        else if(!(cppCopyTypes[[v]] %in% c('copierVector'))) {
            warning(paste0("Note: cppCopyTypes not recognized. Type = ", cppCopyTypes[[v]], "\n"), call. = FALSE)
        }
    }
    ## second pass is for initializations that require everything from first pass be done
    for(v in cppNames) {
        if(is.null(cppCopyTypes[[v]])) next
        if(cppCopyTypes[[v]] == 'copierVector') {
            populateCopierVector(fxnPtr = basePtr, Robject = Robj, vecName = v, dll = dll)
            ## populateCopierVector(fxnPtr = .self$.basePtr, Robject = Robj, vecName = v)
        }
    }
}

buildNimbleFxnInterface <- function(refName,  compiledNodeFun, basePtrCall, where = globalenv()){
    defaults <- list()
    if(inherits(compiledNodeFun, 'symbolTable')) {
        symTab <- compiledNodeFun
        defaults$cnf <- NULL
        warning('No compiled node function provided, so interface will be incomplete')
    } else {
        symTab <- compiledNodeFun$nfProc$setupSymTab
        defaults$cnf <- compiledNodeFun
    }
    ## The following is really equivalent, because it comes *directly* from the place that generates the C++ code
    cppNames <- compiledNodeFun$objectDefs$getSymbolNames() 
    NFBF <-  makeNFBindingFields(symTab, cppNames)
    defaults$cppCT <- makeNimbleFxnCppCopyTypes(symTab, cppNames)
    methodsList <- makeNimbleFxnInterfaceCallMethodCode(compiledNodeFun) ##, compiledNodeFun$nfProc)
    defaults$basePtrCall <- basePtrCall

    # substitute on parsed text string to avoid CRAN issues with .Call registration
    fun <- substitute(function(nfObject, defaults, dll = NULL, project = NULL, ...){		#cModel removed from here
        defaults$cnf$nfProc$evalNewSetupLinesOneInstance(nfObject, check = TRUE) ## in case this instance was not used during nfProc$process
        callSuper(dll = dll, project = project, test = FALSE, ...)
        basePtrCall <- if(is.character(defaults$basePtrCall)) {
            if(inherits(dll, "uninitializedField") | is.null(dll)) stop("Error making a nimbleFxnInterface object: no dll provided")
            lookupSymbol(defaults$basePtrCall)
        } else defaults$basePtrCall
        # avoid R CMD check problem with registration.  basePtrCall is already the result of getNativeSymbolInfo from the dll, if possible from cppDefs_nimbleFunction.R
        .basePtr <<- eval(parse(text = ".Call(basePtrCall)"))
        regLabel <- try(get('name', envir = nfObject), silent = TRUE)
        if(inherits(regLabel, 'try-error') | is.null(regLabel)) regLabel <- environment(nfObject)$className
        eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$register_namedObjects_Finalizer, .basePtr, dll[['handle']], regLabel))
        # .basePtr <<- .Call(basePtrCall)
        cppNames <<- eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getAvailableNames, .basePtr))
        cppCopyTypes <<- defaults$cppCT
        compiledNodeFun <<- defaults$cnf
        vPtrNames <- 	paste0(".", cppNames, "_Ptr")	
        for(vn in seq_along(cppNames) ){
            .self[[vPtrNames[vn]]] <- nimbleInternalFunctions$newObjElementPtr(.basePtr, cppNames[vn], dll = dll)
        }
         if(!missing(nfObject)) {
            setRobject(nfObject)
            ##buildNeededObjects()
            neededObjects <<- nimbleInternalFunctions$buildNeededObjects(Robject, compiledNodeFun, neededObjects, dll, nimbleProject)
            ##copyFromRobject()
            nimbleInternalFunctions$copyFromRobjectViaActiveBindings(Robject, cppNames, cppCopyTypes, .self, dll)
        }
    }, list())

      # if we just have the name of the routine and haven't resolved it, arrange to resolve it when this initialization
      # function is called.  So change the .Call('name') to .Call(lookupSymbol('name')) which will use this objects
      # dll field.
    
    methodsList[[length(methodsList) + 1]] <- fun 
    names(methodsList)[length(methodsList)] <- 'initialize'
    methodsList[[length(methodsList) + 1]] <- quote(function() {
        writeLines(paste0("Derived CnimbleFunctionBase object created by buildNimbleFxnInterface for nimbleFunction with class ", class(Robject)))
    })
    names(methodsList)[length(methodsList)] <- 'show'
    eval(substitute( newClass <-  setRefClass(refName,
                                              fields = FIELDS,
                                              contains = 'CnimbleFunctionBase',
                                              methods = ML,
                                              where = where),
                    list(FIELDS = NFBF, ML = methodsList ) ) )

    ans <- function(nfObject, dll = NULL, project) {
    	wrappedInterfaceBuild <- newClass$new
    	wrappedInterfaceBuild(nfObject, defaults, dll = dll, project = project) ## Only purpose of wrappedInterfaceBuild is to have a helpful name for Rprof that is not "new"
#        newClass$new(nfObject, defaults, dll = dll, project = project)
    }
    return(ans)
}

####
## New class for interfacing multiple compiledNimbleFunctions of the same class

CmultiNimbleFunctionClass <- setRefClass('CmultiNimbleFunctionClass',
                                         fields = list(
                                             nimbleProject = 'ANY',
                                             cppNames = 'ANY',
                                             cppCopyTypes = 'ANY',
                                             basePtrCall = 'ANY',
                                             dll = 'ANY',
                                             ## nfObjectList = 'ANY',
                                             basePtrList = 'ANY',
                                             RobjectList = 'ANY',
                                             neededObjectsList = 'ANY',
                                             compiledNodeFun = 'ANY',    ## a cppNimbleFunctionClass
                                             callEnv = 'ANY'
                                         ),
                                         methods = list(
                                             show = function() {
                                                 cat(paste0('CmultiNimbleFunctionClass object\n'))
                                             },
                                             initialize = function(compiledNodeFun, basePtrCall, project, ...) { ## need to set dll, nimbleProject
                                                 if(missing(project)) stop('Cannot create CmultiNimbleFunctionClass without a project', call. = FALSE)
                                                 if(is.null(project)) stop('Cannot create CmultiNimbleFunctionClass without a NULL project', call. = FALSE)
                                                 nimbleProject <<- project
                                                 neededObjectsList <<- list()
                                                 ##     nfObjectList <<- list()
                                                 basePtrList <<- list()
                                                 RobjectList <<- list()
                                                 dll <<- NULL
                                                 compiledNodeFun <<- compiledNodeFun
                                                 ## basePtrCall is the result of getNativeSymbolInfo with the dll if possible from cppDefs_nimbleFunction.R
                                                 basePtrCall <<- basePtrCall
                                                 callSuper(...)
                                                 symTab <- compiledNodeFun$nfProc$setupSymTab
                                                 cppNames <<- compiledNodeFun$objectDefs$getSymbolNames()
                                                 cppCopyTypes <<- makeNimbleFxnCppCopyTypes(symTab, cppNames)
                                                 ##                                          methodsList <- makeNimbleFxnInterfaceCallMethodCode(compiledNodeFun) ## can do this but need to pass another pointer in
                                                 callCode <- makeNimbleFxnInterfaceCallMethodCode(compiledNodeFun, includeDotSelfAsArg = TRUE, embedInBrackets = TRUE)
                                                 callEnv <<- new.env()
                                                 eval(callCode, envir = callEnv)
                                             },
                                             finalize = function() {
                                                 for(i in seq_along(basePtrList)) {
                                                     if(!is.null(basePtrList[[i]])) {
                                                         nimbleInternalFunctions$nimbleFinalize(basePtrList[[i]])
                                                     }
                                                 }
                                             },
                                             finalizeInstance = function(index) {
                                                 if(!is.null(basePtrList[[index]])) {
                                                     neededObjectsList[[index]] <<- nimbleInternalFunctions$clearNeededObjects(RobjectList[[index]], compiledNodeFun, neededObjectsList[[index]])
                                                     RobjectList[index] <<- list(NULL)
                                                     nimbleInternalFunctions$nimbleFinalize(basePtrList[[index]])
                                                     basePtrList[index] <<- list(NULL)
                                                 }          
                                             },
                                             addInstance = function(nfObject, dll = NULL) { ## role of initialize
                                                 if(!is.null(.self$dll)) {
                                                     if(!identical(dll, .self$dll)) stop('Can not addInstance of a compiled nimbleFunction from different DLLs', call. = FALSE)
                                                 } else {
                                                     if(is.null(dll)) stop('In addInstance, DLL was not set and so must be provided when calling', call. = FALSE)
                                                     dll <<- dll       ## should only occur first time addInstance is called
                                                 }
                                                 compiledNodeFun$nfProc$evalNewSetupLinesOneInstance(nfObject, check = TRUE)
                                                 
                                                 # avoid R CMD check problem with registration:
                                                 newBasePtr <- eval(parse(text = ".Call(basePtrCall)"))
                                                 regLabel <- try(get('name', envir = nfObject), silent = TRUE)
                                                 if(inherits(regLabel, 'try-error') | is.null(regLabel)) regLabel <- environment(nfObject)$className

                                                 eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$register_namedObjects_Finalizer, newBasePtr, dll[['handle']], regLabel))
                                                 
                                                 basePtrList[[length(basePtrList)+1]] <<- newBasePtr
                                                 if(is.nf(nfObject)) newRobject <- nf_getRefClassObject(nfObject)
                                                 else newRobject <- nfObject
                                                 newRobject$.CobjectInterface <- list(.self, length(basePtrList))
                                                 RobjectList[[length(RobjectList)+1]] <<- newRobject
                                                 newNeededObjects <- nimbleInternalFunctions$buildNeededObjects(newRobject, compiledNodeFun, list(), dll, nimbleProject)
                                                 neededObjectsList[[length(neededObjectsList) + 1]] <<- newNeededObjects
                                                 nimble:::copyFromRobject(newRobject, cppNames, cppCopyTypes, newBasePtr, dll)
                                                 if(getNimbleOption('clearNimbleFunctionsAfterCompiling')) compiledNodeFun$nfProc$clearSetupOutputs(newRobject)
                                                 list(.self, length(basePtrList)) ## (this object, index)
                                             },
                                             callMemberFunction = function(index, funName, ...) {
                                                 callEnv[[funName]](..., .basePtr = basePtrList[[index]])
                                             },
                                             memberData = function(index, name, value) { ## value can be missing
                                                 ## This isn't very useful as written for many names because it just gets and sets the external pointers.  It doesn't wrap them in an interface object
                                                 if(!(name %in% cppNames)) stop(paste0('Name ', name, ' is not a valid member variable in the requested object.', call.=FALSE))
                                                 basePtr <- basePtrList[[index]]
                                                 if(!inherits(basePtr, 'externalptr')) stop('Invalid index or basePtr', call. = FALSE)
                                                 ans <- switch(cppCopyTypes[[name]],
                                                               modelVar = {##message('switch modelVar');
                                                                   getSetModelVarPtr(name, value, basePtr, dll = dll)}, ## only makes sense internally
                                                               nimbleFunction ={##message('switch nimbleFunction');
                                                                   getSetNimbleFunction(name, value, basePtr, dll = dll)}, ## ditto
                                                               nimPtrList = {##message('switch nimPtrList');
                                                                   getSetNimPtrList(name, value, basePtr, dll = dll)}, ## ditto
                                                               modelValues = {##message('switch modelValues');
                                                                   getSetModelValues(name, value, basePtr, dll = dll)}, ## ditto
                                                               characterVector = getSetCharacterVector(name, value, basePtr, dll = dll),
                                                               characterScalar = getSetCharacterScalar(name, value, basePtr, dll = dll),
                                                               numericVector ={##message('switch numericVector');
                                                                   getSetNumericVector(name, value, basePtr, dll = dll)},
                                                               doubleScalar = getSetDoubleScalar(name, value, basePtr, dll = dll),
                                                               integerScalar = getSetIntegerScalar(name, value, basePtr, dll = dll),
                                                               logicalScalar = getSetLogicalScalar(name, value, basePtr, dll = dll),
                                                               'Could not get or set a value for this variable')
                                                 ans
                                             }
                                         ))

#' get or set value of member data from a compiled nimbleFunction using a multi-interface
#'
#' Most nimbleFunctions written for direct user interaction allow standard R-object-like access to member data using \code{$} or \code{`[[`}.  However, sometimes compiled nimbleFunctions contained within other compiled nimbleFunctions are interfaced with a light-weight system called a multi-interface.  \code{valueInCompiledNimbleFunction} provides a way to get or set values in such cases.
#'
#' @param cnf Compiled nimbleFunction object
#'
#' @param name Name of the member data
#'
#' @param value If provided, the value to assign to the member data.  If omitted, the value of the member data is returned.
#'
#' @author Perry de Valpine
#'
#' @details The member data of a nimbleFunction are the objects created in \code{setup} code that are used in \code{run} code or other member functions.
#'
#' Whether multi-interfaces are used for nested nimbleFunctions is controlled by the \code{buildInterfacesForCompiledNestedNimbleFunctions} option in \link{nimbleOptions}.
#'
#' To see an example of a multi-interface, see \code{samplerFunctions} in a compiled MCMC interface object.
#'
#' @export
valueInCompiledNimbleFunction <- function(cnf, name, value) { ## value can be missing
    cnf[[1]]$memberData(cnf[[2]], name, value) 
}

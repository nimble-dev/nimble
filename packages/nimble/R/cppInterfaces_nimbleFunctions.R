
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
                                        # error. See initialization function inside buildNimbleFxnInterface
    vNames <- if(missing(cppNames)) names(symTab$symbols) else cppNames
    for(vn in vNames) {
        thisSymbol <- symTab$getSymbolObject(vn)
        if(is.null(thisSymbol)) next
        if(thisSymbol$type == 'model' || thisSymbol$type == 'symbolNodeFunctionVector' || thisSymbol$type == 'symbolModelVariableAccessorVector' ||thisSymbol$type == 'symbolModelValuesAccessorVector' || thisSymbol$type == 'symbolCopierVector' || thisSymbol$type == 'symbolIndexedNodeInfoTable') next ## skip models and NodeFunctionVectors and modelVariableAccessors      
        ptrName = paste0(".", vn, "_Ptr")
        fieldList[[ptrName]] <- "ANY" ## "ANY" 
        ## Model variables:
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
        if(inherits(thisSymbol, 'symbolNimbleList')) { ## copy type 'nimbleList'
            className <- thisSymbol$nlProc$cppDef$name
           nlName <- paste0(".",vn,"_CnimbleList")
           fieldList[[nlName]] <- "ANY" ## This will have the ref class object that interfaces to the C++ nimbleList
           fieldList[[vn]] <-  eval(substitute(
               function(x) {
                   if(missing(x)) {
                       NLNAME
                   } else {
                       if(is.list(x)) { ## can be a list with first element a CmultiNimbleFunction object and second element an index
                           ## Still need to support this in case an extant nimbleList object was part of compilation using multi-interface
                           if(!inherits(x[[1]], 'CmultiNimbleListClass')) stop(paste('Nimble compilation error initializing pointer for nimbleFunction from a CmultiNimbleList object ', NFNAMECHAR, '.'), call. = FALSE)
                           ptrToPtr <- x[[1]]$ptrToPtrList[[ x[[2]] ]]
                           ##message('setting a nimbleFunction for a multiInterface')
                           nimbleInternalFunctions$setSmartPtrFromDoublePtr(VPTR, ptrToPtr, dll = dll)
                       } else {
                           nimbleInternalFunctions$setSmartPtrFromDoublePtr(VPTR, x$.ptrToPtr, dll = dll) 
                       }
                       ## Even if assigned from a multiInterface, we generate a full interface for correct behavior
                       if(inherits(NLNAME, "uninitializedField")) {
                           nestedRgenerator <- nimbleProject$nlCompInfos[[CLASSNAME]]$cppDef$Rgenerator
                           if(is.list(x)) {
                               newNLinterface <- nestedRgenerator( dll = x[[1]]$dll,
                                                                  existingExtPtrs = list(x[[1]]$ptrToSmartPtrList[[ x[[2]] ]],
                                                                                         x[[1]]$ptrToPtrList[[ x[[2]] ]]) )
                           } else {
                               newNLinterface <- nestedRgenerator( dll = x$dll,
                                                                   existingExtPtrs = list(x$.ptrToSmartPtr, x$.ptrToPtr) )
                           }
                           assign(NLNAMECHAR, newNLinterface, inherits = TRUE) ## avoids <<- warnings
                       }
                       
                       NLNAME$resetExtPtrs(VPTR)
                       x
                   }
               }, list(VPTR = as.name(ptrName), CLASSNAME = className, NLNAME = as.name(nlName), NLNAMECHAR = nlName) ) )
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
                    ##nimbleInternalFunctions$getSetCharacterVector(VPTR, VARNAME, x, dll = dll)
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
            return(NULL)
        }
        warning('Warning in trying to build interface: symbol ', vn, ' not understood.', call. = FALSE)
    }
    return(fieldList)
}

makeNimbleListBindingFields <- function(symTab, cppNames, castFunName) {
    fieldList = list(.DUMMY = "ANY")
    vNames <- if(missing(cppNames)) names(symTab$symbols) else cppNames
    for(vn in vNames) {
        thisSymbol <- symTab$getSymbolObject(vn)
        if(is.null(thisSymbol)) next
        if(thisSymbol$type == 'model' || thisSymbol$type == 'symbolNodeFunctionVector' || thisSymbol$type == 'symbolModelVariableAccessorVector' ||thisSymbol$type == 'symbolModelValuesAccessorVector' || thisSymbol$type == 'symbolCopierVector' || thisSymbol$type == 'symbolIndexedNodeInfoTable') next
        ## if(inherits(thisSymbol,'symbolNimArrDoublePtr')) { ## copy type 'modelVar' ##NOT NEEDED
        ## if(inherits(thisSymbol, 'symbolVecNimArrPtr')){ ## copy type 'modelValuesPtr'##NOT NEEDED
        ##if(inherits(thisSymbol, 'symbolNimbleFunction')) { ##NOT NEEDED
        ##if(inherits(thisSymbol, 'symbolModelValues')) {
        ## NOT NEEDEDif(inherits(thisSymbol, 'symbolNimPtrList')) {
        castFunCall <- parse(text = paste0(".Call(dll$", castFunName, ", .ptrToPtr)"), keep.source = FALSE)[[1]]
        if(inherits(thisSymbol, 'symbolNimbleList')) { ## copy type 'nimbleList'
            className <- thisSymbol$nlProc$cppDef$name
            castToPtrPairName <- thisSymbol$nlProc$cppDef$ptrCastToPtrPairFun$name
            eval(substitute( fieldList$VARNAME <- function(x) {
                namedObjectsPtr <- CASTFUNCALL ##.Call(dll$CASTFUN, .ptrToPtr)
                extPtrNL <- nimbleInternalFunctions$newObjElementPtr(namedObjectsPtr, VARNAME, dll = dll)
                ## extPtrNL <- getMemberDataPtr(VARNAME) ## not working
                nimbleInternalFunctions$getSetNimbleList(vptr = extPtrNL, name = VARNAME, value = x, cppDef = symTab$getSymbolObject(VARNAME)$nlProc$cppDef, dll = dll )
                ## if(missing(x)) {
                ##     nestedRgenerator <- nimbleProject$nlCompInfos[[CLASSNAME]]$cppDef$Rgenerator 
                ##     existingExtPtrs <- .Call(dll$CASTTOPTRPAIRNAME, extPtrNL)
                ##     nestedRgenerator( dll = dll, existingExtPtrs = existingExtPtrs )
                   
                ## } else {
                ##     if(is.list(x)) stop("Can't handle multi-interface assigning a list to a list yet")
                ##     nimbleInternalFunctions$setSmartPtrFromDoublePtr(extPtrNL, x$.ptrToPtr, dll = dll)
                ##     x
                ## }
            }, list(VARNAME = vn, CASTFUNCALL = castFunCall, CLASSNAME = className, CASTTOPTRPAIRNAME = castToPtrPairName) ) ) ##CASTFUN = castFunName,
            next
        }
        if(thisSymbol$type == "character") { ## cpp copy type 'character'  : 2 sub-cases (vector and scalar)

            if(thisSymbol$nDim > 0) {   ## character vector (nDim can only be 0 or 1)
                eval(substitute( fieldList$VARNAME <- function(x) {
                    namedObjectsPtr <- CASTFUNCALL ##.Call(dll$CASTFUN, .ptrToPtr)
                    ##vptr <- nimbleInternalFunctions$newObjElementPtr(namedObjectsPtr, name, dll = dll)
                    ##nimbleInternalFunctions$getSetCharacterVector(vptr, VARNAME, value = x, dll = dll)
                    nimbleInternalFunctions$getSetCharacterVector(VARNAME, value = x, namedObjectsPtr, dll = dll)
                }, list(VARNAME = vn, CASTFUNCALL = castFunCall) ) )
                next
            } else {                    ## character scalar
                eval(substitute( fieldList$VARNAME <- function(x) {
                    namedObjectsPtr <- CASTFUNCALL ##.Call(dll$CASTFUN, .ptrToPtr)
                    nimbleInternalFunctions$getSetCharacterScalar(VARNAME, value = x, namedObjectsPtr, dll = dll)
                }, list(VARNAME = vn, CASTFUNCALL = castFunCall) ) )
                next
            }
        }
        if(inherits(thisSymbol, 'symbolBase')) { ## All numeric and logical cases  ## cpp copy type 'numeric': 4 sub-cases
            if(thisSymbol$nDim > 0) {            ## Anything vector
                eval(substitute( fieldList$VARNAME <- function(x) {
                    namedObjectsPtr <- CASTFUNCALL ##.Call(dll$CASTFUN, .ptrToPtr)
                    nimbleInternalFunctions$getSetNumericVector(VARNAME, value = x, namedObjectsPtr, dll = dll)
                }, list(VARNAME = vn, CASTFUNCALL = castFunCall) ) )
                next
            }
            if(thisSymbol$type == "double"){     ## Scalar double
                eval(substitute( fieldList$VARNAME <- function(x) {
                    namedObjectsPtr <- CASTFUNCALL ##.Call(dll$CASTFUN, .ptrToPtr)
                    nimbleInternalFunctions$getSetDoubleScalar(VARNAME, value = x, namedObjectsPtr, dll = dll)
                }, list(VARNAME = vn, CASTFUNCALL = castFunCall) ) )
                next
            }
           if(thisSymbol$type == "integer"){    ## Scalar int
               eval(substitute( fieldList$VARNAME <- function(x) {
                   namedObjectsPtr <- CASTFUNCALL ##.Call(dll$CASTFUN, .ptrToPtr)
                   nimbleInternalFunctions$getSetIntegerScalar(VARNAME, value = x, namedObjectsPtr, dll = dll)
                }, list(VARNAME = vn, CASTFUNCALL = castFunCall) ) )
                next
            }
            if(thisSymbol$type == "logical"){    ## Scalar logical
                eval(substitute( fieldList$VARNAME <- function(x) {
                    namedObjectsPtr <- CASTFUNCALL ##.Call(dll$CASTFUN, .ptrToPtr)
                    nimbleInternalFunctions$getSetLogicalScalar(VARNAME, value = x, namedObjectsPtr, dll = dll)
                }, list(VARNAME = vn, CASTFUNCALL = castFunCall) ) )
                next
            }
            warning("Warning: scalar datatype not current supported for variable ", vn, "\n", call. = FALSE)
            return(NULL)
        }
    }
    return(fieldList)
}

## This simply looks in the dll and then looks everywhere
## It is a band-aid solution for the problem that predefined nimbleLists can have their functions in the sessionSpecificDll
## but what may be at hand is the project dll
nimbleTryGetNativeSymbolInfo <- function(symName, dll) {
    ans <- try(getNativeSymbolInfo(symName, dll), silent = TRUE)
    if(inherits(ans, 'try-error')) {
        ans <- try(getNativeSymbolInfo(symName), silent  = TRUE)
        if(inherits(ans, 'try-error')) stop(paste0('Unable to find compiled function ', symName,'.'), call.=FALSE)
    }
    ans
}

#### functions that are similar to what is created in makeNFBindingFields but are standalone and look up pointers each time
getSetNimbleList <- function(vptr, name, value, cppDef, dll) {
    ## When missing value, we need the cppDef from the symTab of the assignment target
    ## from this we can get the castFun and the catToPtrPairFun
    ## When receiving value, we don't need anything more 
    if(missing(value)) {
        ## This simply looks in the dll and then looks everywhere
        ## It is a band-aid solution for the problem that predefined nimbleLists can have their functions in the sessionSpecificDll
        ## but what may be at hand is the project dll
        nativeSymInfo <- nimbleTryGetNativeSymbolInfo(cppDef$ptrCastToPtrPairFun$name, dll)
        dllToUse <- if(!is.null(nativeSymInfo$package)) nativeSymInfo$package else dll
        
        existingExtPtrs <- eval(call('.Call', nativeSymInfo, vptr) )
        cppDef$Rgenerator( dll = dllToUse, existingExtPtrs = existingExtPtrs )        
    } else {
        if(is.list(value)) {
            ptrToPtr <- value[[1]]$ptrToPtrList[[ value[[2]] ]]
        } else {
            ptrToPtr <- value$.ptrToPtr
        }
        nimbleInternalFunctions$setSmartPtrFromDoublePtr(vptr, ptrToPtr, dll = dll)
        value
    }
}

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

getSetCharacterVector <- function(name, value, basePtr, dll) { ##basePtr, dll) {
    vptr <- nimbleInternalFunctions$newObjElementPtr(basePtr, name, dll = dll)
     if(missing(value)) 
         nimbleInternalFunctions$getCharacterVectorValue(vptr, dll = dll)
     else
         nimbleInternalFunctions$setCharacterVectorValue(vptr, value, dll = dll)
}

getSetCharacterScalar <- function(name, value, basePtr, dll) {
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
                                       .namedObjectsPtr = 'ANY', ## points to same object as .basePtr but cast to C++ base class (confusing b/c .basePtr points to the C++ derived class type)
                                       .finalizationPtr = 'ANY',
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
                                           .namedObjectsPtr <<- NULL
                                           .finalizationPtr <<- NULL
                                           nimbleProject <<- NULL
                                       },
                                       finalize = function() {
                                           nimbleInternalFunctions$nimbleFinalize(.finalizationPtr)
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

CnimbleListBase <- setRefClass('CnimbleListBase',
                                   fields = list(
                                       dll = "ANY",
                                       .basePtr = 'ANY',
                                       .finalizationPtr = 'ANY',
                                       .ptrToSmartPtr = 'ANY',
                                       .ptrToPtr = 'ANY',
                                      ## compiledNodeFun = 'ANY',
                                       Robject = 'ANY', ## this should be the refClassObject, not the function
                                       symTab = 'ANY'
                                      ## cppNames = 'ANY',
                                      ## cppCopyTypes = 'ANY',
                                      ## neededObjects = 'ANY', ## A list of things like modelValues objects, if they don't already exist
                                       ##nimbleProject = 'ANY'
                                       ),
                                   methods = list(
                                       finalizeInternal = function() {
                                           finalize()
                                           .ptrToPtr <<- NULL
                                           .ptrToSmartPtr <<- NULL
                                           .finalizationPtr <<- NULL
                                       },
                                       finalize = function() {
                                           if(!is.null(.finalizationPtr)) nimbleInternalFunctions$nimbleFinalize(.finalizationPtr)
                                       },
                                       lookupSymbol = function(symname) {
                                           if(is.null(dll))
                                               stop("No DLL for this object")
                                           getNativeSymbolInfo(symname, dll)
                                       },
                                       getMemberDataPtr = function(name) {
                                           nimbleInternalFunctions$newObjElementPtr(getNamedObjectsPtr(), name, dll = dll)
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
        if(inherits(thisSymbol, 'symbolNimbleList')) {ans[[thisSymbol$name]] <- 'nimbleList'; next}
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
    ans <- if(!embedInBrackets) quote(list()) else quote({}) 
    numFuns <- length(compiledNodeFun$RCfunDefs)
    if(numFuns == 0) return(ans)
    funNames <- names(compiledNodeFun$RCfunDefs)
    ##funNames[funNames == 'operator()'] <- 'run'
    for(i in seq_along(compiledNodeFun$RCfunDefs)) {
        ## note that the className is really used as a boolean: any non-NULL value triggers treatment as a class, but name isn't needed
        ans[[i+1]] <- compiledNodeFun$RCfunDefs[[i]]$buildRwrapperFunCode(className = compiledNodeFun$nfProc$name, includeLHS = FALSE, returnArgsAsList = FALSE, includeDotSelf = '.basePtr', includeDotSelfAsArg = includeDotSelfAsArg)
        if(embedInBrackets) ans[[i+1]] <- substitute(THISNAME <- FUNDEF, list(THISNAME = as.name(funNames[i]), FUNDEF = ans[[i+1]]))
    }
    if(!embedInBrackets) names(ans) <- c('', funNames)
    ans
}


clearNeededObjects <- function(Robj, compiledNodeFun, neededObjects) {
    if(inherits(compiledNodeFun$nimCompProc$neededObjectNames, 'uninitializedField'))
      return(NULL)
    for(iName in compiledNodeFun$nimCompProc$neededObjectNames) {
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
        if(is.nl(thisObj)) {
          RCO <- thisObj
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
    for(iName in compiledNodeFun$nimCompProc$neededObjectNames) {
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
        if(is.nl(thisObj)) {
          RCO <- thisObj
          if(inherits(RCO$.CobjectInterface, 'uninitializedField') || is.null(RCO$.CobjectInterface)) {
            neededObjects[[iName]] <- nimbleProject$instantiateNimbleList(thisObj, dll, asTopLevel = getNimbleOption('buildInterfacesForCompiledNestedNimbleFunctions'))
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

copyFromRobjectViaActiveBindings <- function(Robj, cppNames, cppCopyTypes, .self, dll) {
  if(is.nl(Robj)) isNL <- TRUE
  else isNL <- FALSE
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
        else if(cppCopyTypes[[v]] == 'nimbleList') {
          modelVar <- Robj[[v]]
          Cnl <- modelVar$.CobjectInterface 
          .self[[v]] <- Cnl
          next
        }
        # else if(cppCopyTypes[[v]] == 'nimbleList' && isNL){
        #   Cnl <- Robj$.CobjectInterface 
        #   .self[[v]] <- Cnl
        #   next
        # }
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

copyFromRobject <- function(Robj, cppNames, cppCopyTypes, basePtr, symTab, dll,
                            useCompiledCopyMethod = FALSE) {
    if(useCompiledCopyMethod) {
        ## Copy some elements from C++ copyFromRobject method
        ## Currently this includes various numeric types as well as and nodeFxnVector
        ## There is a problem creating the C++ method for predefined nimbleLists,
        ##   so currently useCompiledCopyMethod will be FALSE for (all) nimbleLists

        ## use eval() to avoid R CMD check issue with registration 
        eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$copyFromRobject, basePtr, Robj))
        ## .Call(nimbleUserNamespace$sessionSpecificDll$copyFromRobject, basePtr, Robj)
    }
    
    for(v in cppNames) {
        if(is.null(cppCopyTypes[[v]])) next
        if(is.null(Robj[[v]])) {
            warning("Problem in copyFromRobject.  There is an object to be copied that is NULL.  Going to browser.", call. = FALSE)
            browser()
        }
        switch(cppCopyTypes[[v]],
        'modelVar' = {
            modelVar <- Robj[[v]] ## this is a singleVarAccessClass created by replaceModelSingles
            Cmodel <- modelVar$model$CobjectInterface
            varName <- modelVar$var
            getSetModelVarPtr(v, eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getModelObjectPtr, Cmodel$.basePtr, varName)), basePtr, dll = dll)
        },
        'nimbleFunction' = {
            if(useCompiledCopyMethod) {
                NULL
            } else {
            modelVar <- Robj[[v]]
            Cnf <- nf_getRefClassObject(modelVar)$.CobjectInterface ##environment(modelVar)$.CobjectInterface
            if(is.list(Cnf)) {
                valueBasePtr <- Cnf[[1]]$basePtrList[[ Cnf[[2]] ]]
            } else {
                valueBasePtr <- Cnf$.basePtr
            }
            getSetNimbleFunction(v, valueBasePtr, basePtr, dll = dll)
            }
        },
        'nimbleList' = {
            modelVar <- Robj[[v]]
            Cnl <- modelVar$.CobjectInterface
            vptr <- nimbleInternalFunctions$newObjElementPtr(basePtr, v, dll = dll)
            cppDef <- symTab$getSymbolObject(v)$nlProc$cppDef
            getSetNimbleList(vptr, v, Cnl, cppDef, dll = dll)
        },
        'nimPtrList' = {
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
            getSetNimPtrList(v, Cmv, basePtr, dll = dll)
        },
        'modelValues' = { ## somewhat similar to modelVar
            rModelValues <- Robj[[v]]
            Cmv <- rModelValues$CobjectInterface
            if(!Cmv$initialized) {
                k = getsize(rModelValues)
                resize(Cmv, k)
                vNames = rModelValues[['varNames']]
                for(vN in vNames)
                    Cmv[vN,] <- rModelValues[vN,]
                Cmv$symTab <- rModelValues$symTab
                Cmv$initialized <- TRUE
            }
            getSetModelValues(v, Cmv, basePtr, dll = dll)
        },
        'nodeFxnVec' = {
            if(useCompiledCopyMethod) {
                NULL
            } else {
                populateNodeFxnVecNew(fxnPtr = basePtr, Robject = Robj, fxnVecName = v, dll = dll)
            }
        },
        'modelVarAccess' = {
            populateManyModelVarMapAccess(fxnPtr = basePtr, Robject = Robj, manyAccessName = v, dll = dll)
        },
        'modelValuesAccess' = {
            populateManyModelValuesMapAccess(fxnPtr = basePtr, Robject = Robj, manyAccessName = v, dll = dll)
        },
        "modelValuesPtr" = {
            curObj <- Robj[[v]]
            mvPtr = curObj$modelValues$CobjectInterface$componentExtptrs[[curObj$var]]
            getSetModelValuesPtr(v, mvPtr, basePtr, dll = dll)
        },
        'numericList' = {
            stop('numericList is not working\n')
        },
        'indexedNodeInfoTable' = {
            populateIndexedNodeInfoTable(fxnPtr = basePtr, Robject = Robj, indexedNodeInfoTableName = v, dll = dll)
        },
        'characterVector' = {
            getSetCharacterVector(v, Robj[[v]], basePtr, dll = dll)
        },
        'characterScalar' = {
            getSetCharacterScalar(v, Robj[[v]], basePtr, dll = dll)
        },
        'numericVector' = {
        if(useCompiledCopyMethod) {
                NULL
            } else {
                getSetNumericVector(v, Robj[[v]], basePtr, dll = dll)
            }
        },
        'doubleScalar' = {
         if(useCompiledCopyMethod) {
                NULL
            } else {
                getSetDoubleScalar(v, Robj[[v]], basePtr, dll = dll)
            }
        },
        'integerScalar' = {
        if(useCompiledCopyMethod) {
                NULL
            } else {
                getSetIntegerScalar(v, Robj[[v]], basePtr, dll = dll)
            }
        },
        'logicalScalar' = {
            if(useCompiledCopyMethod) {
                NULL
            } else {
                getSetLogicalScalar(v, Robj[[v]], basePtr, dll = dll)
            }
        },
        { ## default:
            if(!(cppCopyTypes[[v]] %in% c('copierVector'))) {
                warning(paste0("Note: cppCopyTypes not recognized. Type = ", cppCopyTypes[[v]], "\n"), call. = FALSE)
            }
        }
        )
    }
    ## second pass is for initializations that require everything from first pass be done
    for(v in cppNames) {
        if(is.null(cppCopyTypes[[v]])) next
        if(cppCopyTypes[[v]] == 'copierVector') {
            populateCopierVector(fxnPtr = basePtr, Robject = Robj, vecName = v, dll = dll)
        }
    }
}

buildNimbleListInterface <- function(refName,  compiledNimbleObj, basePtrCall, where = globalenv()){
    ## This interface is for a "permanent" nimbleList, like one in nimbleFunction member data or simply global environment
    ## But if the element of a nimbleList is another nimbleList, we have to return that interface dynamically, since it may be ephemeral.
    ## It's tempting to fill in interface objects with their pointers still to be filled, but there would be a danger of infinite recursion
    ## 
    defaults <- list()
    if(inherits(compiledNimbleObj, 'symbolTable')) {
        symTab <- compiledNimbleObj
        defaults$cnf <- NULL
        warning('No compiled node function provided, so interface will be incomplete')
        castFunName <- 'dummyCastingFunction'
    } else {
        symTab <- compiledNimbleObj$nimCompProc$getSymbolTable()
        defaults$cnf <- compiledNimbleObj
        castFunName <- compiledNimbleObj$ptrCastFun$name
        castToPtrPairFunName <- compiledNimbleObj$ptrCastToPtrPairFun$name
    }
    isListObj <- inherits(compiledNimbleObj, 'cppNimbleListClass')
    if(!isListObj) stop('compiledNimbleObj must be a nimbleList')
    ## The following is really equivalent, because it comes *directly* from the place that generates the C++ code
    cppNames <- compiledNimbleObj$objectDefs$getSymbolNames()
    NLBF <-  makeNimbleListBindingFields(symTab, cppNames, castFunName)
    defaults$cppCT <- makeNimbleFxnCppCopyTypes(symTab, cppNames)
    defaults$basePtrCall <- basePtrCall
    defaults$extPtrTypeIndex <- compiledNimbleObj$getExtPtrTypeIndex()
    defaults$cppNames <- cppNames
    defaults$nimbleProject <- compiledNimbleObj$nimbleProject
    defaults$symTab <- symTab
    methodsList <- quote(list())
    fun <- substitute(function(nfObject, defaults, dll = NULL, existingExtPtrs = NULL, ...) {
        callSuper(dll = dll, ...)
       ## nimbleProject <<- defaults$nimbleProject ## may become unnnecessary?
        symTab <<- defaults$symTab
        if(is.null(existingExtPtrs)) {
            basePtrCall <- if(is.character(defaults$basePtrCall)) {
                if(inherits(dll, "uninitializedField") | is.null(dll)) stop("Error making a nimbleFxnInterface object: no dll provided")
                lookupSymbol(defaults$basePtrCall)
            } else defaults$basePtrCall
            ## avoid R CMD check problem with registration.  basePtrCall is already the result of getNativeSymbolInfo from the dll, if possible from cppDefs_nimbleFunction.R
            ## .basePtr
            newObjPtrs <- eval(parse(text = ".Call(basePtrCall)"))
            .ptrToSmartPtr <<- newObjPtrs[[1]] ## nimSmartPtrBase* pointing to a smartPtr<derived_nimbleList_class>
                                        #Use .ptrToSmartPtr to get to smartPtr operations. use for finalizer.
            .ptrToPtr <<- newObjPtrs[[2]]      ## void* that is really a **derived_nimbleList_class
                                        #Call the ptrCastFun with .ptrToPtr to get a pointer cast as NamedObjects* 
 
            .finalizationPtr <<- .ptrToSmartPtr
            eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$register_smartPtrBase_Finalizer,
                      .finalizationPtr, ##.basePtr,
                      dll[['handle']], 'CnimbleList'))
        } else {
            ## Unchecked
            .ptrToSmartPtr <<- existingExtPtrs[[1]]
            .ptrToPtr <<- existingExtPtrs[[2]]  
           ## .namedObjectsPtr <<- existingExtPtrs[[  defaults$extPtrTypeIndex['NamedObjects'] ]]
            .finalizationPtr <<- NULL##.ptrToSmartPtr
        }

        if(!missing(nfObject)) { ## for a nimbleList, nfObject could be validly missing
            if(!is.null(existingExtPtrs)) {
              oldCobjectInterface <- nfObject$.CobjectInterface
              if(!is.list(oldCobjectInterface)) stop('Problem promoting nimbleFunction interface from CmultiInterface to full interface')
          }
            ##setRobject(nfObject)
            Robject <<- nfObject
            Robject$.CobjectInterface <<- .self 
          ##buildNeededObjects()
          
          if(!is.null(existingExtPtrs)) {
              ## unchecked
              ##  neededObjects <<- oldCobjectInterface[[1]]$getNeededObjects( oldCobjectInterface[[2]] )
              oldCobjectInterface[[1]]$clearInstance( oldCobjectInterface[[2]] )
          } else {
              nimbleInternalFunctions$copyFromRobjectViaActiveBindings(Robject, defaults$cppNames, defaults$cppCT, .self, dll) 
          }
        }
    }, list(ABC = NULL))
      # if we just have the name of the routine and haven't resolved it, arrange to resolve it when this initialization
      # function is called.  So change the .Call('name') to .Call(lookupSymbol('name')) which will use this objects
      # dll field.

    methodsList[[length(methodsList) + 1]] <- fun
    names(methodsList)[length(methodsList)] <- 'initialize'
    className <- compiledNimbleObj$name
    methodsList[[length(methodsList) + 1]] <- substitute(function() {
        writeLines(paste0("Derived CnimbleListBase object (compiled nimbleList) created by buildNimbleListInterface for nimbleList with class ", 
                          CLASSNAME))
    }, list(CLASSNAME = className))
    names(methodsList)[length(methodsList)] <- 'show'
    resetExtPtrsFun <- substitute(
        function(VPTR) {
            newPtrs <- RESETPTRCALL ##.Call(dll$RESETPTRFUN, VPTR)
            .ptrToSmartPtr <<- newPtrs[[1]]
            .ptrToPtr <<- newPtrs[[2]]
            .finalizationPtr <<- NULL
        }, list(RESETPTRCALL = parse(text = paste0(".Call(dll$",castToPtrPairFunName,", VPTR)"), keep.source=FALSE )[[1]] )##RESETPTRFUN = castToPtrPairFunName)
    )
    methodsList[[length(methodsList) + 1]] <- resetExtPtrsFun
    names(methodsList)[length(methodsList)] <- 'resetExtPtrs'
    getNamedObjectsPtr <- substitute(
        function() {
            CASTFUNCALL ##.Call(dll$CASTFUNNAME, .ptrToPtr)
        }, list(CASTFUNCALL = parse(text = paste0(".Call(dll$", castFunName, ", .ptrToPtr)"), keep.source = FALSE)[[1]] ))##CASTFUNNAME = castFunName))
    methodsList[[length(methodsList) + 1]] <- getNamedObjectsPtr
    names(methodsList)[length(methodsList)] <- 'getNamedObjectsPtr'

    eval(substitute( newClass <-  setRefClass(refName,
                                              fields = FIELDS,
                                              contains = 'CnimbleListBase',
                                              methods = ML,
                                              where = where),
                    list(FIELDS = NLBF, ML = methodsList ) ) )

    ans <- function(nfObject, dll = NULL, project, existingExtPtrs = NULL) {
    	newClass$new(nfObject, defaults, dll = dll, existingExtPtrs = existingExtPtrs) ##get project from defaults now
    }
    return(ans)
}


buildNimbleObjInterface <- function(refName,  compiledNimbleObj, basePtrCall, where = globalenv()){
    defaults <- list()
    if(inherits(compiledNimbleObj, 'symbolTable')) {
        symTab <- compiledNimbleObj
        defaults$cnf <- NULL
        warning('No compiled node function provided, so interface will be incomplete')
    } else {
        symTab <- compiledNimbleObj$nimCompProc$getSymbolTable()
        defaults$cnf <- compiledNimbleObj
    }
    ## The following is really equivalent, because it comes *directly* from the place that generates the C++ code
    cppNames <- compiledNimbleObj$objectDefs$getSymbolNames()
    NFBF <-  makeNFBindingFields(symTab, cppNames)
    defaults$cppCT <- makeNimbleFxnCppCopyTypes(symTab, cppNames)
    defaults$basePtrCall <- basePtrCall
    defaults$extPtrTypeIndex <- compiledNimbleObj$getExtPtrTypeIndex()
    nlClassName <- compiledNimbleObj$name
    
    methodsList <- makeNimbleFxnInterfaceCallMethodCode(compiledNimbleObj) ##, compiledNodeFun$nfProc)
    # substitute on parsed text string to avoid CRAN issues with .Call registration
    fun <- substitute(function(nfObject, defaults, dll = NULL, project = NULL, existingExtPtrs = NULL, ...){		#cModel removed from here
        defaults$cnf$nfProc$evalNewSetupLinesOneInstance(nfObject, check = TRUE)
        callSuper(dll = dll, project = project, test = FALSE, ...)
        
        if(is.null(existingExtPtrs)) {
            basePtrCall <- if(is.character(defaults$basePtrCall)) {
                               if(inherits(dll, "uninitializedField") | is.null(dll)) stop("Error making a nimbleFxnInterface object: no dll provided")
                               lookupSymbol(defaults$basePtrCall)
                           } else defaults$basePtrCall
                                        # avoid R CMD check problem with registration.  basePtrCall is already the result of getNativeSymbolInfo from the dll, if possible from cppDefs_nimbleFunction.R
            ## .basePtr
            regLabel <- try(get('name', envir = nfObject), silent = TRUE)
            if(inherits(regLabel, 'try-error') | is.null(regLabel)) regLabel <- environment(nfObject)$className

            newObjPtrs <- eval(parse(text = ".Call(basePtrCall)"))
            .basePtr <<- newObjPtrs[[1]] ## pointer to *derived* C++ class
            .namedObjectsPtr <<- newObjPtrs[[  defaults$extPtrTypeIndex['NamedObjects'] ]]
            .finalizationPtr <<- .namedObjectsPtr
            eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$register_namedObjects_Finalizer,
                      .finalizationPtr, ##.basePtr,
                      dll[['handle']], regLabel))
        } else {
            .basePtr <<- existingExtPtrs[[1]]
            .namedObjectsPtr <<- existingExtPtrs[[  defaults$extPtrTypeIndex['NamedObjects'] ]]
            if(is.null(.namedObjectsPtr)) stop('Error finding correct pointers')
            .finalizationPtr <<- .namedObjectsPtr
        }
                                        # .basePtr <<- .Call(basePtrCall)
        cppNames <<- eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getAvailableNames, .namedObjectsPtr))##.basePtr))
        cppCopyTypes <<- defaults$cppCT
        compiledNodeFun <<- defaults$cnf
        vPtrNames <- 	paste0(".", cppNames, "_Ptr")
        for(vn in seq_along(cppNames) ){
            .self[[vPtrNames[vn]]] <- nimbleInternalFunctions$newObjElementPtr(.namedObjectsPtr, cppNames[vn], dll = dll) ##.basePtr
        }
        if(!missing(nfObject)) { ## I don't know when nfObject could be missing in a correct usage
            if(!is.null(existingExtPtrs)) {
                oldCobjectInterface <- nfObject$.CobjectInterface
                if(!is.list(oldCobjectInterface)) stop('Problem promoting nimbleFunction interface from CmultiInterface to full interface')
            }
            setRobject(nfObject)
            ##buildNeededObjects()
            
            if(!is.null(existingExtPtrs)) {
                neededObjects <<- oldCobjectInterface[[1]]$getNeededObjects( oldCobjectInterface[[2]] )
                oldCobjectInterface[[1]]$clearInstance( oldCobjectInterface[[2]] )
            } 
            else {
                neededObjects <<- nimbleInternalFunctions$buildNeededObjects(Robject, compiledNodeFun, neededObjects, dll, nimbleProject)
                nimbleInternalFunctions$copyFromRobjectViaActiveBindings(Robject, cppNames, cppCopyTypes, .self, dll)
          }
        }
    }, list())
      # if we just have the name of the routine and haven't resolved it, arrange to resolve it when this initialization
      # function is called.  So change the .Call('name') to .Call(lookupSymbol('name')) which will use this objects
      # dll field.

    methodsList[[length(methodsList) + 1]] <- fun
    names(methodsList)[length(methodsList)] <- 'initialize'
    showTxt <-  "Function"
    methodsList[[length(methodsList) + 1]] <- substitute(function() {
        writeLines(paste0("Derived CnimbleFunctionBase object (compiled nimbleFunction) for nimbleFunction with class ", 
                           CLASSNAME))
    }, list(CLASSNAME = nlClassName)) ## former subs removed, left substitute call for future modifications
    names(methodsList)[length(methodsList)] <- 'show'
    eval(substitute( newClass <-  setRefClass(refName,
                                              fields = FIELDS,
                                              contains = 'CnimbleFunctionBase',
                                              methods = ML,
                                              where = where),
                    list(FIELDS = NFBF, ML = methodsList ) ) )

    ans <- function(nfObject, dll = NULL, project, existingExtPtrs = NULL) {
    	wrappedInterfaceBuild <- newClass$new
    	wrappedInterfaceBuild(nfObject, defaults, dll = dll, project = project, existingExtPtrs = existingExtPtrs) ## Only purpose of wrappedInterfaceBuild is to have a helpful name for Rprof that is not "new"
#        newClass$new(nfObject, defaults, dll = dll, project = project)
    }
    return(ans)
}



####
## New class for interfacing multiple compiledNimbleFunctions of the same class

CmultiNimbleObjClass <- setRefClass('CmultiNimbleObjClass',
                                    fields = list(nimbleProject = 'ANY',
                                                  finalizationPtrList = 'ANY',
                                                  cppNames = 'ANY',
                                                  cppCopyTypes = 'ANY',
                                                  cppNamesOneByOne = 'ANY',
                                                  cppCopyTypesOneByOne = 'ANY',
                                                  basePtrCall = 'ANY',
                                                  dll = 'ANY',
                                                  RobjectList = 'ANY',
                                                  compiledNodeFun = 'ANY'    ## a cppNimbleFunctionClass or cppNimbleListClass
                                                  ),
                                    methods = list(
                                        initialize = function(compiledNodeFun, basePtrCall, ##copyFromRobjectCall,
                                                              project, ...) { ## need to set dll, nimbleProject
                                            nimbleProject <<- project
                                            finalizationPtrList <<- list()
                                            RobjectList <<- list()
                                            dll <<- NULL
                                            compiledNodeFun <<- compiledNodeFun
                                            ## basePtrCall is the result of getNativeSymbolInfo with the dll if possible from cppDefs_nimbleFunction.R
                                            basePtrCall <<- basePtrCall
                                            callSuper(...)
                                            symTab <- compiledNodeFun$nimCompProc$getSymbolTable()
                                            cppNames <<- compiledNodeFun$objectDefs$getSymbolNames()
                                            cppCopyTypes <<- makeNimbleFxnCppCopyTypes(symTab, cppNames)
                                        },
                                        finalize = function() {
                                            for(i in seq_along(finalizationPtrList)) { 
                                                if(!is.null(finalizationPtrList[[i]])) {
                                                    nimbleInternalFunctions$nimbleFinalize(finalizationPtrList[[i]])
                                                }
                                            }
                                        },
                                        memberDataInternal = function(basePtr, index, name, value) { ## value can be missing
                                            ## This isn't very useful as written for many names because it just gets and sets the external pointers.  It doesn't wrap them in an interface objec
                                            ans <- switch(cppCopyTypes[[name]],
                                                          modelVar = {##message('switch modelVar');
                                                              getSetModelVarPtr(name, value, basePtr, dll = dll)}, ## only makes sense internally
                                                          nimbleFunction ={##message('switch nimbleFunction');
                                                              getSetNimbleFunction(name, value, basePtr, dll = dll)}, ## ditto
                                                          nimbleList ={##message('switch nimbleList');
                                                              valueSymbol <- compiledNodeFun$nimCompProc$getSymbolTable()$getSymbolObject(name)
                                                              vptr <- nimbleInternalFunctions$newObjElementPtr(basePtr, name, dll = dll)
                                                              getSetNimbleList(vptr, value, valueSymbol$nlProc$cppDef, dll = dll)
                                                          }, ## ditto
                                                          nimPtrList = {##message('switch nimPtrList');
                                                              getSetNimPtrList(name, value, basePtr, dll = dll)}, ## ditto
                                                          modelValues = {##message('switch modelValues');
                                                              getSetModelValues(name, value, basePtr, dll = dll)}, ## ditto
                                                          characterVector = {
##                                                              vptr <- nimbleInternalFunctions$newObjElementPtr(basePtr, name, dll = dll)
                                                              getSetCharacterVector(name, value, basePtr, dll = dll)
                                                          },
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
                                    

CmultiNimbleFunctionClass <- setRefClass('CmultiNimbleFunctionClass',
                                         contains = 'CmultiNimbleObjClass',
                                         fields = list(
                                             basePtrList = 'ANY',         ## List of pointers cast as derived C++ class
                                             namedObjectsPtrList = 'ANY', ## List of pointers cast as base C++ NamedObjects class
                                             neededObjectsList = 'ANY',
                                             extPtrTypeIndex = 'ANY',
                                             callEnv = 'ANY'
                                         ),
                                         methods = list(
                                             show = function() {
                                                 cat(paste0('CmultiNimbleFunctionClass object\n'))
                                             },
                                             initialize = function(compiledNodeFun, basePtrCall,
                                                                   project, ...) { ## need to set dll, nimbleProject
                                                 callSuper(compiledNodeFun = compiledNodeFun,
                                                           basePtrCall = basePtrCall,
                                                           project = project, ...)
                                                 boolCopyOneByOne <- !(as.character(cppCopyTypes) %in% c('nimbleFunction',
                                                                                                         'nodeFxnVec',
                                                                                                         'numericVector',
                                                                                                         'doubleScalar',
                                                                                                         'integerScalar',
                                                                                                         'logicalScalar'))
                                                 namesForCopyOneByOne <- names(cppCopyTypes)[ boolCopyOneByOne ]
                                                 cppNamesOneByOne <<- cppNames[ cppNames %in% namesForCopyOneByOne ]
                                                 cppCopyTypesOneByOne <<- cppCopyTypes[boolCopyOneByOne]
                                                 
                                                 neededObjectsList <<- list()
                                                 basePtrList <<- list()
                                                 namedObjectsPtrList <<- list()
                                                 extPtrTypeIndex <<- compiledNodeFun$getExtPtrTypeIndex()
                                                 callCode <- makeNimbleFxnInterfaceCallMethodCode(compiledNodeFun, includeDotSelfAsArg = TRUE, embedInBrackets = TRUE)
                                                 callEnv <<- new.env()
                                                 eval(callCode, envir = callEnv)
                                             },
                                             finalizeInstance = function(index) {
                                                 if(!is.null(finalizationPtrList[[index]])) { ## previously basePtrList
                                                     neededObjectsList[[index]] <<- nimbleInternalFunctions$clearNeededObjects(RobjectList[[index]], compiledNodeFun, neededObjectsList[[index]])
                                                     RobjectList[index] <<- list(NULL)
                                                     nimbleInternalFunctions$nimbleFinalize(finalizationPtrList[[index]])
                                                     basePtrList[index] <<- list(NULL)
                                                     namedObjectsPtrList[index] <<- list(NULL)
                                                     finalizationPtrList[index] <<- list(NULL)
                                                 }
                                             },
                                             addInstance = function(nfObject, dll = NULL) { ## role of initialize
                                                 if(!is.null(.self$dll)) {
                                                     if(!identical(dll, .self$dll)) stop('Can not addInstance of a compiled nimbleFunction from different DLLs', call. = FALSE)
                                                 } else {
                                                     if(is.null(dll)) stop('In addInstance, DLL was not set and so must be provided when calling', call. = FALSE)
                                                     dll <<- dll       ## should only occur first time addInstance is called
                                                     if(is.character(basePtrCall)) {
                                                         updatedBasePtrCall <- try( getNativeSymbolInfo(basePtrCall, dll) )
                                                         if(!inherits(updatedBasePtrCall, 'try-error'))
                                                             basePtrCall <<- updatedBasePtrCall
                                                     }
                                                 }
                                                 isNF <- is.nf(nfObject)
                                                 if(isNF) compiledNodeFun$nimCompProc$evalNewSetupLinesOneInstance(nfObject, check = TRUE)
                                                 # avoid R CMD check problem with registration:
                                                 newObjPtrs <- eval(call('.Call', basePtrCall))
                                                 newBasePtr <- newObjPtrs[[1]] ## terminology confusing because this is the derived C++ class
                                                 newNamedObjectsPtr <- newObjPtrs[[ extPtrTypeIndex['NamedObjects'] ]] ## this is the base C++ class
                                                 if(is.null(newNamedObjectsPtr)) stop('Problem: Cannot find right external pointer information')
                                                 regLabel <- try(get('name', envir = nfObject), silent = TRUE)
                                                 if(inherits(regLabel, 'try-error') | is.null(regLabel)) regLabel <- environment(nfObject)$className
                                                 
                                                 basePtrList[[length(basePtrList)+1]] <<- newBasePtr
                                                 namedObjectsPtrList[[length(namedObjectsPtrList)+1]] <<- newNamedObjectsPtr
                                                                                                  
                                                 finalizationPtrList[[length(finalizationPtrList)+1]] <<- newNamedObjectsPtr
                                                 eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$register_namedObjects_Finalizer,
                                                           newNamedObjectsPtr, 
                                                           dll[['handle']], regLabel))
                                                 
                                                 if(isNF) newRobject <- nimble:::nf_getRefClassObject(nfObject)
                                                 else newRobject <- nfObject
                                                 newRobject$.CobjectInterface <- list(.self, length(basePtrList)) ## second element is its index
                                                 RobjectList[[length(RobjectList)+1]] <<- newRobject

                                                 newNeededObjects <- nimbleInternalFunctions$buildNeededObjects(newRobject, compiledNodeFun, list(), dll, nimbleProject)
                                                 neededObjectsList[[length(neededObjectsList) + 1]] <<- newNeededObjects
                                                                                                  
                                                 nimble:::copyFromRobject(newRobject,
                                                                          cppNamesOneByOne,
                                                                          cppCopyTypesOneByOne,
                                                                          newNamedObjectsPtr,
                                                                          symTab = compiledNodeFun$nimCompProc$getSymbolTable(),
                                                                          dll,
                                                                          useCompiledCopyMethod = TRUE) 
                                                 if(getNimbleOption('clearNimbleFunctionsAfterCompiling')) compiledNodeFun$nfProc$clearSetupOutputs(newRobject)
                                                 list(.self, length(basePtrList)) ## (this object, index)
                                             },
                                             clearInstance = function(index) { ## this is called when a Cmulti interface was built and later needs to be replaced with a full interface
                                                 ## But we cannot remove entries from the lists because the indices of remaining objects must not change
                                                 ## therefore we use the mylist[i] <- list(NULL) method
                                                 basePtrList[index] <<- list(NULL)
                                                 namedObjectsPtrList[index] <<- list(NULL)
                                                 RobjectList[index] <<- list(NULL)
                                                 neededObjectsList[index] <<- list(NULL)
                                                 invisible(NULL)
                                             },
                                             getNeededObjects = function(index) {
                                                 neededObjectsList[[index]]
                                             },
                                             getExtPtrs = function(index) {
                                                 list(basePtrList[[index]], namedObjectsPtrList[[index]])
                                             },
                                             callMemberFunction = function(index, funName, ...) {
                                                 callEnv[[funName]](..., .basePtr = basePtrList[[index]])
                                             },
                                             memberData = function(index, name, value) {
                                                 if(!(name %in% cppNames)) stop(paste0('Name ', name, ' is not a valid member variable in the requested object.', call.=FALSE))
                                                 basePtr <- basePtrList[[index]]
                                                 if(!inherits(basePtr, 'externalptr')) stop('Invalid index or basePtr', call. = FALSE)
                                                 memberDataInternal(basePtr, index, name, value)
                                             }
                                         ))

CmultiNimbleListClass <- setRefClass('CmultiNimbleListClass',
                                         contains = 'CmultiNimbleObjClass',
                                         fields = list(
                                             ptrToPtrList = 'ANY',
                                             ptrToSmartPtrList = 'ANY',
                                             finalizationPtrList = 'ANY',
                                             castFunSymbolInfo = 'ANY'
                                         ),
                                         methods = list(
                                             show = function() {
                                                 cat(paste0('CmultiNimbleListlass object\n'))
                                             },
                                             initialize = function(compiledNodeFun, basePtrCall, project, ...) { ## need to set dll, nimbleProject
                                                 callSuper(compiledNodeFun = compiledNodeFun, basePtrCall = basePtrCall, project = project, ...)
                                                 ptrToPtrList <<- list()
                                                 ptrToSmartPtrList <<- list()
                                             },
                                             finalizeInstance = function(index) {
                                                 if(!is.null(finalizationPtrList[[index]])) { ## previously basePtrList
                                                     RobjectList[index] <<- list(NULL)
                                                     nimbleInternalFunctions$nimbleFinalize(finalizationPtrList[[index]])
                                                     ptrToPtrList[index] <<- list(NULL)
                                                     ptrToSmartPtrList[index] <<- list(NULL)
                                                 }
                                             },
                                             addInstance = function(nfObject, dll = NULL) { ## role of initialize
                                                 if(!is.null(.self$dll)) {
                                                     if(!identical(dll, .self$dll)) stop('Can not addInstance of a compiled nimbleFunction from different DLLs', call. = FALSE)
                                                 } else {
                                                     if(is.null(dll)) stop('In addInstance, DLL was not set and so must be provided when calling', call. = FALSE)
                                                     dll <<- dll       ## should only occur first time addInstance is called
                                                     
                                                     if(is.character(basePtrCall)) {
                                                         updatedBasePtrCall <- try( getNativeSymbolInfo(basePtrCall, dll) )
                                                         if(!inherits(updatedBasePtrCall, 'try-error'))
                                                             basePtrCall <<- updatedBasePtrCall
                                                     }
                                                     castFunSymbolInfo <<- getNativeSymbolInfo(compiledNodeFun$ptrCastFun$name, dll)
                                                 }
                                                 newObjPtrs <-  eval(call('.Call', basePtrCall))
                                                 newPtrToSmartPtr <- newObjPtrs[[1]] ## terminology confusing because this is the derived C++ class
                                                 regLabel <- try(get('name', envir = nfObject), silent = TRUE)
                                                 if(inherits(regLabel, 'try-error') | is.null(regLabel)) regLabel <- environment(nfObject)$className
                                                 
                                                 ptrToSmartPtrList[[length(ptrToSmartPtrList)+1]] <<- newPtrToSmartPtr
                                                 ptrToPtrList[[length(ptrToPtrList)+1]] <<- newPtrToPtr <- newObjPtrs[[2]]
                                                 
                                                 newFinalizationPtr <- newPtrToSmartPtr
                                                 if(is.null(newFinalizationPtr)) stop('Problem: Cannot find right external pointer information')
                                                 finalizationPtrList[[length(finalizationPtrList)+1]] <<- newFinalizationPtr
                                                 eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$register_smartPtrBase_Finalizer,
                                                           newFinalizationPtr, 
                                                           dll[['handle']], regLabel))
                                                 newRobject <- nfObject
                                                 newRobject$.CobjectInterface <- list(.self, length(ptrToPtrList)) ## second element is its index
                                                 RobjectList[[length(RobjectList)+1]] <<- newRobject
                                                 namedObjectsPtr <- eval(call('.Call',castFunSymbolInfo, newPtrToPtr))
                                                 nimble:::copyFromRobject(newRobject,
                                                                          cppNames,
                                                                          cppCopyTypes,
                                                                          namedObjectsPtr,
                                                                          symTab = compiledNodeFun$nimCompProc$getSymbolTable(),
                                                                          dll,
                                                                          useCompiledCopyMethod = FALSE)
                                                 list(.self, length(ptrToSmartPtrList)) ## (this object, index)
                                             },
                                             clearInstance = function(index) { ## this is called when a Cmulti interface was built and later needs to be replaced with a full interface
                                                 ## But we cannot remove entries from the lists because the indices of remaining objects must not change
                                                 ## therefore we use the mylist[i] <- list(NULL) method
                                                 ptrToPtrList[index] <<- list(NULL)
                                                 ptrToSmartPtrList[index] <<- list(NULL)
                                                 RobjectList[index] <<- list(NULL)
                                                 invisible(NULL)
                                             },
                                             getNamedObjectsPtr = function(index) {
                                                 eval(call('.Call', castFunSymbolInfo, ptrToPtrList[[index]]))
                                             },
                                             getMemberDataPtr = function(index, name) {
                                                 nimbleInternalFunctions$newObjElementPtr(getNamedObjectsPtr(index), name, dll = dll)
                                             },
                                             memberData = function(index, name, value) {
                                                 if(!(name %in% cppNames)) stop(paste0('Name ', name, ' is not a valid member variable in the requested object.', call.=FALSE))
                                                 ##ptrToPtr <- ptrToPtrList[[index]]
                                                 ##if(!inherits(ptrToPtr, 'externalptr')) stop('Invalid index or ptrToPtr', call. = FALSE)
                                                 ##namedObjectsPtr <- .Call(castFunSymbolInfo, .ptrToPtr)
                                                 vptr <- getMemberDataPtr(index, name) ##nimbleInternalFunctions$newObjElementPtr(namedObjectsPtr, name, dll = dll)
                                                 memberDataInternal(vptr, index, name, value)
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
#' Whether multi-interfaces are used for nested nimbleFunctions is controlled by the \code{buildInterfacesForCompiledNestedNimbleFunctions} option in \code{\link{nimbleOptions}}.
#'
#' To see an example of a multi-interface, see \code{samplerFunctions} in a compiled MCMC interface object.
#'
#' @export
valueInCompiledNimbleFunction <- function(cnf, name, value) { ## value can be missing
    cnf[[1]]$memberData(cnf[[2]], name, value) 
}


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
    fieldList = list(.basePtr = "ANY", .DUMMY = "ANY")  # We use this .DUMMY field to trick R into not mistakenly printing
                                        # an error. See initialization function inside buildNimbleFxnInterface
    vNames <- if(missing(cppNames)) names(symTab$symbols) else cppNames
    for(vn in vNames) {
        thisSymbol <- symTab$getSymbolObject(vn)
        if(is.null(thisSymbol)) next
        if(thisSymbol$type == 'model' || thisSymbol$type == 'symbolNodeFunctionVector' || thisSymbol$type == 'symbolModelVariableAccessorVector' ||thisSymbol$type == 'symbolModelValuesAccessorVector' || thisSymbol$type == 'symbolCopierVector') next ## skip models and NodeFunctionVectors and modelVariableAccessors      
        ptrName = paste0(".", vn, "_Ptr")
        fieldList[[ptrName]] <- "ANY" ## "ANY" 
        ## Model variables:
        if(inherits(thisSymbol, 'symbolOptimReadyFunction'))	next
        if(inherits(thisSymbol,'symbolNimArrDoublePtr')) {
            fieldList[[vn]] <- eval(substitute(
                function(x) {
                    if(missing(x))
                        VPTR
                    else {
                        if(!inherits(x, 'externalptr')) stop(paste('Nimble compilation error initializing ptr ', VPTRname, '.'), call. = FALSE)
                        setDoublePtrFromSinglePtr(VPTR, x)   
                    }
                }, list(VPTR = as.name(ptrName), VPTRname = ptrName ) ) )
            next
        }
        if(inherits(thisSymbol, 'symbolVecNimArrPtr')){
            fieldList[[vn]] <- eval(substitute(
                function(x) {
                    if(missing(x))
                        VPTR
                    else
                        {
                            if(!inherits(x, 'externalptr')) stop(paste('Nimble compilation error initializing ptr ', VPTRname, '.'), call. = FALSE)
                            setDoublePtrFromSinglePtr(VPTR, x)
                        }
                }, list(VPTR = as.name(ptrName), VPTRname = ptrName ) ) )
            next      	
        }
        
        if(inherits(thisSymbol, 'symbolNumericList') ) {
            fieldList[[vn]] <- 'ANY'
            next	
        }
        
        if(inherits(thisSymbol, 'symbolNimbleFunction')) {
            nfName <- paste0(".",vn,"_CnimbleFunction")
            fieldList[[nfName]] <- "ANY" ## This will have the ref class object that interfaces to the C++ nimbleFunction
            fieldList[[vn]] <- eval(substitute(
                function(x) {
                    if(missing(x))
                        NFNAME
                    else {
                        if(!inherits(x, 'CnimbleFunctionBase')) stop(paste('Nimble compilation error initializing nimbleFunction ', NFNAMECHAR, '.'), call. = FALSE)
                        if(!inherits(x$.basePtr, 'externalptr')) stop(paste('Nimble compilation error initializing pointer for nimbleFunction ', NFNAMECHAR, '.'), call. = FALSE)
                        setDoublePtrFromSinglePtr(VPTR, x$.basePtr) ## check field name
                        assign(NFNAMECHAR, x, inherits = TRUE) ## avoids "<<-" warnings
                    }
                }, list(VPTR = as.name(ptrName), NFNAME = as.name(nfName), NFNAMECHAR = nfName) ) )
            next
        }
        if(inherits(thisSymbol, 'symbolModelValues')) { ## similar behavior to symbolNimArrDoublePtr
            mvName <- paste0(".", vn, "_CmodelValues")
            fieldList[[mvName]] <- "ANY" ## This will have the CmodelValues object
            fieldList[[vn]] <- eval(substitute(
                function(x) {
                    if(missing(x))
                        MVNAME
                    else {
                        if(!inherits(x, 'CmodelValues')) stop(paste('Nimble compilation error initializing modelVaues ', MVNAMECHAR, '.'), call. = FALSE)
                        if(!inherits(x$extptr, 'externalptr')) stop(paste('Nimble compilation error initializing pointer for modelValues ', MVNAMECHAR, '.'), call. = FALSE)
                        setDoublePtrFromSinglePtr(VPTR, x$extptr)
                        assign(MVNAMECHAR, x, inherits = TRUE) ## THIS WAY TO AVOID refClass "<<-" warnings 
                    }
                }, list(VPTR = as.name(ptrName), MVNAME = as.name(mvName), MVNAMECHAR = mvName ) ) )
            next
        }
        if(inherits(thisSymbol, 'symbolNimPtrList')) { ## for nimbleFunctionList, set up with some partial generality but not all the way
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
                        if(!checkNimbleFunctionListCpp(x)) stop(paste('Nimble compilation error initializing nimbleFunctionList ', NFLNAMECHAR, '.  Something is not valid in this list.  It may be the contains (base class) value of one or more functions in the list.'), call. = FALSE)
                        setPtrVectorOfPtrs(ACCESSPTR, CONTENTSPTR, length(x$contentsList))
                        for(i in seq_along(x$contentsList)) {
                            if(!inherits(x[[i]]$.basePtr, 'externalptr')) stop(paste('Nimble compilation error initializing pointer ', i, ' of nimPointerList (nimbleFunctionList) ', NFLNAMECHAR, '.'), call. = FALSE)
                            setOnePtrVectorOfPtrs(ACCESSPTR, i, x[[i]]$.basePtr)
                        }
                        assign(NFLNAMECHAR, x, inherits = TRUE)
                    }
                }, list(NFLNAME = as.name(nflName), NFLNAMECHAR = nflName, CONTENTSPTR = as.name(ptrName), ACCESSPTR = as.name(accessorPtrName) ) ) )
            ## getter: return the nimPointerList set up with the list of CinterfaceObjects
            ## setter: call setPtrVectorOfPtrs(accessorExtPtr, contentsExtrPtr.  Then iterate and call setOnePtrVectorOfPtrs(accessorPtr, i, nimPtrList[[i]]$.basePtr)
            next
        }
        if(thisSymbol$type == "character") {
            if(thisSymbol$nDim > 0) {   ## character vector (nDim can only be 0 or 1)
                eval(substitute( fieldList$VARNAME <- function(x){
                    if(missing(x) ) 
                        getCharacterVectorValue(VPTR)
                    else
                        setCharacterVectorValue(VPTR, x)
                }, list(VPTR = as.name(ptrName), VARNAME = vn) ) )
            next
            } else {                    ## character scalar
                eval(substitute( fieldList$VARNAME <- function(x){
                    if(missing(x) ) 
                        getCharacterValue(VPTR)
                    else
                        setCharacterValue(VPTR, x)
                }, list(VPTR = as.name(ptrName), VARNAME = vn) ) )
                next
            }
        }
        if(inherits(thisSymbol, 'symbolBase')) { ## All numeric and logical cases
            if(thisSymbol$nDim > 0) {            ## Anything vector
                eval(substitute( fieldList$VARNAME <- function(x){
                    
                    if(missing(x) ) 
                        getNimValues(VPTR)
                    else 
                        setNimValues(VPTR, x)
                }, list(VPTR = as.name(ptrName), VARNAME = vn) ) )
                next
            }
            
            if(thisSymbol$type == "double"){     ## Scalar double
                eval(substitute( fieldList$VARNAME <- function(x){
                    if(missing(x) ) 
                        getDoubleValue(VPTR)
                    else
                       setDoubleValue(VPTR, x)
                        
                }, list(VPTR = as.name(ptrName), VARNAME = vn) ) )
                next
            }
            if(thisSymbol$type == "integer"){    ## Scalar int
                eval(substitute( fieldList$VARNAME <- function(x){
                    if(missing(x) ) 
                        getIntValue(VPTR)
                    else
                        setIntValue(VPTR, x)                        
                }, list(VPTR = as.name(ptrName), VARNAME = vn) ) )
                next
            }
            if(thisSymbol$type == "logical"){    ## Scalar logical
                eval(substitute( fieldList$VARNAME <- function(x){
                    if(missing(x) ) 
                        getBoolValue(VPTR)
                    else
                        setBoolValue(VPTR, x)
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

CnimbleFunctionBase <- setRefClass('CnimbleFunctionBase',
                                   fields = list(
                                       dll = "ANY",
                                       compiledNodeFun = 'ANY',
                                       Robject = 'ANY', ## this should be the refClassObject, not the function
                                       cppNames = 'ANY',
                                       cppCopyTypes = 'ANY',
                                       neededObjects = 'ANY', ## A list of things like modelValues objects, if they don't already exist
                                       nimbleProject = 'ANY'
                                       ),
                                   methods = list(
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
                                       setRobject = function(Robj) {
                                           if(is.nf(Robj)) Robject <<- nf_getRefClassObject(Robj)
                                           else Robject <<- Robj
                                           Robject$.CobjectInterface <<- .self 
                                       },
                                       ## buildNeededObjects = function(Robj) {
                                       ##     if(missing(Robj)) Robj <- Robject
                                       ##     for(iName in compiledNodeFun$nfProc$neededObjectNames) {
                                       ##         thisObj <- Robj[[iName]]
                                       ##         if(inherits(thisObj, 'modelValuesBaseClass')) {
                                       ##             if(inherits(thisObj$CobjectInterface, 'uninitializedField') || is.null(thisObj$CobjectInterface)) {
                                       ##                 neededObjects[[iName]] <<- nimbleProject$instantiateCmodelValues(thisObj, dll)
                                       ##             }
                                       ##             next
                                       ##         }
                                       ##         if(is.nf(thisObj)) {
                                       ##             RCO <- nf_getRefClassObject(thisObj)
                                       ##             if(inherits(RCO$.CobjectInterface, 'uninitializedField') || is.null(RCO$.CobjectInterface)) {
                                       ##                 neededObjects[[iName]] <<- nimbleProject$instantiateNimbleFunction(thisObj, dll)
                                       ##             }
                                       ##             next
                                       ##         }
                                       ##         if(inherits(thisObj, 'nimbleFunctionList')) {
                                       ##             neededObjects[[iName]] <<- nimPointerList(thisObj$baseClass, length(thisObj$contentsList))
                                       ##             for(i in seq_along(thisObj$contentsList)) {
                                       ##                 RCO <- nf_getRefClassObject(thisObj[[i]])
                                       ##                 if(inherits(RCO$.CobjectInterface, 'uninitializedField') || is.null(RCO$.CobjectInterface)) {
                                       ##                     neededObjects[[iName]][[i]] <<- nimbleProject$instantiateNimbleFunction(thisObj[[i]], dll)
                                       ##                 } else {
                                       ##                     neededObjects[[iName]][[i]] <<- RCO$.CobjectInterface 
                                       ##                 }
                                       ##             }
                                       ##             names(neededObjects[[iName]]$contentsList) <<- names(thisObj$contentsList)
                                       ##             thisObj$CobjectInterface <- neededObjects[[iName]]
                                       ##             next
                                       ##         }
                                       ##         warning('Warning: object ',iName,' not of a type that can be built.', call. = FALSE)
                                       ##     }
                                       ## },
                                       ## copyFromRobject = function(Robj) {
                                       ##     if(missing(Robj)) Robj <- Robject
                                       ##     for(v in cppNames) {
                                       ##         if(is.null(cppCopyTypes[[v]])) next
                                       ##         if(is.null(Robj[[v]])) {
                                       ##             warning("Problem in copyFromRobject.  There is an object to be copied that is NULL.  Going to browser.", call. = FALSE)
                                       ##             browser()
                                       ##         }
                                       ##         if(cppCopyTypes[[v]] == 'modelVar') {
                                       ##             modelVar <- Robj[[v]] ## this is a singleVarAccessClass created by replaceModelSingles
                                       ##             Cmodel <- modelVar$model$CobjectInterface
                                       ##             varName <- modelVar$var
                                       ##             .self[[v]] <<- .Call('getModelObjectPtr', Cmodel$.basePtr, varName)
                                       ##             next
                                       ##         }
                                       ##         else if(cppCopyTypes[[v]] == 'nimbleFunction') {
                                       ##             modelVar <- Robj[[v]]
                                       ##             Cnf <- nf_getRefClassObject(modelVar)$.CobjectInterface ##environment(modelVar)$.CobjectInterface
                                       ##             .self[[v]] <- Cnf
                                       ##             next
                                       ##         }
                                       ##         else if(cppCopyTypes[[v]] == 'nimPtrList') {
                                       ##             if(is.null(Robj[[v]]$contentsList)) {
                                       ##                 warning('Problem in copying a nimPtrList to C++ object. The contentsList is NULL. Going to browser', call. = FALSE)
                                       ##                 browser()
                                       ##             }
                                       ##             if(any(unlist(lapply(Robj[[v]]$contentsList, is.null)))) {
                                       ##                 warning('Problem in copying a nimPtrList to C++ object. The contentsList is NULL')
                                       ##                 browser()
                                       ##             }
                                       ##             modelVar <- Robj[[v]] ## This is a nimPtrList 
                                       ##             Cmv <- modelVar$CobjectInterface ## This was created above in build neededObjects
                                       ##             .self[[v]] <- Cmv
                                       ##         }
                                       ##         else if(cppCopyTypes[[v]] == 'modelValues') { ## somewhat similar to modelVar
                                       ##             rModelValues <- Robj[[v]]
                                       ##             Cmv <- rModelValues$CobjectInterface
                                       ##             k = getsize(rModelValues)
                                       ##             resize(Cmv, k)
                                       ##           	vNames = rModelValues[['varNames']]
                                       ##           	for(vN in vNames)
                                       ##           		Cmv[vN,] <- rModelValues[vN,]
                                       ##           	Cmv$symTab <- rModelValues$symTab	
                                       ##             .self[[v]] <- Cmv
                                       ##             next
                                       ##         }
                                       ##         else if(cppCopyTypes[[v]] == 'nodeFxnVec') {
                                       ##             populateNodeFxnVec(fxnPtr = .basePtr, Robject = Robj, fxnVecName = v) 
                                       ##             next
                                       ##         }
                                       ##         else if(cppCopyTypes[[v]] == 'modelVarAccess'){
                                       ##             populateManyModelVarMapAccess(fxnPtr = .basePtr, Robject = Robj, manyAccessName = v)
                                       ##             next
                                       ##         }
                                       ##         else if(cppCopyTypes[[v]] == 'modelValuesAccess'){
                                       ##             populateManyModelValuesMapAccess(fxnPtr = .basePtr, Robject = Robj, manyAccessName = v)
                                       ##             next
                                       ##         }
                                       ##         else if(cppCopyTypes[[v]] == "modelValuesPtr"){
                                       ##             curObj <- Robj[[v]]
                                       ##             mvPtr = curObj$modelValues$CobjectInterface$componentExtptrs[[curObj$var]]
                                       ##             .self[[v]] <<- mvPtr
                                       ##             next
                                       ##         }
                                       ##         else if(cppCopyTypes[[v]] == 'numericList'){
                                       ##             rawPtr = .Call('getModelObjectPtr', .basePtr, v)
                                       ##             .self[[v]] <<- numericList(buildType = 'C', extPtr = rawPtr)
                                       ##             nRows = Robj[[v]]$nRow
                                       ##             resize(.self[[v]], nRows)
                                       ##             for(i in 1:nRows){
                                       ##                 copyDims = max(c(1, dimOrLength[[v]][[i]]) )
                                       ##                 d1 = copyDims[1]
                                       ##                 d2 = copyDims[2]
                                       ##                 d3 = copyDims[3]
                                       ##                 setSize(.self[[v]], row = i, d1, d2, d3)
                                       ##                 .self[[v]][[i]] <<- Robj[[v]][[i]]
                                       ##             }
                                       ##             next
                                       ##         }
                                       ##         else if(cppCopyTypes[[v]] == 'character') {
                                       ##             .self[[v]] <<- Robj[[v]]
                                       ##         }
                                       ##         else if(cppCopyTypes[[v]] == 'numeric') {
                                       ##             .self[[v]] <<- Robj[[v]]
                                       ##         }
                                       ##         else if(!(cppCopyTypes[[v]] %in% c('copierVector'))) {
                                       ##             warning(paste0("Note: cppCopyTypes not recognized. Type = ", cppCopyTypes[[v]], "\n"), call. = FALSE)
                                       ##         }
                                       ##     }
                                       ##     ## second pass is for initializations that require everything from first pass be done
                                       ##     for(v in cppNames) {
                                       ##         if(is.null(cppCopyTypes[[v]])) next
                                       ##         if(cppCopyTypes[[v]] == 'copierVector') {
                                       ##             populateCopierVector(fxnPtr = .basePtr, Robject = Robj, vecName = v)
                                       ##         }
                                       ##     }
                                       ## },
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
        else if(thisSymbol$type == 'Ronly') next ## skip models
        else if(inherits(thisSymbol,'symbolNimArrDoublePtr')) {ans[[thisSymbol$name]] <- 'modelVar'; next}
        else if(inherits(thisSymbol, 'symbolNodeFunctionVector'))  { ans[[thisSymbol$name]] <- 'nodeFxnVec'; next}
        else if(inherits(thisSymbol, 'symbolModelVariableAccessorVector')) {ans[[thisSymbol$name]] <- 'modelVarAccess';next}
        else if(inherits(thisSymbol, 'symbolModelValuesAccessorVector')) {ans[[thisSymbol$name]] <- 'modelValuesAccess';next}
        else if(inherits(thisSymbol, 'symbolModelValues')) {ans[[thisSymbol$name]] <- 'modelValues'; next}
        else if(inherits(thisSymbol, 'symbolNimbleFunction')) {ans[[thisSymbol$name]] <- 'nimbleFunction'; next}
        else if(inherits(thisSymbol, 'symbolVecNimArrPtr')) {ans[[thisSymbol$name]] <- 'modelValuesPtr'; next}
        else if(inherits(thisSymbol, 'symbolNumericList')) {ans[[thisSymbol$name]] <- 'numericList'; next}
        else if(inherits(thisSymbol, 'symbolNimPtrList')) {ans[[thisSymbol$name]] <- 'nimPtrList'; next}
        else if(inherits(thisSymbol, 'symbolCopierVector')) {ans[[thisSymbol$name]] <- 'copierVector'; next}
        else if(inherits(thisSymbol, 'symbolString')) {ans[[thisSymbol$name]] <- 'character'; next}
        else ans[[thisSymbol$name]] <- 'numeric'
    }
    ans
}

makeNimbleFxnInterfaceCallMethodCode <- function(compiledNodeFun) {
    numFuns <- length(compiledNodeFun$RCfunDefs)
    ans <- quote(list()) 
    if(numFuns == 0) return(ans)
    for(i in seq_along(compiledNodeFun$RCfunDefs)) {
        ## note that the className is really used as a boolean: any non-NULL value triggers treatment as a class, but name isn't needed
        ans[[i+1]] <- compiledNodeFun$RCfunDefs[[i]]$buildRwrapperFunCode(className = compiledNodeFun$nfProc$name, includeLHS = FALSE, returnArgsAsList = FALSE, includeDotSelf = '.basePtr')
    }
    names(ans) <- c('', names(compiledNodeFun$RCfunDefs))
    names(ans)[names(ans) == 'operator()'] <- 'run'
    ans
}

buildNeededObjects = function(Robj, compiledNodeFun, neededObjects, dll) {
    for(iName in compiledNodeFun$nfProc$neededObjectNames) {
        thisObj <- Robj[[iName]]
        if(inherits(thisObj, 'modelValuesBaseClass')) {
            if(inherits(thisObj$CobjectInterface, 'uninitializedField') || is.null(thisObj$CobjectInterface)) {
                neededObjects[[iName]] <<- nimbleProject$instantiateCmodelValues(thisObj, dll)
            }
            next
        }
        if(is.nf(thisObj)) {
            RCO <- nf_getRefClassObject(thisObj)
            if(inherits(RCO$.CobjectInterface, 'uninitializedField') || is.null(RCO$.CobjectInterface)) {
                neededObjects[[iName]] <<- nimbleProject$instantiateNimbleFunction(thisObj, dll)
            }
            next
        }
        if(inherits(thisObj, 'nimbleFunctionList')) {
            neededObjects[[iName]] <<- nimPointerList(thisObj$baseClass, length(thisObj$contentsList))
            for(i in seq_along(thisObj$contentsList)) {
                RCO <- nf_getRefClassObject(thisObj[[i]])
                if(inherits(RCO$.CobjectInterface, 'uninitializedField') || is.null(RCO$.CobjectInterface)) {
                    neededObjects[[iName]][[i]] <<- nimbleProject$instantiateNimbleFunction(thisObj[[i]], dll)
                } else {
                    neededObjects[[iName]][[i]] <<- RCO$.CobjectInterface 
                }
            }
            names(neededObjects[[iName]]$contentsList) <<- names(thisObj$contentsList)
            thisObj$CobjectInterface <- neededObjects[[iName]]
            next
        }
        warning('Warning: object ',iName,' not of a type that can be built.', call. = FALSE)
    }
    neededObjects
}

copyFromRobjectViaActiveBindings = function(Robj, cppNames, cppCopyTypes, .self) {
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
            .self[[v]] <<- .Call('getModelObjectPtr', Cmodel$.basePtr, varName)
            next
        }
        else if(cppCopyTypes[[v]] == 'nimbleFunction') {
            modelVar <- Robj[[v]]
            Cnf <- nf_getRefClassObject(modelVar)$.CobjectInterface ##environment(modelVar)$.CobjectInterface
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
            k = getsize(rModelValues)
            resize(Cmv, k)
            vNames = rModelValues[['varNames']]
            for(vN in vNames)
                Cmv[vN,] <- rModelValues[vN,]
            Cmv$symTab <- rModelValues$symTab	
            .self[[v]] <- Cmv
            next
        }
        else if(cppCopyTypes[[v]] == 'nodeFxnVec') {
            populateNodeFxnVec(fxnPtr = .self$.basePtr, Robject = Robj, fxnVecName = v) 
            next
        }
        else if(cppCopyTypes[[v]] == 'modelVarAccess'){
            populateManyModelVarMapAccess(fxnPtr = .self$.basePtr, Robject = Robj, manyAccessName = v)
            next
        }
        else if(cppCopyTypes[[v]] == 'modelValuesAccess'){
            populateManyModelValuesMapAccess(fxnPtr = .self$.basePtr, Robject = Robj, manyAccessName = v)
            next
        }
        else if(cppCopyTypes[[v]] == "modelValuesPtr"){
            curObj <- Robj[[v]]
            mvPtr = curObj$modelValues$CobjectInterface$componentExtptrs[[curObj$var]]
            .self[[v]] <<- mvPtr
            next
        }
        else if(cppCopyTypes[[v]] == 'numericList'){
            rawPtr = .Call('getModelObjectPtr', .self$.basePtr, v)
            .self[[v]] <<- numericList(buildType = 'C', extPtr = rawPtr)
            nRows = Robj[[v]]$nRow
            resize(.self[[v]], nRows)
            for(i in 1:nRows){
                copyDims = max(c(1, dimOrLength[[v]][[i]]) )
                d1 = copyDims[1]
                d2 = copyDims[2]
                d3 = copyDims[3]
                setSize(.self[[v]], row = i, d1, d2, d3)
                .self[[v]][[i]] <<- Robj[[v]][[i]]
            }
            next
        }
        else if(cppCopyTypes[[v]] == 'character') {
            .self[[v]] <<- Robj[[v]]
        }
        else if(cppCopyTypes[[v]] == 'numeric') {
            .self[[v]] <<- Robj[[v]]
        }
        else if(!(cppCopyTypes[[v]] %in% c('copierVector'))) {
            warning(paste0("Note: cppCopyTypes not recognized. Type = ", cppCopyTypes[[v]], "\n"), call. = FALSE)
        }
    }
    ## second pass is for initializations that require everything from first pass be done
    for(v in cppNames) {
        if(is.null(cppCopyTypes[[v]])) next
        if(cppCopyTypes[[v]] == 'copierVector') {
            populateCopierVector(fxnPtr = .self$.basePtr, Robject = Robj, vecName = v)
        }
    }
}

copyFromRobject = function(Robj, cppNames, cppCopyTypes, .self) {
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
            .self[[v]] <<- .Call('getModelObjectPtr', Cmodel$.basePtr, varName)
            next
        }
        else if(cppCopyTypes[[v]] == 'nimbleFunction') {
            modelVar <- Robj[[v]]
            Cnf <- nf_getRefClassObject(modelVar)$.CobjectInterface ##environment(modelVar)$.CobjectInterface
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
            k = getsize(rModelValues)
            resize(Cmv, k)
            vNames = rModelValues[['varNames']]
            for(vN in vNames)
                Cmv[vN,] <- rModelValues[vN,]
            Cmv$symTab <- rModelValues$symTab	
            .self[[v]] <- Cmv
            next
        }
        else if(cppCopyTypes[[v]] == 'nodeFxnVec') {
            populateNodeFxnVec(fxnPtr = .self$.basePtr, Robject = Robj, fxnVecName = v) 
            next
        }
        else if(cppCopyTypes[[v]] == 'modelVarAccess'){
            populateManyModelVarMapAccess(fxnPtr = .self$.basePtr, Robject = Robj, manyAccessName = v)
            next
        }
        else if(cppCopyTypes[[v]] == 'modelValuesAccess'){
            populateManyModelValuesMapAccess(fxnPtr = .self$.basePtr, Robject = Robj, manyAccessName = v)
            next
        }
        else if(cppCopyTypes[[v]] == "modelValuesPtr"){
            curObj <- Robj[[v]]
            mvPtr = curObj$modelValues$CobjectInterface$componentExtptrs[[curObj$var]]
            .self[[v]] <<- mvPtr
            next
        }
        else if(cppCopyTypes[[v]] == 'numericList'){
            rawPtr = .Call('getModelObjectPtr', .self$.basePtr, v)
            .self[[v]] <<- numericList(buildType = 'C', extPtr = rawPtr)
            nRows = Robj[[v]]$nRow
            resize(.self[[v]], nRows)
            for(i in 1:nRows){
                copyDims = max(c(1, dimOrLength[[v]][[i]]) )
                d1 = copyDims[1]
                d2 = copyDims[2]
                d3 = copyDims[3]
                setSize(.self[[v]], row = i, d1, d2, d3)
                .self[[v]][[i]] <<- Robj[[v]][[i]]
            }
            next
        }
        else if(cppCopyTypes[[v]] == 'character') {
            .self[[v]] <<- Robj[[v]]
        }
        else if(cppCopyTypes[[v]] == 'numeric') {
            .self[[v]] <<- Robj[[v]]
        }
        else if(!(cppCopyTypes[[v]] %in% c('copierVector'))) {
            warning(paste0("Note: cppCopyTypes not recognized. Type = ", cppCopyTypes[[v]], "\n"), call. = FALSE)
        }
    }
    ## second pass is for initializations that require everything from first pass be done
    for(v in cppNames) {
        if(is.null(cppCopyTypes[[v]])) next
        if(cppCopyTypes[[v]] == 'copierVector') {
            populateCopierVector(fxnPtr = .self$.basePtr, Robject = Robj, vecName = v)
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
    
    fun <- substitute(function(nfObject, defaults, dll = NULL, project = NULL, ...){		#cModel removed from here
        defaults$cnf$nfProc$evalNewSetupLinesOneInstance(nfObject, check = TRUE) ## in case this instance was not used during nfProc$process
        callSuper(dll = dll, project = project, test = FALSE, ...)
        basePtrCall <- if(is.character(defaults$basePtrCall)) {
            if(inherits(dll, 'uninitializedField') | is.null(dll)) stop('Error making a nimbleFxnInterface object: no dll provided')
            lookupSymbol(defaults$basePtrCall)
        } else defaults$basePtrCall
        .basePtr <<- .Call(basePtrCall)
        cppNames <<- .Call("getAvailableNames", .basePtr)
        cppCopyTypes <<- defaults$cppCT
        compiledNodeFun <<- defaults$cnf
        vPtrNames <- 	paste0('.', cppNames, '_Ptr')	
        for(vn in seq_along(cppNames) ){
            .self[[vPtrNames[vn]]] <- newObjElementPtr(.basePtr, cppNames[vn])
        }
        if(!missing(nfObject)) {
            setRobject(nfObject)
            ##buildNeededObjects()
            neededObjects <<- buildNeededObjects(Robject, compiledNodeFun, neededObjects, dll)
            ##copyFromRobject()
            copyFromRobjectViaActiveBindings(Robject, cppNames, cppCopyTypes, .self)
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

CmultiNimbleFunction <- setRefClass('CmultiNimbleFunction',
                                    fields = list(
                                        cppNames = 'ANY',
                                        cppCopyTypes = 'ANY',
                                        basePtrCall = 'ANY',
                                        dll = 'ANY',
                                       ## nfObjectList = 'ANY',
                                        basePtrList = 'ANY',
                                        RobjectList = 'ANY',
                                        neededObjectsList = 'ANY',
                                        nimbleProject = 'ANY',
                                        compiledNodeFun = 'ANY'    ## a cppNimbleFunctionClass
                                    ),
                                    methods = list(
                                        initialize = function(compiledNodeFun, basePtrCall, ...) { ## need to set dll, nimbleProject
                                            neededObjectsList <<- list()
                                       ##     nfObjectList <<- list()
                                            basePtrList <<- list()
                                            RobjectList <<- list()
                                            compiledNodeFun <<- compiledNodeFun
                                            basePtrCall <<- basePtrCall
                                            callSuper(...)
                                            symTab <- compiledNodeFun$nfProc$setupSymTab
                                            cppNames <<- compiledNodeFun$objectDefs$getSymbolNames()
                                            cppCopyTypes <<- makeNimbleFxnCppCopyTypes(symTab, cppNames)
##                                          methodsList <- makeNimbleFxnInterfaceCallMethodCode(compiledNodeFun) ## can do this but need to pass another pointer in

                                        },
                                        addInstance = function(nfObject, dll = NULL, project = NULL) { ## role of initialize
                                            if(!is.null(.self$dll)) {
                                                if(!identical(dll, .self$dll)) stop('Can not addInstance of a compiled nimbleFunction from different DLLs', call. = FALSE)
                                            } else {
                                                if(!is.null(dll)) stop('In addInstance, DLL was not set and so must be provided when calling', call. = FALSE)
                                                dll <<- dll
                                            }
                                            compiledNimbleFunction$nfProc$evalNewSetupLinesOneInstance(nfObject, check = TRUE)
                                      ##      nfObjectList <<- c(nfObjectList, nfObject)
                                            newBasePtr <- .Call(basePtrCall)
                                            basePtrList[[length(basePtrList)+1]] <<- newBasePtrxo
                                            if(is.nf(nfObject)) newRobject <- nf_getRefClassObject(nfObject)
                                            else newRobject <- nfObject
                                            newRobject$.CobjectInterface <- .self
                                            RobjectList[[length(RobjectList)+1]] <<- newRobject
                                            newNeededObjects <<- buildNeededObjects(newRobject, compiledNodeFun, list(), dll)
                                            neededObjectsList[[length(neededObjectsList) + 1]] <<- newNeededObjects
                                        },
                                        setRobject = function() {},
                                        buildNeededObjects = function() {},
                                        copyFromRobject = function() {},
                                        lookupSymbol = function() {}
                                    ))

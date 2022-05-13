
populateCopierVector <- function(fxnPtr, Robject, vecName, dll) {
    vecPtr <- eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getModelObjectPtr, fxnPtr, vecName))
    copierVectorObject <- Robject[[vecName]]
    fromPtr <- eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getModelObjectPtr, fxnPtr, copierVectorObject[[1]]))
    toPtr <- eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getModelObjectPtr, fxnPtr, copierVectorObject[[2]]))
    eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$populateCopierVector, vecPtr, fromPtr, toPtr, as.integer(copierVectorObject[[3]]), as.integer(copierVectorObject[[4]])))
}

processModelVarAccess <- function(Robject, manyAccessName) {
    if(class(Robject[[manyAccessName]])[1] != 'processedModelVarAccess') {
        cModel <- Robject[[manyAccessName]][[1]]$CobjectInterface
        if(is(cModel, 'uninitializedField'))
            stop('Compiled C++ model not available; please include the model in your compilation call (or compile it in advance).', call. = FALSE)
        mapInfo <- makeMapInfoFromAccessorVectorFaster(Robject[[manyAccessName]])
        Robject[[manyAccessName]] <- structure(list(mapInfo[[1]], mapInfo[[2]], cModel$.basePtr),
                                               class = 'processedModelVarAccess')
        if(isTRUE(nimbleOptions("enableDerivs")))
            if(isTRUE(nimbleOptions('useADreconfigure')))
                if(is.list(cModel$.ADptrs))
                    Robject[[manyAccessName]][[4]] <- cModel$.ADptrs[[1]]
    }
}

processModelValuesAccess <- function(Robject, manyAccessName) {

}

populateManyModelVarMapAccess <- function(fxnPtr, Robject, manyAccessName, dll) { ## new version
    manyAccessPtr = eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getModelObjectPtr, fxnPtr, manyAccessName))
    processModelVarAccess(Robject, manyAccessName)
    eval(call('.Call',
              nimbleUserNamespace$sessionSpecificDll$populateValueMapAccessorsFromNodeNames,
              manyAccessPtr,
              Robject[[manyAccessName]][[1]],
              Robject[[manyAccessName]][[2]],
              Robject[[manyAccessName]][[3]]))
    ## For AD system, we will copy here also the AD case.
    ## This is klugey, but it avoids re-doing  processModelVarAccess
    ## and anticipates that in the future these objects might be
    ## combined more naturally.
    if(isTRUE(nimbleOptions("enableDerivs")))
        if(isTRUE(nimbleOptions('useADreconfigure'))) {
            ADname <- paste0(manyAccessName, "_AD_")
            ## We don't have a good way to know about the name, so we check from the compiled object
            availableNames <- eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getAvailableNames, fxnPtr))
            if(ADname %in% availableNames) {
                manyAccessPtrAD = eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getModelObjectPtr, fxnPtr, ADname))
                eval(call('.Call',
                          nimbleUserNamespace$sessionSpecificDll$populateValueMapAccessorsFromNodeNames,
                          manyAccessPtrAD,
                          Robject[[manyAccessName]][[1]],
                          Robject[[manyAccessName]][[2]],
                          Robject[[manyAccessName]][[4]]))
            }
        }
}

populateManyModelValuesMapAccess <- function(fxnPtr, Robject, manyAccessName, dll){ ## nearly identical to populateManyModelVarMapAccess
    manyAccessPtr = eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getModelObjectPtr, fxnPtr, manyAccessName))
    cModelValues <- Robject[[manyAccessName]][[1]]$CobjectInterface
    mapInfo <- makeMapInfoFromAccessorVectorFaster(Robject[[manyAccessName]])
    eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$populateValueMapAccessorsFromNodeNames, manyAccessPtr, mapInfo[[1]], mapInfo[[2]], cModelValues$extptr))
}

getNamedObjected <- function(objectPtr, fieldName, dll)
    eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getModelObjectPtr, objectPtr, fieldName))

populateNodeFxnVecNew <- function(fxnPtr, Robject, fxnVecName, dll){
    fxnVecPtr <- getNamedObjected(fxnPtr, fxnVecName, dll = dll)
    indexingInfo <- Robject[[fxnVecName]]$indexingInfo
    declIDs <- indexingInfo$declIDs
    rowIndices <- indexingInfo$unrolledIndicesMatrixRows
    if(is.null(Robject[[fxnVecName]]$model$CobjectInterface) ||
       inherits(Robject[[fxnVecName]]$model$CobjectInterface, 'uninitializedField'))
        stop(paste0("populateNodeFxnVecNew could not access compiled model.\n",
                    "Perhaps you did not compile the model used by your nimbleFunction\n",
                    "along with or before this compilation of the nimbleFunction?"))
    numberedPtrs <- Robject[[fxnVecName]]$model$CobjectInterface$.nodeFxnPointers_byDeclID$.ptr
    
    ## This is not really the most efficient way to do things; eventually 
    ## we want to have nodeFunctionVectors contain just the gids, not nodeNames
    ## gids <- Robject[[fxnVecName]]$model$modelDef$nodeName2GraphIDs(nodes)
    eval(call('.Call',
              nimbleUserNamespace$sessionSpecificDll$populateNodeFxnVectorNew_byDeclID,
              fxnVecPtr,
              as.integer(declIDs),
              numberedPtrs,
              as.integer(rowIndices)))
}

populateNodeFxnVecNew_nimDerivs <- function(fxnPtr, Robject, fxnVecName, dll){
    fxnVecPtr <- getNamedObjected(fxnPtr, fxnVecName, dll = dll)
    indexingInfo <- Robject[[fxnVecName]]$indexingInfo
    declIDs <- indexingInfo$declIDs
    rowIndices <- indexingInfo$unrolledIndicesMatrixRows
    if(is.null(Robject[[fxnVecName]]$model$CobjectInterface) ||
       inherits(Robject[[fxnVecName]]$model$CobjectInterface, 'uninitializedField'))
        stop(paste0("populateNodeFxnVecNew could not access compiled model.\n",
                    "Perhaps you did not compile the model used by your nimbleFunction\n",
                    "along with or before this compilation of the nimbleFunction?"))
    numberedPtrs <- Robject[[fxnVecName]]$model$CobjectInterface$.nodeFxnPointers_byDeclID$.ptr
    
    ## This is not really the most efficient way to do things; eventually 
    ## we want to have nodeFunctionVectors contain just the gids, not nodeNames
    ## gids <- Robject[[fxnVecName]]$model$modelDef$nodeName2GraphIDs(nodes)
    derivsInfo <- Robject[[fxnVecName]]$nimDerivsInfo
    eval(call('.Call',
              nimbleUserNamespace$sessionSpecificDll$populateNodeFxnVectorNew_byDeclID_forDerivs,
              fxnVecPtr, 
              as.integer(declIDs), ##as.integer(derivsInfo$declIDs),
              numberedPtrs,
              as.integer(rowIndices), ##as.integer(derivsInfo$rowIndices),
              derivsInfo))
}


populateIndexedNodeInfoTable <- function(fxnPtr, Robject, indexedNodeInfoTableName, dll) {
    iNITptr <- getNamedObjected(fxnPtr, indexedNodeInfoTableName, dll = dll)
    iNITcontent <- Robject[[indexedNodeInfoTableName]]$unrolledIndicesMatrix
    eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$populateIndexedNodeInfoTable, iNITptr, iNITcontent))
}

# NumberedObjects is a reference class which contains a pointer to a C++ object. This C++ object
# stores void pointers. This pointers are indexed by integers and can be accessesed in R via `[` and `[<-`
# However, the intent is that the pointers will actually be accessed more directly in C++ 
# At this time, used to store pointers to nodeFunctions, which will allow for fast
# population of nodeFunctionVectors. They are indexed by graphID's
numberedObjects <- setRefClass('numberedObjects',
                               fields = c('.ptr' = 'ANY', dll = 'ANY'), 
                               methods = list(
                                   initialize = function(dll){
                                       dll <<- dll
                                       .ptr <<- newNumberedObjects(dll)
                                   },
                                   finalize = function() {
                                       nimbleInternalFunctions$nimbleFinalize(.ptr)
                                   },
                                   getSize = function(){
                                       getSize_NumberedObjects(.ptr, dll)
                                   },
                                   resize = function(size){
                                       resize_NumberedObjects(.ptr, size, dll)
                                   }
                               )
                               )

setMethod('[', 'numberedObjects', function(x, i){
    getNumberedObject(x$.ptr, i, x$dll)
})

setMethod('[<-', 'numberedObjects', function(x, i, value){
    assignNumberedObject(x$.ptr, i, value, x$dll)
    return(x)
})



newNumberedObjects <- function(dll){
    ans <- eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$newNumberedObjects))
    eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$register_numberedObjects_Finalizer, ans, dll[['handle']], "numberedObjects"))
    ans
}

getSize_NumberedObjects <- function(numberedObject, dll){
    eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getSizeNumberedObjects, numberedObject))
}

resize_NumberedObjects <- function(numberedObject, size, dll){
    nil <- eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$resizeNumberedObjects, numberedObject, as.integer(size)) )
}

assignNumberedObject <- function(numberedObject, index, val, dll){
    if(!is(val, 'externalptr'))
        stop('Attempting to assign a val which is not an externalptr to a NumberedObjects')
    if(index < 1 || index > getSize_NumberedObjects(numberedObject, dll) )
        stop('Invalid index')
    nil <- eval(call('.Call', getNativeSymbolInfo('setNumberedObject', nimbleUserNamespace$sessionSpecificDll), numberedObject, as.integer(index), val))
}

getNumberedObject <- function(numberedObject, index, dll){
    if(index < 1 || index > getSize_NumberedObjects(numberedObject, dll) )
        stop('Invalid index')
    eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getNumberedObject, numberedObject, as.integer(index)))	
}

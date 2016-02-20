# A list of the R-wrapped functions in this file is provided, with further explanation below.
# makeCSingleVariableAccessor <- function(rModelPtr, elementName, beginIndex, endIndex)
# resizeCModelAccessors <- function(modelAccessPtr, size)
# populateManyModelVarAccess <- function(fxnPtr, Robject, manyAccessName)
# 
# makeCSingleModelValuesAccessor <- function(rModelValuesPtr, elementName, curRow = 1, beginIndex, endIndex)
# getModelAccessorValues <- function(modelAccessor)
# getModelValuesAccessorValues <- function(modelAccessor)
# setNodeElement <- function(nodePtr, modelElementPtr, nodeElementName)
# newNodeFxnVec <- function(size = 0)
# resizeNodeFxnVec <- function(nodeFxnVecPtr, size)
# addNodeFxn <- function(NVecFxnPtr, NFxnPtr, addAtEnd = TRUE, index = -1)
# removeNodeFxn <- function(NVecFxnPtr, index = 1, removeAll = FALSE)
# newManyVarAccess <- function(size = 0)
# addSingleVarAccess <- function(ManyVarAccessPtr, SingleVarAccessPtr, addAtEnd = TRUE, index = -1)
# removeSingleVarAccess <- function(ManyVarAccessPtr, index, removeAll)
# newManyModelValuesAccess <- function(size)
# addSingleModelValuesAccess <- function(ManyModelValuesAccessPtr, SingleModelValuesAccessPtr, addAtEnd, index)
# removeSingleModelValuesAccess <- function(ManyModelValuesAccessPtr, index, removeAll)

# NumberedObjects 


makeCSingleVariableAccessor <- function(rModelPtr, elementName, beginIndex, endIndex){
    .Call("makeSingleVariableAccessor", rModelPtr, elementName, 
          as.integer(beginIndex), as.integer(endIndex) ) 
}    
        
#  This function make a single variable accessor. You must provide the model in which
#  the variable resides by rModelPtr. elementName is the name of the variable you would
#  like to access. beginIndex and endIndex are R indices (i.e. first element is 1, not 0)
#  for the begining and ending sequence of flat indices this accessor uses

populateCopierVector <- function(fxnPtr, Robject, vecName) {
    vecPtr <- .Call("getModelObjectPtr", fxnPtr, vecName)
    copierVectorObject <- Robject[[vecName]]
    fromPtr <- .Call("getModelObjectPtr", fxnPtr, copierVectorObject[[1]])
    toPtr <- .Call("getModelObjectPtr", fxnPtr, copierVectorObject[[2]])
    .Call('populateCopierVector', vecPtr, fromPtr, toPtr, as.integer(copierVectorObject[[3]]), as.integer(copierVectorObject[[4]]))
}

populateManyModelVarMapAccess <- function(fxnPtr, Robject, manyAccessName) { ## new version
    manyAccessPtr = .Call("getModelObjectPtr", fxnPtr, manyAccessName)
    cModel <- Robject[[manyAccessName]][[1]]$CobjectInterface
    ## cModel <- Robject[[manyAccessName]]$sourceObject$CobjectInterface ## NEW ACCESSORS 
    if(is(cModel, 'uninitializedField'))
        stop('Compiled C++ model not available; please include the model in your compilation call (or compile it in advance).', call. = FALSE)

    ## intermediate version
    ## mapInfo <- makeMapInfoFromAccessorVector(Robject[[manyAccessName]]) ## slower
    ## if(length(mapInfo) > 0) {
    ##     .Call('populateValueMapAccessors', manyAccessPtr, mapInfo, cModel$.basePtr)
    ## }

    ## fastest version
    ## doing it above way is safe, but
    ##   doing it the following way induces the crashing.
    mapInfo <- makeMapInfoFromAccessorVectorFaster(Robject[[manyAccessName]])
    if(length(mapInfo[[1]]) > 0) {
        .Call('populateValueMapAccessorsFromNodeNames', manyAccessPtr, mapInfo[[1]], mapInfo[[2]], cModel$.basePtr)
    }
    
    ## oldest version
    ##if(Robject[[manyAccessName]]$getLength() > 0) { ## NEW ACCESSORS 
        ##.Call('populateValueMapAccessors', manyAccessPtr, Robject[[manyAccessName]]$getMapInfo(), cModel$.basePtr) 
    ##}
}

populateManyModelValuesMapAccess <- function(fxnPtr, Robject, manyAccessName){ ## new version. nearly identical to populateManyModelVarMapAccess
    manyAccessPtr = .Call("getModelObjectPtr", fxnPtr, manyAccessName)
    ##cModelValues <- Robject[[manyAccessName]]$sourceObject$CobjectInterface ## NEW ACCESSORS
    cModelValues <- Robject[[manyAccessName]][[1]]$CobjectInterface

    ## oldest
    ##.Call('populateValueMapAccessors', manyAccessPtr, Robject[[manyAccessName]]$getMapInfo(), cModelValues$extptr) ## NEW ACCESSORS

    ## intermediate
    ##    .Call('populateValueMapAccessors', manyAccessPtr, makeMapInfoFromAccessorVector(Robject[[manyAccessName]]), cModelValues$extptr) ## slower
        
    ##fastest
    mapInfo <- makeMapInfoFromAccessorVectorFaster(Robject[[manyAccessName]]) ##faster
    .Call('populateValueMapAccessorsFromNodeNames', manyAccessPtr, mapInfo[[1]], mapInfo[[2]], cModelValues$extptr)
}

## addNodeFxn_LOOP <- function(x, nodes, fxnVecPtr, countInf){
##     countInf$count <- countInf$count + 1
##     addNodeFxn(fxnVecPtr, nodes[[x]]$.basePtr, addAtEnd = FALSE, index = countInf$count)
## }

getFxnVectorPtr <- function(fxnPtr, fxnVecName)
    .Call('getModelObjectPtr', fxnPtr, fxnVecName)

## populateNodeFxnVec_OLD <- function(fxnPtr, Robject, fxnVecName){
##     fxnVecPtr <- .Call('getModelObjectPtr', fxnPtr, fxnVecName)
##     resizeNodeFxnVec(fxnVecPtr, length(Robject[[fxnVecName]]$nodes))	
##     nodePtrsEnv <- Robject[[fxnVecName]]$model$CobjectInterface$.nodeFxnPointersEnv
##     nil <- .Call('populateNodeFxnVector', fxnVecPtr, Robject[[fxnVecName]]$nodes, nodePtrsEnv)
## }

getNamedObjected <- function(objectPtr, fieldName)
    .Call('getModelObjectPtr', objectPtr, fieldName)

inner_populateNodeFxnVec <- function(fxnVecPtr, gids, numberedPtrs)
    nil <- .Call('populateNodeFxnVector_byGID', fxnVecPtr, as.integer(gids), numberedPtrs)

populateNodeFxnVec <- function(fxnPtr, Robject, fxnVecName){
    fxnVecPtr <- getNamedObjected(fxnPtr, fxnVecName)
    gids <- Robject[[fxnVecName]]$gids
    numberedPtrs <- Robject[[fxnVecName]]$model$CobjectInterface$.nodeFxnPointers_byGID$.ptr
    
    ## This is not really the most efficient way to do things; eventually 
    ## we want to have nodeFunctionVectors contain just the gids, not nodeNames
    ## gids <- Robject[[fxnVecName]]$model$modelDef$nodeName2GraphIDs(nodes)
	
    inner_populateNodeFxnVec(fxnVecPtr, gids, numberedPtrs)
}

populateNodeFxnVecNew <- function(fxnPtr, Robject, fxnVecName){
    fxnVecPtr <- getNamedObjected(fxnPtr, fxnVecName)
    indexingInfo <- Robject[[fxnVecName]]$indexingInfo
    declIDs <- indexingInfo$declIDs
    rowIndices <- indexingInfo$rowIndices
    numberedPtrs <- Robject[[fxnVecName]]$model$CobjectInterface$.nodeFxnPointers_byDeclID$.ptr
    
    ## This is not really the most efficient way to do things; eventually 
    ## we want to have nodeFunctionVectors contain just the gids, not nodeNames
    ## gids <- Robject[[fxnVecName]]$model$modelDef$nodeName2GraphIDs(nodes)
	
    .Call('populateNodeFxnVector_byDeclID', fxnVecPtr, as.integer(declIDs), numberedPtrs, rowIndices)
}


# Currently requires: addSingleModelValuesAccess

makeCSingleModelValuesAccessor <- function(rModelValuesPtr, elementName, curRow = 1, beginIndex, endIndex)
    .Call("makeSingleModelValuesAccessor", rModelValuesPtr, elementName, 
          as.integer(curRow), as.integer(beginIndex), as.integer(endIndex) ) 
#   Same as above but for modelValues instead of variables of a model. Note that the pointer
#   we are passing now must point to a C modelValues object, not a model. To get a modelValues pointer
#   from a model, we would have to first call
#   MVPtr <- getMVPtr(rModelPtr)    (getMVPtr is in BuildInterfaces.R)
#   Then we can pass MVPtr to this function
#   Also, this function requires you to select which row you would like to point to via curRow 
#   (R index, not C++ index)

getModelAccessorValues <- function(modelAccessor)
    .Call("getModelAccessorValues", modelAccessor)
#   This retrieves the values from a modelAccessor. It is very important to note that this is 
#   for a singleVariableAccessor, NOT a singleModelValuesAccessor

getModelValuesAccessorValues <- function(modelAccessor)
    .Call("getMVAccessorValues", modelAccessor)
#   Same as above, but for singleModelValuesAccessors


## newNodeFxnVec <- function(size = 0) 
##     .Call("newNodeFxnVector", as.integer(size)  ) 
    
#   This creates a new NodeFunctionVector. We can declare the size of the vector of nodeFunctions upon building.
#   The default size is 0, which leads to the simplest way to populate the nodeFunctionVector: by just adding on to the
#   end one at a time (see addNodeFun for more details). However, the simpliest way is also slow: this method
#   is of complexity O(n^2). This maybe be a problem if we have to deal with dynamic dependencies, for example. 
#   The more efficient to populate this is to set the size in advance and populate by index. This is of O(n) complexity,
#   but means we must keep track of the index as we populate and know number of indices necessary in advance

## resizeNodeFxnVec <- function(nodeFxnVecPtr, size)
##     nil <- .Call("resizeNodeFxnVector", nodeFxnVecPtr, as.integer(size) ) 
## #	Resizes a nodeFunctionPointer object

## addNodeFxn <- function(NVecFxnPtr, NFxnPtr, addAtEnd = TRUE, index = -1)
##     nil <- .Call("addNodeFun", NVecFxnPtr, NFxnPtr, as.logical(addAtEnd), as.integer(index) ) 
## #   This function adds a nodeFunction pointer (NFxnPtr) to a vector of nodeFunctions (NVecFxnPtr)
## #   This can either add one on to the end of the vector (by setting addAtEnd = TRUE) or 
## #   can insert the a nodeFunction by index (by setting addAtEnd = FALSE and index = R-index of position)

## removeNodeFxn <- function(NVecFxnPtr, index = 1, removeAll = FALSE)
##     nil <- .Call("removeNodeFun", NVecFxnPtr, as.integer(index), as.logical(removeAll) ) 
## #   This function removes either the nodeFunctionPointer at position index (R-index) or removes all if removeAll = TRUE

newManyVarAccess <- function(size = 0)
    .Call("newManyVariableAccessor", as.integer(size) ) 
#   Same as newNodeFxnVec, but a builds a ManyVariableAccessor rather than a nodeVectorFunction

addSingleVarAccess <- function(ManyVarAccessPtr, SingleVarAccessPtr, addAtEnd = TRUE, index = -1)
    nil <- .Call("addSingleVariableAccessor", ManyVarAccessPtr, SingleVarAccessPtr, as.logical(addAtEnd), as.integer(index) )
#   Same as addNodeFxn, but for adding a SingleModelVariableAccessor to a ManyModelVariablesAccessors

  
removeSingleVarAccess <- function(ManyVarAccessPtr, index = 1, removeAll = FALSE)
    nil <- .Call("removeModelVariableAccessor", ManyVarAccessPtr, as.integer(index), as.logical(removeAll) )
#   Same as removeNodeFxn, but for SingleModelVariableAccessors

newManyModelValuesAccess <- function(size)
    .Call("newManyModelValuesAccessor", as.integer(size) ) 
#   Same as newNodeFxnVec, but for ManyModelValuesAccessor

addSingleModelValuesAccess <- function(ManyModelValuesAccessPtr, SingleModelValuesAccessPtr, addAtEnd, index = - 1)
    nil <- .Call("addSingleModelValuesAccessor", ManyModelValuesAccessPtr, SingleModelValuesAccessPtr, as.logical(addAtEnd), as.integer(index) ) 
#   Same as addNodeFxn, but for adding singleModelValuesAccess to ManyModelValuesAccessor

removeSingleModelValuesAccess <- function(ManyModelValuesAccessPtr, index, removeAll = FALSE)
    nil <- .Call("removeModelValuesAccessor", ManyModelValuesAccessPtr, as.integer(index), as.logical(removeAll) ) 
#   Same as removeNodeFxn, but for ManyModelValuessAccessor




# NumberedObjects is a reference class which contains a pointer to a C++ object. This C++ object
# stores void pointers. This pointers are indexed by integers and can be accessesed in R via `[` and `[<-`
# However, the intent is that the pointers will actually be accessed more directly in C++ 
# At this time, used to store pointers to nodeFunctions, which will allow for fast
# population of nodeFunctionVectors. They are indexed by graphID's
numberedObjects <- setRefClass('numberedObjects', fields = c('.ptr' = 'ANY'), 
                               methods = list(
                                   initialize = function(){
                                       .ptr <<- newNumberedObjects()
                                   },
                                   getSize = function(){
                                       getSize_NumberedObjects(.ptr)
                                   },
                                   resize = function(size){
                                       resize_NumberedObjects(.ptr, size)
                                   }
                               )
                               )

setMethod('[', 'numberedObjects', function(x, i){
    getNumberedObject(x$.ptr, i)
})

setMethod('[<-', 'numberedObjects', function(x, i, value){
    assignNumberedObject(x$.ptr, i, value)
    return(x)
})



newNumberedObjects <- function(){
    .Call('newNumberedObjects')
}

getSize_NumberedObjects <- function(numberedObject){
    .Call('getSizeNumberedObjects', numberedObject)
}

resize_NumberedObjects <- function(numberedObject, size){
    nil <- .Call('resizeNumberedObjects', numberedObject, as.integer(size) )
}

assignNumberedObject <- function(numberedObject, index, val){
    if(!is(val, 'externalptr'))
        stop('Attempting to assign a val which is not an externalptr to a NumberedObjects')
    if(index < 1 || index > getSize_NumberedObjects(numberedObject) )
        stop('Invalid index')
    nil <- .Call('setNumberedObject', numberedObject, as.integer(index), val)
}

getNumberedObject <- function(numberedObject, index){
    if(index < 1 || index > getSize_NumberedObjects(numberedObject) )
        stop('Invalid index')
    .Call('getNumberedObject', numberedObject, as.integer(index))	
}


numberedModelValuesAccessors <- setRefClass('numberedModelValuesAccessors',
                                            fields = c('.ptr' = 'ANY'),
                                            methods = list(
                                                initialize = function(){ 
                                                    .ptr  <<- .Call('new_SingleModelValuesAccessor_NumberedObjects')
                                                },
                                                getSize = function(){
                                                    getSize_NumberedObjects(.ptr)
                                                },
                                                resize = function(size){
                                                    resize_NumberedObjects(.ptr, size)
                                                }))

setMethod('[', 'numberedModelValuesAccessors', function(x, i){
    getNumberedObject(x$.ptr, i)
})

setMethod('[<-', 'numberedModelValuesAccessors', function(x, i, value){
    assignNumberedObject(x$.ptr, i, value)
    return(x)
})
		
		
		
numberedModelVariableAccessors <- setRefClass('numberedModelVariableAccessors',
                                              fields = c('.ptr' = 'ANY'),
                                              methods = list(
                                                  initialize = function(){ 
                                                      .ptr  <<- .Call('new_SingleModelVariablesAccessor_NumberedObjects')
                                                  },
                                                  getSize = function(){
                                                      getSize_NumberedObjects(.ptr)
                                                  },
                                                  resize = function(size){
                                                      resize_NumberedObjects(.ptr, size)
                                                  }))

setMethod('[', 'numberedModelValuesAccessors', function(x, i){
    getNumberedObject(x$.ptr, i)
})

setMethod('[<-', 'numberedModelValuesAccessors', function(x, i, value){
    assignNumberedObject(x$.ptr, i, value)
    return(x)
})

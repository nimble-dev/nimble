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





makeCSingleVariableAccessor <- function(rModelPtr, elementName, beginIndex, endIndex){
   .Call("makeSingleVariableAccessor", rModelPtr, elementName, 
        as.integer(beginIndex), as.integer(endIndex) ) 
    }    
        
#  This function make a single variable accessor. You must provide the model in which
#  the variable resides by rModelPtr. elementName is the name of the variable you would
#  like to access. beginIndex and endIndex are R indices (i.e. first element is 1, not 0)
#  for the begining and ending sequence of flat indices this accessor uses

resizeCModelAccessors <- function(modelAccessPtr, size)
	.Call('resizeManyModelVarAccessor', modelAccessPtr, as.integer(size) ) 
#	Resizes a manyModelVariablesAccessor


populateManyModelVarAccess <- function(fxnPtr, Robject, manyAccessName)
	{
	manyAccessPtr = .Call("getModelObjectPtr", fxnPtr, manyAccessName)
	resizeCModelAccessors(manyAccessPtr, length(Robject[[manyAccessName]]$modelVariableAccessors) )
	rModelPtr <- Robject[[manyAccessName]]$model$CobjectInterface$.basePtr
	for(i in seq_along(Robject[[manyAccessName]]$modelVariableAccessors) ){
		acc = Robject[[manyAccessName]]$modelVariableAccessors[[i]]
		singleAccPtr = makeCSingleVariableAccessor(rModelPtr = rModelPtr, elementName = acc$var, beginIndex = acc$first, endIndex = acc$last)
		addSingleVarAccess(manyAccessPtr, singleAccPtr, addAtEnd = FALSE, index = i)
	}
}
# Populates a manyModelVariablesAccessor, as called from copyFromRobject() in the CnimbleFunctionBase

resizeCModelValuesAccessors <- function(modelValuesAccessPtr, size)
	.Call('resizeManyModelValuesAccessor', modelValuesAccessPtr, as.integer(size) ) 
#	Resizes a manyModelvaluesAccessor

populateManyModelValuesAccess <- function(fxnPtr, Robject, manyAccessName){
	manyAccessPtr = .Call("getModelObjectPtr", fxnPtr, manyAccessName)
	resizeCModelValuesAccessors(manyAccessPtr, length(Robject[[manyAccessName]]$modelValuesAccessors) )
	rModelValuesPtr <- Robject[[manyAccessName]]$modelValues$CobjectInterface$extptr
	for(i in seq_along(Robject[[manyAccessName]]$modelValuesAccessors) ){
		acc = Robject[[manyAccessName]]$modelValuesAccessors[[i]]
		singleAccPtr = makeCSingleModelValuesAccessor(rModelValuesPtr = rModelValuesPtr, elementName = acc$var, beginIndex = acc$first, endIndex = acc$last)
		addSingleModelValuesAccess(manyAccessPtr, singleAccPtr, addAtEnd = FALSE, index = i)
	}
}

populateNodeFxnVec <- function(fxnPtr, Robject, fxnVecName){
	fxnVecPtr = .Call("getModelObjectPtr", fxnPtr, fxnVecName)
	resizeNodeFxnVec(fxnVecPtr, length(Robject[[fxnVecName]]$nodes) )
	cModel <- Robject[[fxnVecName]]$model$CobjectInterface	#$.basePtr
##	allNodePtrs = getNodeFxnPtrs(cModel)
	count = 0	
	for(node in Robject[[fxnVecName]]$nodes){
		count = count + 1
		nodeP = cModel$nodes[[node]]$.basePtr ##allNodePtrs[[node]]
		addNodeFxn(fxnVecPtr, nodeP, addAtEnd = FALSE, index = count)
		}
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

# setNodeElement <- function(nodePtr, modelElementPtr, nodeElementName)
# 	.Call("setNodeModelPtr", nodePtr, modelElementPtr, as.character(nodeElementName) ) 
# #	This function initializes a pointer to a model element for a node function. The nodeElementName must
# #	match the namedObjects name for the element in the nodeFunction.

	
newNodeFxnVec <- function(size = 0) 
    .Call("newNodeFxnVector", as.integer(size)  ) 
    
#   This creates a new NodeFunctionVector. We can declare the size of the vector of nodeFunctions upon building.
#   The default size is 0, which leads to the simplest way to populate the nodeFunctionVector: by just adding on to the
#   end one at a time (see addNodeFun for more details). However, the simpliest way is also slow: this method
#   is of complexity O(n^2). This maybe be a problem if we have to deal with dynamic dependencies, for example. 
#   The more efficient to populate this is to set the size in advance and populate by index. This is of O(n) complexity,
#   but means we must keep track of the index as we populate and know number of indices necessary in advance

resizeNodeFxnVec <- function(nodeFxnVecPtr, size)
	nil <- .Call("resizeNodeFxnVector", nodeFxnVecPtr, as.integer(size) ) 
#	Resizes a nodeFunctionPointer object

addNodeFxn <- function(NVecFxnPtr, NFxnPtr, addAtEnd = TRUE, index = -1)
  nil <- .Call("addNodeFun", NVecFxnPtr, NFxnPtr, as.logical(addAtEnd), as.integer(index) ) 
#   This function adds a nodeFunction pointer (NFxnPtr) to a vector of nodeFunctions (NVecFxnPtr)
#   This can either add one on to the end of the vector (by setting addAtEnd = TRUE) or 
#   can insert the a nodeFunction by index (by setting addAtEnd = FALSE and index = R-index of position)

removeNodeFxn <- function(NVecFxnPtr, index = 1, removeAll = FALSE)
  nil <- .Call("removeNodeFun", NVecFxnPtr, as.integer(index), as.logical(removeAll) ) 
#   This function removes either the nodeFunctionPointer at position index (R-index) or removes all if removeAll = TRUE

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

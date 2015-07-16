cGetNRow <- function(cMV, compIndex = 1)
{
  nRow = .Call(getNativeSymbolInfo("getNRow"), cMV$componentExtptrs[[compIndex]])
  return(nRow)
}

cAddBlank <- function(cMV, addNum)
{
  addNum = as.integer(addNum)
  for(i in 1:length( cMV$componentExtptrs) )
    status = .Call(getNativeSymbolInfo('addBlankModelValueRows'), cMV$componentExtptrs[[i]], addNum)
  if(status == FALSE)
    stop("Error adding rows to ModelValues")
}

cCopyVariableRows <- function(cMVFrom, cMVTo, varIndex, rowsFrom = 1:cGetNRow(cMVFrom), rowsTo = 1:cGetNRow(cMVFrom) ) 
{
  if(length(varIndex) > 1)
    stop("cCopyVariableRows only takes on varIndex at a time")
  rowsFrom = as.integer(rowsFrom )
  rowsTo = as.integer(rowsTo )
  if(cMVFrom$extptrCall != cMVTo$extptrCall)
    stop("ModelValues are not of the same type")
  fromPtr <- cMVFrom$componentExtptrs[[varIndex]]
  toPtr <- cMVTo$componentExtptrs[[varIndex]]
  status = .Call(getNativeSymbolInfo('copyModelValuesElements'), fromPtr, toPtr, rowsFrom, rowsTo)
  if(status == FALSE)
    stop("Did not correctly copy from one ModelValues to another")
}

newObjElementPtr = function(rPtr, name){
  .Call("getModelObjectPtr", rPtr, name)
} 

getNimValues <- function(elementPtr, pointDepth = 1){
  if(!inherits(elementPtr, "externalptr"))    
    return(NULL)
  .Call("Nim_2_SEXP", elementPtr, as.integer(pointDepth) ) 
}

setNimValues <- function(elementPtr, values, pointDepth = 1, allowResize = TRUE){
  ptrExp <- substitute(elementPtr)
  storage.mode(values) <- 'numeric'
  if(!inherits(elementPtr, "externalptr"))    
    return(NULL)
      jnk = .Call("SEXP_2_Nim", elementPtr, as.integer(pointDepth), values, allowResize)
  values
}

setPtrVectorOfPtrs <- function(accessorPtr, contentsPtr, length) {
    if(!inherits(accessorPtr, 'externalptr')) return(NULL)
    if(!inherits(contentsPtr, 'externalptr')) return(NULL)
    if(!is.numeric(length)) return(NULL)
    .Call('setPtrVectorOfPtrs', accessorPtr, contentsPtr, as.integer(length))
    contentsPtr
}

setOnePtrVectorOfPtrs <- function(accessorPtr, i, contentsPtr) {
    if(!inherits(accessorPtr, 'externalptr')) return(NULL)
    if(!is.numeric(i)) return(NULL)
    if(!inherits(contentsPtr, 'externalptr')) return(NULL)
    .Call('setOnePtrVectorOfPtrs', accessorPtr, as.integer(i-1), contentsPtr)
    contentsPtr
}

setDoublePtrFromSinglePtr <- function(elementPtr, value) {
    if(!inherits(elementPtr, 'externalptr')) return(NULL)
    if(!inherits(value, 'externalptr')) return(NULL)
    .Call('setDoublePtrFromSinglePtr', elementPtr, value)
    value
}

getDoubleValue <- function(elementPtr, pointDepth = 1){
  if(!inherits(elementPtr, "externalptr") ) 
    return(NULL)
  .Call("double_2_SEXP", elementPtr, as.integer(pointDepth) ) 
}

setDoubleValue <- function(elementPtr, value,  pointDepth = 1){
  if(!inherits(elementPtr, "externalptr"))    
    return(NULL)
  jnk = .Call("SEXP_2_double", elementPtr, as.integer(pointDepth), value)
  value
}

getIntValue <- function(elementPtr, pointDepth = 1){
  if(!inherits(elementPtr, "externalptr") ) 
    return(NULL)
  .Call("int_2_SEXP", elementPtr, as.integer(pointDepth) ) 
}

setIntValue <- function(elementPtr, value,  pointDepth = 1){
  if(!inherits(elementPtr, "externalptr"))    
    return(NULL)
  jnk = .Call("SEXP_2_int", elementPtr, as.integer(pointDepth), value )
}

getBoolValue <- function(elementPtr, pointDepth = 1){
    if(!inherits(elementPtr, "externalptr") ) 
        return(NULL)
    .Call("bool_2_SEXP", elementPtr, as.integer(pointDepth) ) 
}

setBoolValue <- function(elementPtr, value,  pointDepth = 1){
    if(!inherits(elementPtr, "externalptr"))    
        return(NULL)
    jnk = .Call("SEXP_2_bool", elementPtr, as.integer(pointDepth), value )
}

getBoolValue <- function(elementPtr, pointDepth = 1){
    if(!inherits(elementPtr, "externalptr") ) 
        return(NULL)
    .Call("bool_2_SEXP", elementPtr, as.integer(pointDepth) ) 
}

setCharacterValue <- function(elementPtr, value){
    if(!inherits(elementPtr, "externalptr"))    
        return(NULL)
    jnk = .Call("SEXP_2_string", elementPtr, value )
}

getCharacterValue <- function(elementPtr){
    if(!inherits(elementPtr, "externalptr") ) 
        return(NULL)
    .Call("string_2_SEXP", elementPtr ) 
}

setCharacterVectorValue <- function(elementPtr, value){
    if(!inherits(elementPtr, "externalptr"))    
        return(NULL)
    jnk = .Call("SEXP_2_stringVector", elementPtr, value )
}

getCharacterVectorValue <- function(elementPtr){
    if(!inherits(elementPtr, "externalptr") ) 
        return(NULL)
    .Call("stringVector_2_SEXP", elementPtr) 
}

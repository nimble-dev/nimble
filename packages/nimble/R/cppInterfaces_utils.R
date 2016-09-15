cGetNRow <- function(cMV, compIndex = 1)
{
  nRow = .Call(getNativeSymbolInfo("getNRow", cMV$dll), cMV$componentExtptrs[[compIndex]])
  return(nRow)
}

## cAddBlank <- function(cMV, addNum)
## {
##   addNum = as.integer(addNum)
##   for(i in 1:length( cMV$componentExtptrs) )
##     status = .Call(getNativeSymbolInfo('addBlankModelValueRows', cMV$dll), cMV$componentExtptrs[[i]], addNum)
##   if(status == FALSE)
##     stop("Error adding rows to ModelValues")
## }

## cCopyVariableRows <- function(cMVFrom, cMVTo, varIndex, rowsFrom = 1:cGetNRow(cMVFrom), rowsTo = 1:cGetNRow(cMVFrom), dll ) 
## {
##   if(length(varIndex) > 1)
##     stop("cCopyVariableRows only takes on varIndex at a time")
##   rowsFrom = as.integer(rowsFrom )
##   rowsTo = as.integer(rowsTo )
##   if(cMVFrom$extptrCall != cMVTo$extptrCall)
##     stop("ModelValues are not of the same type")
##   fromPtr <- cMVFrom$componentExtptrs[[varIndex]]
##   toPtr <- cMVTo$componentExtptrs[[varIndex]]
##   status = .Call(getNativeSymbolInfo('copyModelValuesElements', dll), fromPtr, toPtr, rowsFrom, rowsTo)
##   if(status == FALSE)
##     stop("Did not correctly copy from one ModelValues to another")
## }

newObjElementPtr = function(rPtr, name, dll){
  .Call(getNativeSymbolInfo("getModelObjectPtr", dll), rPtr, name)
} 

getNimValues <- function(elementPtr, pointDepth = 1, dll){
  if(!inherits(elementPtr, "externalptr"))    
    return(NULL)
  .Call(getNativeSymbolInfo("Nim_2_SEXP", dll), elementPtr, as.integer(pointDepth) ) 
}

setNimValues <- function(elementPtr, values, pointDepth = 1, allowResize = TRUE, dll){
  ptrExp <- substitute(elementPtr)
  storage.mode(values) <- 'numeric'
  if(!inherits(elementPtr, "externalptr"))    
    return(NULL)
      jnk = .Call(getNativeSymbolInfo("SEXP_2_Nim", dll), elementPtr, as.integer(pointDepth), values, allowResize)
  values
}

setPtrVectorOfPtrs <- function(accessorPtr, contentsPtr, length, dll) {
    if(!inherits(accessorPtr, 'externalptr')) return(NULL)
    if(!inherits(contentsPtr, 'externalptr')) return(NULL)
    if(!is.numeric(length)) return(NULL)
    .Call(getNativeSymbolInfo('setPtrVectorOfPtrs', dll), accessorPtr, contentsPtr, as.integer(length))
    contentsPtr
}

setOnePtrVectorOfPtrs <- function(accessorPtr, i, contentsPtr, dll) {
    if(!inherits(accessorPtr, 'externalptr')) return(NULL)
    if(!is.numeric(i)) return(NULL)
    if(!inherits(contentsPtr, 'externalptr')) return(NULL)
    .Call(getNativeSymbolInfo('setOnePtrVectorOfPtrs', dll), accessorPtr, as.integer(i-1), contentsPtr)
    contentsPtr
}

setDoublePtrFromSinglePtr <- function(elementPtr, value, dll) {
    if(!inherits(elementPtr, 'externalptr')) return(NULL)
    if(!inherits(value, 'externalptr')) return(NULL)
    .Call(getNativeSymbolInfo('setDoublePtrFromSinglePtr', dll), elementPtr, value)
    value
}

getDoubleValue <- function(elementPtr, pointDepth = 1, dll){
  if(!inherits(elementPtr, "externalptr") ) 
    return(NULL)
  .Call(getNativeSymbolInfo("double_2_SEXP", dll), elementPtr, as.integer(pointDepth) ) 
}

setDoubleValue <- function(elementPtr, value,  pointDepth = 1, dll){
  if(!inherits(elementPtr, "externalptr"))    
    return(NULL)
  jnk = .Call(getNativeSymbolInfo("SEXP_2_double", dll), elementPtr, as.integer(pointDepth), value)
  value
}

getIntValue <- function(elementPtr, pointDepth = 1, dll){
  if(!inherits(elementPtr, "externalptr") ) 
    return(NULL)
  .Call(getNativeSymbolInfo("int_2_SEXP", dll), elementPtr, as.integer(pointDepth) ) 
}

setIntValue <- function(elementPtr, value,  pointDepth = 1, dll){
  if(!inherits(elementPtr, "externalptr"))    
    return(NULL)
  jnk = .Call(getNativeSymbolInfo("SEXP_2_int", dll), elementPtr, as.integer(pointDepth), value )
}

getBoolValue <- function(elementPtr, pointDepth = 1, dll){
    if(!inherits(elementPtr, "externalptr") ) 
        return(NULL)
    .Call(getNativeSymbolInfo("bool_2_SEXP", dll), elementPtr, as.integer(pointDepth) ) 
}

setBoolValue <- function(elementPtr, value,  pointDepth = 1, dll){
    if(!inherits(elementPtr, "externalptr"))    
        return(NULL)
    jnk = .Call(getNativeSymbolInfo("SEXP_2_bool", dll), elementPtr, as.integer(pointDepth), value )
}

setCharacterValue <- function(elementPtr, value, dll){
    if(!inherits(elementPtr, "externalptr"))    
        return(NULL)
    jnk = .Call(getNativeSymbolInfo("SEXP_2_string", dll), elementPtr, value )
}

getCharacterValue <- function(elementPtr, dll){
    if(!inherits(elementPtr, "externalptr") ) 
        return(NULL)
    .Call(getNativeSymbolInfo("string_2_SEXP", dll), elementPtr ) 
}

setCharacterVectorValue <- function(elementPtr, value, dll){
    if(!inherits(elementPtr, "externalptr"))    
        return(NULL)
    jnk = .Call(getNativeSymbolInfo("SEXP_2_stringVector", dll), elementPtr, value )
}

getCharacterVectorValue <- function(elementPtr, dll){
    if(!inherits(elementPtr, "externalptr") ) 
        return(NULL)
    .Call(getNativeSymbolInfo("stringVector_2_SEXP", dll), elementPtr) 
}

nimbleFinalize <- function(extptr) {
    eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$RNimble_Ptr_ManualFinalizer, extptr))
}

cGetNRow <- function(cMV, compIndex = 1)
{
  nRow = eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getNRow, cMV$componentExtptrs[[compIndex]]))
  return(nRow)
}

newObjElementPtr <- function(rPtr, name, dll){
  eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getModelObjectPtr, rPtr, name))
}

getNimValues <- function(elementPtr, pointDepth = 1, dll){
  if(!inherits(elementPtr, "externalptr"))
    return(NULL)
  eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$Nim_2_SEXP, elementPtr, as.integer(pointDepth)) )
}

setNimValues <- function(elementPtr, values, pointDepth = 1, allowResize = TRUE, dll){
  ptrExp <- substitute(elementPtr)
  storage.mode(values) <- 'numeric'
  if(!inherits(elementPtr, "externalptr"))
    return(NULL)
      jnk = eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$SEXP_2_Nim, elementPtr, as.integer(pointDepth), values, allowResize))
  values
}

setPtrVectorOfPtrs <- function(accessorPtr, contentsPtr, length, dll) {
    if(!inherits(accessorPtr, 'externalptr')) return(NULL)
    if(!inherits(contentsPtr, 'externalptr')) return(NULL)
    if(!is.numeric(length)) return(NULL)
    eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$setPtrVectorOfPtrs, accessorPtr, contentsPtr, as.integer(length)))
    contentsPtr
}

setOnePtrVectorOfPtrs <- function(accessorPtr, i, contentsPtr, dll) {
    if(!inherits(accessorPtr, 'externalptr')) return(NULL)
    if(!is.numeric(i)) return(NULL)
    if(!inherits(contentsPtr, 'externalptr')) return(NULL)
    eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$setOnePtrVectorOfPtrs, accessorPtr, as.integer(i-1), contentsPtr))
    contentsPtr
}

setDoublePtrFromSinglePtr <- function(elementPtr, value, dll) {
    if(!inherits(elementPtr, 'externalptr')) return(NULL)
    if(!inherits(value, 'externalptr')) return(NULL)
    eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$setDoublePtrFromSinglePtr, elementPtr, value))
    value
}

setSmartPtrFromSinglePtr <- function(elementPtr, value, dll) {
  if(!inherits(elementPtr, 'externalptr')) return(NULL)
  if(!inherits(value, 'externalptr')) return(NULL)
  eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$setSmartPtrFromSinglePtr, elementPtr, value))
  value
}

setSmartPtrFromDoublePtr <- function(elementPtr, value, dll) {
  if(!inherits(elementPtr, 'externalptr')) return(NULL)
  if(!inherits(value, 'externalptr')) return(NULL)
  eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$setSmartPtrFromDoublePtr, elementPtr, value))
  value
}

getDoubleValue <- function(elementPtr, pointDepth = 1, dll){
  if(!inherits(elementPtr, "externalptr") )
      return(NULL)
  eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$extract_double_2_SEXP, elementPtr, as.integer(pointDepth)))
}

setDoubleValue <- function(elementPtr, value,  pointDepth = 1, dll){
  if(!inherits(elementPtr, "externalptr"))
      return(NULL)
  eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$populate_SEXP_2_double, elementPtr, as.integer(pointDepth), value))
  value
}

getIntValue <- function(elementPtr, pointDepth = 1, dll){
  if(!inherits(elementPtr, "externalptr") )
      return(NULL)
  eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$extract_int_2_SEXP, elementPtr, as.integer(pointDepth)))
}

setIntValue <- function(elementPtr, value,  pointDepth = 1, dll){
  if(!inherits(elementPtr, "externalptr"))
      return(NULL)
  eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$populate_SEXP_2_int, elementPtr, as.integer(pointDepth), value))
  value
}

getBoolValue <- function(elementPtr, pointDepth = 1, dll){
    if(!inherits(elementPtr, "externalptr") )
        return(NULL)
     eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$extract_bool_2_SEXP, elementPtr, as.integer(pointDepth)))
}

setBoolValue <- function(elementPtr, value,  pointDepth = 1, dll){
    if(!inherits(elementPtr, "externalptr"))
        return(NULL)
    eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$populate_SEXP_2_bool, elementPtr, as.integer(pointDepth), value))
    value
}

setCharacterValue <- function(elementPtr, value, dll){
    if(!inherits(elementPtr, "externalptr")) return(NULL)
    eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$populate_SEXP_2_string, elementPtr, value))
    value
}

getCharacterValue <- function(elementPtr, dll){
    if(!inherits(elementPtr, "externalptr") )
        return(NULL)
    eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$extract_string_2_SEXP, elementPtr))
}

setCharacterVectorValue <- function(elementPtr, value, dll){
    if(!inherits(elementPtr, "externalptr"))
        return(NULL)
    eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$populate_SEXP_2_stringVector, elementPtr, value))
    value
}

getCharacterVectorValue <- function(elementPtr, dll){
    if(!inherits(elementPtr, "externalptr") )
        return(NULL)
    eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$extract_stringVector_2_SEXP, elementPtr))
}

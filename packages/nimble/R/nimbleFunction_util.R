existsFunctionEnvVar <- function(f, var) {
  return(exists(x = var, envir = environment(f), inherits = FALSE))
}

getFunctionEnvVar <- function(f, var) {
  return(environment(f)[[var]])
}

nfGetDefVar <- function(f, var) {
    return(environment(nf_getGeneratorFunction(f))[[var]])
}

is.nf <- function(f) {
    if(inherits(f, 'nimbleFunctionBase')) return(TRUE)
    #	$runRelated
    return(is.function(f) && 
               existsFunctionEnvVar(f, 'nfRefClassObject') ) 	
}

is.nfGenerator <- function(f) {
    return(is.function(f) && 
               existsFunctionEnvVar(f, 'generatorFunction') &&
               existsFunctionEnvVar(f, 'nfRefClassDef') &&
               existsFunctionEnvVar(f, 'nfRefClass'))
}

nf_getRefClassObject <- function(f) {
    if(is.nfGenerator(f))     stop('trying to access RefClassObject from nimbleFunction generator.\nError: need to use the specialized nimbleFunction')
    
    #	$runRelated
    if(inherits(f, 'nimbleFunctionBase')) 	return(f)
    if(inherits(f, 'CnimbleFunctionBase'))	return(f)
    
    if(!is.nf(f))             stop('invalid nimbleFunction argument\n')
    return(getFunctionEnvVar(f, 'nfRefClassObject'))
}

nf_getGeneratorFunction <- function(f) {
    if(is.nfGenerator(f))    return(f)
    if(is.nf(f))             return(nf_getRefClassObject(f)$.generatorFunction)
    
    #if(is.refObject(f))		return(f$.generatorFunction)
    
    stop('invalid nimbleFunction argument\n')
}

nf_getInstances <- function(f) {
    if(is.nfGenerator(f))    return(getFunctionEnvVar(f, 'instances'))
    if(is.nf(f))             return(getFunctionEnvVar(nf_getGeneratorFunction(f), 'instances'))
    stop('invalid nimbleFunction argument\n')
}

nf_getMethodList <- function(f) {
    if(is.nfGenerator(f))    return(getFunctionEnvVar(f, 'methodList'))
    if(is.nf(f))             return(getFunctionEnvVar(nf_getGeneratorFunction(f), 'methodList'))
    stop('invalid nimbleFunction argument\n')
}

nf_getSetupOutputNames <- function(f, hidden = FALSE) {
    nameFunction <- if(hidden)     function(x) x     else     nf_namesNotHidden
    if(is.nfGenerator(f))    return(nameFunction(names(getFunctionEnvVar(f, 'nfRefClass')$fields())))
    if(is.nf(f))             return(nameFunction(names(getFunctionEnvVar(nf_getGeneratorFunction(f), 'nfRefClass')$fields())))
    stop('invalid nimbleFunction argument\n')
}





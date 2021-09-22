existsFunctionEnvVar <- function(f, var) {
  return(exists(x = var, envir = environment(f), inherits = FALSE))
}

getFunctionEnvVar <- function(f, var) {
  return(environment(f)[[var]])
}

nfGetDefVar <- function(f, var) {
    return(environment(nf_getGeneratorFunction(f))[[var]])
}

#' check if a nimbleFunction
#'
#' Checks an object to determine if it is a nimbleFunction (i.e., a function created by \code{nimbleFunction} using setup code).
#'
#' @param f object to be tested
#'
#' @param inputIsName logical indicating whether the function is provided as the character name of the function or the function object itself
#'
#' @param where Optional argument needed due to R package namespace issues but which should not need to be provided by a user.
#' 
#' @seealso \code{\link{nimbleFunction}} for how to create a nimbleFunction
#' @export
is.nf <- function(f, inputIsName = FALSE, where = -1) {
    if(inputIsName) f <- get(f, pos = where)
    if(inherits(f, 'nimbleFunctionBase')) return(TRUE)
    return(is.function(f) && !is.null(environment(f)) &&  
               existsFunctionEnvVar(f, 'nfRefClassObject') ) 	
}

#' @export
is.Cnf <- function(f, inputIsName = FALSE, where = -1) {
    if(inputIsName) f <- get(f, pos = where)
    if(inherits(f, 'CnimbleFunctionBase')) return(TRUE)
    return(FALSE)
}

is.nfGenerator <- function(f, inputIsName = FALSE, where = -1) {
    if(inputIsName) f <- get(f, pos = where)
    return(is.function(f) && 
               existsFunctionEnvVar(f, 'generatorFunction') &&
               existsFunctionEnvVar(f, 'nfRefClassDef') &&
               existsFunctionEnvVar(f, 'nfRefClass'))
}

nf_getRefClassObject <- function(f) {
    if(is.nfGenerator(f))     stop('trying to access RefClassObject from nimbleFunction generator.\nError: need to use the specialized nimbleFunction')
    if(inherits(f, 'nimbleFunctionBase')) 	return(f)
    if(inherits(f, 'CnimbleFunctionBase'))	return(f)
    
    if(!is.nf(f))             stop('invalid nimbleFunction argument\n')
    return(getFunctionEnvVar(f, 'nfRefClassObject'))
}

nf_getGeneratorFunction <- function(f) {
    if(is.nfGenerator(f))    return(f)
    if(is.nf(f))             return(nf_getRefClassObject(f)$.generatorFunction)
    if(is.Cnf(f))            return(nf_getGeneratorFunction(f$Robject))
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

#'
#' Get nimbleFunction definition
#'
#' Returns a list containing the nimbleFunction definition components (setup function, run function, and other member methods) for the supplied nimbleFunction generator or specialized instance.
#'
#' @param nf A nimbleFunction generator, or a compiled or un-compiled specialized nimbleFunction.
#'
#' @export
#' @author Daniel Turek
getDefinition <- function(nf) {
    nfGen <- nf_getGeneratorFunction(nf)
    nfEnv <- environment(nfGen)
    defList <- c(list(setup=nfEnv$setup, run=nfEnv$run), nfEnv$methods)
    defList
}

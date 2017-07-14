#' Create a nimbleFunction that wraps a call to external compiled code
#'
#' Given C header information, a function that takes scalars or pointers can be called from a compiled nimbleFunction.  If non-scalar return values are needed, an argument can be selected to behave as the return value in nimble.
#'
#' @param prototype Argument type information.  This can be provided as an R function using \code{nimbleFunction} type declarations or as a list of \code{nimbleType} objects.  See examples below.
#' @param returnType (Optional) Return object type information.  This can be provided similarly to \code{prototype} as either a \code{nimbleFunction} type declaration or as a \code{nimbleType} object.  In the latter case, the name will be ignored.  If omitted, the return type will be \code{void}.
#' @param Cfun Name of the external function (character).
#' @param headerfile Name (possibly including file path) of the header file where Cfun is declared.
#' @param compiledfile Name (possibly including path) of the .o file where Cfun has been compiled.
#'
#' @author Perry de Valpine
#' @export
#' @details This is a beta feature with limited functionality.  The only argument types allowed in Cfun are \code{double}, \code{int}, and \code{bool}, corresponding to \code{nimbleFunction} types \code{double}, \code{integer}, and \code{logical}, respectively.
#'
#' If the dimensionality is greater than zero, the arguments in \code{Cfun} should be pointers.  This means it will typically be necessary to pass additional integer arguments telling \code{Cfun} the size(s) of non-scalar arguments.
#'
#' The return argument can only be a scalar or void.  Since non-scalar arguments are passed by pointer, you can use an argument to return results from \code{Cfun}.  If you wish to have a \code{nimbleFunction} that uses one argument as a return object, you can wrap the result of \code{nimbleExternalCall} in another \code{nimbleFunction} that allocates the return object.
#'
#' @return A \code{nimbleFunction} that takes the indicated input arguments, calls \code{Cfun}, and returns the result.
#'
#' @examples
#' ## not run
#' sink('add1.h')
#' cat('
#' extern "C" {
#' void my_internal_function(double *p, double*ans, int n);
#' }
#' ')
#' sink()
#' sink('add1.cpp') 
#' cat('
#' #include <cstdio>
#' #include "add1.h"
#' void my_internal_function(double *p, double *ans, int n) {
#'   printf("In my_internal_function\\n"); /* cat reduces the double slash to single slash */ 
#'   for(int i = 0; i < n; i++) 
#'     ans[i] = p[i] + 1.0;
#' }
#' ')
#' sink()
#' system('g++ add1.cpp -c -o add1.o')
#' Radd1 <- nimbleExternalCall(function(x = double(1), ans = double(1), n = integer()){}, Cfun = #' 'my_internal_function', headerFile = 'add1.h', oFile = 'add1.o')
#' ## If you need to use a function with non-scalar return object in model code, you can wrap it #' in another nimbleFunction like this:
#' model_add1 <- nimbleFunction(
#'     run = function(x = double(1)) {
#'         ans <- numeric(length(x))
#'         Radd1(x, ans, length(x))
#'         return(ans)
#'         returnType(double(1))
#'     })
#' demoCode <- nimbleCode({
#'     for(i in 1:4) {x[i] ~ dnorm(0,1)} ## just to get a vector
#'     y[1:4] <- model_add1(x[1:4])
#' })
#' demoModel <- nimbleModel(demoCode, inits = list(x = rnorm(4)))
#' CdemoModel <- compileNimble(demoModel)
nimbleExternalCall <- function(prototype, returnType = void(), Cfun, headerFile, oFile) {
    ## construct a nimbleFunction to wrap a call to Cfun
    returnType <- substitute(returnType)
    if(is.list(prototype)) {
        args <- nimbleTypeList2argTypeList(prototype)
        fun <- function(){}
        formals(fun) <- args
    } else {
        if(!is.function(prototype)) stop("Invalid prototype argument")
        fun <- prototype
        args <- formals(fun)
    }
    argsSymTab <- nimble:::argTypeList2symbolTable(args)
    argNames <- names(args)
    replacedArgNames <- argNames
    ## Populate a set of lines to convert from NimArr<>s to C pointers...
    convertLines <- list()
    ## ... and back again
    unconvertLines <- list()
    for(i in seq_along(args)) {
        thisSymbol <- argsSymTab$getSymbolObject( argNames[i] )
        thisType <- thisSymbol$type
        thisNdim <- thisSymbol$nDim
        if(thisNdim > 0) {
            ## for argument "x", make an x_internalPOINTER_
            ptrName <- paste0(argNames[i], '_internalPOINTER_')
            replacedArgNames[i] <- ptrName
            ## Make line "x_internalPOINTER_ <- nimbleConvert(x)"
            newConvertLine <- substitute( A <- nimbleConvert(B), list(A = as.name(ptrName), B = as.name(argNames[i])))
            convertLines[[ length(convertLines) + 1]] <- newConvertLine
            ## Make line "nimbleUnconvert(x_internalPOINTER_, x)"
            newUnconvertLine <- substitute( nimbleUnconvert(A, B), list(A = as.name(ptrName), B = as.name(argNames[i])))
            unconvertLines[[ length(unconvertLines) + 1]] <- newUnconvertLine
        }
    }
    if(inherits(returnType, 'nimbleType'))
        returnType <- nimbleType2argType(returnType)[[1]]
    returnSymbol <- nimble:::argType2symbol(returnType)
    ## Make the call to Cfun, without yet a return value
    externalCallLine <- as.call(c(list(as.name(Cfun)), lapply(replacedArgNames, as.name)))
    returnLines <- list()
    if(returnSymbol$type != 'void') {
        ## Insert the return value as a variable "RETURNVALUE"
        externalCallLine <- substitute(RETURNVALUE <- EXTERNALCALL, list(EXTERNALCALL = externalCallLine))
        ## Insert "return(RETURNVALUE)" and "returnType(<correct type>)"
        returnLines <- list(
            quote(return(RETURNVALUE)),
            substitute(returnType(RT), list(TR = returnType))
        )
    }
    ## put all the lines together
    allLines <- c(list(as.name("{")), convertLines, list(externalCallLine), unconvertLines, returnLines)
    body(fun) <- as.call(allLines)
    ans <- quote(RCfunction(fun, check = FALSE))
    ans <- eval(ans)
    ## Stick header information into the nfMethodRC
    environment(ans)$nfMethodRCobject$externalHincludes <- paste0('\"',headerFile,'\"')
    ## If the user provided the ofile with a .c or .o ending, remove and replace with .cpp. That is a kluge because later .cpp is replaced with .o.
    oFile <- gsub('\\.c', '', oFile)
    oFile <- gsub('\\.o', '', oFile)
    oFile <- paste0(oFile, '.cpp') ## This later gets stripped and replaced with .o but it must be .cpp to be stripped (until we someday generalize that handling)
    ## Stick the "oFile", with .cpp extension as a kluge, into the nfMethodRC
    environment(ans)$nfMethodRCobject$externalCPPincludes <- oFile
    return(ans)
}

#' Make an R function callable from compiled nimbleFunctions.
#'
#' Normally compiled nimbleFunctions call other compiled nimbleFunctions.  nimbleRcall enables any R function (with viable argument types and return values) to be called (and evaluated in R) from compiled nimbleFunctions.
#'
#' @param fun
#' @param returnType
#' @param Rfun
#' @param envir
#'
#' @details
#'
#' @return
#'
#' @seealso
#'
#' @export
#' @author Perry de Valpine
nimbleRcall <- function(fun, returnType = void(), Rfun, envir = .GlobalEnv) {
    returnType <- substitute(returnType)
    args <- formals(fun)
    argsSymTab <- nimble:::argTypeList2symbolTable(args)
    argNames <- names(args)
    replacedArgNames <- argNames
    convertLines <- list()
    SEXPsetupLines <- list()
    SEXPsetupLines[[1]] <- substitute(SLANGpiece <- SLANG <- PROTECT(allocVector(LANGSXP, length )), list(length = length(args) + 1))
    SEXPsetupLines[[2]] <- substitute(SET_NEXT_LANG_ARG(SLANGpiece, install(FUNNAME)), list(FUNNAME = as.character(Rfun)))
    for(i in seq_along(args)) {
        thisSymbol <- argsSymTab$getSymbolObject( argNames[i] )
        thisType <- thisSymbol$type
        thisNdim <- thisSymbol$nDim
        SEXPargName <- paste0(argNames[i], '_SEXP_')
        replacedArgNames[i] <- SEXPargName
        newConvertLine <- substitute( A <- NimArr_2_SEXP(B), list(A = as.name(SEXPargName), B = as.name(argNames[i]))) ## MakeSEXPcopyLineFrom almost does this but it is at cppDef level
        convertLines[[ length(convertLines) + 1]] <- newConvertLine
        SEXPsetupLines[[ length(SEXPsetupLines) + 1 ]] <- substitute(nimVerbatim(SET_TAG(SLANGpiece, install(ARGNAME))), list(ARGNAME = as.character(argNames[i])))
        SEXPsetupLines[[ length(SEXPsetupLines) + 1 ]] <- substitute(nimVerbatim(SET_NEXT_LANG_ARG(SLANGpiece, SEXPARGNAME)), list(SEXPARGNAME = as.name(SEXPargName)))
    }
    returnSymbol <- nimble:::argType2symbol(returnType)
    externalCallLine <- quote(SANS <- Reval(SLANG, R_GlobalEnv))
    if(returnSymbol$type != 'void') {
        resultLines <- list(substitute(declare(ans, RT), list(RT = returnType)),
                            quote(nimVerbatim(SEXP_2_NimArr(SANS, ans))),
                            quote(nimVerbatim(UNPROTECT(1))),
                            quote(return(ans)),
                            substitute(returnType(RT), list(RT = returnType))
                            )
    } else {
        resultLines <- list(quote(UNPROTECT(1)),
                            substitute(returnType(RT), list(RT = returnType))
                            )
    }
    allLines <- c(list(as.name("{")), convertLines, SEXPsetupLines, list(externalCallLine), resultLines)
    body(fun) <- as.call(allLines)
    ans <- quote(nimbleFunction(run = fun))
    ans <- eval(ans)
}


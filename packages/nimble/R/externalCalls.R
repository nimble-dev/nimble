#' Create a nimbleFunction that wraps a call to external compiled code
#'
#' Given C header information, a function that takes scalars or pointers can be called from a compiled nimbleFunction.  If non-scalar return values are needed, an argument can be selected to behave as the return value in nimble.
#'
#' @param prototype Argument type information.  This can be provided as an R function using \code{nimbleFunction} type declarations or as a list of \code{nimbleType} objects.
#' @param returnType Return object type information.  This can be provided similarly to \code{prototype} as either a \code{nimbleFunction} type declaration or as a \code{nimbleType} object.  In the latter case, the name will be ignored. If there is no return value, this should be \code{void()}.
#' @param Cfun Name of the external function (character).
#' @param headerFile Name (possibly including file path) of the header file where Cfun is declared.
#' @param oFile Name (possibly including path) of the .o file where Cfun has been compiled.  Spaces in the path may cause problems.
#' @param where An optional \code{where} argument passed to \code{setRefClass} for where the reference class definition generated for this nimbleFunction will be stored.  This is needed due to R package namespace issues but should never need to be provided by a user.
#'
#' @author Perry de Valpine
#' @export
#' @details The only argument types allowed in Cfun are \code{double}, \code{int}, and \code{bool}, corresponding to \code{nimbleFunction} types \code{double}, \code{integer}, and \code{logical}, respectively.
#'
#' If the dimensionality is greater than zero, the arguments in \code{Cfun} should be pointers.  This means it will typically be necessary to pass additional integer arguments telling \code{Cfun} the size(s) of non-scalar arguments.
#'
#' The return argument can only be a scalar or void.  Since non-scalar arguments are passed by pointer, you can use an argument to return results from \code{Cfun}.  If you wish to have a \code{nimbleFunction} that uses one argument of \code{Cfun} as a return object, you can wrap the result of \code{nimbleExternalCall} in another \code{nimbleFunction} that allocates the return object.  This is useful for using \code{Cfun} in a \code{nimbleModel}.  See example below.
#'
#' Note that a \code{nimbleExternalCall} can only be executed in a compiled \code{nimbleFunction}, not an uncompiled one.
#' 
#' If you have problems with spaces in file paths (e.g. for \code{oFile}), try compiling everything locally by including \code{dirName = "."} as an argument to \code{compileNimble}.
#' 
#' @return A \code{nimbleFunction} that takes the indicated input arguments, calls \code{Cfun}, and returns the result.
#'
#' @seealso \code{\link{nimbleRcall}} for calling arbitrary R code from compiled \code{nimbleFunction}s.
#' 
#' @examples
#' \dontrun{
#' sink('add1.h')
#' cat('
#'  extern "C" {
#'  void my_internal_function(double *p, double*ans, int n);
#'  }
#' ')
#' sink()
#' sink('add1.cpp') 
#' cat('
#'  #include <cstdio>
#'  #include "add1.h"
#'  void my_internal_function(double *p, double *ans, int n) {
#'    printf("In my_internal_function\\n");
#'      /* cat reduces the double slash to single slash */ 
#'    for(int i = 0; i < n; i++) 
#'      ans[i] = p[i] + 1.0;
#'  }
#' ')
#' sink()
#' system('g++ add1.cpp -c -o add1.o')
#' Radd1 <- nimbleExternalCall(function(x = double(1), ans = double(1),
#' n = integer()){}, Cfun =  'my_internal_function',
#' headerFile = file.path(getwd(), 'add1.h'), returnType = void(),
#' oFile = file.path(getwd(), 'add1.o'))
#' ## If you need to use a function with non-scalar return object in model code,
#' ## you can wrap it  in another nimbleFunction like this:
#' model_add1 <- nimbleFunction(
#'      run = function(x = double(1)) {
#'          ans <- numeric(length(x))
#'          Radd1(x, ans, length(x))
#'          return(ans)
#'          returnType(double(1))
#'      })
#' demoCode <- nimbleCode({
#'      for(i in 1:4) {x[i] ~ dnorm(0,1)} ## just to get a vector
#'      y[1:4] <- model_add1(x[1:4])
#' })
#' demoModel <- nimbleModel(demoCode, inits = list(x = rnorm(4)),
#' check = FALSE, calculate = FALSE)
#' CdemoModel <- compileNimble(demoModel, showCompilerOutput = TRUE)
#' }
nimbleExternalCall <- function(prototype, returnType, Cfun, headerFile, oFile, where = getNimbleFunctionEnvironment()) {
    ## construct a nimbleFunction to wrap a call to Cfun
    force(where)
    returnTypeExpr <- substitute(returnType)
    if(is.list(prototype)) {
        args <- nimbleTypeList2argTypeList(prototype)
        fun <- function(){}
        formals(fun) <- args
    } else {
        if(!is.function(prototype)) stop("Invalid prototype argument")
        fun <- prototype
        args <- formals(fun)
    }
    argsSymTab <- argTypeList2symbolTable(args)
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
    if(inherits(try(returnType, silent = TRUE), 'nimbleType'))
        returnType <- nimbleType2argType(returnType)[[1]]
    else
        returnType <- returnTypeExpr
    returnSymbol <- argType2symbol(returnType)
    ## Make the call to Cfun, without yet a return value
    externalCallExpr <- as.call(c(list(as.name(Cfun)), lapply(replacedArgNames, as.name)))
    externalCallLine <- substitute(asReturnSymbol(ECE, type = RStype, nDim = RSnDim), list(ECE = externalCallExpr, RStype = returnSymbol$type, RSnDim = returnSymbol$nDim))
    returnLines <- list()
    if(returnSymbol$type != 'void') {
        ## Insert the return value as a variable "RETURNVALUE"
        externalCallLine <- substitute(RETURNVALUE <- EXTERNALCALL, list(EXTERNALCALL = externalCallLine))
        ## Insert "return(RETURNVALUE)" and "returnType(<correct type>)"
        returnLines <- list(
            quote(return(RETURNVALUE)),
            substitute(returnType(RT), list(RT = returnType))
        )
    }
    ## put all the lines together
    allLines <- c(list(as.name("{")), convertLines, list(externalCallLine), unconvertLines, returnLines)
    body(fun) <- as.call(allLines)
    ans <- RCfunction(fun, check = FALSE, where = where)
    ## Stick header information into the nfMethodRC
    environment(ans)$nfMethodRCobject$externalHincludes <- paste0('\"',headerFile,'\"')
    ## If the user provided the ofile with a .c or .o ending, remove and replace with .cpp. That is a kluge because later .cpp is replaced with .o.
    oFile <- gsub('\\.c', '', oFile)
    oFile <- gsub('\\.o', '', oFile)
    oFile <- paste0(oFile, '.cpp') ## This later gets stripped and replaced with .o but it must be .cpp to be stripped (until we someday generalize that handling)
    ## Stick the "oFile", with .cpp extension as a kluge, into the nfMethodRC
    if(grepl(" ", oFile)) warning("The space in the oFile name may cause a problem.")
    environment(ans)$nfMethodRCobject$externalCPPincludes <- oFile
    return(ans)
}

#' Make an R function callable from compiled nimbleFunctions (including nimbleModels).
#'
#' Normally compiled nimbleFunctions call other compiled nimbleFunctions.  nimbleRcall enables any R function (with viable argument types and return values) to be called (and evaluated in R) from compiled nimbleFunctions.
#'
#' @param prototype Argument type information for Rfun.  This can be provided as an R function using \code{nimbleFunction} type declarations or as a list of \code{nimbleType} objects.
#' @param returnType Return object type information.  This can be provided similarly to \code{prototype} as either a \code{nimbleFunction} type declaration or as a \code{nimbleType} object.  In the latter case, the name will be ignored. If there is no return value this should be \code{void()}.
#' @param Rfun The name of an R function to be called from compiled nimbleFunctions.
#' @param where An optional \code{where} argument passed to \code{setRefClass} for where the reference class definition generated for this nimbleFunction will be stored.  This is needed due to R package namespace issues but should never need to be provided by a user.
#'
#' @details The \code{nimbleFunction} returned by \code{nimbleRcall} can be used in other \code{nimbleFunction}s.  When called from a compiled \code{nimbleFunction} (including from a model), arguments will be copied according to the declared types, the function named by \code{Rfun} will be called, and the returned object will be copied if necessary.  The example below shows use of an R function in a compiled \code{nimbleModel}.
#'
#' A \code{nimbleFunction} returned by \code{nimbleRcall} can only be used in a compiled \code{nimbleFunction}.  \code{Rfun} itself should work in an uncompiled \code{nimbleFunction}. 
#' 
#' @return A \code{nimbleFunction} that wraps a call to \code{Rfun} with type-declared arguments and return object.
#'
#' @seealso \code{\link{nimbleExternalCall}} for calling externally provided C (or other) compiled code.
#'
#' @export
#' 
#' @author Perry de Valpine
#'
#' @examples
#' \dontrun{
#' ## Say we want an R function that adds 2 to every value in a vector
#' add2 <- function(x) {
#'    x + 2 
#' }
#' Radd2 <- nimbleRcall(function(x = double(1)){}, Rfun = 'add2',
#' returnType = double(1))
#' demoCode <- nimbleCode({
#'     for(i in 1:4) {x[i] ~ dnorm(0,1)} 
#'     z[1:4] <- Radd2(x[1:4])
#' })
#' demoModel <- nimbleModel(demoCode, inits = list(x = rnorm(4)),
#' check = FALSE, calculate = FALSE)
#' CdemoModel <- compileNimble(demoModel)
#' }
nimbleRcall <- function(prototype, returnType, Rfun, where = getNimbleFunctionEnvironment()) {
    force(where)
    returnTypeExpr <- substitute(returnType)
    if(is.list(prototype)) {
        args <- nimbleTypeList2argTypeList(prototype)
        fun <- function(){}
        formals(fun) <- args
    } else {
        if(!is.function(prototype)) stop("Invalid prototype argument")
        fun <- prototype
        args <- formals(fun)
    }
    argsSymTab <- argTypeList2symbolTable(args)
    argNames <- names(args)
    replacedArgNames <- argNames
    convertLines <- list()
    SEXPsetupLines <- list()
    SEXPsetupLines[[1]] <- substitute(SLANGpiece <- SLANG <- PROTECT(Rf_allocVector(LANGSXP, length )), list(length = length(args) + 1))
    SEXPsetupLines[[2]] <- substitute(SET_NEXT_LANG_ARG(SLANGpiece, Rf_install(FUNNAME)), list(FUNNAME = as.character(Rfun)))
    for(i in seq_along(args)) {
        thisSymbol <- argsSymTab$getSymbolObject( argNames[i] )
        thisType <- thisSymbol$type
        thisNdim <- thisSymbol$nDim
        SEXPargName <- paste0(argNames[i], '_SEXP_')
        replacedArgNames[i] <- SEXPargName
        newConvertLine <- substitute( A <- PROTECT(NimArr_2_SEXP(B)), list(A = as.name(SEXPargName), B = as.name(argNames[i]))) ## MakeSEXPcopyLineFrom almost does this but it is at cppDef level
        convertLines[[ length(convertLines) + 1]] <- newConvertLine
        SEXPsetupLines[[ length(SEXPsetupLines) + 1 ]] <- substitute(nimVerbatim(SET_TAG(SLANGpiece, Rf_install(ARGNAME))), list(ARGNAME = as.character(argNames[i])))
        SEXPsetupLines[[ length(SEXPsetupLines) + 1 ]] <- substitute(nimVerbatim(SET_NEXT_LANG_ARG(SLANGpiece, SEXPARGNAME)), list(SEXPARGNAME = as.name(SEXPargName)))
    }
    if(inherits(try(returnType, silent = TRUE), 'nimbleType'))
        returnType <- nimbleType2argType(returnType)[[1]]
    else
        returnType <- returnTypeExpr
    returnSymbol <- argType2symbol(returnType)
    externalCallLines <- list(quote(nimVerbatim(PutRNGstate())),
                              quote(SANS <- PROTECT(Rf_eval(SLANG, R_GlobalEnv))),
                              quote(nimVerbatim(GetRNGstate())))
    if(returnSymbol$type != 'void') {
        resultLines <- list(substitute(declare(ans, RT), list(RT = returnType)),
                            quote(nimVerbatim(SEXP_2_NimArr(SANS, ans))),
                            substitute(nimVerbatim(UNPROTECT(2 + NARGS)), list(NARGS = length(args))),
                            quote(return(ans)),
                            substitute(returnType(RT), list(RT = returnType))
                            )
    } else {
        resultLines <- list(quote(UNPROTECT(1)),
                            substitute(returnType(RT), list(RT = returnType))
                            )
    }
    allLines <- c(list(as.name("{")), convertLines, SEXPsetupLines, externalCallLines, resultLines)
    body(fun) <- as.call(allLines)
    ans <- quote(RCfunction(fun, check = FALSE, where = where))
    ans <- eval(ans)
    ## By calling RCfunction, there is an nfMethodRC with the necessary body
    ## for subsequent compiler processing.
    ## Now we can replace the body of the R function with something that
    ## will actually execute in R.
    Rcall <- as.call(c(list(as.name(as.character(Rfun))),
                       lapply(argNames, as.name)))
    body(ans) <- substitute({ RCALL },
                            list(RCALL = Rcall))
    ans
}

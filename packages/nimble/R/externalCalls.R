#' Create a nimbleFunction that wraps a call to external compiled code
#'
#' Given C header information, a function that takes scalars or pointers can be called from a compiled nimbleFunction.  If non-scalar return values are needed, an argument can be selected to behave as the return value in nimble.
#'
#' @param fun
#' @param returnType
#' @param Cfun
#' @param headerfile
#' @param oFile
#'
#' @author Perry de Valpine
#' @export
#' @details
#'
#' @return
#'
#' @seealso
#'
#' @examples
#' ## not run
nimbleExternalCall <- function(fun, returnType = void(), Cfun, headerFile, oFile) {
    ## construct a nimbleFunction to wrap a call to Cfun
    returnType <- substitute(returnType)
    args <- formals(fun)
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
    ans <- quote(nimbleFunction(run = fun))
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



nimbleExternalCall <- function(fun, returnType = void(), Cfun, headerFile, oFile) {
    ## construct a nimbleFunction
    returnType <- substitute(returnType)
    args <- formals(fun)
    argsSymTab <- nimble:::argTypeList2symbolTable(args)
    argNames <- names(args)
    replacedArgNames <- argNames
    convertLines <- list()
    unconvertLines <- list()
    for(i in seq_along(args)) {
        thisSymbol <- argsSymTab$getSymbolObject( argNames[i] )
        thisType <- thisSymbol$type
        thisNdim <- thisSymbol$nDim
        if(thisNdim > 0) {
            ptrName <- paste0(argNames[i], '_internalPOINTER_')
            replacedArgNames[i] <- ptrName
            newConvertLine <- substitute( A <- nimbleConvert(B), list(A = as.name(ptrName), B = as.name(argNames[i])))
            convertLines[[ length(convertLines) + 1]] <- newConvertLine
            newUnconvertLine <- substitute( nimbleUnconvert(A, B), list(A = as.name(ptrName), B = as.name(argNames[i])))
            unconvertLines[[ length(unconvertLines) + 1]] <- newUnconvertLine
        }
    }
    returnSymbol <- nimble:::argType2symbol(returnType)
    externalCallLine <- as.call(c(list(as.name(Cfun)), lapply(replacedArgNames, as.name)))
    returnLines <- list()
    if(returnSymbol$type != 'void') {
        externalCallLine <- substitute(RETURNVALUE <- EXTERNALCALL, list(EXTERNALCALL = externalCallLine))
        returnLines <- list(
            quote(return(RETURNVALUE)),
            substitute(returnType(RT), list(TR = returnType))
        )
    }
    allLines <- c(list(as.name("{")), convertLines, list(externalCallLine), unconvertLines, returnLines)
    body(fun) <- as.call(allLines)
    ans <- quote(nimbleFunction(run = fun))
    ans <- eval(ans)
    ## right here stick header information into the nfMethodRC
    environment(ans)$nfMethodRCobject$externalHincludes <- paste0('\"',headerFile,'\"')
    oFile <- gsub('\\.c', '', oFile)
    oFile <- gsub('\\.o', '', oFile)
    oFile <- paste0(oFile, '.cpp') ## This later gets stripped and replaced with .o but it must be .cpp to be stripped (until we someday generalize that handling)
    environment(ans)$nfMethodRCobject$externalCPPincludes <- oFile
    return(ans)
}

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


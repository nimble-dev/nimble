
RCfunctionDef <- setRefClass('RCfunctionDef',
                             contains = 'cppFunctionDef',
                             fields = list(
                                 SEXPinterfaceFun = 'ANY',
                                 SEXPinterfaceCname = 'ANY',	#character
                                 RCfunProc = 'ANY' ## a RCfunProcessing object
                                 ), 
                             methods = list(
                                 initialize = function(...) {
                                 	 SEXPinterfaceCname <<- character()
                                     Hincludes <<- c(Hincludes,
                                                     nimbleIncludeFile("NimArr.h"),
                                                     "<Rinternals.h>",
                                                     nimbleIncludeFile("accessorClasses.h"),
                                                     nimbleIncludeFile("nimDists.h"))
                                     CPPincludes <<- c(CPPincludes,
                                                       '<Rmath.h>',
                                                       '<math.h>',
                                                       nimbleIncludeFile("EigenTypedefs.h"),
                                                       nimbleIncludeFile("Utils.h"),
                                                       nimbleIncludeFile("accessorClasses.h"))
                                     CPPusings <<- c(CPPusings) 
                                     callSuper(...)
                                 },
                                 getDefs = function() {
                                     list(.self,
                                          if(!inherits(SEXPinterfaceFun, 'uninitializedField')) SEXPinterfaceFun)
                                 },
                                 getHincludes = function() {
                                     Hinc <- c(Hincludes,
                                               if(!inherits(SEXPinterfaceFun, 'uninitializedField')) SEXPinterfaceFun$getHincludes())
                                     Hinc
                                 },
                                 getCPPincludes = function() {
                                     CPPinc <- c(CPPincludes,
                                                 if(!inherits(SEXPinterfaceFun, 'uninitializedField')) SEXPinterfaceFun$getCPPincludes())
                                     CPPinc
                                 },
                                 getCPPusings = function() {
                                     CPPuse <- unique(c(CPPusings,
                                                        if(!inherits(SEXPinterfaceFun, 'uninitializedField')) SEXPinterfaceFun$getCPPusings()))
                                     CPPuse
                                 },
                                 genNeededTypes = function() {
                                     for(i in seq_along(RCfunProc$neededRCfuns)) {
                                         neededType<- RCfunProc$neededRCfuns[[i]]
                                         if(inherits(neededType, 'nfMethodRC')) {
                                             thisCppDef <- nimbleProject$getRCfunCppDef(neededType, NULLok = TRUE)
                                             if(is.null(thisCppDef)) {
                                                 thisCppDef <- nimbleProject$needRCfunCppClass(neededType, genNeededTypes = TRUE)
                                                 neededTypeDefs[[neededType$uniqueName]] <<- thisCppDef
                                             } else {
                                                 Hincludes <<- c(Hincludes, thisCppDef)
                                                 CPPincludes <<- c(CPPincludes, thisCppDef)
                                             }
                                             next
                                         } else {
                                             stop("There is a neededType for an RCfun that is not another RCfun. That is not allowed.", call. = FALSE)
                                         }
                                     }
                                 },
                                 buildFunction = function(RCfun, parentST = NULL) {
                                     RCfunProc <<- RCfun
                                     name <<- RCfunProc$name
                                     argNames <- RCfunProc$compileInfo$origLocalSymTab$getSymbolNames() ## this has only the original arguments
                                     args <<- symbolTable2cppVars(RCfunProc$compileInfo$newLocalSymTab, argNames, include = argNames, parentST = parentST)
                                     allNames <- RCfunProc$compileInfo$newLocalSymTab$getSymbolNames() ## this has had local variables added
                                     localArgs <- symbolTable2cppVars(RCfunProc$compileInfo$newLocalSymTab, argNames, include = allNames[!(allNames %in% argNames)], parentST = args)
                                     code <<- cppCodeBlock(code = RCfunProc$compileInfo$nimExpr,
                                                           objectDefs = localArgs)
                                     returnType <<- RCfunProc$compileInfo$returnSymbol$genCppVar()
                                     invisible(NULL)
                                 },
                                 buildRwrapperFunCode = function(className = NULL, eval = FALSE, includeLHS = TRUE, returnArgsAsList = TRUE, includeDotSelf = '.self', env = globalenv(), dll = NULL, includeDotSelfAsArg = FALSE) {
                                     returnVoid <- returnType$baseType == 'void'
                                     asMember <- !is.null(className)
                                     argsCode = RCfunProc$RCfun$arguments
                                     argNames <- names(argsCode)
                                     
                                     if(is.character(SEXPinterfaceCname) && is.null(dll) && eval) {
                                         warning("creating a .Call() expression with no DLL information")
                                         browser()                                     
                                     }
                                     dotCall <- substitute(.Call(SEXPname), list(SEXPname = SEXPinterfaceCname))
                                     for(i in seq_along(argNames)) dotCall[[i+2]] <- as.name(argNames[i])
                                     if(asMember & is.character(includeDotSelf)) dotCall[[length(argNames) + 3]] <- as.name(includeDotSelf)
                                     if(returnArgsAsList) {
                                         argNamesAssign <- if(length(argNames) > 0) paste0('\"',argNames, '\"') else character(0)
                                         if(!returnVoid) argNamesAssign <- c(argNamesAssign, '\"return\"')
                                         if(length(argNamesAssign) > 0)
                                             namesAssign <- parse(text = paste0('names(ans) <- c(', paste(argNamesAssign, collapse = ', '), ')'), keep.source = FALSE)[[1]]
                                         else
                                             namesAssign <- quote(ans <- NULL)
                                     } else {
                                         if(length(argNames)+!returnVoid > 0 & !returnVoid)
                                             namesAssign <- parse(text = paste0('ans <- ans[[',length(argNames)+!returnVoid,']]'), keep.source = FALSE)[[1]]
                                         else
                                             namesAssign <- quote(ans <- invisible(NULL))
                                     }
                                     
                                     argNamesCall = argNames
                                     for(i in seq_along(argNamesCall) ){
                                     	if(argsCode[i] != '')
                                     		argNamesCall[i] = paste(argNames[i], " = ", as.character(argsCode[[i]]) )
                                     }
                                     if(includeDotSelfAsArg) argNamesCall <- c(argNamesCall, includeDotSelf)
                                     
                                     funCode <- parse(text = paste0('function(', paste0(argNamesCall, collapse = ','),') A'), keep.source = FALSE)[[1]]
                                     bodyCode <- substitute({ans <- DOTCALL; NAMESASSIGN; ans}, list(DOTCALL = dotCall, NAMESASSIGN = namesAssign))
                                     funCode[[3]] <- bodyCode
                                     funCode[[4]] <- NULL
                                     if(includeLHS) funCode <- substitute(FUNNAME <- FUNCODE, list(FUNNAME = as.name(paste0('R',name)), FUNCODE = funCode))
                                     if(eval) {
                                         fun = eval(funCode) 
                                         environment(fun) = env #??? may want this to be environment() or the default value for env to be environment()
                                         if(!is.null(dll))   {
                                        # replace the name of the symbol in the .Call() with the resolved symbol.
					     body(fun)[[2]][[3]][[2]] = getNativeSymbolInfo(SEXPinterfaceCname, dll)
					 }

                                         fun
                                     } else 
                                         funCode
                                 },
                                 buildSEXPinterfaceFun = function(className = NULL) {
                                     asMember <- !is.null(className)
                                     objects <- symbolTable2cppVars(RCfunProc$compileInfo$origLocalSymTab)
                                     argNames <- RCfunProc$compileInfo$origLocalSymTab$getSymbolNames()
                                     Snames <- character(length(argNames))
                                     copyLines <- list()
                                     interfaceArgs <- symbolTable()
                                     objects$setParentST(interfaceArgs)
                                     returnVoid <- returnType$baseType == 'void'
                                     
                                     for(i in seq_along(argNames)) {
                                         Snames[i] <- Rname2CppName(paste0('S_', argNames[i]))
                                         ## For each argument to the RCfunction we need a corresponding SEXP argument to the interface function
                                         interfaceArgs$addSymbol(cppSEXP(name = Snames[i]))
                                         
                                         ## and we need a line to copy from the SEXP to the local variable
                                         ## The to argument uses the origLocalSymbolObject rather than the objects (which has cppVars) because that has the nDim
                                         ## The name of that and the new one in objects must match
                                         copyLines[[i]] <- buildCopyLineFromSEXP(interfaceArgs$getSymbolObject(Snames[i]),
                                                                                 RCfunProc$compileInfo$origLocalSymTab$getSymbolObject(argNames[i]))
                                     }

                                     RHScall <- as.call(c(list(as.name(name)),
                                                              lapply(argNames, as.name)))
                                     if(asMember) {
                                         ## Add a final argument for the extptr
                                         interfaceArgs$addSymbol(cppSEXP(name = 'SextPtrToObject'))
                                         ## And make the RHScall
                                         RHScall <- substitute(cppMemberDereference(
                                             template(static_cast, cppPtrType(CN))(R_ExternalPtrAddr(SextPtrToObject)), RHS),
                                                               list(CN = as.name(className), RHS = RHScall))
                                     }
                                     
                                     if(returnVoid) {
                                         fullCall <- RHScall
                                     } else {
                                         objects$addSymbol(cppSEXP(name = 'S_returnValue_1234')) ## Object for the return statement: "return(S_returnValue_1234)"
                                         LHSvar <- RCfunProc$compileInfo$returnSymbol$genCppVar()
                                         LHSvar$name <- "LHSvar_1234"
                                         objects$addSymbol(LHSvar)
                                         fullCall <- substitute(LHS <- RHS, list(LHS = as.name(LHSvar$name), RHS = RHScall))
                                     }
                                     ## Put GetRNGstate() and PutRNGstate() around the call.    
                                     fullCall <- substitute({GetRNGstate(); FULLCALL; PutRNGstate()}, list(FULLCALL = fullCall))
                                     
                                     returnAllArgs <- TRUE
                                     ## Pack up all inputs and the return value in a list.
                                     if(returnAllArgs) {
                                         numArgs <- length(argNames)
                                         if(numArgs + !returnVoid > 0) {
                                             objects$addSymbol(cppSEXP(name = 'S_returnValue_LIST_1234'))
                                             returnListLines <- returnCopyLines <- vector('list', length = numArgs+!returnVoid)
                                             allocVectorLine <- substitute(PROTECT(S_returnValue_LIST_1234 <- allocVector(VECSXP, nAp1)), list(nAp1 = numArgs + !returnVoid))
                                             if(numArgs > 0) {
                                                 for(i in 1:numArgs) {
                                                     returnCopyLines[[i]] <- buildCopyLineToSEXP(RCfunProc$compileInfo$origLocalSymTab$getSymbolObject(argNames[i]),
                                                                                                 interfaceArgs$getSymbolObject(Snames[i]))
                                                     returnListLines[[i]] <- substitute(SET_VECTOR_ELT(S_returnValue_LIST_1234, Im1, THISSEXP),
                                                                                        list(Im1 = i-1, THISSEXP = as.name(Snames[i])))
                                                 }
                                             }
                                             if(!returnVoid) {
                                                 rsName <- RCfunProc$compileInfo$returnSymbol$name
                                                 RCfunProc$compileInfo$returnSymbol$name <<- LHSvar$name
                                                 returnCopyLines[[numArgs+1]] <- buildCopyLineToSEXP(RCfunProc$compileInfo$returnSymbol,
                                                                                                     objects$getSymbolObject('S_returnValue_1234'))
                                                 RCfunProc$compileInfo$returnSymbol$name <<- rsName
                                                 returnListLines[[numArgs+1]] <- substitute(SET_VECTOR_ELT(S_returnValue_LIST_1234, I, THISSEXP),
                                                                                            list(I = numArgs, THISSEXP = as.name('S_returnValue_1234')))
                                             }
                                             returnLine <- quote(return(S_returnValue_LIST_1234))
                                             unprotectLine <- substitute(UNPROTECT(N), list(N = numArgs + 1 + !returnVoid))
                                             allCode <- embedListInRbracket(c(copyLines, list(fullCall), list(allocVectorLine),
                                                                              returnCopyLines, returnListLines, list(unprotectLine), list(returnLine)))
                                  
                                         } else { ## No input or return objects
                                             returnLine <- quote(return(R_NilValue))
                                             allCode <- embedListInRbracket(c(copyLines, list(fullCall),
                                                                              list(returnLine)))
                                         }
                                     } else {
                                         writeLines("Haven't written the single return case yet")        
                                     }
                                     SEXPinterfaceCname <<- paste0('CALL_',Rname2CppName(paste0(if(!is.null(className)) paste0(className,'_') else NULL, name))) ##Rname2CppName needed for operator()
                                     SEXPinterfaceFun <<- cppFunctionDef(name = SEXPinterfaceCname,
                                                                         args = interfaceArgs,
                                                                         code = cppCodeBlock(code = RparseTree2ExprClasses(allCode), objectDefs = objects),
                                                                         returnType = cppSEXP(),
                                                                         externC = TRUE,
                                                                         CPPincludes = list(nimbleIncludeFile("RcppUtils.h")))
                                     invisible(NULL)
                                 }
                                 
                                 )
                             )

SEXPscalarConvertFunctions <- list(double  = 'SEXP_2_double',
                                   integer = 'SEXP_2_int',
                                   logical = 'SEXP_2_bool',
                                   character = 'STRSEXP_2_string')


toSEXPscalarConvertFunctions <- list(double  = 'double_2_SEXP',
                                     integer = 'int_2_SEXP',
                                     logical = 'bool_2_SEXP',
                                     character = 'string_2_STRSEXP')

buildCopyLineFromSEXP <- function(fromSym, toSym) {
    if(inherits(toSym, 'symbolBasic')) {
        if(toSym$nDim == 0) {
            ans <- substitute(TO <- CONVERT(FROM), list(TO = as.name(toSym$name),
                                                        FROM = as.name(fromSym$name),
                                                        CONVERT = as.name(SEXPscalarConvertFunctions[[toSym$type]]) ) )
        } else {
            if(toSym$type == 'character') {
                ans <- substitute(STRSEXP_2_vectorString(FROM, TO), list(TO = as.name(toSym$name),
                                                                           FROM = as.name(fromSym$name)))
            } else {
                ans <- substitute(template(SEXP_2_NimArr,NDIM)(FROM, TO), list(TO = as.name(toSym$name),
                                                                               FROM = as.name(fromSym$name),
                                                                               NDIM = toSym$nDim))
            }
        }
        return(ans)
    }
    stop(paste("Error, don't know how to make a SEXP copy line for something of class", class(toSym)))
}

buildCopyLineToSEXP <- function(fromSym, toSym) {
    if(inherits(fromSym, 'symbolBasic')) {
        if(fromSym$nDim == 0) {
            ans <- substitute(PROTECT(TO <- CONVERT(FROM)), list(TO = as.name(toSym$name),
                                                        FROM = as.name(fromSym$name),
                                                        CONVERT = as.name(toSEXPscalarConvertFunctions[[fromSym$type]] ) ) )
        } else {
            if(fromSym$type == 'character') {
                ans <- substitute(PROTECT(TO <- vectorString_2_STRSEXP(FROM)), list(TO = as.name(toSym$name),
                                                                           FROM = as.name(fromSym$name)))
            } else {
                ans <- substitute(PROTECT(TO <- template(NimArr_2_SEXP, NDIM)(FROM)), list(TO = as.name(toSym$name),
                                                                                           FROM = as.name(fromSym$name),
                                                                                           NDIM = fromSym$nDim))
            }
        }
        return(ans)
    }
    stop(paste("Error, don't know how to make a copy line to SEXP for something of class", class(fromSym)))
}


## Fun can be an nfMethodRC (returned by RCfunction()), an RCfunProcessing, or a cppProject (each representing steps of processing)

compileRCfunction <- function(fun, name, fileName, dirName, writeFiles = TRUE, compileCpp = TRUE, loadSO = TRUE, debug = FALSE, debugCpp = FALSE, returnInternals = FALSE, .useLib = UseLibraryMakevars) {
    if(missing(name)) name <- deparse(substitute(fun))
    Cname <- Rname2CppName(name)
    if(missing(fileName)) fileName <- Cname
    givenDirName <- if(missing(dirName)) {
                         dirName <- fileName
                         FALSE
                      } else
                         TRUE

    if(inherits(fun, 'nfMethodRC')) {
        RCFP <- RCfunProcessing$new(fun, Cname)
        RCFP$process(debug = debug, debugCpp = debugCpp)
        RCF <- RCfunctionDef$new()
        RCF$buildFunction(RCFP)
        RCF$buildSEXPinterfaceFun()
    }
    if(inherits(fun, 'RCfunProcessing')) {
        RCF <- fun
    }

    if(!inherits(fun, 'cppProjectClass')) {
        cppProj <- cppProjectClass$new(dirName = dirName)
        cppProj$addFunction(RCF, Cname)
    } else {
        cppProj <- fun
        Cname <- names(cppProj$cppDefs)[1] ## Assume for now there is only 1
          # set the dirName only if the caller explicitly specified it.
          # Otherwise, it should have been set in a previous call that created the cppProj.
        if(givenDirName || (length(cppProj$dirName) == 0 || cppProj$dirName == ""))
           cppProj$dirName = dirName
    }

   
    if(writeFiles) {
        cppProj$writeFiles(Cname)
    }
    if(compileCpp) {
        cppProj$compileFile(Cname, .useLib = .useLib)
    }
    if(loadSO) {
        cppProj$loadSO(Cname)
    }
    if(!loadSO) returnInternals <- TRUE
    if(returnInternals) {
        return(cppProj)
    } else {
        cppProj$cppDefs[[1]]$buildRwrapperFunCode(includeLHS = FALSE, eval = TRUE, dll = cppProj$dll)
    }
}

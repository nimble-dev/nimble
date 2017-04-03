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
                                     const <<- RCfunProc$const
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

                                     # attempt to mimic cxxfunction in inline pkg
                                     # dotCall <- quote( .Call( EXTERNALNAME ))
                                     # dotCall[[1L]] <- .Call
                                     # dotCall[[2L]] <- getNativeSymbolInfo(SEXPinterfaceCname, dll )$address
                                     
                                     # avoid R CMD check problem with registration
                                     # ok not to use getNativeSymbolInfo with a dll argument because SEXPinterfaceCname can't possible be in nimble.so, so it is unique to the project dll.
                                     txt <- ".Call(SEXPname)"
                                     dotCall <- eval(substitute(substitute(txt1, list(SEXPname = SEXPinterfaceCname)), list(txt1 = parse(text = txt)[[1]])))
                                     
                                     # dotCall <- substitute(.Call(SEXPname), list(SEXPname = SEXPinterfaceCname))

                                     for(i in seq_along(argNames)) dotCall[[i+2]] <- as.name(argNames[i])
                                     if(asMember & is.character(includeDotSelf)) dotCall[[length(argNames) + 3]] <- as.name(includeDotSelf)
                                     if(returnArgsAsList) {
                                       ansReturnName <- substitute(ans$return, list())
                                         argNamesAssign <- if(length(argNames) > 0) paste0('\"',argNames, '\"') else character(0)
                                         if(!returnVoid) argNamesAssign <- c(argNamesAssign, '\"return\"')
                                         if(length(argNamesAssign) > 0)
                                             namesAssign <- parse(text = paste0('names(ans) <- c(', paste(argNamesAssign, collapse = ', '), ')'), keep.source = FALSE)[[1]]
                                         else
                                             namesAssign <- quote(ans <- NULL)
                                     } else {
                                       ansReturnName <- substitute(ans, list())
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
                                     if(inherits(RCfunProc$compileInfo$returnSymbol, 'symbolNimbleList')){
                                         ##returnListObj <- RCfunProc$compileInfo$returnSymbol$nlProc$instances[[1]]
                                       ##returnListDefs <- returnListObj$nimbleListDef
                                       ##returnListNestedLists <- returnListObj$nestedListGenList
                                         ## refClassDef <- returnListObj$.refClassDef
                                         ##returnListNew <- new(refClassDef, NLDEFCLASSOBJECT = returnListDefs, 
                                         ##                   NESTEDGENLIST = returnListNestedLists)

                                         nlGenerator <- RCfunProc$compileInfo$returnSymbol$nlProc$nlGenerator
                                         returnListDefs <- nl.getListDef(nlGenerator)
                                         returnListNestedLists <- nl.getNestedGens(nlGenerator)
                                         refClassGen <- nl.getDefinitionContent(nlGenerator, 'nlRefClass')
                                         returnListNew <- refClassGen$new(NLDEFCLASSOBJECT = returnListDefs, 
                                                                          NESTEDGENLIST = returnListNestedLists)
                                         returnNimListGen <- list(new = function(){return(returnListNew)})
                                         
                                       getGenList <- function(nimListGen, genList = list()){
                                         nimList <- nimListGen$new()
                                         if(is.null(genList[[nimList$nimbleListDef$className]]))
                                           genList[[nimList$nimbleListDef$className]] <- nimListGen
                                         nestedNimLists <- nimList$nestedListGenList
                                         for(i in seq_along(nestedNimLists)){
                                           tempGenList <- nestedNimLists[[i]]
                                           genList <- getGenList(tempGenList, genList)
                                         }
                                         return(genList)
                                       }
                                       
                                       genList <- getGenList(returnNimListGen)

                                       listCode <-  substitute({eval(makeNewNimListFromC <- function(name){
                                         genList <- GENLIST
                                         returnList <- genList[[name]]$new();
                                         return(returnList)}, envir = globalenv());
                                         assign("makeNewNimListSEXPRESSIONFromC", makeNewNimListFromC, envir = globalenv())},
                                                               list(GENLIST = genList))
                                     }
                                     else{
                                       listCode <- NULL
                                     }
                                       
                                     funCode <- parse(text = paste0('function(', paste0(argNamesCall, collapse = ','),') A'), keep.source = FALSE)[[1]]
                                     ## the first warning may be removed later if there is no CnativeSymbolInfo_ to be created or if eval is FALSE (as for a nimbleFunction member
                                     if(asMember & is.character(includeDotSelf))
                                         bodyCode <- substitute({
                                           if(is.null(CnativeSymbolInfo_)) {warning("Trying to call compiled nimbleFunction that does not exist (may have been cleared)."); return(NULL)};
                                           if(is.null(DOTSELFNAME)) stop('Object for calling this function is NULL (may have been cleared)');
                                            LISTCODE; ans <- DOTCALL; NAMESASSIGN; ans}, list(DOTCALL = dotCall, NAMESASSIGN = namesAssign,
                                                                                             DOTSELFNAME = includeDotSelf, LISTCODE = listCode))
                                         
                                      else
                                         bodyCode <- substitute({
                                             if(is.null(CnativeSymbolInfo_)) {warning("Trying to call compiled nimbleFunction that does not exist (may have been cleared)."); return(NULL)};
                                              LISTCODE; ans <- DOTCALL; NAMESASSIGN;ans}, 
                                             list(LISTCODE = listCode, DOTCALL = dotCall,  NAMESASSIGN = namesAssign))
                                     funCode[[3]] <- bodyCode
                                     funCode[[4]] <- NULL
                                     if(includeLHS) funCode <- substitute(FUNNAME <- FUNCODE, list(FUNNAME = as.name(paste0('R',name)), FUNCODE = funCode))
                                     if(eval) {
                                         fun = eval(funCode)
                                         newenv <- eval(quote(new.env()), envir = env)
                                         environment(fun) = newenv #??? may want this to be environment() or the default value for env to be environment()
                                         ##environment(fun) = env #??? may want this to be environment() or the default value for env to be environment()
                                         if(!is.null(dll))   {
                                        # replace the name of the symbol in the .Call() with the resolved symbol.
					                              ##body(fun)[[2]][[3]][[2]] = getNativeSymbolInfo(SEXPinterfaceCname, dll)
                                             body(fun)[[4]][[3]][[2]] = quote(CnativeSymbolInfo_)
                                             assign('CnativeSymbolInfo_', getNativeSymbolInfo(SEXPinterfaceCname, dll), envir = newenv)
					 } else {
                                             body(fun)[[2]] <- NULL ## remove the check for valid CnativeSymbolInfo_
                                         }

                                         fun
                                     } else {
                                         funCode[[3]][[2]] <- NULL
                                         funCode
                                     }
                                 },
					 buildSEXPinterfaceFun = function(className = NULL) {
					   asMember <- !is.null(className)
					   objects <- symbolTable2cppVars(RCfunProc$compileInfo$origLocalSymTab)
					   argNames <- RCfunProc$compileInfo$origLocalSymTab$getSymbolNames()
					   Snames <- character(length(argNames))
					   copyLines <- list()
					   conditionalLineList <- list()
					   interfaceArgs <- symbolTable()
					   objects$setParentST(interfaceArgs)
					   returnVoid <- returnType$baseType == 'void'
					   copyLineCounter <- 1
					   
					   for(i in seq_along(argNames)) {
					     Snames[i] <- Rname2CppName(paste0('S_', argNames[i]))
					     ## For each argument to the RCfunction we need a corresponding SEXP argument to the interface function
					     interfaceArgs$addSymbol(cppSEXP(name = Snames[i]))
					     
					     ## and we need a line to copy from the SEXP to the local variable
					     ## The to argument uses the origLocalSymbolObject rather than the objects (which has cppVars) because that has the nDim
					     ## The name of that and the new one in objects must match
					     tempLines <- buildCopyLineFromSEXP(interfaceArgs$getSymbolObject(Snames[i]),
					                                        RCfunProc$compileInfo$origLocalSymTab$getSymbolObject(argNames[i]))
					     if(inherits(RCfunProc$compileInfo$origLocalSymTab$getSymbolObject(argNames[i]), 'symbolNimbleList')){
					       CPPincludes <<- c(CPPincludes, nimbleIncludeFile("smartPtrs.h"))
					     }
					     copyLines <- c(copyLines, tempLines)
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
					       allocVectorLine <- substitute(PROTECT(S_returnValue_LIST_1234 <- allocVector(VECSXP, nAp1)), list(nAp1 = numArgs + !returnVoid))
					       conditionalLineList <- list()
					       returnListLines <- list()
					       if(numArgs > 0) {
					         for(i in 1:numArgs) {
					           argSymTab <- RCfunProc$compileInfo$origLocalSymTab$getSymbolObject(argNames[i])
					           ## generateConditionalLines() is the same as a call to buildCopyLineToSEXP() with one exception:
					           ## if the from (c++) argument is a nimbleList, we first check whether that nimbleList contains a pointer to a sexp object.
					           ## If so, we copy the nl to that object.  If not, we create a new nl sexp object, and copy the c++ nl into that .
					           conditionalLineList <- c(conditionalLineList, 
					                                    generateConditionalLines(argSymTab, interfaceArgs$getSymbolObject(Snames[i])))
					           returnListLines <- c(returnListLines,
					                                substitute(SET_VECTOR_ELT(S_returnValue_LIST_1234, Im1, THISSEXP),
					                                           list(Im1 = i-1, THISSEXP = as.name(Snames[i]))))
					           isList <- inherits(argSymTab, 'symbolNimbleList')
					           if(isList){
					             resetRCopiedFlag  <- paste0(argSymTab$name,"->resetFlags();")
					             resetRCopiedFlagLine <- substitute(cppLiteral(resetText), list(resetText = resetRCopiedFlag))
					             returnListLines <- c(returnListLines,
					                                  resetRCopiedFlagLine)
					           }
					         }
					       }
					       if(!returnVoid) {
					         RCfunProc$compileInfo$returnSymbol$name <<- LHSvar$name
					         returnSymTab <- RCfunProc$compileInfo$returnSymbol
					         
					         conditionalLineList <- c(generateConditionalLines(returnSymTab,
					                                                           objects$getSymbolObject('S_returnValue_1234')), conditionalLineList)   
					         
					         returnListLines <- c(returnListLines, 
					                              substitute(SET_VECTOR_ELT(S_returnValue_LIST_1234, I, THISSEXP),
					                                         list(I = numArgs, THISSEXP = as.name('S_returnValue_1234'))))
					         isList <- inherits(returnSymTab, 'symbolNimbleList')
					         if(isList){
					           resetRCopiedFlag  <- paste0(returnSymTab$name,"->resetFlags();")
					           resetRCopiedFlagLine <- substitute(cppLiteral(resetText), list(resetText = resetRCopiedFlag))
					           returnListLines <- c(returnListLines,
					                                resetRCopiedFlagLine)
					         }
					       }
					       returnLine <- quote(return(S_returnValue_LIST_1234))
					       unprotectLine <- substitute(UNPROTECT(N), list(N = numArgs + 1 + !returnVoid))
					       allCode <- embedListInRbracket(c(copyLines, list(fullCall), list(allocVectorLine),
					                                        conditionalLineList, 
					                                        returnListLines,
					                                        list(unprotectLine),
					                                        list(returnLine)))
					       
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
					 
                             ))


SEXPscalarConvertFunctions <- list(double  = 'SEXP_2_double',
                                   integer = 'SEXP_2_int',
                                   logical = 'SEXP_2_bool',
                                   character = 'STRSEXP_2_string')


toSEXPscalarConvertFunctions <- list(double  = 'double_2_SEXP',
                                     integer = 'int_2_SEXP',
                                     logical = 'bool_2_SEXP',
                                     character = 'string_2_STRSEXP')

buildCopyLineFromSEXP <- function(fromSym, toSym) {
    if(inherits(toSym, c('symbolNimbleList', 'symbolNimbleListGenerator'))){
      ans <- list()
      ansText  <- paste0(toSym$name, " = new ",toSym$nlProc$name, ";")
      ans[[1]] <- substitute(cppLiteral(answerText), list(answerText = ansText))
      ansText  <- paste0(toSym$name, "->copyFromSEXP(", fromSym$name, ");")
      ans[[2]] <- substitute(cppLiteral(answerText), list(answerText = ansText))
      return(ans)
    }
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
    if(inherits(toSym, 'symbolInternalType')) {
        thisInternalType <- as.character(toSym[['argList']][[1]])
        if(thisInternalType == 'indexedNodeInfoClass') {            
            ans <- substitute(TO <- indexedNodeInfo(SEXP_2_vectorDouble(FROM)), list(TO = as.name(toSym$name),
                                                                                     FROM = as.name(fromSym$name)))
            return(ans)
        } else{
            stop(paste("Error, don't know how to make a SEXP copy line for something of class internal type, case", thisInternalType))
        }
    }
    stop(paste("Error, don't know how to make a SEXP copy line for something of class", class(toSym)))
}

buildCopyLineToSEXP <- function(fromSym, toSym, writeCall = FALSE, conditionalText = "") {
    if(inherits(fromSym, c('symbolNimbleList', 'symbolNimbleListGenerator'))){
      if(writeCall == TRUE)  ansText  <- paste0(conditionalText, fromSym$name, "->createNewSEXP();")
      else ansText  <- paste0(conditionalText,  'PROTECT(', toSym$name, ' = ', fromSym$name, "->copyToSEXP());")
      ans <- substitute(cppLiteral(answerText), list(answerText = ansText))
      return(ans)
    }
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
    if(inherits(fromSym, 'symbolInternalType')) {
        thisInternalType <- as.character(fromSym[['argList']][[1]])
        if(thisInternalType == 'indexedNodeInfoClass') {
            ans <- substitute(PROTECT(TO <- (vectorDouble_2_SEXP(FROM))), list(TO = as.name(toSym$name),
                                                                            FROM = as.name(fromSym$name) ) )
            return(ans)
        } else {
            stop(paste("Error, don't know how to make a SEXP copy line for something of class internal type, case", thisInternalType))
        }
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


generateConditionalLines <- function(LHSSymTab, 
                                     RHSSymTab){
  conditionalLines <- list()
  isList <- inherits(LHSSymTab, 'symbolNimbleList')
  if(isList){
    conditionalText <- paste0('if (!(*',  LHSSymTab$name, ').RObjectPointer) ')
    conditionalLines <-  c(conditionalLines, buildCopyLineToSEXP(LHSSymTab, RHSSymTab,
                                            writeCall = TRUE, conditionalText = conditionalText))
  }
  conditionalLines <-  c(conditionalLines, buildCopyLineToSEXP(LHSSymTab, RHSSymTab,
                                                writeCall = FALSE))
  return(conditionalLines)
}

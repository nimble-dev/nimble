cppNimbleListClass <- setRefClass('cppNimbleListClass',
                                  contains = 'cppNimbleClassClass',
                                  fields = list(
                                  ),
                                  methods = list(
                                      initialize = function(nimCompProc, debugCpp = FALSE, fromModel = FALSE, ...) {
                                          callSuper(nimCompProc, debugCpp, fromModel, ...)
                                          inheritance <<- c(inheritance, 'pointedToBase')
                                      },
                                      buildCmultiInterface = function(dll = NULL) {
                                          sym <- if(!is.null(dll))
                                                      getNativeSymbolInfo(SEXPgeneratorFun$name, dll)
                                                  else
                                                      SEXPgeneratorFun$name
                                          # message('CmultiInterface for nimbleList does not exist')
                                          CmultiInterface <<- CmultiNimbleFunctionClass(compiledNodeFun = .self, basePtrCall = sym, project = nimbleProject)
                                      },
                                      buildRgenerator = function(where = globalenv(), dll = NULL) {
                                          sym <- if(!is.null(dll))
                                                     getNativeSymbolInfo(SEXPgeneratorFun$name, dll)
                                                 else
                                                     SEXPgeneratorFun$name
                                           Rgenerator <<- buildNimbleObjInterface(paste0(name,'_refClass') , .self, sym, where = where)
                                          # message('Rgenerator for nimbleList does not exist')
                                      },
                                      buildAll = function(where = where) {
                                        buildCopyFromSexp()
                                        buildCopyToSexp()
                                        buildWriteToSexp()
                                        callSuper(where)
                                      },
                                      buildWriteToSexp = function(){
                                        interfaceArgs <- symbolTable()
                                        newListLine <- list()
                                        returnType <- "SEXP"
                                        listElementTable <- symbolTable()
                                        listElementTable$addSymbol(cppSEXP(name = "S_newNimList"))
                                        newListLine[[1]] <- substitute(PROTECT(S_newNimList <- makeNewNimbleList()),
                                                                           list())
                                        newListLine[[2]] <- substitute(copyToSEXP(S_newNimList), list())
                                        newListLine[[3]] <-   substitute(UNPROTECT(1), list())
                                        returnLine <- list(quote(cppLiteral("return(S_newNimList);")))
                                        allCode <- embedListInRbracket(c(newListLine, returnLine))
                                        functionDefs[[paste0(name, "_writeTo")]] <<- cppFunctionDef(name = "writeToSEXP",
                                                                                                    args = interfaceArgs,
                                                                                                    code = cppCodeBlock(code = RparseTree2ExprClasses(allCode), objectDefs = listElementTable),
                                                                                                    returnType = cppSEXP(),
                                                                                                    externC = FALSE,
                                                                                                    CPPincludes = list(nimbleIncludeFile("RcppUtils.h"),
                                                                                                                       nimbleIncludeFile("smartPtrs.h")))

                                      },
                                      buildCopyToSexp = function(){
                                        functionArgName <-  Rname2CppName(paste0('S_nimList_'))
                                        interfaceArgs <- symbolTable()
                                        interfaceArgs$addSymbol(cppSEXP(name = functionArgName))
                                        argNames <- nimCompProc$symTab$getSymbolNames()
                                        protectLines <- list()
                                        writeToSexpLines <- list()
                                        copyToListLines <- list()
                                        nameLines <- list()
                                        Snames <- character(length(argNames))
                                        Rnames <- character(length(argNames))
                                        returnType <- "void"
                                        listElementTable <- symbolTable()
                                        numArgs <- length(argNames)
                                        
                                        environmentCPPName <- Rname2CppName('.xData')
                                        listElementTable$addSymbol(cppSEXP(name = environmentCPPName))
                                        envLine <- substitute({PROTECT(ENVNAME <- allocVector(STRSXP, 1));
                                                              SET_STRING_ELT(ENVNAME, 0, mkChar(".xData"));}, 
                                                              list(ENVNAME = as.name(environmentCPPName)))
                                                           
                                        printLine <-  substitute({Rprintf("copyto: doub is : %f\n", doub);})

                                             
                                        for(i in seq_along(argNames)) {
                                          Snames[i] <- Rname2CppName(paste0('S_', argNames[i]))
                                          listElementTable$addSymbol(cppSEXP(name = Snames[i]))
                                          Rnames[i] <- paste0('R_', argNames)
                                          listElementTable$addSymbol(cppSEXP(name = Rnames[i]))
                                          nameLines[[i]] <- substitute({PROTECT(RNAME <- allocVector(STRSXP, 1));
                                            SET_STRING_ELT(RNAME, 0, mkChar(ARGNAME));}, 
                                            list(RNAME = as.name(Rnames[i]), ARGNAME = argNames[i]))
                                          
                                          writeToSexpLines[[i]] <- buildCopyLineToSEXP(nimCompProc$symTab$getSymbolObject(argNames[i]),
                                                                                       listElementTable$getSymbolObject(Snames[i]))
                                          copyToListLines[[i]] <- substitute(defineVar(install(ARGNAME), VALUE, GET_SLOT(S_nimList_, XDATA)),
                                                                             list(ARGNAME = argNames[i], VALUE = as.name(Snames[i]),
                                                                                  XDATA = as.name(environmentCPPName)))
                                            
                                        }
                                        unprotectLineNoReturn <- list(substitute(UNPROTECT(N), list(N = 2*numArgs + 1)))
                                        allCode <- embedListInRbracket(c(list(printLine), list(envLine), nameLines, writeToSexpLines,
                                                                         copyToListLines, unprotectLineNoReturn))
                                        
                                        functionDefs[[paste0(name, "_copyTo")]] <<- cppFunctionDef(name = "copyToSEXP",
                                                                               args = interfaceArgs,
                                                                               code = cppCodeBlock(code = RparseTree2ExprClasses(allCode), objectDefs = listElementTable),
                                                                               returnType = cppVoid(),
                                                                               externC = FALSE,
                                                                               CPPincludes = list(nimbleIncludeFile("RcppUtils.h"),
                                                                                                  nimbleIncludeFile("smartPtrs.h")))
                                      },
                                      buildCopyFromSexp = function(){
                                        functionArgName <-  Rname2CppName(paste0('S_nimList_'))
                                        interfaceArgs <- symbolTable()
                                        interfaceArgs$addSymbol(cppSEXP(name = functionArgName))
                                        
                                        objects <- symbolTable2cppVars(nimCompProc$symTab)
                                        argNames <- nimCompProc$symTab$getSymbolNames() ##not actually arg names any more, but names of nimbleListElements
                                        nameLines <- list()
                                        copyFromListLines <- list()
                                        copyLines <- list()
                                        Snames <- character(length(argNames))
                                        Rnames <- character(length(argNames))
                                        returnType <- "void"
                                        listElementTable <- symbolTable()
                                        
                                        environmentCPPName <- Rname2CppName('.xData')
                                        listElementTable$addSymbol(cppSEXP(name = environmentCPPName))
                                        envLine <- substitute({PROTECT(ENVNAME <- allocVector(STRSXP, 1));
                                          SET_STRING_ELT(ENVNAME, 0, mkChar(".xData"));}, 
                                          list(ENVNAME = as.name(environmentCPPName)))
                                        environmentCPP <- Rname2CppName("nimListEnvironment")
                                        listElementTable$addSymbol(cppSEXP(name = environmentCPP))
                                        envLine2 <- substitute({
                                          ENVNAME <- GET_SLOT(S_nimList_, XDATA);}, 
                                          list(ENVNAME = as.name(environmentCPP),
                                               XDATA = as.name(environmentCPPName)))
                                        
                                        printLine <-  substitute({Rprintf("copyfrom: doub is : %f\n", doub);})
                                        
                                        for(i in seq_along(argNames)) {
                                          Snames[i] <- Rname2CppName(paste0('S_', argNames[i]))
                                          listElementTable$addSymbol(cppSEXP(name = Snames[i]))
                                          
                                          Rnames[i] <- paste0('R_', argNames)
                                          listElementTable$addSymbol(cppSEXP(name = Rnames[i]))
                                          nameLines[[i]] <- substitute({PROTECT(RNAME <- allocVector(STRSXP, 1));
                                            SET_STRING_ELT(RNAME, 0, mkChar(ARGNAME));}, 
                                            list(RNAME = as.name(Rnames[i]), ARGNAME = argNames[i]))
                                            
                                            copyFromListLines[[i]] <- substitute(PROTECT(SVAR <- findVarInFrame(ENVNAME, install(RNAME))),
                                                                               list(RNAME = argNames[i], ENVNAME = as.name(environmentCPP),
                                                                                    SVAR = as.name(Snames[i])))
                                            copyLines[[i]] <- buildCopyLineFromSEXP(listElementTable$getSymbolObject(Snames[i]),
                                                                                    nimCompProc$symTab$getSymbolObject(argNames[i]))
                                        }
                                        
                                        numArgs <- length(argNames)
                                        unprotectLine <- substitute(UNPROTECT(N), list(N = 2*numArgs + 1))
                                        allCode <- embedListInRbracket(c(envLine, envLine2, nameLines, copyFromListLines, copyLines,
                                                                         list(printLine), list(unprotectLine)))
                                        functionDefs[[paste0(name, "_copyFrom")]] <<- cppFunctionDef(name = "copyFromSEXP",
                                                                                args = interfaceArgs,
                                                                                code = cppCodeBlock(code = RparseTree2ExprClasses(allCode), objectDefs = listElementTable),
                                                                                returnType = cppVoid(),
                                                                                externC = FALSE,
                                                                                CPPincludes = list(nimbleIncludeFile("RcppUtils.h"),
                                                                                                   nimbleIncludeFile("smartPtrs.h")))
                                      }
                                      )
                                  )


buildCopyLineIntoSEXP <- function(fromSym, toSym) {
  if(inherits(fromSym, c('symbolNimbleList', 'symbolNimbleListGenerator'))){
    ansText  <- paste0(fromSym$name, "->copyToSEXP(", toSym$name, ");")
    ans <- substitute(cppLiteral(ANSWERTEXT), list(ANSWERTEXT = ansText))
    return(ans)
  }
  if(inherits(fromSym, 'symbolBasic')) {
    if(fromSym$nDim == 0) {
      ans <- substitute(PROTECT(CONVERT(FROM, TO)), list(TO = as.name(toSym$name),
                                                           FROM = as.name(fromSym$name),
                                                           CONVERT = as.name(copyToSEXPscalarConvertFunctions[[fromSym$type]] ) ) )
    } else {
      if(fromSym$type == 'character') {
        ans <- substitute(PROTECT(copyVectorString_2_STRSEXP(FROM, TO)), list(TO = as.name(toSym$name),
                                                                            FROM = as.name(fromSym$name)))
      } else {
        ans <- substitute(PROTECT(template(copyNimArr_2_SEXP, NDIM)(FROM, TO)), list(TO = as.name(toSym$name),
                                                                                   FROM = as.name(fromSym$name),
                                                                                   NDIM = fromSym$nDim))
      }
    }
    return(ans)
  }
  if(inherits(fromSym, 'symbolInternalType')) {
    thisInternalType <- as.character(fromSym[['argList']][[1]])
    if(thisInternalType == 'indexedNodeInfoClass') {
      ans <- substitute(PROTECT(TO <- (copyVectorDouble_2_SEXP(FROM))), list(TO = as.name(toSym$name),
                                                                         FROM = as.name(fromSym$name) ) )
      return(ans)
    } else {
      stop(paste("Error, don't know how to make a SEXP copy line for something of class internal type, case", thisInternalType))
    }
  }
  stop(paste("Error, don't know how to make a copy line to SEXP for something of class", class(fromSym)))
}

copyToSEXPscalarConvertFunctions <- list(double  = 'copyDouble_2_SEXP',
                                     integer = 'copyInt_2_SEXP',
                                     logical = 'copyBool_2_SEXP',
                                     character = 'copyString_2_STRSEXP')
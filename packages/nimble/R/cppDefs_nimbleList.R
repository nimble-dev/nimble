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
                                        listElementTable$addSymbol(cppSEXP(name = "S_listName"))
                                      
                                        newListLine[[1]] <- substitute({PROTECT(S_listName <- allocVector(STRSXP, 1));
                                          SET_STRING_ELT(S_listName, 0, mkChar(LISTNAME));}, 
                                          list(LISTNAME = nimCompProc$nimbleListObj$className))
                                        newListLine[[2]] <- substitute(PROTECT(S_newNimList <- makeNewNimbleList(S_listName)),
                                                                           list())
                                        newListLine[[3]] <- substitute(printf("made new nl"))
                                        
                                        newListLine[[4]] <- substitute(copyToSEXP(S_newNimList), list())
                                        newListLine[[5]] <-   substitute(UNPROTECT(2), list())
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
                                        Snames <- character(length(argNames))
                                        returnType <- "void"
                                        listElementTable <- symbolTable()
                                        numArgs <- length(argNames)
                                        
                                        environmentCPPName <- Rname2CppName('S_.xData')  ## create SEXP for ref class environment 
                                        listElementTable$addSymbol(cppSEXP(name = environmentCPPName))
                                        envLine <- substitute({PROTECT(ENVNAME <- allocVector(STRSXP, 1));
                                                              SET_STRING_ELT(ENVNAME, 0, mkChar(".xData"));}, 
                                                              list(ENVNAME = as.name(environmentCPPName)))
                                                           browser()
                                        for(i in seq_along(argNames)) {
                                          Snames[i] <- Rname2CppName(paste0('S_', argNames[i]))
                                          listElementTable$addSymbol(cppSEXP(name = Snames[i]))
                                          writeToSexpLines[[i]] <- buildCopyLineToSEXP(nimCompProc$symTab$getSymbolObject(argNames[i]),
                                                                                       listElementTable$getSymbolObject(Snames[i]), returnCall = TRUE)
                                          copyToListLines[[i]] <- substitute(defineVar(install(ARGNAME), VALUE, GET_SLOT(S_nimList_, XDATA)),
                                                                             list(ARGNAME = argNames[i], VALUE = as.name(Snames[i]),
                                                                                  XDATA = as.name(environmentCPPName)))
                                        }
                                        unprotectLineNoReturn <- list(substitute(UNPROTECT(N), list(N = numArgs + 1)))
                                        allCode <- embedListInRbracket(c(list(envLine), writeToSexpLines,
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
                                        copyFromListLines <- list()
                                        copyLines <- list()
                                        Snames <- character(length(argNames))
                                        returnType <- "void"
                                        listElementTable <- symbolTable()
                                        
                                        environmentCPPName <- Rname2CppName('S_.xData')  ## create SEXP for ref class environment 
                                        listElementTable$addSymbol(cppSEXP(name = environmentCPPName))
                                        envLine <- substitute({PROTECT(ENVNAME <- allocVector(STRSXP, 1));
                                          SET_STRING_ELT(ENVNAME, 0, mkChar(".xData"));}, 
                                          list(ENVNAME = as.name(environmentCPPName)))

                                        for(i in seq_along(argNames)) {
                                          Snames[i] <- Rname2CppName(paste0('S_', argNames[i]))
                                          listElementTable$addSymbol(cppSEXP(name = Snames[i]))
                                            copyFromListLines[[i]] <- substitute(PROTECT(SVAR <- findVarInFrame(GET_SLOT(S_nimList_, XDATA), install(ARGNAME))),
                                                                               list(ARGNAME = argNames[i], 
                                                                                    SVAR = as.name(Snames[i]),
                                                                                    XDATA = as.name(environmentCPPName)))
                                            copyLine <-  buildCopyLineFromSEXP(listElementTable$getSymbolObject(Snames[i]),
                                                                               nimCompProc$symTab$getSymbolObject(argNames[i]))
                                            copyLines <- c(copyLines, copyLine)
                                        }
                                        numArgs <- length(argNames)
                                        unprotectLine <- substitute(UNPROTECT(N), list(N = numArgs + 1))
                                        allCode <- embedListInRbracket(c(envLine, copyFromListLines, copyLines,
                                                                         list(unprotectLine)))
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

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
                                        callSuper(where)
                                      },
                                      buildCopyToSexp = function(){
                                        functionArgName <-  Rname2CppName(paste0('S_nimList_'))
                                        interfaceArgs <- symbolTable()
                                        interfaceArgs$addSymbol(cppSEXP(name = functionArgName))
                                        argNames <- nimCompProc$symTab$getSymbolNames()
                                        
                                        objectsST <- nimCompProc$symTab$copy()
                                        for(i in seq_along(argNames)){
                                          objectsST$setSymbolField(argNames[[i]], "name", paste0("temp_", argNames[[i]], "_"))
                                        }
                                        objects <- symbolTable2cppVars(objectsST)
                                        
                                         ##not actually arg names any more, but names of nimbleListElements
                                        protectLines <- list()
                                        copyLines <- list()
                                        rewriteLines <- list()
                                        Snames <- character(length(argNames))
                                        returnType <- "void"
                                        printLines <- list()
                                        printLines[[1]] <- substitute(Rprintf("Running copyTo func"), list())
                                        # objects <= symbolTable()
                                        # objects$setParentST(listElementTable)
                                        for(i in seq_along(argNames)) {
                                          Snames[i] <- Rname2CppName(paste0('S_', argNames[i]))
                                          objects$addSymbol(cppSEXP(name = Snames[i]))
                                          # protectLines[[i]] <- cppLiteral(paste0(
                                          #   "PROTECT(",Snames[i]," = getClassElement(S_nimList_, ",argNames[i],"))"))
                                          
                                          copyLines[[i]] <- buildCopyLineFromSEXP(objects$getSymbolObject(Snames[i]),
                                                                                  objectsST$getSymbolObjects()[[i]])
                                          protectLines[[i]] <- substitute(PROTECT(Svar <- getClassElement(S_nimList_, argName)),
                                                                          list(Svar = as.name(Snames[i]), argName = argNames[i]))
                                          rewriteLines[[i]] <- substitute(tempNm <- argNm, list(argNm = as.name(argNames[i]),
                                                                                                tempNm = as.name(objects$getSymbolNames()[i])))
                                        }
                                        
                                        numArgs <- length(argNames)
                                        unprotectLine <- substitute(UNPROTECT(N), list(N = numArgs))
                                        allCode <- embedListInRbracket(c(printLines, protectLines, copyLines, rewriteLines, list(unprotectLine)))
                                       
                                        # SEXPinterfaceCname <<- paste0('CALL_',Rname2CppName(paste0(if(!is.null(className)) paste0(className,'_') else NULL, name))) ##Rname2CppName needed for operator()
                                              
                                        #may not want [[name]] below but instead something else ie 'copier', 'listCopier', etc.
                                        functionDefs[[paste0(name, "_copyTo")]] <<- cppFunctionDef(name = "copyToSEXP",
                                                                               args = interfaceArgs,
                                                                               code = cppCodeBlock(code = RparseTree2ExprClasses(allCode), objectDefs = objects),
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
                                        protectLines <- list()
                                        copyLines <- list()
                                        Snames <- character(length(argNames))
                                        returnType <- "void"
                                        listElementTable <- symbolTable()
                                        # objects <= symbolTable()
                                        # objects$setParentST(listElementTable)
                                        printLines <- list()
                                        printLines[[1]] <- substitute(Rprintf("Running copyFrom func\n"), list())
                                        
                                        printLinesLater <- list()
                                        printLinesLater[[1]] <- as.name(paste0("int nn = LENGTH(S_a);"))
                                        printLinesLater[[2]] <- substitute(Rprintf("got length a= %d\n", nn), list()) 
                                        
                                        
                                        # printLinesLater[[1]] <- as.name(paste0("double A;"))
                                        # printLinesLater[[2]] <- substitute(A <- REAL(S_a)[1], list())
                                        # printLinesLater[[3]] <- substitute(Rprintf("got a= %g\n", A), list()) 
                                        
                                        for(i in seq_along(argNames)) {
                                          Snames[i] <- Rname2CppName(paste0('S_', argNames[i]))
                                          listElementTable$addSymbol(cppSEXP(name = Snames[i]))
                                          # protectLines[[i]] <- cppLiteral(paste0(
                                          #   "PROTECT(",Snames[i]," = getClassElement(S_nimList_, ",argNames[i],"))"))
                                          
                                          protectLines[[i]] <- substitute(PROTECT(Svar <- getClassElement(S_nimList_, argName)),
                                                                          list(Svar = as.name(Snames[i]), argName = argNames[i]))
                                          copyLines[[i]] <- buildCopyLineFromSEXP(listElementTable$getSymbolObject(Snames[i]),
                                                                                  nimCompProc$symTab$getSymbolObject(argNames[i]))
                                        }
                                        
                                        numArgs <- length(argNames)
                                        unprotectLine <- substitute(UNPROTECT(N), list(N = numArgs))
                                        allCode <- embedListInRbracket(c(printLines, protectLines, printLinesLater, copyLines, list(unprotectLine)))
                                        
                                        # SEXPinterfaceCname <<- paste0('CALL_',Rname2CppName(paste0(if(!is.null(className)) paste0(className,'_') else NULL, name))) ##Rname2CppName needed for operator()
                                        
                                        #may not want [[name]] below but instead something else ie 'copier', 'listCopier', etc.
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

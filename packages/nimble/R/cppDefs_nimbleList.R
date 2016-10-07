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
                                        argNames <- nimCompProc$symTab$getSymbolNames()
                                        protectLines <- list()
                                        copyLinesReturn <- list()
                                        writeLinesReturn <- list()
                                        nameLinesReturn <- list()
                                        
                                        Snames <- character(length(argNames))
                                        returnType <- "SEXP"
                                        listElementTable <- symbolTable()
                                        numArgs <- length(argNames)
                                        
                                        nameLineInitText <- paste0("SEXP nms = PROTECT(allocVector(STRSXP, ", numArgs, "));")
                                        nameLinesReturn[[1]] <- substitute(cppLiteral(nameLineText), 
                                                                           list(nameLineText = nameLineInitText))
                                        for(i in seq_along(argNames)) {
                                          Snames[i] <- Rname2CppName(paste0('S_', argNames[i]))
                                          listElementTable$addSymbol(cppSEXP(name = Snames[i]))
                                          copyLinesReturn[[i]] <- buildCopyLineToSEXP(nimCompProc$symTab$getSymbolObject(argNames[i]),
                                                                                    listElementTable$getSymbolObject(Snames[i]))
                                          writeLinesReturn[[i]] <- substitute(SET_VECTOR_ELT(S_returnList, index, Sname),
                                                                              list(index = i-1, Sname = as.name(Snames[i])))
                                          
                                          nameLineReturnText <- paste0("SET_STRING_ELT(nms, ", i-1, ', mkChar("', argNames[i], '"));')
                                          nameLinesReturn[[i+1]] <-  substitute(cppLiteral(nameLineText), 
                                                                                list(nameLineText = nameLineReturnText))
                                        }
                                        
                                        nameLinesReturn[[numArgs+2]] <-quote(cppLiteral("setAttrib(S_returnList, R_NamesSymbol, nms);"))
      
                                        copyLinesReturn[[numArgs+1]] <- quote(cppLiteral("SEXP S_returnList;"))
                                        copyLinesReturn[[numArgs+2]] <- substitute(S_returnList <- PROTECT(allocVector(VECSXP, numArgs)),
                                                                                   list(numArgs = numArgs))
                                        
                                        unprotectLineReturn <- list(substitute(UNPROTECT(N), list(N = numArgs+1)))
                                        returnLine <- list(quote(cppLiteral("return(S_returnList);")))

                                        allCode <- embedListInRbracket(c(copyLinesReturn, writeLinesReturn, nameLinesReturn,
                                                                         unprotectLineReturn, returnLine))
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
                                        copyLinesNoReturn <- list()
                                        copyLinesReturn <- list()
                                        rewriteLines <- list()
                                        writeLinesReturn <- list()
                                        Snames <- character(length(argNames))
                                        returnType <- "void"
                                        listElementTable <- symbolTable()
                                        numArgs <- length(argNames)
                                        
                                        printLines <- list()
                                        printLines[[1]] <- as.name("double aVal;")
                                        printLines[[2]] <- as.name("aVal = a(0);")
                                        printLines[[3]] <- substitute(Rprintf(" a = %g\n", aVal), list())
                                        # printLines[[1]] <- as.name("double aVal;")
                                        # printLines[[2]] <- as.name("aVal = REAL(S_a)[0];")
                                        # printLines[[3]] <- substitute(Rprintf("got a = %g\n", aVal), list())
                                        # objects$setParentST(listElementTable)
                                        for(i in seq_along(argNames)) {
                                          Snames[i] <- Rname2CppName(paste0('S_', argNames[i]))
                                          listElementTable$addSymbol(cppSEXP(name = Snames[i]))
                                          
                                          # ##test here
                                          # protectLines[[i]] <- buildCopyLineToSEXP(nimCompProc$symTab$getSymbolObject(argNames[i]),
                                          #                                             listElementTable$getSymbolObject(Snames[i]))
                                          # copyLinesNoReturn[[i]] <- substitute(setClassElement(S_nimList_, Svar, argName),
                                          #                                     list(Svar = as.name(Snames[i]), argName = argNames[i]))
                                          # 
                                          # ## end here
                                          protectLines[[i]] <- substitute(PROTECT(Svar <- getClassElement(S_nimList_, argName)),
                                                                          list(Svar = as.name(Snames[i]), argName = argNames[i]))

                                          copyText  <- paste0("copyNimArr_2_SEXP<1>( ", argNames[i], ", &", Snames[i], ");")
                                          copyLinesNoReturn[[i]] <- substitute(cppLiteral(copText), list(copText = copyText))
                                        }
                                        unprotectLineNoReturn <- list(substitute(UNPROTECT(N), list(N = numArgs)))
                                        allCode <- embedListInRbracket(c(protectLines, copyLinesNoReturn, unprotectLineNoReturn))

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
                                        protectLines <- list()
                                        copyLines <- list()
                                        Snames <- character(length(argNames))
                                        returnType <- "void"
                                        listElementTable <- symbolTable()
                                        # objects <= symbolTable()
                                        # objects$setParentST(listElementTable)
                                       
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
                                        allCode <- embedListInRbracket(c(protectLines, copyLines, list(unprotectLine)))
                                        
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

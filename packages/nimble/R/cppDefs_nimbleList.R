cppNimbleListClass <- setRefClass('cppNimbleListClass',
                                  contains = 'cppNimbleClassClass',
                                  fields = list(
                                    eigenList = 'ANY'
                                  ),
                                  methods = list(
                                      initialize = function(nimCompProc, debugCpp = FALSE, fromModel = FALSE, ...) {
                                        callSuper(nimCompProc, debugCpp, fromModel, ...)
                                        inheritance <<- c(inheritance, 'pointedToBase')
                                      },
                                      getDefs = function() {
                                        if(eigenList){
                                          if(inherits(SEXPfinalizerFun, 'uninitializedField'))
                                            list(SEXPgeneratorFun)
                                          else
                                            list(SEXPgeneratorFun, SEXPfinalizerFun)
                                        }
                                        else callSuper()
                                        # callSuper()
                                      },
                                      buildCmultiInterface = function(dll = NULL) {
                                          sym <- if(!is.null(dll))
                                                      getNativeSymbolInfo(SEXPgeneratorFun$name, dll)
                                                  else
                                                      SEXPgeneratorFun$name
                                          CmultiInterface <<- CmultiNimbleFunctionClass(compiledNodeFun = .self, basePtrCall = sym, project = nimbleProject)
                                      },
                                      buildRgenerator = function(where = globalenv(), dll = NULL) {
                                          sym <- if(!is.null(dll))
                                                     getNativeSymbolInfo(SEXPgeneratorFun$name, dll)
                                                 else
                                                     SEXPgeneratorFun$name
                                           Rgenerator <<- buildNimbleObjInterface(paste0(name,'_refClass') , .self, sym, where = where)
                                      },
                                      addListSymbols = function(){
                                        RPointerSymbol <- symbolSEXP(name = 'RObjectPointer')
                                        objectDefs$addSymbol(RPointerSymbol$genCppVar())
                                        RCopiedSymbol <- symbolBasic(name = 'RCopiedFlag', type = 'logical',
                                                                     nDim = 0)
                                        objectDefs$addSymbol(RCopiedSymbol$genCppVar())
                                      },
                                      buildAll = function(where = where) {
                                        addListSymbols()
                                        if(eigenList == TRUE){
                                          # makeCppNames()
                                          # buildConstructorFunctionDef()
                                          # buildSEXPgenerator(finalizer = "namedObjects_Finalizer")
                                          # buildRgenerator(where = where)
                                          # buildCmultiInterface()
                                          callSuper(where)
                                        }
                                        else{
                                          buildCopyFromSexp()
                                          buildCopyToSexp()
                                          buildCreateNewSexp()
                                          callSuper(where)
                                        }
                                      },
                                      buildCreateNewSexp = function(){
                                        interfaceArgs <- symbolTable()
                                        newListLine <- list()
                                        listElementTable <- symbolTable()
                                        listElementTable$addSymbol(cppSEXP(name = "S_newNimList"))
                                        listElementTable$addSymbol(cppSEXP(name = "S_listName"))
                                      
                                        newListLine[[1]] <- substitute({PROTECT(S_listName <- allocVector(STRSXP, 1));
                                          SET_STRING_ELT(S_listName, 0, mkChar(LISTNAME));}, 
                                          list(LISTNAME = nimCompProc$nimbleListObj$className))
                                        newListLine[[2]] <- substitute(PROTECT(S_newNimList <- makeNewNimbleList(S_listName)),
                                                                           list())
                                        newListLine[[3]] <- quote(cppLiteral('RObjectPointer = S_newNimList;'))
                                        newListLine[[4]] <-   substitute(UNPROTECT(2), list())
                                        allCode <- embedListInRbracket(c(newListLine))
                                        functionDefs[[paste0(name, "_createNewSEXP")]] <<- cppFunctionDef(name = "createNewSEXP",
                                                                                                    args = interfaceArgs,
                                                                                                    returnType = cppVoid(),
                                                                                                    code = cppCodeBlock(code = RparseTree2ExprClasses(allCode), objectDefs = listElementTable),
                                                                                                    externC = FALSE,
                                                                                                    CPPincludes = list(nimbleIncludeFile("RcppUtils.h"),
                                                                                                                       nimbleIncludeFile("smartPtrs.h")))

                                      },
                                      buildCopyToSexp = function(){
                                        interfaceArgs <- symbolTable()
                                        elementNames <- nimCompProc$symTab$getSymbolNames()
                                        protectLines <- list()
                                        writeToSexpLines <- list()
                                        copyToListLines <- list()
                                        Snames <- character(length(elementNames))
                                        listElementTable <- symbolTable()
                                        numElements <- length(elementNames)
                                        conditionalLineList <- list()
                                        
                                        conditionalClauseStart <- list(quote(cppLiteral('if (!RCopiedFlag){')))
                                        conditionalClauseEnd <- list(quote(cppLiteral('}')))
                                        environmentCPPName <- Rname2CppName('S_.xData')  ## create SEXP for ref class environment 
                                        listElementTable$addSymbol(cppSEXP(name = environmentCPPName))
                                        envLine <- substitute({PROTECT(ENVNAME <- allocVector(STRSXP, 1));
                                          SET_STRING_ELT(ENVNAME, 0, mkChar(".xData"));}, 
                                          list(ENVNAME = as.name(environmentCPPName)))
                                        
                                        for(i in seq_along(elementNames)){
                                          Snames[i] <- Rname2CppName(paste0('S_', elementNames[i]))
                                          listElementTable$addSymbol(cppSEXP(name = Snames[i]))
                                          elementSymTab <- nimCompProc$symTab$getSymbolObject(elementNames[i])
                                          conditionalLineList <- c(conditionalLineList, generateConditionalLines(nimCompProc$symTab$getSymbolObject(elementNames[i]),
                                                                                                                 listElementTable$getSymbolObject(Snames[i])))
                                          
                                          copyToListLines[[i]] <- substitute(defineVar(install(ELEMENTNAME), VALUE, GET_SLOT(ROBJ, XDATA)),
                                                                             list(ELEMENTNAME = elementNames[i], VALUE = as.name(Snames[i]),
                                                                                  ROBJ = as.name('RObjectPointer'),
                                                                                  XDATA = as.name(environmentCPPName)))
                                        }
                                        
                                        setFlagLine <- list(substitute(RCopiedFlag <- true,
                                                                       list()))
                                        returnLine <- list(substitute(return(ROBJ),
                                                           list(ROBJ = as.name('RObjectPointer'))))
                                        unprotectLine <- list(substitute(UNPROTECT(N), list(N = numElements + 1)))
                                        allCode <- embedListInRbracket(c(conditionalClauseStart, list(envLine),  conditionalLineList,
                                                                         copyToListLines, setFlagLine, unprotectLine,
                                                                         conditionalClauseEnd, returnLine))
                                        functionDefs[[paste0(name, "_copyTo")]] <<- cppFunctionDef(name = "copyToSEXP",
                                                                               args = interfaceArgs,
                                                                               code = cppCodeBlock(code = RparseTree2ExprClasses(allCode), objectDefs = listElementTable),
                                                                               returnType = cppSEXP(),
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
                                        storeSexpLine <- list(quote(cppLiteral('RObjectPointer =  S_nimList_;')))

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
                                        allCode <- embedListInRbracket(c(storeSexpLine, envLine, 
                                                                         copyFromListLines, copyLines,
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

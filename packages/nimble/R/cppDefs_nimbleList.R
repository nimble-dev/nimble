cppNimbleListClass <- setRefClass('cppNimbleListClass',
                                  contains = 'cppNimbleClassClass',
                                  fields = list(
                                      eigenList = 'ANY',
                                      ptrCastFun = 'ANY'
                                  ),
                                  methods = list(
                                      initialize = function(nimCompProc, debugCpp = FALSE, fromModel = FALSE, ...) {
                                        callSuper(nimCompProc, debugCpp, fromModel, ...)
                                        inheritance <<- c(inheritance, 'pointedToBase')
                                      },
                                      getDefs = function() {
                                        ans <- if(eigenList){ ## prevents redefinition of EIGEN_EIGENCLASS or EIGEN_SVDCLASS
                                          if(inherits(SEXPfinalizerFun, 'uninitializedField'))
                                            list(SEXPgeneratorFun)
                                          else
                                            list(SEXPgeneratorFun, SEXPfinalizerFun)
                                        }
                                               else callSuper()
                                        c(ans, ptrCastFun)
                                      },
                                      buildSEXPgenerator = function(finalizer = NULL) { ## build a function that will provide a new object and return an external pointer
                                          ## This differs from general cppClassDef case
                                          ## because we will return a pointer to a smartPtr to a nimbleList object
                                          ## the pointer to the smartPtr will be case to base class nimSmartPtrBase, which has a virtual casting
                                          ## member function
                                          ##
                                          ## and also we'll return a ptr to a the realPtr (a double pointer to the object)
                                          
                                          extPtrTypes <<- c('ptrToSmartPtr','ptrToPtrToObject')
                                          CBobjectDefs <- list(cppVarFull(name = 'ptrToSmartPtr', baseType = 'nimSmartPtrBase', ptr = 1),
                                                               cppVarFull(name = 'smartPtr', baseType = 'nimSmartPtr', templateArgs = name),
                                                               cppVarFull(name = 'ptrToPtr', baseType = 'void', ptr = 1),
                                                               SptrToSmartPtr = cppSEXP(name = 'SptrToSmartPtr'),
                                                               SptrToPtr = cppSEXP(name = 'SptrToPtr'),
                                                               Sans = cppSEXP(name = 'Sans'))
                                          newCodeLine <- cppLiteral(c(paste0('smartPtr = new ', name,';'),
                                                                      'ptrToSmartPtr = dynamic_cast<nimSmartPtrBase*>(&smartPtr);',
                                                                      'ptrToPtr = smartPtr.getVoidPtrToRealPtr();',
                                                                      '(smartPtr)->NO_hw();',
                                                                      paste0('reinterpret_cast<',name,'*>(*static_cast<void**>(ptrToPtr))->NO_hw();'),
                                                                      'PROTECT(SptrToSmartPtr = R_MakeExternalPtr(ptrToSmartPtr, R_NilValue, R_NilValue));',
                                                                      'PROTECT(SptrToPtr = R_MakeExternalPtr(ptrToPtr, R_NilValue, R_NilValue));'))
                                          allocVectorLine <- cppLiteral(paste0('PROTECT(Sans = allocVector(VECSXP,', 2, '));'))
                                          
                                          packListLines <- cppLiteral(c('SET_VECTOR_ELT(Sans,0,SptrToSmartPtr);',
                                                                        'SET_VECTOR_ELT(Sans,1,SptrToPtr);'
                                                                        ))
                                          
                                          notificationLine <- if(nimbleOptions()$messagesWhenBuildingOrFinalizingCppObjects)
                                                                  paste0('std::cout<< \"In generator for ', name, '. Created at pointer \" << newObj << \"\\n\";')
                                                              else character(0)
                                          codeLines <- substitute({
                                              ## Finalizer registration now happens through nimble's finalizer mapping system.
                                              UNPROTECT(3)
                                              return(Sans)
                                          }, list())
                                          allCodeList <- list(newCodeLine, cppLiteral(notificationLine), allocVectorLine, packListLines, codeLines)
                                          allCode <- putCodeLinesInBrackets(allCodeList)
                                          SEXPgeneratorFun <<- cppFunctionDef(name = paste0('new_',name),
                                                                              args = list(),
                                                                              code = cppCodeBlock(code = allCode, objectDefs = CBobjectDefs, skipBrackets = TRUE),
                                                                              returnType = cppSEXP(),
                                                                              externC = TRUE)
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
                                           Rgenerator <<- buildNimbleListInterface(paste0(name,'_refClass') , .self, sym, where = where)
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
                                          buildCastPtrToNamedObjectsPtrFun()
                                          if(eigenList){
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
                                      buildCastPtrToNamedObjectsPtrFun = function() {
                                          ## SEXP castPointer_XYZ(SEXP input) {
                                          ##   return( R_MakeExternalPtr(dynamic_cast<NamedObjects*>(static_cast<CLASSNAME*>(R_ExternalPtrAddr(input)))) );
                                          ##}
                                          newName <- paste0(name, "_castPtrPtrToNamedObjectsPtrSEXP")
                                          args = list(cppSEXP(name = 'input'))
                                          allCodeList <- list(cppLiteral(paste0("return( R_MakeExternalPtr(dynamic_cast<NamedObjects*>(reinterpret_cast<",name,"*>(*static_cast<void**>(R_ExternalPtrAddr(input)))), R_NilValue, R_NilValue));")))
                                          allCode <- putCodeLinesInBrackets(allCodeList)
                                          ptrCastFun <<- cppFunctionDef(name = newName,
                                                                        args = args,
                                                                        code = cppCodeBlock(code = allCode, objectDefs = list(), skipBrackets = TRUE),
                                                                        returnType = cppSEXP(),
                                                                        externC = TRUE)
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

cppNimbleListClass <- setRefClass('cppNimbleListClass',
                                  contains = 'cppNimbleClassClass',
                                  fields = list(
                                      ptrCastFun = 'ANY',
                                      ptrCastToPtrPairFun = 'ANY',
                                      predefined = 'ANY'
                                  ),
                                  methods = list(
                                      initialize = function(nimCompProc, debugCpp = FALSE, fromModel = FALSE, ...) {
                                        callSuper(nimCompProc, debugCpp, fromModel, ...)
                                        predefined <<- nimCompProc$nimbleListObj$predefined
                                        inheritance <<- c(inheritance, 'pointedToBase')
                                      },
                                      getDefs = function() {
                                        if(predefined) {
                                            ## Prevent redefinition of nimbleList classes that live
                                            ## in predefinedNimbleLists.h and predefinedNumbleLists.cpp.
                                            return(list())
                                        } else if(name == 'EIGEN_SVDCLASS' || name == 'EIGEN_EIGENCLASS') {
                                            ## Handle these two EIGEN_ classes specially. They are implemented as an inheritance hierarchy,
                                            ## so we cannot generate the class definitions (Nick and Perry understand why).
                                            ## Note that this section of code is only ever run under the control of generateStaticCode.R to
                                            ## regenerate C++ code in predefinedNimbleLists.h and predefinedNimbleLists.cpp.
                                            defs <- c(SEXPgeneratorFun, ptrCastFun, ptrCastToPtrPairFun)
                                            if(!inherits(SEXPfinalizerFun, 'uninitializedField')) {
                                                defs <- c(defs, SEXPfinalizerFun)
                                            }
                                            return(defs)
                                        } else {
                                            return(c(callSuper(), ptrCastFun, ptrCastToPtrPairFun))
                                        }
                                      },
                                      buildCastPtrToNamedObjectsPtrFun = function() {
                                          ## Creating code like:
                                          ## SEXP  nfRefClass_84_castPtrPtrToNamedObjectsPtrSEXP ( SEXP input )  {
                                          ##  return( R_MakeExternalPtr(dynamic_cast<NamedObjects*>(reinterpret_cast<nfRefClass_84*>(*static_cast<void**>(R_ExternalPtrAddr(input)))), R_NilValue, R_NilValue));
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
                                      buildCastPtrToPtrPairFun = function() {
## example of code being generated:
## SEXP  nfRefClass_84_castDerivedPtrPtrToPairOfPtrsSEXP ( SEXP input )  {
##   nimSmartPtrBase * ptrToSmartPtrBase;
##   nimSmartPtr<nfRefClass_84> * ptrToSmartPtr;
##   void * ptrToPtr;
##   SEXP SptrToSmartPtrBase;
##   SEXP SptrToPtr;
##   SEXP Sans;

##   ptrToSmartPtr = static_cast<nimSmartPtr<nfRefClass_84> *>(R_ExternalPtrAddr(input));
##   ptrToSmartPtrBase = dynamic_cast<nimSmartPtrBase*>(ptrToSmartPtr);
##   ptrToPtr = ptrToSmartPtr->getVoidPtrToRealPtr();
##   reinterpret_cast<nfRefClass_84*>(*static_cast<void**>(ptrToPtr))->NO_hw();
 
##   PROTECT(SptrToSmartPtrBase = R_MakeExternalPtr(ptrToSmartPtrBase, R_NilValue, R_NilValue));
##   PROTECT(SptrToPtr = R_MakeExternalPtr(ptrToPtr, R_NilValue, R_NilValue));
##   PROTECT(Sans = Rf_allocVector(VECSXP,2));
##   SET_VECTOR_ELT(Sans,0,SptrToSmartPtrBase);
##   SET_VECTOR_ELT(Sans,1,SptrToPtr);
##   UNPROTECT(3);
##   return(Sans);
## }
                                          newName <- paste0(name, "_castDerivedPtrPtrToPairOfPtrsSEXP")

                                          args = list(input = cppSEXP(name = 'input'))
                                          CBobjectDefs <- list(cppVarFull(name = 'ptrToSmartPtrBase', baseType = 'nimSmartPtrBase', ptr = 1),
                                                               cppVarFull(name = 'ptrToSmartPtr', baseType = 'nimSmartPtr', templateArgs = name, ptr = 1),
                                                               cppVarFull(name = 'ptrToPtr', baseType = 'void', ptr = 1),
                                                               SptrToSmartPtrBase = cppSEXP(name = 'SptrToSmartPtrBase'),
                                                               SptrToPtr = cppSEXP(name = 'SptrToPtr'),
                                                               Sans = cppSEXP(name = 'Sans'))
                                          newCodeLine <- cppLiteral(c(paste0('ptrToSmartPtr = static_cast<nimSmartPtr<',name,'>* >(R_ExternalPtrAddr(input));'),
                                                                      'ptrToSmartPtrBase = dynamic_cast<nimSmartPtrBase*>(ptrToSmartPtr);',
                                                                      'ptrToPtr = ptrToSmartPtr->getVoidPtrToRealPtr();',
                                                                      'PROTECT(SptrToSmartPtrBase = R_MakeExternalPtr(ptrToSmartPtrBase, R_NilValue, R_NilValue));',
                                                                      'PROTECT(SptrToPtr = R_MakeExternalPtr(ptrToPtr, R_NilValue, R_NilValue));'))
                                          allocVectorLine <- cppLiteral(paste0('PROTECT(Sans = Rf_allocVector(VECSXP,', 2, '));'))
                                          
                                          packListLines <- cppLiteral(c('SET_VECTOR_ELT(Sans,0,SptrToSmartPtrBase);',
                                                                        'SET_VECTOR_ELT(Sans,1,SptrToPtr);'
                                                                        ))
                                          
                                          codeLines <- substitute({
                                              ## Finalizer registration now happens through nimble's finalizer mapping system.
                                              UNPROTECT(3)
                                              return(Sans)
                                          }, list())
                                          allCodeList <- list(newCodeLine, allocVectorLine, packListLines, codeLines)
                                          allCode <- putCodeLinesInBrackets(allCodeList)
                                          ptrCastToPtrPairFun <<- cppFunctionDef(name = newName,
                                                                              args = args,
                                                                              code = cppCodeBlock(code = allCode, objectDefs = CBobjectDefs, skipBrackets = TRUE),
                                                                              returnType = cppSEXP(),
                                                                              externC = TRUE)                                          
                                      },
                                      buildSEXPgenerator = function(finalizer = NULL) {
                                          ## Build a function that will provide a new object and return an external pointer.
                                          ##
                                          ## This differs from general cppClassDef case
                                          ## because we will return a pointer to a smartPtr to a nimbleList object.
                                          ## The pointer to the smartPtr will be case to base class nimSmartPtrBase, which has a virtual casting
                                          ## member function.
                                          ##
                                          ## We'll also return a ptr to a the realPtr (a double pointer to the object).
                                          ##
                                          ## Example output:
                                          ##
                                          ## SEXP  new_nfRefClass_84 (  )  {
                                          ##   nimSmartPtr<nfRefClass_84> * ptrToSmartPtr;
                                          ##   nfRefClass_84 * newObj;
                                          ##   SEXP SptrToSmartPtr;
                                          ##
                                          ##   newObj = new  nfRefClass_84 ;
                                          ##   ptrToSmartPtr = new nimSmartPtr<nfRefClass_84>;
                                          ##   ptrToSmartPtr->setPtrFromT(newObj);
                                          ##   (*ptrToSmartPtr)->NO_hw();
                                          ##
                                          ##   PROTECT(SptrToSmartPtr = R_MakeExternalPtr(ptrToSmartPtr, R_NilValue, R_NilValue));
                                          ##   UNPROTECT(1);
                                          ##   return(nfRefClass_84_castDerivedPtrPtrToPairOfPtrsSEXP(SptrToSmartPtr));
                                          ## }

                                          extPtrTypes <<- c('ptrToSmartPtr','ptrToPtrToObject')

                                          CBobjectDefs <- list(cppVarFull(name = 'ptrToSmartPtr', baseType = 'nimSmartPtr', templateArgs = name, ptr = 1),
                                                               cppVarFull(name = 'newObj', baseType = name, ptr = 1),
                                                               SptrToSmartPtr = cppSEXP(name = 'SptrToSmartPtr'),
                                                               Sans = cppSEXP(name = "Sans"))
                                          newCodeLine <- cppLiteral(c(paste('newObj = new ',name,';'),
                                                                      paste0('ptrToSmartPtr = new nimSmartPtr<',name,'>;'),
                                                                      'ptrToSmartPtr->setPtrFromT(newObj);',
                                                                      'PROTECT(SptrToSmartPtr = R_MakeExternalPtr(ptrToSmartPtr, R_NilValue, R_NilValue));'))
                                        # String pasting required to get the format PROTECT(Y = f(X)) instead of Y = PROTECT(f(X)).
                                        # It might not matter, but we are being careful to make the PROTECT's look just like in Tomas Kalibera's blog post (and CRAN's code analysis)
                                          codeLines <- cppLiteral(c(paste0('PROTECT(Sans = ', ptrCastToPtrPairFun$name, '(SptrToSmartPtr));'),
                                                                    'UNPROTECT(2);',
                                                                    'return(Sans);'))
                                          ## codeLines <- substitute({
                                          ##     ## Finalizer registration now happens through nimble's finalizer mapping system.
                                          ##     Sans = PROTECT(ptrCastToPtrPairFun(SptrToSmartPtr))
                                          ##     UNPROTECT(2)
                                          ##     return(Sans)
                                          ## }, list(ptrCastToPtrPairFun = as.name(ptrCastToPtrPairFun$name)))
                                          allCodeList <- list(newCodeLine, codeLines)
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
                                          CmultiInterface <<- CmultiNimbleListClass(compiledNodeFun = .self, basePtrCall = sym, project = nimbleProject)
                                      },
                                      buildRgenerator = function(where = globalenv(), dll = NULL) {
                                          sym <- if(!is.null(dll))
                                                     getNativeSymbolInfo(SEXPgeneratorFun$name, dll)
                                                 else
                                                     SEXPgeneratorFun$name
                                           Rgenerator <<- buildNimbleListInterface(paste0(name,'_refClass'), .self, sym, where = where)
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
                                          buildCastPtrToPtrPairFun()
                                          if(!predefined){
                                              buildCopyFromSexp()
                                              buildCopyToSexp()
                                              buildCreateNewSexp()
                                              buildResetFlags()
                                              addCopyFromRobject()
                                          }
                                          callSuper(where)
                                      },
                                      buildCreateNewSexp = function(){
                                        interfaceArgs <- symbolTable()
                                        newListLine <- list()
                                        listElementTable <- symbolTable()
                                        listElementTable$addSymbol(cppSEXP(name = "S_newNimList"))
                                        listElementTable$addSymbol(cppSEXP(name = "S_listName"))
                                      
                                        newListLine[[1]] <- substitute({PROTECT(S_listName <- Rf_allocVector(STRSXP, 1));
                                          SET_STRING_ELT(S_listName, 0, PROTECT(Rf_mkChar(LISTNAME)));},
                                          list(LISTNAME = nimCompProc$nimbleListObj$className))
                                        newListLine[[2]] <- substitute(PROTECT(S_newNimList <- makeNewNimbleList(S_listName)),
                                                                           list())
                                        newListLine[[3]] <- quote(cppLiteral('R_PreserveObject(RObjectPointer = S_newNimList);'))
                                        newListLine[[4]] <-   substitute(UNPROTECT(2+1), list())
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
                                        envLine <- substitute({PROTECT(ENVNAME <- Rf_allocVector(STRSXP, 1));
                                          SET_STRING_ELT(ENVNAME, 0, PROTECT(Rf_mkChar(".xData")));},
                                          list(ENVNAME = as.name(environmentCPPName)))
                                        
                                        for(i in seq_along(elementNames)){
                                          Snames[i] <- Rname2CppName(paste0('S_', elementNames[i]))
                                          listElementTable$addSymbol(cppSEXP(name = Snames[i]))
                                          elementSymTab <- nimCompProc$symTab$getSymbolObject(elementNames[i])
                                          conditionalLineList <- c(conditionalLineList, generateConditionalLines(nimCompProc$symTab$getSymbolObject(elementNames[i]),
                                                                                                                 listElementTable$getSymbolObject(Snames[i])))
                                          
                                          copyToListLines[[i]] <- substitute(Rf_defineVar(Rf_install(ELEMENTNAME), VALUE, PROTECT(GET_SLOT(ROBJ, XDATA))),
                                                                             list(ELEMENTNAME = elementNames[i], VALUE = as.name(Snames[i]),
                                                                                  ROBJ = as.name('RObjectPointer'),
                                                                                  XDATA = as.name(environmentCPPName)))
                                        }
                                        
                                        setFlagLine <- list(substitute(RCopiedFlag <- true,
                                                                       list()))
                                        returnLine <- list(substitute(return(ROBJ),
                                                           list(ROBJ = as.name('RObjectPointer'))))
                                        unprotectLine <- list(substitute(UNPROTECT(N), list(N = 2 * numElements + 1 + 1)))
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
                                        storeSexpLine <- list(quote(cppLiteral('R_PreserveObject(RObjectPointer =  S_nimList_);')))

                                        environmentCPPName <- Rname2CppName('S_.xData')  ## create SEXP for ref class environment 
                                        listElementTable$addSymbol(cppSEXP(name = environmentCPPName))
                                        envLine <- substitute({PROTECT(ENVNAME <- Rf_allocVector(STRSXP, 1));
                                          SET_STRING_ELT(ENVNAME, 0, PROTECT(Rf_mkChar(".xData")));},
                                          list(ENVNAME = as.name(environmentCPPName)))

                                        for(i in seq_along(argNames)) {
                                          Snames[i] <- Rname2CppName(paste0('S_', argNames[i]))
                                          listElementTable$addSymbol(cppSEXP(name = Snames[i]))
                                            copyFromListLines[[i]] <- substitute(PROTECT(SVAR <- Rf_findVarInFrame(PROTECT(GET_SLOT(S_nimList_, XDATA)), Rf_install(ARGNAME))),
                                                                               list(ARGNAME = argNames[i], 
                                                                                    SVAR = as.name(Snames[i]),
                                                                                    XDATA = as.name(environmentCPPName)))
                                            copyLine <-  buildCopyLineFromSEXP(listElementTable$getSymbolObject(Snames[i]),
                                                                               nimCompProc$symTab$getSymbolObject(argNames[i]))
                                            copyLines <- c(copyLines, copyLine)
                                        }
                                        numArgs <- length(argNames)
                                        unprotectLine <- substitute(UNPROTECT(N), list(N = 2 * numArgs + 1 + 1))
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
                                      },
                                      buildResetFlags = function(){
                                        interfaceArgs <- symbolTable()
                                        elementNames <- nimCompProc$symTab$getSymbolNames()
                                        resetNestedFlagLines <- list()
                                        listElementTable <- symbolTable()
                                        
                                        resetFlagLine <- list(substitute(RCopiedFlag <- false,
                                                                         list()))

                                        for(i in seq_along(elementNames)){
                                          elementSymTab <- nimCompProc$symTab$getSymbolObject(elementNames[i])
                                          isList <- inherits(elementSymTab, 'symbolNimbleList')
                                          if(isList){
                                            resetRCopiedFlag  <- paste0(elementSymTab$name,"->resetFlags();")
                                            resetRCopiedFlagLine <- substitute(cppLiteral(resetText), list(resetText = resetRCopiedFlag))
                                            resetNestedFlagLines <- c(resetNestedFlagLines,
                                                                 resetRCopiedFlagLine)
                                          }
                                        }
                                        allCode <- embedListInRbracket(c(resetFlagLine, resetNestedFlagLines))
                                        functionDefs[[paste0(name, "_resetFlags")]] <<- cppFunctionDef(name = "resetFlags",
                                                                                                   args = interfaceArgs,
                                                                                                   code = cppCodeBlock(code = RparseTree2ExprClasses(allCode), objectDefs = listElementTable),
                                                                                                   returnType = cppVoid(),
                                                                                                   externC = FALSE,
                                                                                                   CPPincludes = list(nimbleIncludeFile("RcppUtils.h"),
                                                                                                                      nimbleIncludeFile("smartPtrs.h")))
                                      }

                                      )
                                  )

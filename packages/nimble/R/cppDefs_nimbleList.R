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
                                        buildSEXPCopier()
                                        callSuper(where)
                                      },
                                      buildSEXPCopier = function(){
                                        browser()
                                        objects <- symbolTable2cppVars(nimCompProc$symTab)
                                        argNames <- nimCompProc$symTab$getSymbolNames()
                                        copyLines <- list()
                                        Snames <- character(length(argNames))
                                        returnType <- "void"
                                        interfaceArgs <- symbolTable()
                                        objects$setParentST(interfaceArgs)
                      
                                        for(i in seq_along(argNames)) {
                                          Snames[i] <- Rname2CppName(paste0('S_', argNames[i]))
                                          interfaceArgs$addSymbol(cppSEXP(name = Snames[i]))
                                          copyLines[[i]] <- buildCopyLineFromSEXP(interfaceArgs$getSymbolObject(Snames[i]),
                                                                                  nimCompProc$symTab$getSymbolObject(argNames[i]))
                                        }
                                        
                                        numArgs <- length(argNames)
                                            if(numArgs> 0) {
                                              objects$addSymbol(cppSEXP(name = 'S_returnValue_LIST_1234'))
                                              returnListLines <- returnCopyLines <- vector('list', length = numArgs)
                                              allocVectorLine <- substitute(PROTECT(S_returnValue_LIST_1234 <- allocVector(VECSXP, nAp1)), list(nAp1 = numArgs))
                                                for(i in 1:numArgs) {
                                                  returnCopyLines[[i]] <- buildCopyLineToSEXP(nimCompProc$symTab$getSymbolObject(argNames[i]),
                                                                                                                       interfaceArgs$getSymbolObject(Snames[i]))
                                                  returnListLines[[i]] <- substitute(SET_VECTOR_ELT(S_returnValue_LIST_1234, Im1, THISSEXP),
                                                                                     list(Im1 = i-1, THISSEXP = as.name(Snames[i])))
                                                }
                                              }
                                              unprotectLine <- substitute(UNPROTECT(N), list(N = numArgs + 1))
                                              allCode <- embedListInRbracket(c(copyLines, list(allocVectorLine),
                                                                               returnCopyLines, returnListLines, list(unprotectLine)))
                                       
                                        # SEXPinterfaceCname <<- paste0('CALL_',Rname2CppName(paste0(if(!is.null(className)) paste0(className,'_') else NULL, name))) ##Rname2CppName needed for operator()
                                              
                                        #may not want [[name]] below but instead something else ie 'copier', 'listCopier', etc.
                                        functionDefs[[name]] <<- cppFunctionDef(name = "copyFromSEXP",
                                                                               args = interfaceArgs,
                                                                               code = cppCodeBlock(code = RparseTree2ExprClasses(allCode), objectDefs = objects),
                                                                               returnType = cppVoid(),
                                                                               externC = FALSE)
                                      }
                                      )
                                  )



# 
# buildSEXPinterfaceFun = function(className = NULL) {
#   asMember <- !is.null(className)
#   objects <- symbolTable2cppVars(RCfunProc$compileInfo$origLocalSymTab)
#   argNames <- RCfunProc$compileInfo$origLocalSymTab$getSymbolNames()
#   browser()
#   argNames <- c()
#   copyLines <- list()
#   Snames <- character(length(argNames))
#   interfaceArgs <- symbolTable()
#   objects$setParentST(interfaceArgs)
#   returnVoid <- returnType$baseType == 'void'
#   for(i in seq_along(argNames)) {
#     Snames[i] <- Rname2CppName(paste0('S_', argNames[i]))
#     ## For each argument to the RCfunction we need a corresponding SEXP argument to the interface function
#     interfaceArgs$addSymbol(cppSEXP(name = Snames[i]))
#     
#     ## and we need a line to copy from the SEXP to the local variable
#     ## The to argument uses the origLocalSymbolObject rather than the objects (which has cppVars) because that has the nDim
#     ## The name of that and the new one in objects must match
#     copyLines[[i]] <- buildCopyLineFromSEXP(interfaceArgs$getSymbolObject(Snames[i]),
#                                             RCfunProc$compileInfo$origLocalSymTab$getSymbolObject(argNames[i]))
#   }
#   
#   RHScall <- as.call(c(list(as.name(name)),
#                        lapply(argNames, as.name)))
#   if(asMember) {
#     ## Add a final argument for the extptr
#     interfaceArgs$addSymbol(cppSEXP(name = 'SextPtrToObject'))
#     ## And make the RHScall
#     RHScall <- substitute(cppMemberDereference(
#       template(static_cast, cppPtrType(CN))(R_ExternalPtrAddr(SextPtrToObject)), RHS),
#       list(CN = as.name(className), RHS = RHScall))
#   }
#   
#   if(returnVoid) {
#     fullCall <- RHScall
#   } else {
#     objects$addSymbol(cppSEXP(name = 'S_returnValue_1234')) ## Object for the return statement: "return(S_returnValue_1234)"
#     LHSvar <- RCfunProc$compileInfo$returnSymbol$genCppVar()
#     LHSvar$name <- "LHSvar_1234"
#     objects$addSymbol(LHSvar)
#     fullCall <- substitute(LHS <- RHS, list(LHS = as.name(LHSvar$name), RHS = RHScall))
#   }
#   ## Put GetRNGstate() and PutRNGstate() around the call.    
#   fullCall <- substitute({GetRNGstate(); FULLCALL; PutRNGstate()}, list(FULLCALL = fullCall))
#   
#   returnAllArgs <- TRUE
#   ## Pack up all inputs and the return value in a list.
#   if(returnAllArgs) {
#     numArgs <- length(argNames)
#     if(numArgs + !returnVoid > 0) {
#       objects$addSymbol(cppSEXP(name = 'S_returnValue_LIST_1234'))
#       returnListLines <- returnCopyLines <- vector('list', length = numArgs+!returnVoid)
#       allocVectorLine <- substitute(PROTECT(S_returnValue_LIST_1234 <- allocVector(VECSXP, nAp1)), list(nAp1 = numArgs + !returnVoid))
#       if(numArgs > 0) {
#         for(i in 1:numArgs) {
#           returnCopyLines[[i]] <- buildCopyLineToSEXP(RCfunProc$compileInfo$origLocalSymTab$getSymbolObject(argNames[i]),
#                                                       interfaceArgs$getSymbolObject(Snames[i]))
#           returnListLines[[i]] <- substitute(SET_VECTOR_ELT(S_returnValue_LIST_1234, Im1, THISSEXP),
#                                              list(Im1 = i-1, THISSEXP = as.name(Snames[i])))
#         }
#       }
#       if(!returnVoid) {
#         rsName <- RCfunProc$compileInfo$returnSymbol$name
#         RCfunProc$compileInfo$returnSymbol$name <<- LHSvar$name
#         returnCopyLines[[numArgs+1]] <- buildCopyLineToSEXP(RCfunProc$compileInfo$returnSymbol,
#                                                             objects$getSymbolObject('S_returnValue_1234'))
#         RCfunProc$compileInfo$returnSymbol$name <<- rsName
#         returnListLines[[numArgs+1]] <- substitute(SET_VECTOR_ELT(S_returnValue_LIST_1234, I, THISSEXP),
#                                                    list(I = numArgs, THISSEXP = as.name('S_returnValue_1234')))
#       }
#       returnLine <- quote(return(S_returnValue_LIST_1234))
#       unprotectLine <- substitute(UNPROTECT(N), list(N = numArgs + 1 + !returnVoid))
#       allCode <- embedListInRbracket(c(copyLines, list(fullCall), list(allocVectorLine),
#                                        returnCopyLines, returnListLines, list(unprotectLine), list(returnLine)))
#     } else { ## No input or return objects
#       returnLine <- quote(return(R_NilValue))
#       allCode <- embedListInRbracket(c(copyLines, list(fullCall),
#                                        list(returnLine)))
#     }
#   } else {
#     writeLines("Haven't written the single return case yet")        
#   }
#   SEXPinterfaceCname <<- paste0('CALL_',Rname2CppName(paste0(if(!is.null(className)) paste0(className,'_') else NULL, name))) ##Rname2CppName needed for operator()
#   SEXPinterfaceFun <<- cppFunctionDef(name = SEXPinterfaceCname,
#                                       args = interfaceArgs,
#                                       code = cppCodeBlock(code = RparseTree2ExprClasses(allCode), objectDefs = objects),
#                                       returnType = cppSEXP(),
#                                       externC = TRUE,
#                                       CPPincludes = list(nimbleIncludeFile("RcppUtils.h")))
#   invisible(NULL)
# }

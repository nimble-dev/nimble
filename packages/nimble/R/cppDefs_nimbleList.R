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
                                      genNeededTypes = function(debugCpp = FALSE, fromModel = FALSE){
                                        CPPincludes <<- c(CPPincludes, nimbleIncludeFile("smartPtrs.h"))
                                        callSuper(debugCpp, fromModel)
                                      },
                                      buildAll = function(where = where) {
                                        # buildSEXPCopier()
                                        callSuper(where)
                                      },
                                      buildSEXPCopier = function(){
                                        browser()
                                        # argNames <- RCfunProc$compileInfo$origLocalSymTab$getSymbolNames()
                                        # 
                                        # inputArgs <- list(cppSEXP(name = 'Sinput'))
                                        # code <- putCodeLinesInBrackets(list(cppLiteral(c(notificationLine, castLine, deleteLine))))
                                        # 
                                        # 
                                        # #may not want [[name]] below but instead something else ie 'copier', 'listCopier', etc.
                                        # functionDefs[[name]] <- cppFunctionDef(name = paste0(name,'_Copier'),
                                        #                                        args = interfaceArgs,
                                        #                                        code = cppCodeBlock(code = RparseTree2ExprClasses(allCode), objectDefs = objects),
                                        #                                        returnType = cppSEXP(),
                                        #                                        externC = TRUE)
                                      }
                                      )
                                  )

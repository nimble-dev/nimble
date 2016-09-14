cppNimbleListClass <- setRefClass('cppNimbleListClass',
                                  contains = 'cppNimbleClassClass',
                                  fields = list(
                                  ),
                                  methods = list(
                                      initialize = function(nimCompProc, debugCpp = FALSE, fromModel = FALSE, ...) {
                                          callSuper(nimCompProc, debugCpp, fromModel, ...)
                                          # inheritance <<- c(inheritance, 'pointedToBase')
                                      },
                                      buildCmultiInterface = function(dll = NULL) {
                                          sym <- if(!is.null(dll))
                                                      getNativeSymbolInfo(SEXPgeneratorFun$name, dll)
                                                  else
                                                      SEXPgeneratorFun$name
                                          message('CmultiInterface for nimbleList does not exist')
                                          ##CmultiInterface <<- CmultiNimbleFunctionClass(compiledNodeFun = .self, basePtrCall = sym, project = nimbleProject)
                                      },
                                      buildRgenerator = function(where = globalenv(), dll = NULL) {
                                          sym <- if(!is.null(dll))
                                                     getNativeSymbolInfo(SEXPgeneratorFun$name, dll)
                                                 else
                                                     SEXPgeneratorFun$name
                                          Rgenerator <<- buildNimbleFxnInterface(paste0(name,'_refClass') , .self, sym, where = where)
                  
                                          # message('Rgenerator for nimbleList does not exist')
                                      })
                                  )

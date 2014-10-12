## Small class for information on compilation of each nf method
RCfunctionCompileClass <- setRefClass('RCfunctionCompileClass',
                                      fields = list(
                                          origRcode = 'ANY',
                                          origLocalSymTab = 'ANY',
                                          nimExpr = 'ANY',
                                          newLocalSymTab = 'ANY',
                                          returnSymbol = 'ANY',
                                          newRcode = 'ANY',
                                          typeEnv = 'ANY'	#environment
                                          ),
                                          methods = list(initialize <- function(...){typeEnv <<- new.env(); callSuper(...)}
                                          ))

RCvirtualFunProcessing <- setRefClass('RCvirtualFunProcessing',
                                      fields = list(
                                          name = 'ANY',		#character
                                          RCfun = 'ANY', ##nfMethodRC
                                          compileInfo = 'ANY' ## RCfunctionCompileClass``
                                          ),
                                      methods = list(
                                          initialize = function(f = NULL, funName) {
                                              if(!is.null(f)) {
                                                  if(missing(funName)) {
                                                      sf <- substitute(f)
                                                      name <<- Rname2CppName(deparse(sf))
                                                  } else {
                                                      name <<- funName
                                                  }
                                                  if(is.function(f)) {
                                                      RCfun <<- nfMethodRC$new(f)
                                                  } else {
                                                      if(!inherits(f, 'nfMethodRC')) {
                                                          stop('Error: f must be a function or an RCfunClass')
                                                      }
                                                      RCfun <<- f
                                                  }
                                                  compileInfo <<- RCfunctionCompileClass$new(origRcode = RCfun$code, newRcode = RCfun$code)
                                              }
                                          },
                                          showCpp = function() {
                                              writeCode(nimGenerateCpp(compileInfo$nimExpr, compileInfo$newLocalSymTab))
                                          },
                                          setupSymbolTables = function(parentST = NULL) {
                                              compileInfo$origLocalSymTab <<- argTypeList2symbolTable(RCfun$argInfo) ## will be used for function args.  must be a better way.
                                              compileInfo$newLocalSymTab <<- argTypeList2symbolTable(RCfun$argInfo)
                                              if(!is.null(parentST)) {
                                                  compileInfo$origLocalSymTab$setParentST(parentST)
                                                  compileInfo$newLocalSymTab$setParentST(parentST)
                                              }
                                              compileInfo$returnSymbol <<- argType2symbol(RCfun$returnType, "returnValues")
                                          },
                                          process = function(...) {
                                              if(inherits(compileInfo$origLocalSymTab, 'uninitializedField')) {
                                                  setupSymbolTables()
                                              }
                                          }))

RCfunction <- function(f, name = NA, returnCallable = TRUE) {
    if(is.na(name)) name <- rcFunLabelMaker()
    nfm <- nfMethodRC$new(f, name)
    if(returnCallable) nfm$generateFunctionObject(keep.nfMethodRC = TRUE) else nfm
}

is.rcf <- function(x) {
    if(inherits(x, 'nfMethodRC')) return(TRUE)
    if(is.function(x)) {
        if(is.null(environment(x))) return(FALSE)
        if(exists('nfMethodRCobject', envir = environment(x), inherits = FALSE)) return(TRUE)
    }
    FALSE
}

rcFunLabelMaker <- labelFunctionCreator('rcFun')

RCfunProcessing <- setRefClass('RCfunProcessing',
                               contains = 'RCvirtualFunProcessing',
                               fields = list(
                                   neededRCfuns = 'list' ## nfMethodRC objects
                                   ),
                               methods = list(
                                   process = function(debug = FALSE, debugCpp = FALSE, debugCppLabel = character()) {
                                       
                                       if(!is.null(nimbleOptions$debugRCfunProcessing)) {
                                           if(nimbleOptions$debugRCfunProcessing) {
                                               debug <- TRUE
                                               writeLines('Debugging RCfunProcessing (nimbleOptions$debugRCfunProcessing is set to TRUE)') 
                                           }
                                       }
                                   	
                                       if(inherits(compileInfo$origLocalSymTab, 'uninitializedField')) {
                                           setupSymbolTables()
                                       }
                                       
                                       if(debug) {
                                           writeLines('**** READY FOR makeExprClassObjects *****')
                                           browser()
                                       }
                                                                              
                                       ## set up exprClass object
                                       compileInfo$nimExpr <<- RparseTree2ExprClasses(compileInfo$newRcode)
                                       
                                       if(debug) {
                                           print('nimDeparse(compileInfo$nimExpr)')
                                           writeCode(nimDeparse(compileInfo$nimExpr))
                                           writeLines('***** READY FOR processSpecificCalls *****')
                                           browser()
                                       }

                                       exprClasses_processSpecificCalls(compileInfo$nimExpr, compileInfo$newLocalSymTab)

                                       if(debug) {
                                           print('nimDeparse(compileInfo$nimExpr)')
                                           writeCode(nimDeparse(compileInfo$nimExpr))
                                           writeLines('***** READY FOR buildInterms *****')
                                           browser()
                                       }
                                       
                                       ## build intermediate variables
                                       exprClasses_buildInterms(compileInfo$nimExpr)
                                       if(debug) {
                                           print('nimDeparse(compileInfo$nimExpr)')
                                           writeCode(nimDeparse(compileInfo$nimExpr))
                                           print('compileInfo$newLocalSymTab')
                                           print(compileInfo$newLocalSymTab)
                                           writeLines('***** READY FOR initSizes *****')
                                           browser()
                                       }
                                       
                                       compileInfo$typeEnv <<- exprClasses_initSizes(compileInfo$nimExpr, compileInfo$newLocalSymTab)
                                       if(debug) {
                                           print('ls(compileInfo$typeEnv)')
                                           print(ls(compileInfo$typeEnv))
                                           print('lapply(compileInfo$typeEnv, function(x) x$show())')
                                           lapply(compileInfo$typeEnv, function(x) x$show())
                                           writeLines('***** READY FOR setSizes *****')
                                      browser()
                                       }

                                       compileInfo$typeEnv[['neededRCfuns']] <<- list()

                                       passedArgNames <- as.list(names(RCfun$argInfo))
                                       names(passedArgNames) <- names(RCfun$argInfo)
                                       compileInfo$typeEnv[['passedArgumentNames']] <<- passedArgNames ## only the names are used.  

                                       exprClasses_setSizes(compileInfo$nimExpr, compileInfo$newLocalSymTab, compileInfo$typeEnv)
                                       neededRCfuns <<- compileInfo$typeEnv[['neededRCfuns']]
                                       
                                       if(debug) {
                                           print('compileInfo$nimExpr$show(showType = TRUE)')
                                           print(compileInfo$nimExpr$show(showType = TRUE))
                                           print('compileInfo$nimExpr$show(showAssertions = TRUE)')
                                           print(compileInfo$nimExpr$show(showAssertions = TRUE))
                                           writeLines('***** READY FOR insertAssertions *****')
                                           browser()
                                       }
                                       
                                       exprClasses_insertAssertions(compileInfo$nimExpr)
                                       if(debug) {
                                           print('compileInfo$nimExpr$show(showAssertions = TRUE)')
                                           compileInfo$nimExpr$show(showAssertions = TRUE)
                                           print('compileInfo$nimExpr$show(showToEigenize = TRUE)')
                                           compileInfo$nimExpr$show(showToEigenize = TRUE)
                                           print('nimDeparse(compileInfo$nimExpr)')
                                           writeCode(nimDeparse(compileInfo$nimExpr))
                                           writeLines('***** READY FOR labelForEigenization *****')
                                           browser()
                                       }
                                       
                                       exprClasses_labelForEigenization(compileInfo$nimExpr)
                                       
                                       if(debug) {
                                           print('nimDeparse(compileInfo$nimExpr)')
                                           writeCode(nimDeparse(compileInfo$nimExpr))
                                           writeLines('***** READY FOR eigenize *****')
                                           browser()
                                       }

                                       exprClasses_eigenize(compileInfo$nimExpr, compileInfo$newLocalSymTab, compileInfo$typeEnv)
                                       if(debug) {
                                           print('nimDeparse(compileInfo$nimExpr)')
                                           writeCode(nimDeparse(compileInfo$nimExpr))
                                           print('compileInfo$newLocalSymTab')
                                           print(compileInfo$newLocalSymTab)
                                           print('ls(compileInfo$typeEnv)')
                                           print(ls(compileInfo$typeEnv))
                                           
                                           writeLines('***** READY FOR liftMaps*****')
                                           browser()
                                       }

                                       exprClasses_liftMaps(compileInfo$nimExpr, compileInfo$newLocalSymTab, compileInfo$typeEnv)
                                       if(debug) {
                                           print('nimDeparse(compileInfo$nimExpr)')
                                           writeCode(nimDeparse(compileInfo$nimExpr))
                                           writeLines('***** READY FOR cppOutput*****')
                                           browser()
                                       }
                                       
                                       if(debugCpp) {
                                           if(debug) writeLines('*** Inserting debugging')
                                           exprClasses_addDebugMarks(compileInfo$nimExpr, paste(debugCppLabel, name))
                                           if(debug) {
                                               print('nimDeparse(compileInfo$nimExpr)')
                                               writeCode(nimDeparse(compileInfo$nimExpr))
                                               writeLines('***** READY FOR cppOutput*****')
                                               browser()
                                           }
                                       }
                                       if(debug & debugCpp) {
                                           print('writeCode(nimGenerateCpp(compileInfo$nimExpr, newMethods$run$newLocalSymTab))')
                                           writeCode(nimGenerateCpp(compileInfo$nimExpr, compileInfo$newLocalSymTab))
                                       }
                                   }
                                   )
                               )

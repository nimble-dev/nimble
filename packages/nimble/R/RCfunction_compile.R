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
                                          nameSubList = 'ANY',
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
                                              argInfoWithMangledNames <- RCfun$argInfo
                                              numArgs <- length(argInfoWithMangledNames)
                                              if(numArgs>0) names(argInfoWithMangledNames) <- paste0("ARG", 1:numArgs, "_", Rname2CppName(names(argInfoWithMangledNames)),"_")
                                              nameSubList <<- lapply(names(argInfoWithMangledNames), as.name)
                                              names(nameSubList) <<- names(RCfun$argInfo)
                                              compileInfo$origLocalSymTab <<- argTypeList2symbolTable(argInfoWithMangledNames) ## will be used for function args.  must be a better way.
                                              compileInfo$newLocalSymTab <<- argTypeList2symbolTable(argInfoWithMangledNames)
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

nf_substituteExceptFunctionsAndDollarSigns <- function(code, subList) {
    ## this is almost like doing substitute with code and subList, but it doesn't traverse the RHS of a '$' operator
    ## and it doesn't replace and function names
    if(is.character(code)) return(code)
    if(is.numeric(code)) return(code)
    if(is.logical(code)) return(code)
    if(is.name(code)) {
        maybeAns <- subList[[as.character(code)]]
        return( if(is.null(maybeAns)) code else maybeAns )
    }
    if(is.call(code)) {
        if(length(code) == 1) return(code)
        else
            if(is.call(code[[1]])) indexRange <- 1:length(code)
            else
                if(as.character(code[[1]])=='$') indexRange <- 2
                else 
                    indexRange <- 2:length(code)
        for(i in indexRange) code[[i]] <- nf_substituteExceptFunctionsAndDollarSigns(code[[i]], subList)
        return(code)
    }
    if(is.list(code)) { ## Keyword processing stage of compilation may have stuck lists into the argument list of a call (for maps)
        code <- lapply(code, nf_substituteExceptFunctionsAndDollarSigns, subList)
    }
    stop(paste("Error doing replacement for code ", deparse(code)))
}

RCfunProcessing <- setRefClass('RCfunProcessing',
                               contains = 'RCvirtualFunProcessing',
                               fields = list(
                                   neededRCfuns = 'list' ## nfMethodRC objects
                                   ),
                               methods = list(
                                   process = function(debug = FALSE, debugCpp = FALSE, debugCppLabel = character(), doKeywords = TRUE) {
                                       
                                       if(!is.null(nimbleOptions()$debugRCfunProcessing)) {
                                           if(nimbleOptions()$debugRCfunProcessing) {
                                               debug <- TRUE
                                               writeLines('Debugging RCfunProcessing (nimbleOptions()$debugRCfunProcessing is set to TRUE)') 
                                           }
                                       }

                                       if(!is.null(nimbleOptions()$debugCppLineByLine)) {
                                           if(nimbleOptions()$debugCppLineByLine) {
                                               debugCpp <- TRUE
                                               if(length(debugCppLabel) == 0) debugCppLabel <- name
                                           }
                                       }

                                       if(doKeywords) {
                                           matchKeywords()
                                           processKeywords()
                                       }
                                                                              
                                       if(inherits(compileInfo$origLocalSymTab, 'uninitializedField')) {
                                           setupSymbolTables()
                                       }

                                      
                                       if(debug) {
                                           writeLines('**** READY FOR makeExprClassObjects *****')
                                           browser()
                                       }
                                       if(length(nameSubList) > 0)
                                           compileInfo$newRcode <<- nf_substituteExceptFunctionsAndDollarSigns(compileInfo$newRcode, nameSubList)
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

                                       passedArgNames <- as.list(compileInfo$origLocalSymTab$getSymbolNames()) 
                                       names(passedArgNames) <- compileInfo$origLocalSymTab$getSymbolNames() 
                                       compileInfo$typeEnv[['passedArgumentNames']] <<- passedArgNames ## only the names are used.  
  
                                       exprClasses_setSizes(compileInfo$nimExpr, compileInfo$newLocalSymTab, compileInfo$typeEnv)
                                       neededRCfuns <<- compileInfo$typeEnv[['neededRCfuns']]
                                       
                                       if(debug) {
                                           print('compileInfo$nimExpr$show(showType = TRUE) -- broken')
                                          ## print(compileInfo$nimExpr$show(showType = TRUE))
                                           print('compileInfo$nimExpr$show(showAssertions = TRUE) -- possible broken')
                                           ## print(compileInfo$nimExpr$show(showAssertions = TRUE))
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
                                           #print(compileInfo$newLocalSymTab)
                                           #print('ls(compileInfo$typeEnv)')
                                           #print(ls(compileInfo$typeEnv))
                                           
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
                                   },
                                   processKeywords = function(nfProc = NULL) {
                                       compileInfo$newRcode <<- processKeywords_recurse(compileInfo$origRcode, nfProc)
                                   },
                                   processKeywords_recurse = function(code, nfProc = NULL){
                                       cl = length(code)
                                       if(cl == 1) {
                                           if(is.call(code)) {
                                               if(length(code[[1]]) > 1) code[[1]] <- processKeywords_recurse(code[[1]], nfProc)
                                           }
                                           return(code)
                                       }
                                       
                                       if(length(code[[1]]) == 1) {
                                           code <- processKeyword(code, nfProc)
                                       }
    
                                       cl = length(code)
    
                                       if(is.call(code)) {
                                           if(length(code[[1]]) > 1) code[[1]] <- processKeywords_recurse(code[[1]], nfProc)
                                           if(cl >= 2) {
                                               for(i in 2:cl) {
                                                   code[[i]] <- processKeywords_recurse(code[[i]], nfProc)
                                               }
                                           }
                                       }
                                       return(code)
                                   },
                                   matchKeywords = function(nfProc = NULL) {
                                       compileInfo$origRcode <<- matchKeywords_recurse(compileInfo$origRcode, nfProc) ## nfProc needed for member functions of nf objects
                                   },
                                   matchKeywords_recurse = function(code, nfProc = NULL) {
                                       cl = length(code)
                                       if(cl == 1){ ## There are no arguments
                                           if(is.call(code)){  
                                               if(length(code[[1]]) > 1) code[[1]] <- matchKeywords_recurse(code[[1]], nfProc) ## recurse on the "a$b" part of a$b() (or the "a(b)" part of a(b)()), etc
                                           }
                                           return(code)
                                       }
                                       if(length(code[[1]]) == 1) ## a simple call like a(b,c), not a$b(c)
                                           code <- matchKeywordCode(code, nfProc)
                                       
                                       if(is.call(code)) {
                                           if(length(code[[1]]) > 1) {
                                               if(deparse(code[[1]][[1]]) == '$') code <- matchKeywordCodeMemberFun(code, nfProc) ## handle a$b(c) as one unit
                                               else code[[1]] <- matchKeywords_recurse(code[[1]], nfProc) ## handle "a(b)" part of a(b)(c), which is probably *never* triggered
                                           }
                                           if(cl >= 2) { ## recurse through arguments
                                               for(i in 2:cl) {
                                                   code[[i]] <- matchKeywords_recurse(code[[i]], nfProc)
                                               }
                                           }
                                       }
                                       return(code)
                                   }

                                   )
                               )

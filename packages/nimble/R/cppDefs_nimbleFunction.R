##

cppVirtualNimbleFunctionClass <- setRefClass('cppVirtualNimbleFunctionClass',
    contains = 'cppClassDef',
    fields = list(
        nfProc = 'ANY'
    ),
    methods = list(
        initialize = function(nfProc, ...) {
            callSuper(...)
            if(!missing(nfProc)) processNFproc(nfProc)
            useGenerator <<-  FALSE
            baseClassObj <- environment(nfProc$nfGenerator)$contains
            if(is.null(baseClassObj)) {
                addInheritance("NamedObjects")
            } else {
                if(is.character(baseClassObj)) addInheritance(baseClassObj)
                else {
                    virtual <- environment(baseClassObj)$virtual
                    if(isTRUE(virtual))
                        baseClassName <- environment(baseClassObj)$className
                    else
                        baseClassName <- environment(baseClassObj)$CclassName
                    addInheritance(baseClassName)
                }
            }
        },
        processNFproc = function(nfp) {
            nfProc <<- nfp
            assign('cppDef', .self, envir = environment(nfProc$nfGenerator))
            for(i in names(nfp$RCfunProcs)) { ## This is what we should do for cppNimbleFunctions too
                functionDefs[[i]] <<- RCfunctionDef(virtual = TRUE, abstract = TRUE)
                functionDefs[[i]]$buildFunction(nfp$RCfunProcs[[i]])
            }
        }
    )
)

## cppNimbleClassClass defines commonalities between cppNimbleFunctionClass and cppNimbleListClass, both of which are classes in nimble
cppNimbleClassClass <- setRefClass('cppNimbleClassClass',
                                   contains = 'cppNamedObjectsClass',
                                   fields = list(
                                       ## Inherits a functionDefs list for member functions
                                       ## Inherits an objectDefs list for member data
                                       SEXPmemberInterfaceFuns = 'ANY', ## List of SEXP interface functions, one for each member function
                                       nimCompProc = 'ANY', ## nfProcessing or nlProcessing object, needed to get the member data symbol table post-compilation

                                       Rgenerator = 'ANY' , ## function to generate and wrap a new object from an R object
                                       CmultiInterface = 'ANY', ## object for interfacing multiple C instances when a top-level interface is not needed
                                       built = 'ANY',
                                       loaded = 'ANY',
                                       Cwritten = 'ANY',
                                       RCfunDefs = 'ANY'
                                       
                                   ),
                                   methods = list(
                                       getDefs = function() {
                                           callSuper()
                                       },
                                       getHincludes = function() {
                                           callSuper()
                                       },
                                       getCPPincludes = function() {
                                           callSuper()
                                       },
                                       getCPPusings = function() {
                                           callSuper()
                                       },
                                       genNeededTypes = function(debugCpp = FALSE, fromModel = FALSE) {
                                           for(i in seq_along(nimCompProc$neededTypes)) {
                                               neededType<- nimCompProc$neededTypes[[i]]
                                               if(inherits(neededType, 'nfMethodRC')) {
                                                   thisCppDef <- nimbleProject$getRCfunCppDef(neededType, NULLok = TRUE)
                                                   if(is.null(thisCppDef)) {
                                                       thisCppDef <- nimbleProject$needRCfunCppClass(neededType, genNeededTypes = TRUE, fromModel = fromModel)
                                                       neededTypeDefs[[neededType$uniqueName]] <<- thisCppDef
                                                   } else {
                                                       Hincludes <<- c(Hincludes, thisCppDef)
                                                       CPPincludes <<- c(CPPincludes, thisCppDef)
                                                   }
                                                   next
                                               }
                                               if(inherits(neededType, 'symbolModelValues')) {
                                                   thisCppDef <- nimbleProject$getModelValuesCppDef(neededType$mvConf, NULLok = TRUE)
                                                   if(is.null(thisCppDef)) {
                                                       thisCppDef <- nimbleProject$needModelValuesCppClass(neededType$mvConf, fromModel = fromModel)
                                                       mvClassName <- environment(neededType$mvConf)$className
                                                       neededTypeDefs[[mvClassName]] <<- thisCppDef
                                                   } else {
                                                       Hincludes <<- c(Hincludes, thisCppDef)
                                                       CPPincludes <<- c(CPPincludes, thisCppDef)
                                                   }
                                                   next
                                               }
                                               if(inherits(neededType, 'symbolNimbleFunction')) {
                                                   generatorName <- environment(neededType$nfProc$nfGenerator)$name
                                                   thisCppDef <- nimbleProject$getNimbleFunctionCppDef(generatorName = generatorName)
                                                   if(is.null(thisCppDef)) {
                                                       className <- names(nimCompProc$neededTypes)[i]
                                                       if(neededType$type == 'nimbleFunction')
                                                           thisCppDef <- nimbleProject$buildNimbleFunctionCompilationInfo(generatorName = generatorName, fromModel = fromModel)
                                                       else if(neededType$type == 'nimbleFunctionVirtual')
                                                           thisCppDef <- nimbleProject$buildVirtualNimbleFunctionCompilationInfo(vfun = generatorName)
                                                       else stop('symbolNimbleFunction does not have type nimbleFunction or nimbleFunctionVirtual')
                                                       neededTypeDefs[[ className ]] <<- thisCppDef
                                                   } else {
                                                       Hincludes <<- c(Hincludes, thisCppDef)
                                                       CPPincludes <<- c(CPPincludes, thisCppDef)
                                                   }
                                                   next
                                               }
                                               if(inherits(neededType, 'symbolNimbleList')) {
                                                 CPPincludes <<- c(CPPincludes, nimbleIncludeFile("smartPtrs.h"))
                                                 generatorName <- neededType$nlProc$name
                                                 if(generatorName == 'NIMBLE_ADCLASS') Hincludes <<- c(Hincludes, nimbleIncludeFile("nimbleCppAD.h"))
                                                 thisCppDef <- nimbleProject$getNimbleListCppDef(generatorName = generatorName)
                                                 if(is.null(thisCppDef)){
                                                      className <- names(nimCompProc$neededTypes)[i]
                                                      thisCppDef <- nimbleProject$buildNimbleListCompilationInfo(className = generatorName, fromModel = fromModel)
                                                      neededTypeDefs[[ className ]] <<- thisCppDef
                                                      Hincludes <<- c(Hincludes, thisCppDef)
                                                      CPPincludes <<- c(CPPincludes, thisCppDef)
                                                 }
                                                 next
                                               }
                                               if(inherits(neededType, 'symbolNimbleFunctionList')) {
                                                   baseClassName <- environment(neededType$baseClass)$name
                                                   thisCppDef <- nimbleProject$getNimbleFunctionCppDef(generatorName = baseClassName)
                                                   if(is.null(thisCppDef)) {
                                                       thisCppDef <- nimbleProject$buildVirtualNimbleFunctionCompilationInfo(vfun = neededType$baseClass)
                                                       neededTypeDefs[[baseClassName]] <<- thisCppDef
                                                   } else {
                                                       Hincludes <<- c(Hincludes, thisCppDef)
                                                       CPPincludes <<- c(CPPincludes, thisCppDef)
                                                   }
                                               }
                                           }
                                       },
                                       initialize = function(nimCompProc, debugCpp = FALSE, fromModel = FALSE, ...) {
                                           callSuper(...) ## must call this first because it sets objectDefs to list()
                                           RCfunDefs <<- list()
                                           if(!missing(nimCompProc)) processNimCompProc(nimCompProc, debugCpp = debugCpp, fromModel = fromModel)
                                           built <<- FALSE
                                           loaded <<- FALSE
                                           Cwritten <<- FALSE
                                       },
                                       processNimCompProc = function(ncp, debugCpp = FALSE, fromModel = FALSE) {
                                           ncp$cppDef <- .self
                                           nimCompProc <<- ncp
                                           genNeededTypes(debugCpp = debugCpp, fromModel = fromModel)
                                           objectDefs <<- symbolTable2cppVars(ncp$getSymbolTable())
                                       },
                                       addCopyFromRobject = function() {
                                           ## The next line creates the cppCopyTypes exactly the same way as in buildNimbleObjInterface
                                           ## and CmultiNimbleObjClass::initialize.
                                           cppCopyTypes <- makeNimbleFxnCppCopyTypes(nimCompProc$getSymbolTable(), objectDefs$getSymbolNames())
                                           copyFromRobjectDefs <- makeCopyFromRobjectDef(className = nfProc$name, cppCopyTypes, .self$nfProc$instances[[1]])
                                           functionDefs[['copyFromRobject']] <<- copyFromRobjectDefs$copyFromRobjectDef
                                          ## SEXPmemberInterfaceFuns[['copyFromRobject']] <<- copyFromRobjectDefs$copyFromRobjectInterfaceDef
                                       },
                                       buildAll = function(where = where) {
                                           makeCppNames()
                                           buildConstructorFunctionDef()
                                           buildSEXPgenerator(finalizer = "namedObjects_Finalizer")
                                           buildRgenerator(where = where)
                                           buildCmultiInterface()
                                       },
                                       buildRgenerator = function() {message('whoops, base class version of buildRgenerator')},
                                       buildCmultiInterface = function() {message('whoops, base class version of buildCmultiInterface')},
                                       makeCppNames = function() {
                                           Rnames2CppNames <<- as.list(Rname2CppName(objectDefs$getSymbolNames()))
                                           names(Rnames2CppNames) <<- objectDefs$getSymbolNames()
                                       },
                                       buildConstructorFunctionDef = function() {
                                         newNestedListLines <- list()
                                         flagLine <- list()
                                         pointerLine <- list()
                                         if(!(is.null(nimCompProc[['nimbleListObj']]))){
                                           flagText <- paste0('RCopiedFlag = false')
                                           flagLine[[1]] <- substitute(FLAGTEXT, 
                                                            list(FLAGTEXT = as.name(flagText)))
                                           pointerText <- paste0('RObjectPointer = NULL')
                                           pointerLine[[1]] <- substitute(POINTERTEXT,
                                                                          list(POINTERTEXT = as.name(pointerText)))
                                           for(i in seq_along(nimCompProc$neededTypes)){
                                             newListText <- paste0(nimCompProc$neededTypes[[i]]$name, " = new ", as.name(names(nimCompProc$neededTypes)[i]))
                                             newNestedListLines[[i]] <- substitute(NEWLISTTEXT,
                                                                                   list(NEWLISTTEXT = as.name(newListText)))
                                           }
                                         }

                                         code <- putCodeLinesInBrackets(c(flagLine, pointerLine, newNestedListLines,
                                                                          list(namedObjectsConstructorCodeBlock())))
                                         conFunDef <- cppFunctionDef(name = name,
                                                                     returnType = emptyTypeInfo(),
                                                                     code = cppCodeBlock(code = code, skipBrackets = TRUE))
                                         functionDefs[['constructor']] <<- conFunDef
                                       }
                                   ),
                                   )

cppNimbleFunctionClass <- setRefClass('cppNimbleFunctionClass',
                                      contains = 'cppNimbleClassClass',
                                      fields = list(
                                          nfProc = 'ANY', ## an nfProcessing class, needed to get the member data symbol table post-compilation
                                          parentsSizeAndDims = 'ANY',
                                          ADconstantsInfo = 'ANY'
                                          ),
                                          methods = list(
                                              getDefs = function() {
                                                  c(callSuper(), SEXPmemberInterfaceFuns) 
                                              },
                                              getHincludes = function() {
                                                  c(callSuper(), unlist(lapply(SEXPmemberInterfaceFuns, function(x) x$getHincludes()), recursive = FALSE))
                                              },
                                              getCPPincludes = function() {
                                                  c(callSuper(), unlist(lapply(SEXPmemberInterfaceFuns, function(x) x$getCPPincludes()), recursive = FALSE))
                                              },
                                              getCPPusings = function() {
                                                  CPPuse <- unique(c(callSuper(), unlist(lapply(SEXPmemberInterfaceFuns, function(x) x$getCPPusings()))))
                                                  CPPuse
                                              },
                                              initialize = function(nfProc, isNode, debugCpp = FALSE, fromModel = FALSE, ...) {
                                                  callSuper(nfProc, debugCpp, fromModel, ...)
                                                  if(!missing(nfProc)) processNFproc(nfProc, debugCpp = debugCpp, fromModel = fromModel)
                                                  if(isNode) {
                                                      inheritance <<- inheritance[inheritance != 'NamedObjects']
                                                      baseClassObj <- environment(nfProc$nfGenerator)$contains
                                                      if(is.null(baseClassObj)) {
                                                          inheritance <<- c(inheritance, 'nodeFun')
                                                          parentsSizeAndDims <<- environment(nfProc$nfGenerator)$parentsSizeAndDims
                                                          ADconstantsInfo <<- environment(nfProc$nfGenerator)$ADconstantsInfo
                                                          
                                                      }
                                                  }
                                              },
                                              processNFproc = function(nfp, debugCpp = FALSE, fromModel = FALSE) {
                                                  nfProc <<- nimCompProc
                                                  buildFunctionDefs()
                                                  ## This is slightly klugey
                                                  ## The objectDefs here are for the member data
                                                  ## We need them to be the parentST for each member function
                                                  ## However the building of the cpp objects is slightly out of order, with the
                                                  ## member functions already having been built during nfProcessing.
                                                  for(i in seq_along(functionDefs)) {
                                                      functionDefs[[i]]$args$parentST <<- objectDefs
                                                  }
                                                  SEXPmemberInterfaceFuns <<- lapply(functionDefs,
                                                                                     function(x)
                                                                                         if(inherits(x$SEXPinterfaceFun, "cppFunctionDef"))
                                                                                             x$SEXPinterfaceFun
                                                                                         else
                                                                                             NULL
                                                                                     )
                                                  SEXPmemberInterfaceFuns <<-
                                                      SEXPmemberInterfaceFuns[ !unlist(lapply(SEXPmemberInterfaceFuns, is.null)) ]
                                                      
                                                  nimCompProc <<- nfProc
                                              },
                                              buildFunctionDefs = function() {
                                                  for(i in seq_along(nfProc$RCfunProcs)) {
                                                      RCname <- names(nfProc$RCfunProcs)[i]
                                                      functionDefs[[RCname]] <<- RCfunctionDef$new() ## all nodeFunction members are const f
                                                      functionDefs[[RCname]]$buildFunction(nfProc$RCfunProcs[[RCname]])
                                                      if(!grepl("ADproxyModel", RCname)) ## don't build "CALL_" function for AD functions
                                                          functionDefs[[RCname]]$buildSEXPinterfaceFun(className = nfProc$name)
                                                      RCfunDefs[[RCname]] <<- functionDefs[[RCname]]
                                                  }
                                              },
                                              addTypeTemplateFunction = function( funName,
                                                                                 derivControl = list()) {
                                                  static <- isTRUE(derivControl[['static']])
                                                  regularFun <- RCfunDefs[[funName]]
                                                  useModelInfo <- list()
                                                  if(static) {
                                                      newFunName <- paste0(funName, '_AD_')
                                                      ## We need the template<TYPE_> function and also
                                                      ## whether a calculate was detected and if so its nodeFxnVector_name
                                                      res <-  makeTypeTemplateFunction(newFunName, regularFun)
                                                      functionDefs[[newFunName]] <<- res$fun
                                                      useModelInfo$nodeFxnVector_name <- res[['nodeFxnVector_name']]
                                                  } else {
                                                      ## Wee need a separate one for the "reconfigure" system that is non-static 
                                                      if(isTRUE(nimbleOptions('useADreconfigure'))) {
                                                          newFunName <- paste0(funName, '_AD2_')
                                                          res <- makeTypeTemplateFunction(newFunName,
                                                                                          regularFun,
                                                                                          static = FALSE)
                                                          functionDefs[[newFunName]] <<- res$fun
                                                          useModelInfo$nodeFxnVector_name <- res[['nodeFxnVector_name']]
                                                      }
                                                  }
                                                  useModelInfo
                                              },
                                              addADtapingFunction = function( funName,
                                                                             independentVarNames,
                                                                             dependentVarNames,
                                                                             useModelInfo = NULL,
                                                                             derivControl = list()) {
                                                  regularFun <- RCfunDefs[[funName]]
                                                  static <- isTRUE(derivControl[['static']])
                                                  if(static) {
                                                      ADfunName <- paste0(funName, '_AD_')
                                                      newFunName <- paste0(funName, '_callForADtaping_')
                                                      functionDefs[[newFunName]] <<- makeADtapingFunction(newFunName,
                                                                                                          regularFun,
                                                                                                          ADfunName,
                                                                                                          independentVarNames,
                                                                                                          dependentVarNames,
                                                                                                          nfProc$isNode,
                                                                                                          className = name)
                                                  }
                                                  if(!static) {
                                                      if(isTRUE(nimbleOptions("useADreconfigure"))) {
                                                          ADfunName2 <- paste0(funName, '_AD2_')
                                                          newFunName2 <- paste0(funName, '_callForADtaping2_')
                                                          functionDefs[[newFunName2]] <<- makeADtapingFunction2(newFunName2,
                                                                                                                regularFun,
                                                                                                                ADfunName2,
                                                                                                                independentVarNames,
                                                                                                                dependentVarNames,
                                                                                                                nfProc$isNode,
                                                                                                                className = name,
                                                                                                                useModelInfo = useModelInfo)
                                                      }
                                                  }
                                                  invisible(NULL)
                                              },
                                              addADargumentTransferFunction = function(funName,
                                                                                       independentVarNames,
                                                                                       funIndex,
                                                                                       derivControl = list()) {
                                                  static <- isTRUE(derivControl[['static']])
                                                  regularFun <- RCfunDefs[[funName]]
                                                  ## funIndex <- which(names(environment(nfProc$nfGenerator)$enableDerivs == funName)) ## needed for correct index for allADtapePtrs_
                                                  if(static) {
                                                      newFunName <- paste0(funName, '_ADargumentTransfer_')
                                                      functionDefs[[newFunName]] <<- makeADargumentTransferFunction(newFunName,
                                                                                                                    regularFun,
                                                                                                                    independentVarNames,
                                                                                                                    funIndex,
                                                                                                                    parentsSizeAndDims,
                                                                                                                    ADconstantsInfo)
                                                  } else {
                                                  ## Version 2 is for the "reconfig" version
                                                      if(isTRUE(nimbleOptions("useADreconfigure"))) {
                                                          newFunName2 <- paste0(funName, '_ADargumentTransfer2_')
                                                          ## For now, use same funIndex since this will use the non-static tape vector, but we may want to be more careful
                                                          functionDefs[[newFunName2]] <<- makeADargumentTransferFunction2(newFunName2,
                                                                                                                          regularFun,
                                                                                                                          independentVarNames,
                                                                                                                          funIndex,
                                                                                                                          parentsSizeAndDims,
                                                                                                                          ADconstantsInfo)
                                                      }
                                                  }
                                                  invisible(NULL)
                                              },
                                              addStaticInitClass = function(staticMethods) {
                                                  neededTypeDefs[['staticInitClass']] <<- makeStaticInitClass(.self, staticMethods) ##
                                                  invisible(NULL)
                                              },
                                              addADclassContentOneFun = function(funName,
                                                                                 derivControl,
                                                                                 funIndex) {
                                                  outSym <- nfProc$RCfunProcs[[funName]]$compileInfo$returnSymbol
                                                  checkADargument(funName, outSym, returnType = TRUE)
                                                  if(length(nfProc$RCfunProcs[[funName]]$nameSubList) == 0) stop(paste0('Derivatives cannot be enabled for method ', funName, ', since this method has no arguments.'))
                                                  if(!nfProc$isNode){
                                                      for(iArg in seq_along(functionDefs[[funName]]$args$symbols)){
                                                          arg <- functionDefs[[funName]]$args$symbols[[iArg]]
                                                          argSym <- nfProc$RCfunProcs[[funName]]$compileInfo$origLocalSymTab$getSymbolObject(arg$name)
                                                          argName <- names(nfProc$RCfunProcs[[funName]]$nameSubList)[iArg]
                                                          checkADargument(funName, argSym, argName = argName)
                                                      }
                                                  }
                                                  useModelInfo <- addTypeTemplateFunction(funName,
                                                                                          derivControl = derivControl) ## returns either NULL (if there is no model$calculate) or a nodeFxnVector name (if there is)
                                                  independentVarNames <- names(functionDefs[[funName]]$args$symbols)
                                                  if(nfProc$isNode) independentVarNames <- independentVarNames[-1]  ## remove ARG1_INDEXEDNODEINFO__ from independentVars
                                                  addADtapingFunction(funName,
                                                                      independentVarNames = independentVarNames,
                                                                      dependentVarNames = 'ANS_',
                                                                      useModelInfo,
                                                                      derivControl = derivControl )
                                                  addADargumentTransferFunction(funName,
                                                                                independentVarNames = independentVarNames,
                                                                                funIndex = funIndex,
                                                                                derivControl = derivControl)
                                              },
                                              checkADargument = function(funName, argSym, argName = NULL, returnType = FALSE){
                                                  argTypeText <- if(returnType) 'returnType' else 'argument'
                                                  if(argSym$type != 'double')
                                                      stop(paste0('The ', argName, ' ', argTypeText, ' of the ', funName, ' method is not a double, this method cannot have derivatives enabled.'))
                                        # if(!(argSym$nDim %in% c(0,1)))
                                        #    stop(paste0('The ', argName, ' ', argTypeText, ' of the ', funName, ' method must be a double scalar or double vector for derivatives to be enabled.'))
                                                  if((argSym$nDim != 0) && is.na(argSym$size))
                                                      stop(paste0('To enable derivatives, size must be given for the ',
                                                                  argName, ' ', argTypeText, ' of the ', funName,
                                                                  ' method,  e.g. double(1, 3) for a length 3 vector.' ))
                                              },
                                              addADclassContent = function() {
                                        # CPPincludes <<- c("<TMB/distributions_R.hpp>", CPPincludes)
                                                  enableDerivs <- environment(nfProc$nfGenerator)$enableDerivs
                                                  numEnabledFuns <- length(enableDerivs) 
                                                  boolStaticAD <- unlist(lapply(enableDerivs,
                                                                         function(x) isTRUE(x[['static']])))
                                                  ## would maybe have been easier to split on boolStaticAD
                                                  includeStatic <- any(boolStaticAD)
                                                  includeNonStatic <- any(!boolStaticAD)
                                                  numNonStatic <- sum(!boolStaticAD)
                                                  enableNames <- names(enableDerivs)
                                                  staticNames <- enableNames[boolStaticAD]
                                                  nonStaticNames <- enableNames[!boolStaticAD]
                                                  
                                                  if(includeNonStatic)
                                                      Hincludes <<- c("<array>", Hincludes)
                                                  ## These are moved to CPPincludes so they are #included before Rmath.h
                                                  CPPincludes <<- c(nimbleIncludeFile("nimbleCppAD.h"),
                                                                  nimbleIncludeFile("nimDerivs_TMB.h"), CPPincludes)
                                                  addInheritance("nimbleFunctionCppADbase")
                                                  if(includeStatic) {
                                                      objectDefs$addSymbol(cppVarFull(baseType = 'vector',
                                                                                      templateArgs = list(cppVarFull(baseType = 'CppAD::ADFun',
                                                                                                                     templateArgs = list('double'),
                                                                                                                     ptr = 1)),
                                                                                      static = TRUE,
                                                                                      name = 'allADtapePtrs_'))
                                                      whichBoolStaticAD <- which(boolStaticAD)
                                                      for(iiAD in seq_along(whichBoolStaticAD)) {
                                                          iAD <- whichBoolStaticAD[iiAD]
                                                          addADclassContentOneFun(enableNames[iAD],
                                                                                  enableDerivs[[iAD]],
                                                                                  funIndex = iiAD)
                                                      }
                                                  }
                                                  if(includeNonStatic) {
                                                      if(isTRUE(nimbleOptions('useADreconfigure'))) {
                                                          ## 1. Add myADtapePtrs_.  allADtapePtrs_ is static. myADtapePtrs_ is not static.
                                                          objectDefs$addSymbol(cppVarFull(baseType = 'std::array',
                                                                                          templateArgs = list(cppVarFull(baseType = 'CppAD::ADFun',
                                                                                                                         templateArgs = list('double'),
                                                                                                                         ptr = 1),
                                                                                                              numNonStatic),
                                                                                          static = FALSE,
                                                                                          name = 'myADtapePtrs_'))
                                                          ## 2. Add a destructor
                                                          ## Note that i becomes (i)-1 in cppOutput processing
                                                          ## but for the cppLiteral we need to code it directly.
                                                          destForCode <- substitute(for(i in 1:N) {
                                                                                        if(myADtapePtrs_[i])
                                                                                            cppLiteral("delete myADtapePtrs_[(i)-1];")
                                                                                    },
                                                                                    list(N = numNonStatic))
                                                          destSymTab <- symbolTable$new()
                                                          destSymTab$addSymbol( nimble:::cppInt(name = "i") )
                                                          destCode <- putCodeLinesInBrackets(destForCode) ## placeholder
                                                          nimDestCode <- RparseTree2ExprClasses(destCode)
                                                          destName <- paste0("~", name)
                                                          destFunDef <- cppFunctionDef(name = destName,
                                                                                       returnType = emptyTypeInfo(),
                                                                                       code = cppCodeBlock(code = nimDestCode,
                                                                                                           objectDefs = destSymTab,
                                                                                                           skipBrackets = TRUE))
                                                          functionDefs[["destructor"]] <<- destFunDef
                                                          
                                                          whichBoolNonStaticAD <- which(!boolStaticAD)
                                                          for(iiAD in seq_along(whichBoolNonStaticAD)) {
                                                              iAD <- whichBoolNonStaticAD[iiAD]
                                                              addADclassContentOneFun(enableNames[iAD],
                                                                                      enableDerivs[[iAD]],
                                                                                      funIndex = iiAD)
                                                          }
                                                      } else {
                                                          warning(paste0('non-static AD case requested for ',
                                                                         paste(nonStaticNames, sep=" ", collapse = " "),
                                                                         ' but nimbleOption useADreconfigure is not TRUE'))
                                                      }
                                                  }
                                                  objectDefs$addSymbol(cppVarFull(name = 'ADtapeSetup',
                                                                                  baseType = 'nimbleCppADinfoClass'))
                                                  ## for(iAD in seq_along(enableDerivs)) {
                                                  ##     addADclassContentOneFun(enableNames[iAD],
                                                  ##                             enableDerivs[[iAD]],
                                                  ##                             funIndex = iAD)
                                                  ## }
                                                  ## for(adEnabledFun in environment(nfProc$nfGenerator)$enableDerivs){
                                                  ##     addADclassContentOneFun(adEnabledFun)
                                                  ## }
                                                  ## static declaration in the class definition
                                                  ## globals to hold the global static definition
                                                  if(includeStatic) {
                                                      globals <- cppGlobalObjects(name = paste0('staticGlobals_', name),
                                                                                  staticMembers = TRUE)
                                                      globals$objectDefs[['staticGlobalTape']] <-
                                                          cppVarFull(baseType = 'vector',
                                                                     templateArgs = list(cppVarFull(baseType = 'CppAD::ADFun',
                                                                                                    templateArgs = list('double'),
                                                                                                    ptr = 1)),
                                                                     name = paste0(name,'::allADtapePtrs_'))
                                                      neededTypeDefs[['allADtapePtrs_']] <<- globals
                                                      
                                                      addStaticInitClass(staticNames)
                                                  }
                                                  invisible(NULL)
                                              },
                                              buildCmultiInterface = function(dll = NULL) {
                                                  sym <- if(!is.null(dll))
                                                             getNativeSymbolInfo(SEXPgeneratorFun$name, dll)
                                                         else
                                                             SEXPgeneratorFun$name
                                                  ## copyFromRobject_sym <- if(!is.null(dll))
                                                  ##                            getNativeSymbolInfo(SEXPmemberInterfaceFuns[['copyFromRobject']]$name, dll)
                                                  ##                        else
                                                  ##                            SEXPmemberInterfaceFuns[['copyFromRobject']]$name
                                                  CmultiInterface <<- CmultiNimbleFunctionClass(compiledNodeFun = .self,
                                                                                                basePtrCall = sym,
                                                                                                ##copyFromRobjectCall = copyFromRobject_sym,
                                                                                                project = nimbleProject)
                                              },
                                              buildRgenerator = function(where = globalenv(), dll = NULL) {
                                                  sym <- if(!is.null(dll))
                                                             getNativeSymbolInfo(SEXPgeneratorFun$name, dll)
                                                         else
                                                            SEXPgeneratorFun$name

                                                  Rgenerator <<- buildNimbleObjInterface(paste0(name,'_refClass') , .self, sym, where = where)
                                              },
                                              promoteCallable = function(R_NimbleFxn, asTopLevel = TRUE){
                                                  ## see comment in nimbleProjectClass::compileNimbleFunction
                                                  oldCobjectInterface <- nf_getRefClassObject(R_NimbleFxn)$.CobjectInterface
                                                  if(!is.list(oldCobjectInterface)) return(oldCobjectInterface)
                                                  existingExtPtrs <- oldCobjectInterface[[1]]$getExtPtrs(oldCobjectInterface[[2]])
                                                  thisDll <- oldCobjectInterface[[1]]$dll
                                                  newCobjectInterface <- Rgenerator(R_NimbleFxn, thisDll, project = nimbleProject, existingExtPtrs = existingExtPtrs)
                                                  newCobjectInterface
                                              },
                                              buildCallable = function(R_NimbleFxn, dll = NULL, asTopLevel = TRUE){
                                                  cppInterfaceObject <- NULL
                                                  if(asTopLevel) {
                                                      if(!is.null(Rgenerator))
                                                          cppInterfaceObject <- Rgenerator(R_NimbleFxn, dll, project = nimbleProject)
                                                      } else { ## actually this particular pathway should never be taken. asTopLevel = FALSE will occur only for nimbleProject$instantiateNimbleFunction
                                                          if(!is.null(CmultiInterface))
                                                              cppInterfaceObject <- CmultiInterface$addInstance(R_NimbleFxn, dll)
                                                      }
                                                  cppInterfaceObject
                                              },
                                              buildAll = function(where = where) {
                                                  baseClassObj <- environment(nfProc$nfGenerator)$contains
                                                  
                                                  if(!is.null(baseClassObj)) {
                                                      inheritance <<- inheritance[inheritance != 'NamedObjects']
                                                      baseClassName <- environment(baseClassObj)$name
                                                      addInheritance(baseClassName)
                                                      addAncestors('NamedObjects')
                                                  }
                                                  if(nimbleOptions('experimentalEnableDerivs') &&
                                                     length(environment(nfProc$nfGenerator)$enableDerivs) > 0) {
                                                      addADclassContent()
                                                  }
                                                  if('nodeFun' %in% .self$inheritance) {
                                                      updateADproxyModelMethods(.self)
                                                  }
                                                     
                                                  addCopyFromRobject()
                                                  
                                                  callSuper(where)
                                              }
                                          ),
                                      )

## The next block of code has the initial setup for an AST processing stage
## to make modifications for AD based on context etc.
modifyForAD_handlers <- list(eigenBlock = 'modifyForAD_eigenBlock',
                             calculate = 'modifyForAD_calculate',
                             getValues = 'modifyForAD_getSetValues',
                             setValues = 'modifyForAD_getSetValues')

exprClasses_modifyForAD <- function(code, symTab,
                                    workEnv = list2env(list(wrap_in_value = FALSE))) {
  if(code$isName) {
    if(!isTRUE(workEnv[['wrap_in_value']]))
      return(invisible())
    else {
      symObj <- symTab$getSymbolObject(code$name)
      if(!is.null(symObj)) {
        if(identical(symObj$baseType, "CppAD::AD")) {
          insertExprClassLayer(code$caller, code$callerArgID, 'Value') ## It is ok to leave some fields (type, sizeExpr) unpopulated at this late processing stage
        }
        ##code$name <- paste0("Value(", code$name, ")") ## For any more generality, Value should be in the AST, but for now paste it.
      }
    }
  }
  if(code$isCall) {
    recurse_modifyForAD(code, symTab, workEnv)
    handler <- modifyForAD_handlers[[code$name]]
    if(!is.null(handler))
      eval(call(handler, code, symTab, workEnv))
  }
  if(is.null(code$caller)) { ## It is the top level `{`
    if(identical(code$name, "{")) {
      if(isTRUE(workEnv[['hasCalculate']])) {
        if(is.null( workEnv[['nodeFxnVector_name']] ) )
          stop("While setting up C++ AD code: calculate found but nodeFxnVector_name not set.")
        ## insert setup_extraInput_step(nfv)
        set_model_tape_info_line <- RparseTree2ExprClasses(cppLiteral("set_CppAD_tape_info_for_model my_tape_info_RAII_(model_nodes_nodeFxnVector_includeData_TRUE_SU__derivs_, TYPE_::get_tape_id_nimble(), TYPE_::get_tape_handle_nimble());"))
        
        setupExtraInputLine <- RparseTree2ExprClasses(
          substitute(setup_extraInput_step(NFV),
                     list(NFV = as.name(workEnv$nodeFxnVector_name)))
        )
        ## This inserts a single line at the beginning by
        ## creating a new `{` block.
        ## If anything more general is needed later, we could use
        ## the exprClasses_insertAssertions system
        newExpr <- newBracketExpr(args = c(list(set_model_tape_info_line,
                                                setupExtraInputLine),
                                           code$args))
        code$args <- list()
        setArg(code, 1, newExpr)
      }
    }
  }
  
  invisible(NULL)
}

recurse_modifyForAD <- function(code, symTab, workEnv) {
  for(i in seq_along(code$args)) {
    if(inherits(code$args[[i]], 'exprClass')) {
      exprClasses_modifyForAD(code$args[[i]], symTab, workEnv)
    }
  }
  invisible(NULL)
}

modifyForAD_getSetValues <- function(code, symTab, workEnv) {
  if(code$name == 'setValues') {
    code$name <- 'setValues_AD_AD'
  }
  accessorName <- code$args[[2]]$name
  already_AD_name <- grepl("_AD_$", accessorName) ## _AD_ is at the end. I'm not sure it would ever already be there, so this is defensive.
  if(!already_AD_name) {
    classSymTab <- symTab$getParentST()$getParentST()
    newSymbol <- classSymTab$getSymbolObject(accessorName)$copy()
    newSymbol$name <- paste0(newSymbol$name, '_AD_')
    classSymTab$addSymbol(newSymbol)
    code$args[[2]]$name <- newSymbol$name
  }
  invisible(NULL)
}

modifyForAD_eigenBlock <- function(code, symTab, workEnv) {
  orig_wrap_in_value <- workEnv[['wrap_in_value']]
  workEnv[['wrap_in_value']] <- TRUE
  recurse_modifyForAD(code, symTab, workEnv)
  workEnv[['wrap_in_value']] <- orig_wrap_in_value
  invisible(NULL)
}

modifyForAD_calculate <- function(code, symTab, workEnv) {
  workEnv$hasCalculate <- TRUE
  workEnv$nodeFxnVector_name <- code$args[[1]]$name
  recurse_modifyForAD(code, symTab, workEnv)
  code$name <- "calculate_ADproxyModel"
  code$args[[2]] <- 1 ## this sets includeExtraOutputStep = true
  invisible(NULL)
}

updateADproxyModelMethods <- function(.self) {
    ## Update return type and names of functions like dnorm -> nimDerivs_dnorm
    functionNames <- names(.self$functionDefs)
    ADproxyModel_functionNames <- functionNames[ grepl("_ADproxyModel", functionNames ) ]
    if(length(ADproxyModel_functionNames) > 0) {
        ## The following two headers were added to CPPincludes because nimDerives_TMB must
        ## come before Rmath.h, so that functions like pbeta do not get re-defined Rf_pbeta,
        ## which causes problems in TMB headers.
        .self$CPPincludes <- c(nimbleIncludeFile("nimbleCppAD.h"),
                               nimbleIncludeFile("nimDerivs_TMB.h"), .self$CPPincludes)
    }
    for(fn in ADproxyModel_functionNames) {
        thisDef <- .self$functionDefs[[fn]]
        thisDef$returnType <- cppVarSym2templateTypeCppVarSym( thisDef$returnType,
                                                              replacementBaseType = "CppAD::AD",
                                                              replacementTemplateArgs = "double" )
        parentST <- thisDef$code$objectDefs$getParentST()
        thisDef$code$objectDefs <-
            symbolTable2templateTypeSymbolTable(thisDef$code$objectDefs,
                                                replacementBaseType = "CppAD::AD",
                                                replacementTemplateArgs = "double" )
        thisDef$code$objectDefs$setParentST(parentST)
        thisDef$code$cppADCode <- 2L 
        ADtypeDefs <- symbolTable()
        ADtypeDefs$addSymbol(cppVarFull(baseType = "typedef Matrix<CppAD::AD<double>, Dynamic, Dynamic>", name = "MatrixXd") )
        thisDef$code$typeDefs <- ADtypeDefs
##        exprClasses_modifyForAD(thisDef$code$code, thisDef$code$objectDefs)
    }
    classST <- .self$objectDefs
    classSymNames <- classST$getSymbolNames()
    ADproxySymNames <- classSymNames[ grepl("ADproxyModel_", classSymNames ) ]
    for(sn in ADproxySymNames) {
        newSym <-
            cppVarSym2templateTypeCppVarSym(classST$getSymbolObject(sn),
                                            replacementBaseType = "CppAD::AD",
                                            replacementTemplateArgs = "double")
        classST$addSymbol(newSym, allowReplace = TRUE)
    }
  
    NULL
}

makeSingleCopyCall <- function(varName, cppCopyType) {
    switch(cppCopyType,
           'nimbleFunction' = {
               cppLiteral(paste0("COPY_NIMBLE_FXN_FROM_R_OBJECT(\"", varName, "\");")) 
           },
           'numericVector' = {
               cppLiteral(paste0("COPY_NUMERIC_VECTOR_FROM_R_OBJECT(\"", varName, "\");"))
           },
           'doubleScalar' = {
               cppLiteral(paste0("COPY_DOUBLE_SCALAR_FROM_R_OBJECT(\"", varName, "\");"))
           },
           'integerScalar' = {
               cppLiteral(paste0("COPY_INTEGER_SCALAR_FROM_R_OBJECT(\"", varName, "\");"))
           },
           'logicalScalar' = {
               cppLiteral(paste0("COPY_LOGICAL_SCALAR_FROM_R_OBJECT(\"", varName, "\");"))
           },
           'nodeFxnVec' = {
               cppLiteral(paste0("COPY_NODE_FXN_VECTOR_FROM_R_OBJECT(\"", varName, "\");"))
           },
           'nodeFxnVec_nimDerivs' = {
               cppLiteral(paste0("COPY_NODE_FXN_VECTOR_DERIVS_FROM_R_OBJECT(\"", varName, "\");"))
           },
           'modelVarAccess' = {
               cppLiteral(paste0("COPY_VALUE_MAP_ACCESSORS_FROM_NODE_NAMES(\"", varName, "\");"))
           },
           
           NULL)
}

makeCopyFromRobjectDef <- function(className, cppCopyTypes, Robj) {
    ## Make method for copying from R object
    copyFromRobjectDef <- RCfunctionDef()
    copyFromRobjectDef$name <- 'copyFromRobject'
    copyFromRobjectDef$returnType <- cppVoid()
    copyFromRobjectDef$args <- symbolTable()
    localVars <- symbolTable()
    RobjectVarSym <- cppSEXP(name = 'Robject')
    copyFromRobjectDef$args$addSymbol(RobjectVarSym)
    copyCalls <- list()
    varNames <- names(cppCopyTypes)
    for(i in seq_along(cppCopyTypes)) {
      ## if(cppCopyTypes[[i]] == "nodeFxnVec"){
      ##   if(!is.null(Robj[[varNames[i]]]$nimDerivsInfo)){
      ##     cppCopyTypes[[i]] = "nodeFxnVec_derivs"
      ##   }
      ## }
      copyCalls[[varNames[i]]] <- makeSingleCopyCall(varNames[i], cppCopyTypes[[i]])
    }

    if(length(copyCalls) == 0) {
        ## create an empty function
        allRcode <- do.call('call',
                            list('{'),
                            quote = TRUE
                            )
    } else {
        unprotectCount <- 2 + length(copyCalls) ## 2 from SETUP_S_xData
        ## copyCalls <- list(cppLiteral("SETUP_S_xData;"),
        ##                   cppLiteral("COPY_NUMERIC_VECTOR_FROM_R_OBJECT(Robject, \"v\");"),
        ##                   cppLiteral(paste0("UNPROTECT(",unprotectCount,");"))
        ##                   )
        allRcode <- do.call('call',
                            c(list('{'),
                              list(cppLiteral("SETUP_S_xData;")),
                              copyCalls,
                              list(cppLiteral(paste0("UNPROTECT(",unprotectCount,");")))),
                            quote = TRUE
                            )
    }
    allCode <- RparseTree2ExprClasses(allRcode)
    copyFromRobjectDef$code <- cppCodeBlock(code = allCode,
                                            objectDefs = localVars)

    ## Make SEXP interface function to call from R:
    ## SEXPinterfaceCname <- paste0("CALL_", className, "_copyFromRobject")
    ## interfaceArgs <- symbolTable()
    ## interfaceArgs$addSymbol(RobjectVarSym)
    ## interfaceArgs$addSymbol(cppSEXP(name = 'SextPtrToObject'))

    ## RHScall <- as.call(list(as.name('copyFromRobject'),
    ##                         as.name('Robject')))
    ## castedCall <- substitute(cppMemberDereference(
    ##     template(static_cast, cppPtrType(CN))(R_ExternalPtrAddr(SextPtrToObject)), RHS),
    ##     list(CN = as.name(className), RHS = RHScall))
    ## returnLine <- quote(return(R_NilValue))
    
    ## interfaceCode <- embedListInRbracket(list(castedCall, returnLine))
    ## copyFromRobjectInterfaceDef <- cppFunctionDef(name = SEXPinterfaceCname,
    ##                                               args = interfaceArgs,
    ##                                               code = cppCodeBlock(code = RparseTree2ExprClasses(interfaceCode),
    ##                                                                   objectDefs = localVars), ## also empty local vars
    ##                                               returnType = cppSEXP(),
    ##                                               externC = TRUE)
    copyFromRobjectInterfaceDef <- NULL
    
    list(copyFromRobjectDef = copyFromRobjectDef,
         copyFromRobjectInterfaceDef = copyFromRobjectInterfaceDef)
}

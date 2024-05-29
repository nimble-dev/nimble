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
                abstract <- !isFALSE( environment(nfProc$nfGenerator)$methodControl[[i]]$required ) # default required = TRUE.  required is a synonym for abstract.  If it's abstract, it's required in derived classes.
                functionDefs[[i]] <<- RCfunctionDef(virtual = TRUE, abstract = abstract)
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
                                         buildDerivsInCopier <- FALSE
                                         if(inherits(nimCompProc, 'virtualNFprocessing')) {
                                           buildDerivs <- environment(nimCompProc$nfGenerator)[['buildDerivs']]
                                           buildDerivsInCopier <- (length(buildDerivs) > 0) & (!isFALSE(buildDerivs))
                                         }
                                         copyFromRobjectDefs <- makeCopyFromRobjectDef(cppCopyTypes,
                                                                                       buildDerivs = buildDerivsInCopier)
                                         functionDefs[['copyFromRobject']] <<- copyFromRobjectDefs$copyFromRobjectDef
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
                                              addTypeTemplateFunction = function(funName,
                                                                                 derivControl = list()) {
                                                  regularFun <- RCfunDefs[[funName]]
                                                  useModelInfo <- list()
                                                  ## if(static) {
                                                  ##     newFunName <- paste0(funName, '_AD_')
                                                  ##     ## We need the template<TYPE_> function and also
                                                  ##     ## whether a calculate was detected and if so its nodeFxnVector_name
                                                  ##     res <-  makeTypeTemplateFunction(newName = newFunName,
                                                  ##                                      .self = regularFun,
                                                  ##                                      static = TRUE,
                                                  ##                                      useRecordingInfo = TRUE,
                                                  ##                                      derivControl = derivControl)
                                                  ##     functionDefs[[newFunName]] <<- res$fun
                                                  ##     useModelInfo$nodeFxnVector_name <- res[['nodeFxnVector_name']]
                                                  ## } else {
                                                      ## We need a separate one for the "reconfigure" system that is non-static 
##                                                      if(isTRUE(nimbleOptions('useADreconfigure'))) {
                                                  newFunName <- paste0(funName, '_AD2_')
                                                  res <- makeTypeTemplateFunction(newName = newFunName,
                                                                                  .self = regularFun,
                                                                                  ##static = FALSE,
                                                                                  useRecordingInfo = TRUE,
                                                                                  derivControl = derivControl)
                                                  functionDefs[[newFunName]] <<- res$fun
                                                  useModelInfo$nodeFxnVector_name <- res[['nodeFxnVector_name']]
                                                  
##                                                      }
##                                                  }
                                                  useModelInfo
                                              },
                                              add_deriv_function = function(funName,
                                                                            independentVarNames,
                                                                            derivControl = list()) {
                                                  regularFun <- RCfunDefs[[funName]]
                                                  ## static <- isTRUE(derivControl[['static']])
                                                  ## if(static) {
                                                  ##     message("add_deriv_function support for static is not written.")
                                                  ## } else {
                                                  if(isTRUE(getNimbleOption("useADreconfigure"))) {
                                                      newFunName <- paste0(funName, '_deriv_')
                                                      argTransName <- paste0(funName, '_ADargumentTransfer2_')
                                                      functionDefs[[newFunName]] <<- make_deriv_function(regularFun,
                                                                                                         newFunName,
                                                                                                         independentVarNames,
                                                                                                         argTransName)

                                                      ## Add meta-taping type template function of the new function
                                                      newMetaFunName <- paste0(funName, '_AD2__deriv_')
                                                      argTransName <- paste0(funName, '_ADargumentTransfer2__AD2_')
                                                      
                                                      functionDefs[[newMetaFunName]] <<- make_deriv_function(regularFun,
                                                                                                             newMetaFunName,
                                                                                                             independentVarNames,
                                                                                                             argTransName,
                                                                                                             meta = TRUE)
                                                  } else {
                                                      message("Don't know how to use static=FALSE without nimbleOptions('useADreconfigure')." )
                                                  }
                                              ## }
                                              },
                                              addADtapingFunction = function( funName,
                                                                             independentVarNames,
                                                                             dependentVarNames,
                                                                             useModelInfo = NULL,
                                                                             derivControl = list()) {
                                                regularFun <- RCfunDefs[[funName]]
                                                if(is.null(useModelInfo))
                                                  useModelInfo <- list()
##                                                  static <- isTRUE(derivControl[['static']])
                                                  ## if(static) {
                                                  ##     ADfunName <- paste0(funName, '_AD_')
                                                  ##     newFunName <- paste0(funName, '_callForADtaping2_')
                                                  ##     functionDefs[[newFunName]] <<- makeADtapingFunction2(newFunName,
                                                  ##                                                          regularFun,
                                                  ##                                                          ADfunName,
                                                  ##                                                          independentVarNames,
                                                  ##                                                          dependentVarNames,
                                                  ##                                                          nfProc$isNode,
                                                  ##                                                          className = name)
                                                  ## }
                                                  ## if(!static) {
                                                  ##     if(isTRUE(nimbleOptions("useADreconfigure"))) {
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
##                                                      }
##                                                  }
                                                  invisible(NULL)
                                              },
                                              addADargumentTransferFunction = function(funName,
                                                                                       independentVarNames,
                                                                                       funIndex,
                                                                                       useModelInfo = NULL,
                                                                                       derivControl = list()) {
                                                  ## static <- isTRUE(derivControl[['static']])
                                                  regularFun <- RCfunDefs[[funName]]
                                                  ## funIndex <- which(names(environment(nfProc$nfGenerator)$enableDerivs == funName)) ## needed for correct index for allADtapePtrs_
                                                  ## if(static) {
                                                  ##     newFunName <- paste0(funName, '_ADargumentTransfer_')
                                                  ##     functionDefs[[newFunName]] <<- makeADargumentTransferFunction(newFunName,
                                                  ##                                                                   regularFun,
                                                  ##                                                                   independentVarNames,
                                                  ##                                                                   funIndex,
                                                  ##                                                                   parentsSizeAndDims,
                                                  ##                                                                   ADconstantsInfo)
                                                  ## } else {
                                                  ## Version 2 is for the "reconfig" version
##                                                      if(isTRUE(nimbleOptions("useADreconfigure"))) {
                                                          newFunName2 <- paste0(funName, '_ADargumentTransfer2_')
                                                          ## For now, use same funIndex since this will use the non-static tape vector, but we may want to be more careful
                                                          callForTapingName <- paste0(funName, '_callForADtaping2_')
                                                          functionDefs[[newFunName2]] <<- makeADargumentTransferFunction2(newFunName2,
                                                                                                                          regularFun,
                                                                                                                          callForTapingName,
                                                                                                                          independentVarNames,
                                                                                                                          funIndex,
                                                                                                                          parentsSizeAndDims,
                                                                                                                          ADconstantsInfo,
                                                                                                                          useModelInfo,
                                                                                                                          derivControl)
                                                          newFunName2meta <- paste0(funName, '_ADargumentTransfer2__AD2_')
                                                          ## For now, use same funIndex since this will use the non-static tape vector, but we may want to be more careful
                                                          callForTapingName <- paste0(funName, '_callForADtaping2_')
                                                          functionDefs[[newFunName2meta]] <<- makeADargumentTransferFunction2(newFunName2meta,
                                                                                                                              regularFun,
                                                                                                                              callForTapingName,
                                                                                                                              independentVarNames,
                                                                                                                              funIndex,
                                                                                                                              parentsSizeAndDims,
                                                                                                                              ADconstantsInfo,
                                                                                                                              useModelInfo,
                                                                                                                              derivControl,
                                                                                                                              metaTape = TRUE)
##                                                      }
                                                  ## }
                                                  invisible(NULL)
                                              },
                                              ## addStaticInitClass = function(staticMethods) {
                                              ##     neededTypeDefs[['staticInitClass']] <<- makeStaticInitClass(.self, staticMethods) ##
                                              ##     invisible(NULL)
                                              ## },
                                              addADclassContentOneFun = function(funName,
                                                                                 derivControl,
                                                                                 funIndex) {
                                                ##outSym <- nfProc$RCfunProcs[[funName]]$compileInfo$returnSymbol
                                                ## TO-DO: revisit checking arguments for deriv compability
                                                  ## isStatic <- isTRUE(derivControl[['static']])
                                                  ## if(!isTRUE(derivControl[['meta']]) & isStatic)
                                                  ##   checkADargument(funName, outSym, returnType = TRUE)
                                                  RCfunProc <- nfProc$RCfunProcs[[funName]]
                                                  nameSubList <- RCfunProc$nameSubList
                                                  compileInfo <- nfProc$compileInfos[[funName]]
                                                  if(length(nameSubList) == 0)
                                                    stop(paste0('Derivatives cannot be built for method ',
                                                                funName,
                                                                ', since this method has no arguments.'))
                                                  ## if(!nfProc$isNode & !isTRUE(derivControl[['meta']]) & isStatic){
                                                  ##     for(iArg in seq_along(functionDefs[[funName]]$args$symbols)){
                                                  ##         arg <- functionDefs[[funName]]$args$symbols[[iArg]]
                                                  ##         argSym <- nfProc$RCfunProcs[[funName]]$compileInfo$origLocalSymTab$getSymbolObject(arg$name)
                                                  ##         argName <- names(nfProc$RCfunProcs[[funName]]$nameSubList)[iArg]
                                                  ##         checkADargument(funName, argSym, argName = argName)
                                                  ##     }
                                                  ## }
                                                  useModelInfo <- addTypeTemplateFunction(funName,
                                                                                          derivControl = derivControl) ## returns either NULL (if there is no model$calculate) or a nodeFxnVector name (if there is)
                                                  ## independentVarNames <- names(functionDefs[[funName]]$args$symbols)
                                                  independentVarNames <- as.character(nameSubList)
                                               ## ivBaseTypes <- unlist(lapply(functionDefs[[funName]]$args$symbols, `[[`, "baseType"))
                                                  ivBaseTypes <- lapply(compileInfo$newLocalSymTab$symbols[independentVarNames], `[[`, "type")
                                                  ## independentVarNames <- names(ivBaseTypes)
                                                includeIV <- ivBaseTypes == "double" ## previously ivBaseTypes != "bool"
                                                ignore <- derivControl[["ignore"]]
                                                if(!is.null(ignore)) {
                                                    noDeriv_mangledArgNames <- as.character(unlist(nameSubList[ignore]))
                                                    includeIV <- includeIV & !(independentVarNames %in% noDeriv_mangledArgNames)
                                                }
                                                independentVarNames <- as.character(independentVarNames[includeIV])

                                                  if(nfProc$isNode) independentVarNames <- independentVarNames[-1]  ## remove ARG1_INDEXEDNODEINFO__ from independentVars
                                                  if(!isTRUE(derivControl[['meta']])) {
                                                      addADtapingFunction(funName,
                                                                          independentVarNames = independentVarNames,
                                                                          dependentVarNames = 'ANS_',
                                                                          useModelInfo,
                                                                          derivControl = derivControl )
                                                      addADargumentTransferFunction(funName,
                                                                                    independentVarNames = independentVarNames,
                                                                                    funIndex = funIndex,
                                                                                    useModelInfo = useModelInfo,
                                                                                    derivControl = derivControl)
                                                      add_deriv_function(funName,
                                                                         independentVarNames,
                                                                         derivControl)
                                                      
                                                      funIndex + 1 ## function return value increments by one in non-meta case
                                                  } else {
                                                      funIndex ## but does not increment in meta case
                                                  }
                                              },
                                              ## addNimDerivsCalculateContentOneFun = function(funName,
                                              ##                                               derivControl,
                                              ##                                               funIndex) {
                                              ##     ## We only need to add one line here, which can go at the beginning
                                              ##     ## set_tape_ptr(ADtapeSetup, myADtapePtrs_[<funIndex>-1], true);
                                              ##     thisFunctionDef <- functionDefs[[funName]]
                                              ##     thisCode <- thisFunctionDef$code$code ## lines of code
                                              ##     newFunIndex <- exprClasses_updateNimDerivsCalculate(thisCode, funIndex)
                                              ##     newFunIndex
                                              ## },
                                              checkADargument = function(funName, argSym, argName = NULL, returnType = FALSE){
                                                  argTypeText <- if(returnType) 'returnType' else 'argument'
                                                  if(argSym$type != 'double')
                                                      stop(paste0('The ', argName, ' ', argTypeText, ' of the ', funName, ' method is not a double, this method cannot have derivatives built.'))
                                        # if(!(argSym$nDim %in% c(0,1)))
                                        #    stop(paste0('The ', argName, ' ', argTypeText, ' of the ', funName, ' method must be a double scalar or double vector for derivatives to be enabled.'))
                                                  if((argSym$nDim != 0) && is.na(argSym$size))
                                                      stop(paste0('To build derivatives, size must be given for the ',
                                                                  argName, ' ', argTypeText, ' of the ', funName,
                                                                  ' method,  e.g. double(1, 3) for a length 3 vector.' ))
                                              },
                                              addADclassContent = function() {
                                                constructorCode <- NULL ## default return object
                                                buildDerivs <- environment(nfProc$nfGenerator)$buildDerivs

                                                buildNames <- names(buildDerivs)
                                                ADinUse <- FALSE
                                                for(iED in seq_along(buildDerivs)) {
                                                  if(!isTRUE(buildDerivs[[iED]][['isNode']])) {
                                                    ADinUse <- TRUE
                                                    if(!isTRUE(buildDerivs[[iED]][['calculate']])) {
                                                      addADclassContentOneFun(buildNames[iED],
                                                                              buildDerivs[[iED]],
                                                                              funIndex = iED) ## funIndex might be deprecated
                                                    }
                                                  }
                                                }

                                                if(ADinUse) {
                                                  Hincludes <<- c(nimbleIncludeFile("nimbleCppADbaseClass.h"),
                                                                  Hincludes)
                                                  CPPincludes <<- c(nimbleIncludeFile("nimbleCppAD.h"),
                                                                    nimbleIncludeFile("nimDerivs_dists.h"),
                                                                    ## The ".cpp" extension make it end up linked as a .o
                                                                    ## in cppProjectClass$writeFiles. Prior to that, it is
                                                                    ## copied and locally compiled by cppProjectClass
                                                                    "nimbleCppADbaseClass.cpp", CPPincludes)
                                                  addInheritance("nimbleFunctionCppADbase")
                                                  addADinfoObjects(.self)
                                                }
                                                constructorCode                                                      
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
                                                  handleDerivs <- isTRUE(getNimbleOption("enableDerivs")) &&
                                                      length(environment(nfProc$nfGenerator)$buildDerivs) > 0
                                                  if(handleDerivs) {
                                                      constructorCode <- addADclassContent() ## Might generate code to insert into constructor, which is built later
                                                      if('nodeFun' %in% .self$inheritance) {
                                                        updateADproxyModelMethods(.self)
                                                      }
                                                  }
                                                  
                                                  addCopyFromRobject()
                                                  
                                                  callSuper(where)
                                                  if(handleDerivs) {
                                                      if(!is.null(constructorCode)) { ## insert AD-related code to constructor
                                                          conDef <- functionDefs[['constructor']] ## Very hacky way to insert this code.
                                                          curlen <- length(conDef$code$code[[2]]$code)
                                                          conDef$code$code[[2]]$code[[curlen + 1 ]] <- constructorCode
                                                      }
                                                  }
                                              }
                                          ),
                                      )


## updateNimDerivsCalculateHandlers <- list(nimDerivs_calculate = "updateNimDerivsCalculate_update")

## exprClasses_updateNimDerivsCalculate <- function(code, funIndex) {
##   if(code$isName) return(funIndex)
##   if(code$isCall) {
##     for(i in seq_along(code$args)) {
##       if(inherits(code$args[[i]], 'exprClass')) {
##         funIndex <- exprClasses_updateNimDerivsCalculate(code$args[[i]], funIndex)
##       }
##     }
##     handler <- updateNimDerivsCalculateHandlers[[code$name]]
##     if(!is.null(handler)) funIndex <- eval(call(handler, code, funIndex))
##   }
##   funIndex
## }

## updateNimDerivsCalculate_update <- function(code, funIndex) {
##   ## newExpr <- substitute(set_tape_ptr(ADtapeSetup, myADtapePtrs_[NUM], 1),
##   ##                       list(NUM = funIndex))
##   ## newExpr <- RparseTree2ExprClasses(newExpr)
##   newArg1 <- RparseTree2ExprClasses(quote(ADtapeSetup))
##   newArg2 <- substitute(myADtapePtrs_[NUM],
##                         list(NUM = funIndex))
##   newArg2 <- RparseTree2ExprClasses(newArg2)
##   oldArgs <- code$args
##   for(i in rev(seq_along(code$args))) {
##     setArg(code, i+2, oldArgs[[i]])
##   }
##   setArg(code, 1, newArg1)
##   setArg(code, 2, newArg2)
##   rm(oldArgs)
##   funIndex + 1
## }

## The next block of code has the initial setup for an AST processing stage
## to make modifications for AD based on context etc.
modifyForAD_handlers <- c(list(
    `[` = 'modifyForAD_indexingBracket',
    # pow = 'modifyForAD_issuePowWarning',
    eigenBlock = 'modifyForAD_eigenBlock',
    calculate = 'modifyForAD_calculate',
    getValues = 'modifyForAD_getSetValues',
    setValues = 'modifyForAD_getSetValues',
    nfMethod = 'modifyForAD_nfMethod',
    chainedCall = 'modifyForAD_chainedCall',
    getDerivs_wrapper = 'modifyForAD_getDerivs_wrapper',
    AssignEigenMap = 'modifyForAD_AssignEigenMap',
    ADbreak = 'modifyForAD_ADbreak',
    `*` = 'modifyForAD_matmult',
    eigInverse = 'modifyForAD_matinverse',
    MAKE_FIXED_VECTOR = 'modifyForAD_MAKE_FIXED_VECTOR'),
  #  makeCallList(assignmentOperators, 'modifyForAD_assign'), # Handled manually inside modifyForAD_recurseForAD
    makeCallList(recyclingRuleOperatorsAD, 'modifyForAD_RecyclingRule'),
    makeCallList(c('EIGEN_FS', 'EIGEN_BS', 'EIGEN_SOLVE', 'EIGEN_CHOL', 'nimNewMatrixD', 'nimCd',
                   'nimRepd','nimDiagonalD', 'nimNonseqIndexedd'),
                 'modifyForAD_prependNimDerivs'))

exprClasses_modifyForAD <- function(code, symTab,
                                    workEnv = list2env(list(wrap_in_value = FALSE,
                                                            .anyAD = FALSE))) {
  if(code$isName) {
    symObj <- symTab$getSymbolObject(code$name, inherits = TRUE)
    # Check if the symbol is TYPE_ or CppAD::AD or a NimArr of TYPE_ or CppAD::AD
    symIsAD <- !is.null(symObj) && (identical(symObj$baseType, "TYPE_") ||
                                    identical(symObj$baseType, "CppAD::AD") ||
                                    identical(symObj$baseType, "NimArr") && (identical(symObj$templateArgs[[2]], "TYPE_") ||
                                                                             (inherits(symObj$templateArgs[[2]], "cppVarFull") &&
                                                                              identical(symObj$templateArgs[[2]]$baseType, "CppAD::AD"))))
    if(isTRUE(symIsAD)) workEnv$.anyAD <- TRUE

    # Is this "wrap_in_value" step ever used or is it cruft?
    if(!isTRUE(workEnv[['wrap_in_value']]))
      return(invisible())
    else {
      if(isTRUE(symIsAD)) ## N.B. Previously this came from symTab$getSymbolObject(code$name, inherits = FALSE)
        if(identical(symObj$baseType, "CppAD::AD")) # What objects get this type, rather than TYPE_
          insertExprClassLayer(code$caller, code$callerArgID, 'Value') ## It is ok to leave some fields (type, sizeExpr) unpopulated at this late processing stage
    }
  }
  if(code$isCall) {
    recurse_modifyForAD(code, symTab, workEnv)
    handler <- modifyForAD_handlers[[code$name]]
    if(!is.null(handler))
      eval(call(handler, code, symTab, workEnv))
    if(workEnv$RsymTab$symbolExists(code$name, TRUE)) ## Could be a nimbleFunction class method
      modifyForAD_recordable(code, symTab, workEnv)
    else {
      modifyForAD_RCfunction(code, symTab, workEnv)
    }
  }
  if(is.null(code$caller)) { ## It is the top level `{`
    if(identical(code$name, "{")) {
      if(!isTRUE(workEnv$.inModel)) {
        if(isTRUE(workEnv$.illegalADcall)) {
          allMsgs <- paste(workEnv$.illegalADmsgs, collapse="\n")
          stop(paste("Problem(s) creating derivative code.\n", allMsgs))
        }
      }
      if(isTRUE(workEnv[['hasCalculate']])) {
        if(is.null( workEnv[['nodeFxnVector_name']] ) )
          stop("While setting up C++ AD code: calculate found but nodeFxnVector_name not set.")
        ## insert setup_extraInput_step(nfv)
        set_model_tape_info_line <- RparseTree2ExprClasses(cppLiteral(
            paste0("set_CppAD_tape_info_for_model my_tape_info_RAII_(",
                   workEnv[['nodeFxnVector_name']],
                   ", recordingInfo_.tape_id(), recordingInfo_.tape_handle());")))
                  ## ", TYPE_::get_tape_id_nimble(), TYPE_::get_tape_handle_nimble());")))
        

            set_atomic_info_line <- RparseTree2ExprClasses(cppLiteral(
                paste0("set_CppAD_atomic_info_for_model(",
                       workEnv[['nodeFxnVector_name']],
                       ", recordingInfo_.atomic_vec_ptr());")))
                      ## ", CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage());")))

        ## This inserts a single line at the beginning by
        ## creating a new `{` block.
        ## If anything more general is needed later, we could use
        ## the exprClasses_insertAssertions system
        newExpr <- newBracketExpr(args = c(list(set_model_tape_info_line,
                                                set_atomic_info_line),
                                           code$args))
        code$args <- list()
        setArg(code, 1, newExpr)
      }
    }
  }  
  invisible(NULL)
}

recurse_modifyForAD <- function(code, symTab, workEnv) {
  if(is.list(code$aux))
    if(isTRUE(code$aux$.avoidAD))
      return(invisible(NULL))
  anyAD_on_input <- workEnv$.anyAD
  anyAD_thisCall <- FALSE
  assignment <- isTRUE(any(code$name == assignmentOperators))
  for(i in seq_along(code$args)) {
    if(inherits(code$args[[i]], 'exprClass')) {
      if(assignment) workEnv$onLHS <- i == 1
      workEnv$.anyAD <- FALSE
      exprClasses_modifyForAD(code$args[[i]], symTab, workEnv)
      anyAD_thisCall <- anyAD_thisCall || isTRUE(workEnv$.anyAD)
      if(assignment) workEnv$onLHS <- FALSE
    }
  }
  workEnv$.anyAD <- anyAD_thisCall
  invisible(NULL)
}

modifyForAD_indexingBracket <- function(code, symTab, workEnv) {
  if(!isTRUE(workEnv$.anyAD)) return(invisible(NULL))
  if(is.null(code$type)) return(invisible(NULL)) # arises from constructions like setMap that lack type annotation
  if(length(code$args) > 4) return(invisible(NULL)) # 4D or higher is not handled, assumed to be static indexing
  newName <- "stoch_ind_get"
  if(isTRUE(workEnv$onLHS))
    if(code$caller$name %in% assignmentOperators)
      newName <- "stoch_ind_set"
  code$name <- newName
  invisible(NULL)
}

## modifyForAD_issuePowWarning <- function(code, symTab, workEnv) {
##   message("   [Note] Operator `pow` may cause derivative problems with negative arguments. If the exponent is guaranteed to be an integer, use `pow_int` instead.")
##   invisible(NULL)
## }

modifyForAD_matmult <- function(code, symTab, workEnv) {
  if(isTRUE(getNimbleOption("useADmatMultAtomic")))
    if(length(code$args)==2)  ## something is wrong if there are not 2 args
      if(!(isEigScalar(code$args[[1]]) || isEigScalar(code$args[[2]]))) ## are both non-scalar?
        code$name <- 'nimDerivs_matmult'
  invisible(NULL)
}

modifyForAD_matinverse <- function(code, symTab, workEnv) {
  if(isTRUE(getNimbleOption("useADmatInverseAtomic")))
    code$name <- 'nimDerivs_matinverse'
  invisible(NULL)
}

modifyForAD_MAKE_FIXED_VECTOR <- function(code, symTab, workEnv) {
  if(length(code$args) >= 5) {
    go <- if(length(code$args) >= 6)
            isTRUE(code$args[[6]])
    else TRUE
    if(go)
      if(identical(code$args[[5]], "double"))
        code$args[[5]] <- "CppAD::AD<double>"
  }
}

modifyForAD_ADbreak <- function(code, symTab, workEnv) {
  code$name <- 'CppAD::Value'
  invisible(NULL)
}

modifyForAD_RecyclingRule <- function(code, symTab, workEnv) {
  ## The conversion to nimDerivs_really has to do with the C++ return type
  ## and that should always convert to CppAD::AD<double>, which is achieved
  ## by the following changes.  In other words, there is no need to parse
  ## through the arguments and see which are AD types and which not.
  code$name <- paste0('nimDerivs_', code$name)
  code$name <- gsub('::', '::nimDerivs_', code$name)
  code$name <- gsub("<MatrixXd>", "<MatrixXd_TYPE_>", code$name)
  invisible(NULL)
}

modifyForAD_AssignEigenMap <- function(code, symTab, workEnv) {
  EigMapName <- code$args[[1]]$name
  ## We extract the variable name by undoing the construction of the Eigen Map variable name.
  ## This would be nicer to do with generic functions that undo the label creation.
  ## At this time of this writing, hacking the solution here was more practical.
  ## 
  ## This method of extracting the name of the variable we are making a map into
  ## attempted to invert the name-generation and name-mangling steps that led to it.
  ## However, it is not fully general.  In addition, it is not clear if this would
  ## ever really make sense.  Although noDerivs_vars can make a variable enter
  ## a deriv-enabled function with its raw types, we should never really be doing
  ## any vectorized operation with it.  The reason is that the result will be CppAD::AD<double> (or TYPE_),
  ## since we do not attempt to somehow track AD and non-AD types through the AST processing.
  ## And Eigen requires same-type operations, I believe, because we manually cast when necessary among
  ## basic types, but we do not cast between AD and non-AD types.

  ## However, we try this for cases where it succeeds.  Notably it is useful for the results of
  ## lifting "order = nimC(0, 1)" from nimDerivs, for example.
  ## The pieces created from that are now protected from AD modification.

  go <- TRUE
  varSym <- NULL
  ## We try to look for the relevant base symbol.  If we can't find it, by default we do the AD modification
  # unpack the getPtr step, which might have an offset added via "+"
  ptrExpr <- code$args[[2]]
  if(ptrExpr$name == "+")
    ptrExpr <- ptrExpr$args[[1]]
  if(ptrExpr$name == "getPtr") {
    varName <- ptrExpr$args[[1]]$name
    varSym <- symTab$getSymbolObject(varName, inherits = TRUE)
  }
  if(!isTRUE(workEnv$.inModel)) {
      if(!is.null(varSym)) {
          scalarType <- varSym$templateArgs[[2]]
          go <- identical(scalarType, "TYPE_")
      }
  }
  if(go) {
    EigMapSym <- symTab$getSymbolObject(EigMapName)
    EigMapSym_baseType <- EigMapSym$baseType
    if(EigMapSym_baseType == "EigenMapStrd") {
      EigMapSym$baseType <- "EigenMapStrd_TYPE_"
      #MatrixXd_label <- code$args[[3]]
      code$args[[3]]$name <- "MatrixXd_TYPE_"
    }
    if(EigMapSym_baseType == "Map") {
      EigMapSym$templateArgs[[1]] <- "MatrixXd_TYPE_"
      code$args[[3]]$name <- "MatrixXd_TYPE_"
    }
  }
  invisible(NULL)
}

modifyForAD_prependNimDerivs <- function(code, symTab, workEnv) {
  origName <- code$name
  atomic <- TRUE
  if(origName == 'EIGEN_FS' | origName == 'EIGEN_BS')
    if(!isTRUE(getNimbleOption("useADsolveAtomic")))
        atomic <- FALSE
  if(origName == "EIGEN_CHOL")
    if(!isTRUE(getNimbleOption("useADcholAtomic")))
        atomic <- FALSE
  if(atomic)
    code$name <- paste0("nimDerivs_", origName)
  else
    code$name <- paste0("nimDerivs_", origName, '_no_atomic')
  invisible(NULL)
}

modifyForAD_RCfunction <- function(code, symTab, workEnv) {
  # If there is no AD happening, skip everything
  if(isTRUE(workEnv$.anyAD)) {
    lacksADsupport <- any(code$name == callsNotAllowedInAD)
    # If there is lack of AD support due to disallowed keyword, skip to error handling
    if(!lacksADsupport) {
      isRCwithoutAD <- FALSE
      # If there is code$aux, this is a potential AD call
      if(!is.null(code$aux)) {
        buildDerivs <- code$aux[['buildDerivs']]
        # If code$aux$buildDerivs is set, this is a potential AD call
        if(!is.null(buildDerivs)) {
          # If code$aux$buildDerivs is TRUE, insert the template argument
          if(isTRUE(buildDerivs)) {
            code$name <- paste0(code$name, "< CppAD::AD<double> >")
          } else {
            # If code$aux$buildDerivs is not TRUE, this is an RC function without AD support
            isRCwithoutAD <- TRUE
          }
        }
      }
      # So far lacksADsupport must be FALSE. Now we update it to include case isRCwithoutAD
      lacksADsupport <- isRCwithoutAD
    }
    if(lacksADsupport) {
      workEnv$.illegalADcall <- TRUE
      if(!is.character(workEnv$.illegalADmsgs)) workEnv$.illegalADmsgs <- character()
      workEnv$.illegalADmsgs <- c(workEnv$.illegalADmsgs,
                                  exprClassProcessingErrorMsg(
                                    code,
                                    "A function that might be intended for derivative tracking does not have buildDerivs set to enable that."
                                  ))
    }
  }
  ## if(isTRUE(workEnv$.anyAD)) {
  ##   isRCwithoutAD <- (!is.null(code$aux)) && (!isTRUE(code$aux[['buildDerivs']]))
  ##   isRCwithAD <- (!is.null(code$aux)) && (isTRUE(code$aux[['buildDerivs']]))
  ##   lacksADsupport <- isRCwithoutAD ||
  ##                      any(code$name == callsNotAllowedInAD)
  ##   supportsAD <- !lacksADsupport
  ##   if(supportsAD) {
  ##     if(isRCwithAD)
  ##       code$name <- paste0(code$name, "< CppAD::AD<double> >")
  ##   } else {
  ##     browser()
  ##     workEnv$.illegalADcall <- TRUE
  ##     if(!is.character(workEnv$.illegalADmsgs)) workEnv$.illegalADmsgs <- character()
  ##     workEnv$.illegalADmsgs <- c(workEnv$.illegalADmsgs,
  ##                                 exprClassProcessingErrorMsg(
  ##                                   code,
  ##                                   "A function that might be intended for derivative tracking does not have buildDerivs set to enable that."
  ##                                 ))
  ##   }
  ## }
  invisible(NULL)
}

modifyForAD_getDerivs_wrapper <- function(code, symTab, workEnv) {
  ## This really modifies the "argumentTransfer" call in the first argument of getDerivs_wrapper
  arg1 <- code$args[[1]]
  arg1$name <- paste0(arg1$name, "_AD2_")
  code$name <- "getDerivs_wrapper_meta"
  nArgs <- length(code$args)
  newExpr <- RparseTree2ExprClasses(quote(recordingInfo_))
  setArg(code, nArgs + 1, newExpr)
  invisible(NULL)
}

modifyForAD_nfMethod <- function(code, symTab, workEnv) {
  if(code$args[[1]]$name != "cppPointerDereference")
    message("   [Note] In modifyForAD_nfMethod, was expecting cppPointerDereference.  There must be another case that needs implementation.")
  objName <- code$args[[1]]$args[[1]]$name
  NFsymObj <- workEnv$RsymTab$getSymbolObject(objName, TRUE)
  methodName <- code$args[[2]]
  if(!is.character(methodName))
    message("   [Note] In modifyForAD_nfMethod, was expecting method name to be character.  There must be another case that needs implementation.")
  methodSymObj <- NFsymObj$nfProc$compileInfos[[methodName]]$newLocalSymTab$getSymbolObject(methodName, TRUE)
  buildDerivs <- methodSymObj$nfMethodRCobj$buildDerivs
  if(!is.null(buildDerivs))
    if(!isFALSE(buildDerivs)) {
      code$args[[2]] <- paste0( code$args[[2]], "_AD2_")
      code$name <- "nfMethodAD" ## This is a tag for modifyForAD_chainedCall
    }
  invisible(NULL)
}

modifyForAD_chainedCall <- function(code, symTab, workEnv) {
  arg1name <- code$args[[1]]$name
  if(arg1name == "nfMethodAD") { ## A tag set in modifyForAD_nfMethod
    setArg(code, length(code$args) + 1, RparseTree2ExprClasses(quote(recordingInfo_)))
    code$args[[1]]$name <- "nfMethod" ## tag is no longer needed. revert to regular nfMethodAD
  }
  invisible(NULL)
}

# This is when the current method (which has buildDerivs set, else we wouldn't be here)
# has a call to another method. We assume that method should also have buildDerivs set
# any thus we should append the "_AD2_" on its name and throw an error if it doesn't
# have buildDerivs set.
# However, if no arguments seem to involve AD, we leave the call unmodified.
modifyForAD_recordable <- function(code, symTab, workEnv) {
  symObj <- workEnv$RsymTab$getSymbolObject(code$name, TRUE)
  buildDerivs <- symObj$nfMethodRCobj$buildDerivs
  anyAD <- workEnv$.anyAD
  callSupportsAD <- (!is.null(buildDerivs)) && (!isFALSE(buildDerivs))
  if(anyAD) {
    if(callSupportsAD) {
      code$name <- paste0(code$name, "_AD2_")
      setArg(code, length(code$args) + 1, RparseTree2ExprClasses(quote(recordingInfo_)))
    } else {
      workEnv$.illegalADcall <- TRUE
      if(!is.character(workEnv$.illegalADmsgs)) workEnv$.illegalADmsgs <- character()
      workEnv$.illegalADmsgs <- c(workEnv$.illegalADmsgs,
                                  exprClassProcessingErrorMsg(
                                    code,
                                    "A function that might be intended for derivative tracking does not have buildDerivs set to enable that."
                                  ))
    }
  }
  ## if(!is.null(buildDerivs))
  ##   if(!isFALSE(buildDerivs)) {
  ##     code$name <- paste0(code$name, "_AD2_")
  ##     setArg(code, length(code$args) + 1, RparseTree2ExprClasses(quote(recordingInfo_)))
  ##   }
  invisible(NULL)
}

modifyForAD_getSetValues <- function(code, symTab, workEnv) {
  if(code$name == 'setValues') {
    code$name <- 'setValues_AD_AD_taping'
  }
  accessorName <- code$args[[2]]$name
  already_AD_name <- grepl("_AD_$", accessorName) ## _AD_ is at the end. I'm not sure it would ever already be there, so this is defensive.
  if(!already_AD_name) {
    classSymTab <- symTab$getParentST()$getParentST()
    newSymbolName <- paste0(accessorName, '_AD_')
    if(!classSymTab$symbolExists(newSymbolName)) {
      newSymbol <- classSymTab$getSymbolObject(accessorName)$copy()
      if(newSymbol$name != accessorName)
        warning("Something is wrong with internal processing of names in managing use of values() for nimDerivs.")
      newSymbol$name <- newSymbolName
      classSymTab$addSymbol(newSymbol)
    }
    code$args[[2]]$name <- newSymbolName
  }
  thirdArg <- RparseTree2ExprClasses(as.name(accessorName))
  setArg(code, 3, thirdArg)
  fourthArg <- RparseTree2ExprClasses(quote(recordingInfo_))
  setArg(code, 4, fourthArg)
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
  setArg(code, length(code$args) + 1, RparseTree2ExprClasses(quote(cppLiteral("recordingInfo_"))))
  invisible(NULL)
}

updateADproxyModelMethods <- function(.self) {
  ## Update return type and names of functions like dnorm -> nimDerivs_dnorm
  functionNames <- names(.self$functionDefs)
  ADproxyModel_functionNames <- functionNames[ grepl("_ADproxyModel", functionNames ) ]
  if(length(ADproxyModel_functionNames) > 0) {
    ## The following two headers were added to CPPincludes because nimDerives_dists must
    ## come before Rmath.h
    .self$CPPincludes <- c(nimbleIncludeFile("nimbleCppAD.h"),
                           nimbleIncludeFile("nimDerivs_dists.h"), .self$CPPincludes)
  }

  classST <- .self$objectDefs
  classSymNames <- classST$getSymbolNames()
  ADproxySymNames <- classSymNames[ grepl("ADproxyModel_", classSymNames ) ]
  for(sn in ADproxySymNames) {
    newSym <-
      cppVarSym2templateTypeCppVarSym(classST$getSymbolObject(sn),
                                      replacementBaseType = "CppAD::AD",
                                      replacementTemplateArgs = "double") ## These end up the header file so they should not use TYPE_ as that is not typedef'd there
    classST$addSymbol(newSym, allowReplace = TRUE)
  }

  for(fn in ADproxyModel_functionNames) {
    thisDef <- .self$functionDefs[[fn]]
    thisDef$returnType <- cppVarSym2templateTypeCppVarSym( thisDef$returnType,
                                                          replacementBaseType = "CppAD::AD",
                                                          replacementTemplateArgs = "double" )
    parentST <- thisDef$code$objectDefs$getParentST()
    thisDef$code$objectDefs <-
      symbolTable2templateTypeSymbolTable(thisDef$code$objectDefs) #,
    ##                                                replacementBaseType = "CppAD::AD",
    ##                                                replacementTemplateArgs = "double" )
    thisDef$code$objectDefs$setParentST(parentST)
    thisDef$code$cppADCode <- 2L
    ADtypeDefs <- symbolTable()
    ADtypeDefs$addSymbol(cppVarFull(baseType = "typedef typename EigenTemplateTypes<CppAD::AD<double> >::typeEigenMapStrd", name = "EigenMapStrd_TYPE_") )
    ADtypeDefs$addSymbol(cppVarFull(baseType = "typedef Matrix<CppAD::AD<double>, Dynamic, Dynamic>", name = "MatrixXd_TYPE_") )
    ADtypeDefs$addSymbol(cppVarFull(baseType = "typedef CppAD::AD<double>", name = "TYPE_") )
    thisDef$code$typeDefs <- ADtypeDefs
    workEnv <- new.env()
    workEnv$RsymTab <- thisDef$RCfunProc$compileInfo$newLocalSymTab
    workEnv$.inModel <- TRUE
    exprClasses_modifyForAD(thisDef$code$code, thisDef$code$objectDefs, workEnv)
    if(isTRUE(workEnv$.illegalADcall)) {
      handlingOption <- getNimbleOption("unsupportedDerivativeHandling")
      if(identical(handlingOption, "ignore")) {
        # do nothing
      } else if(identical(handlingOption, "warn")) {
        firstSym <- cppVarFull(name = "first", static = TRUE, baseType = "bool", constructor = "{true}")
        thisDef$code$objectDefs <- symbolTable$new()
        thisDef$code$typeDefs <- symbolTable$new()
        thisDef$code$objectDefs$setParentST(parentST)
        thisDef$code$objectDefs$addSymbol(firstSym)
        errorTrappingCode <- quote({
          if(first) {
            nimPrint("\n  [Warning] In C++, a model node is being included in derivatives that does not support them.\n",
                     "            This message will only be shown once. Results will not be valid.\n",
                     "            To make this an error, set `nimbleOptions(unsupportedDerivativeHandling='error')`.\n",
                     "            To skip this error handling entirely (and allow possible malformed code to be used),\n",
                     "            set `nimbleOptions(unsupportedDerivativeHandling='ignore')`.")
            first <- false
          }
          return(`CppAD::AD<double>`(0))
        })
        errorTrappingExpr <- RparseTree2ExprClasses(errorTrappingCode)
        thisDef$code$code <- errorTrappingExpr
      } else { #default to error
        thisDef$code$objectDefs <- symbolTable$new()
        thisDef$code$typeDefs <- symbolTable$new()
        thisDef$code$objectDefs$setParentST(parentST)
        msg <- paste0("\n  [Error] In C++, a model node is being included in derivatives that does not support them.\n",
                     "          To make this a warning, set `nimbleOptions(unsupportedDerivativeHandling='warn')`.\n",
                     "          To skip this error handling entirely (and allow possible malformed code to be used),\n",
                     "          set `nimbleOptions(unsupportedDerivativeHandling='ignore')`.")
        errorTrappingCode <- substitute(
        {
          nimStop(MSG)
          return(`CppAD::AD<double>`(0))
        }, list(MSG = msg))
        errorTrappingExpr <- RparseTree2ExprClasses(errorTrappingCode)
        thisDef$code$code <- errorTrappingExpr
      }
    }
  }
  ## classST <- .self$objectDefs
  ## classSymNames <- classST$getSymbolNames()
  ## ADproxySymNames <- classSymNames[ grepl("ADproxyModel_", classSymNames ) ]
  ## for(sn in ADproxySymNames) {
  ##   newSym <-
  ##     cppVarSym2templateTypeCppVarSym(classST$getSymbolObject(sn),
  ##                                     replacementBaseType = "CppAD::AD",
  ##                                     replacementTemplateArgs = "double") ## These end up the header file so they should not use TYPE_ as that is not typedef'd there
  ##   classST$addSymbol(newSym, allowReplace = TRUE)
  ## }
  NULL
}

makeSingleCopyCall <- function(varName, cppCopyType,
                               buildDerivs = FALSE) {
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
               cppLiteral(paste0("COPY_VALUE_MAP_ACCESSORS_FROM_NODE_NAMES(\"", varName, "\"",
                                 if(isTRUE(buildDerivs)) ", 1" else ", 0", ");"))
           },
           NULL)
}

makeCopyFromRobjectDef <- function(cppCopyTypes,
                                   buildDerivs = FALSE) {
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
      copyCalls[[varNames[i]]] <- makeSingleCopyCall(varNames[i], cppCopyTypes[[i]],
                                                     buildDerivs = buildDerivs)
    }

    if(length(copyCalls) == 0) {
        ## create an empty function
        allRcode <- do.call('call',
                            list('{'),
                            quote = TRUE
                            )
    } else {
        unprotectCount <- 2 + length(copyCalls) ## 2 from SETUP_S_xData
        allRcode <- do.call('call',
                            c(list('{'),
                              list(cppLiteral("SETUP_S_xData;")),
                              copyCalls,
                              list(cppLiteral(paste0("UNPROTECT(",unprotectCount,"+1);")))),
                            quote = TRUE
                            )
    }
    allCode <- RparseTree2ExprClasses(allRcode)
    copyFromRobjectDef$code <- cppCodeBlock(code = allCode,
                                            objectDefs = localVars)

    copyFromRobjectInterfaceDef <- NULL
    
    list(copyFromRobjectDef = copyFromRobjectDef,
         copyFromRobjectInterfaceDef = copyFromRobjectInterfaceDef)
}

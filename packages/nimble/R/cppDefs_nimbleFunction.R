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
                                                 ))

## cppNimbleClassClass defines commonalities between cppNimbleFunctionClass and cppNimbleListClass, both of which are classes in nimble
cppNimbleClassClass <- setRefClass('cppNimbleClassClass',
                                   contains = 'cppNamedObjectsClass',
                                   fields = list(
                                       ## Inherits a functionDefs list for member functions
                                       ## Inherits an objectDefs list for member data
                                       SEXPmemberInterfaceFuns = 'ANY', ## List of SEXP interface functions, one for each member function
                                       nimCompProc = 'ANY', ## an nfProcessing or nlProcessing class, needed to get the member data symbol table post-compilation

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
                                                 thisCppDef <- nimbleProject$getNimbleFunctionCppDef(generatorName = generatorName)
                                                 eigListClassNames <- sapply(nlEigenReferenceList, function(x){return(x$className)})
                                                ##EIGEN_EIGENCLASS and EIGEN_SVDCLASS only need partial cpp defs, as some is in permanent c++ code
                                                if(is.null(thisCppDef) && (generatorName %in% eigListClassNames)){
                                                    className <- names(nimCompProc$neededTypes)[i]
                                                    nimbleProject$buildNimbleListCompilationInfo(className = generatorName, fromModel = fromModel, eigenList = TRUE)
                                                 }
                                                 else if(is.null(thisCppDef)) {
                                                   className <- names(nimCompProc$neededTypes)[i]
                                                   thisCppDef <- nimbleProject$buildNimbleListCompilationInfo(className = generatorName, fromModel = fromModel)
                                                   neededTypeDefs[[ className ]] <<- thisCppDef
                                                   Hincludes <<- c(Hincludes, thisCppDef)
                                                   CPPincludes <<- c(CPPincludes, thisCppDef)
                                                 }
                                                 next
                                               }
                                               if(inherits(neededType, 'symbolOptimReadyFunction')){
                                                   typeName <- neededType$genName
                                                   thisCppDef <- nimbleProject$getNimbleFunctionCppDef(generatorName = typeName)                                                     	 
                                                   if(is.null(thisCppDef)){
                                                       thisCppDef <- cppOptimObject(name = typeName, nfSym = neededType)
                                                       neededTypeDefs[[ typeName ]] <<- thisCppDef
                                                   } else {
                                                       Hincludes <<- c(Hincludes, thisCppDef)
                                                       CPPincludes <<- c(CPPincludes, thisCppDef)
                                                   }
                                                   
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
                                          nfProc = 'ANY' ## an nfProcessing class, needed to get the member data symbol table post-compilation
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
                                                  SEXPmemberInterfaceFuns <<- lapply(functionDefs, function(x) x$SEXPinterfaceFun)
                                                  nimCompProc <<- nfProc
                                              },
                                              buildFunctionDefs = function() {
                                                  for(i in seq_along(nfProc$RCfunProcs)) {
                                                      RCname <- names(nfProc$RCfunProcs)[i]
                                                      functionDefs[[RCname]] <<- RCfunctionDef$new() ## all nodeFunction members are const f
                                                      functionDefs[[RCname]]$buildFunction(nfProc$RCfunProcs[[RCname]])
                                                      functionDefs[[RCname]]$buildSEXPinterfaceFun(className = nfProc$name)
                                                      RCfunDefs[[RCname]] <<- functionDefs[[RCname]]
                                                  }
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
                                              promoteCallable = function(R_NimbleFxn, asTopLevel = TRUE){
                                                  ## see comment in nimbleProjectClass::compileNimbleFunction
                                                  oldCobjectInterface <- nf_getRefClassObject(R_NimbleFxn)$.CobjectInterface
                                                  if(!is.list(oldCobjectInterface)) return(oldCobjectInterface)
                                                  existingBasePtr <- oldCobjectInterface[[1]]$getBasePtr(oldCobjectInterface[[2]])
                                                  thisDll <- oldCobjectInterface[[1]]$dll
                                                  newCobjectInterface <- Rgenerator(R_NimbleFxn, thisDll, project = nimbleProject, existingBasePtr = existingBasePtr)
                                                  newCobjectInterface
                                              },
                                              buildCallable = function(R_NimbleFxn, dll = NULL, asTopLevel = TRUE){
                                             ##     cat('buildCallable\n')
                                                  if(asTopLevel) {
                                                      cppInterfaceObject <- Rgenerator(R_NimbleFxn, dll, project = nimbleProject)
                                                  } else { ## actually this particular pathway should never be taken. asTopLevel = FALSE will occur only for nimbleProject$instantiateNimbleFunction 
                                                      cppInterfaceObject <- CmultiInterface$addInstance(R_NimbleFxn, dll)
                                                  }
                                                  return(cppInterfaceObject)
                                              },
                                              buildAll = function(where = where) {
                                                  baseClassObj <- environment(nfProc$nfGenerator)$contains
                                                  if(!is.null(baseClassObj)) {
                                                      inheritance <<- inheritance[inheritance != 'NamedObjects']
                                                      baseClassName <- environment(baseClassObj)$name
                                                      addInheritance(baseClassName)
                                                  }
                                                  callSuper(where)
                                              }
                                          ),
                                      )

## This and the functions below were original stand-along compilation functions
## Now they are superceded by compileNimble.
## A stand-alone call could still work but these may have faller out of date.
compileNimbleFunction <- function(fun, dirName, all.instances = TRUE, individual.instances = NA, name = deparse(substitute(fun) ),  fileName = Rname2CppName(name), writeFiles = !(environment(fun)$Cwritten), compileCpp = !(environment(fun)$compiled), loadSO = !(environment(fun)$loadedSO), debug = FALSE, debugCpp = FALSE, returnInternals = FALSE ) {
    if(missing(dirName))    dirName <- makeDefaultDirName()
    
        if(!is.nfGenerator(fun)) stop('Error in compileNimbleFunction: fun should be a nimbleFunction generator')
	cppProj <- nfWriteCompileAndLoadSO(RFun = fun, dirName = dirName, name = name, fileName = fileName,
                                           writeFiles = writeFiles, compileCpp = compileCpp, loadSO = loadSO, debug = debug, debugCpp = debugCpp)

        if(returnInternals) return(cppProj)
        
	if(all.instances == FALSE){
            if(!is.list(individual.instances))
                return(nfBuildCInterface(cppProj, instance = individual.instances, name = name) ) 
            else{
                functionList = list()
                for(i in 1:length(individual.instances) ) 
                    functionList[[i]] <- nfBuildCInterface(cppProj, instance = individual.instances[[i]], name = name)
            }
            return(functionList)
	}

	functionList = list()
	k = length(environment(fun)$instances)
	if(k == 1)
            {
		nfRefClassObject <- environment(fun)$instances[[1]]
		wrappedFunction <- function(...) { nfRefClassObject$run() }
		environment(wrappedFunction) <- new.env(parent = parent.frame())
		environment(wrappedFunction)$nfRefClassObject <- nfRefClassObject
		return(nfBuildCInterface(cppProj, instance = wrappedFunction, name = name) )
            }
	
	for(i in 1:k) {
            nfRefClassObject <- environment(fun)$instances[[i]]
            wrappedFunction <- function(...) { nfRefClassObject$run() }
            environment(wrappedFunction) <- new.env(parent = parent.frame())
            environment(wrappedFunction)$nfRefClassObject <- nfRefClassObject
            functionList[[i]] <- nfBuildCInterface(cppProj, instance = wrappedFunction, name = name)
	}
	return(functionList)		
    }

## nf should be a nfProcessing object OR a nimbleFunction (generator or function).
makeCppNIMBLEfunction <- function(nf, name, debug = FALSE, debugCpp = FALSE, isNode = FALSE, where = globalenv()) { 
    if(inherits(nf, 'nfProcessing') | inherits(nf, 'virtualNFprocessing')) {
        if(!missing(name)) {
            writeLines(paste('Warning, name', name, 'will be overwritten by the name of the nfProcessing object', nfp$name))
        }
        name <- nf$name
        nfp <- nf
    } else {
        if(missing(name)) name <- deparse(substitute(nf))
    }
    Cname <- Rname2CppName(name)
    if(is.nf(nf) | is.nfGenerator(nf)) {
        virtual <- FALSE
        if(is.nfGenerator(nf)) {
            if(environment(nf)$virtual) virtual <- TRUE
        }
        if(virtual)
            nfp <- virtualNFprocessing$new(nf, Cname)
        else {
            genFun <- nf_getGeneratorFunction(nf)
            nfp <- environment(genFun)$nfProc
            if(is.null(nfp)) {
                nfp <- nfProcessing$new(nf, Cname)
                environment(genFun)$nfProc <- nfp
            } else {
                nfp$setClassName(Cname)
            }
        }
        nfp$process(debug = debug, debugCpp = debugCpp)
    }

    if(inherits(nfp, 'nfProcessing')) {
        ans <- cppNimbleFunctionClass$new(name = Cname, nfProc = nfp, isNode = isNode, debugCpp = debugCpp) 
        ans$buildAll(where = where)
    } else {
        if(inherits(nfp, 'virtualNFprocessing'))
            ans <- cppVirtualNimbleFunctionClass(name = Cname, nfProc = nfp)
        else
            stop('Some problem in makeCppNIMBLEfunction')
    }
    
    ans
}

nfWriteCompileAndLoadSO <- function(RFun, dirName, name = deparse(substitute(RFun)), fileName = Rname2CppName(name),
                                    writeFiles = TRUE, compileCpp = TRUE, loadSO = TRUE, debug = FALSE, debugCpp = FALSE){
    Cname = Rname2CppName(name)
    NFC <- makeCppNIMBLEfunction(RFun, Cname, debug = debug, debugCpp = debugCpp)
    cppProj <- cppProjectClass$new(dirName = dirName)
    cppProj$addClass(NFC, Cname)

    if(writeFiles){
        cppProj$writeFiles(fileName)
        environment(RFun)$Cwritten <- TRUE
    }
    if(compileCpp) {
        cppProj$compileFile(fileName)
        environment(RFun)$compiled <- TRUE
    }
    if(loadSO) {
        cppProj$loadSO(fileName)
        environment(RFun)$loadedSO <- TRUE
    }
    return(cppProj)
}

nfBuildCInterface <- function(cppProj, instance, name = deparse(substitute(RFun) )){
	return( cppProj$cppDefs[[name]]$buildCallable(instance, cppProj$dll ) )
}


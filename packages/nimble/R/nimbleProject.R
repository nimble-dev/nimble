##.nimbleProjectClassMasterList <- new.env(, emptyenv())

projectNameCreator <- labelFunctionCreator('P')

nfCompilationInfoClass <- setRefClass('nfCompilationInfoClass',
                                      fields = list(
                                          nfProc = 'ANY',      ## an nfProcessing object 
                                          nfGenerator = 'ANY', ## a nfGenerator, which is a function with special stuff in its environment
                                          cppDef = 'ANY',       ## a cppNimbleFunctionClass object
                                          labelMaker = 'ANY',    ## a label maker function
                                          virtual =  'ANY',		#'logical',
                                          RinitTypesProcessed = 'ANY',		# 'logical', ## setupTypesForUsingFunction() 
                                          Rcompiled =  'ANY',		#'logical',
                                          written =  'ANY',		#'logical',
                                          cppCompiled =  'ANY',		#'logical',
                                          loaded =  'ANY',		#'logical',
                                          fromModel =  'ANY',		#'logical',
                                          Rinstances =  'ANY'		#'list'
                                          ),
                                      methods = list(
                                      		initialize = function(...){Rinstances <<- list(); callSuper(...)},
                                          addRinstance = function(nfi) {Rinstances[[ length(Rinstances)+1 ]] <<- nfi}
                                          ))

mvInfoClass <- setRefClass('mvInfoClass',
                           fields = list(
                               mvSpec = 'ANY', ## a custom modelValues class
                               cppClassName =  'ANY',		#'character',
                               cppClass = 'ANY', ## a cppModelValuesClass object,
                               fromModel =  'ANY',		#'logical',
                               RmvObjs =  'ANY'		#'list'
                               ),
                           methods = list(
                           		initialize = function(...){RmvObjs <<- list(); callSuper(...)},
                               addRmv = function(Rmv) RmvObjs[[length(RmvObjs)+1]] <<- Rmv))

RCfunInfoClass <- setRefClass('RCfunInfoClass',
                              fields = list(
                                  nfMethodRCobj = 'ANY', ## an mfMethodRC
                                  RCfunProc     = 'ANY', ## an RCfunProcessing or NULL
                                  cppClass      = 'ANY',  ## an RCfunctionDef or NULL
                                  fromModel     =  'ANY'		#'logical'
                                  ))

modelDefInfoClass <- setRefClass('modelDefInfoClass',
                                 fields = list(
                                     labelMaker = 'ANY'
                                     ))

nimbleProjectClass <- setRefClass('nimbleProjectClass',
                             fields = list(
                                 RCfunInfos         =  'ANY',		#'list', ## a list of RCfunInfoClass objects
                                 RCfunCppInterfaces =  'ANY',		#'list', 
                                 mvInfos            =  'ANY',		#'list', ## a list of mvInfoClass objects
                                 modelDefInfos      =  'ANY',		#'list',
                                 modelCppInterfaces =  'ANY',		#'list',
                                 models             =  'ANY',		#'list',
                                 nimbleFunctions    =  'ANY',		#'list',
                                 nimbleFunctionCppInterfaces =  'ANY',		#'list',
                                 nfCompInfos        =  'ANY',		#'list', ## list of nfCompilationInfoClass objects
                                 cppProjects        =  'ANY',		#'list', ## list of cppProjectClass objects, 1 for each dll to be produced
                                 dirName            =  'ANY',		#'character',
                                 nimbleLabel        =  'ANY',		#'character',
                                 refClassDefsEnv    =  'ANY',		#'environment',
                                 projectName        =  'ANY'		#'character'
                                 ),
                             methods = list(
                                 show = function() {
                                     writeLines(paste0('nimbleProject object'))
                                 },
                                 initialize = function(dir = NULL, name = '') {
                                 	RCfunInfos <<- list()
                                 	RCfunCppInterfaces <<- list()
                                 	mvInfos <<- list()
                                 	modelDefInfos <<- list()
                                 	modelCppInterfaces <<- list()
                                 	models <<- list()
                                 	nimbleFunctions <<- list()
                                 	nimbleFunctionCppInterfaces <<- list()
                                 	nfCompInfos <<- list()
                                 	cppProjects <<- list()
                                 	refClassDefsEnv <<- new.env()
                                     dirName <<- if(is.null(dir)) makeDefaultDirName() else dir
                                     if(name == '') projectName <<- projectNameCreator() else projectName <<- name
                                  },
                                 resetFunctions = function() {
                                     ## clear everything except models and nimbleFunctions from models
                                     for(i in names(mvInfos)) {
                                         if(length(mvInfos[[i]]$fromModel) > 0) {
                                             if(!mvInfos[[i]]$fromModel) {
                                                 mvInfos[[i]]$cppClass <<- NULL
                                                 for(j in seq_along(mvInfos[[i]]$RmvObjs))
                                                     mvInfos[[i]]$RmvObjs[[j]]$CobjectInterface <<- NULL
                                                 mvInfos[[i]] <<- NULL
                                             }
                                         }
                                     }
                                     for(i in names(RCfunInfos)) {
                                         if(length(RCfunInfos[[i]]$fromModel) > 0) {
                                             if(!RCfunInfos[[i]]$fromModel) {
                                                 thisName <- RCfunInfos[[i]]$uniqueName
                                                 cppProjects[[thisName]] <<- NULL
                                                 RCfunInfos[[i]] <<- NULL
                                                 RCfunCppInterfaces[[i]] <<- NULL
                                             }
                                         }
                                     }
                                     for(i in names(nfCompInfos)) {
                                         if(length(nfCompInfos[[i]]$fromModel) > 0) {
                                             if(!nfCompInfos[[i]]$fromModel) {
                                                 for(j in seq_along(nfCompInfos[[i]]$Rinstances)) {
                                                     thisEnv <- environment(nfCompInfos[[i]]$Rinstances[[j]])
                                                     thisRCO <- nf_getRefClassObject(nfCompInfos[[i]]$Rinstances[[j]])
                                                     if(exists('name', envir = thisRCO, inherits = FALSE)) {
                                                         thisname <- thisRCO$name
                                                         nimbleFunctions[[ thisname ]] <<- NULL
                                                         nimbleFunctionCppInterfaces[ thisname ] <<- NULL
                                                         rm('name', envir = thisRCO)
                                                     }
                                                     thisRCO[['nimbleProject']] <- NULL
                                                     thisRCO$.CobjectInterface <- NULL
##                                                     if(exists('.CobjectInterface', envir = thisEnv, inherits = FALSE)) {
##                                                         rm('.CobjectInterface', envir = thisEnv)
##                                                     }
                                                 }
                                                 nfCompInfos[[i]] <<- NULL
                                                 cppProjects[[i]] <<- NULL
                                             }
                                         }
                                     }
                                 },
                                     
                                 addModelValuesClass = function(mvSpec, fromModel = FALSE) {
                                     mvClassName <- environment(mvSpec)$className
                                     if(!is.null(mvInfos[[mvClassName]])) stop('Trying to add a modelValues class with the same name as one already in this project', call. = FALSE)
                                     mvInfos[[mvClassName]] <<-  mvInfoClass(cppClassName = mvClassName, cppClass = NULL, mvSpec = mvSpec, fromModel = fromModel)
                                 },
                                 getModelValuesCppDef = function(mvSpec, NULLok = FALSE) {
                                     mvClassName <- environment(mvSpec)$className
                                     if(is.null(mvInfos[[mvClassName]])) {
                                         if(!NULLok) stop('Project does not know about this modelValues class but the cppDef is being requested', call. = FALSE)
                                         else return(NULL)
                                     }
                                     ans <- mvInfos[[mvClassName]]$cppClass
                                     if(is.null(ans)) {
                                         if(!NULLok) stop('cppDef for a requested modelValues class has not been generated but was requested.', call. = FALSE)
                                         else return(NULL)
                                     }
                                     ans
                                 },
                                 addModel = function(model) {
                                     if(!inherits(model, 'RModelBaseClass')) stop('model provided to project is not an RModelBaseClass', call. = FALSE)
                                     modelDefName <- model$modelDef$name
                                     if(is.null(modelDefInfos[[modelDefName]])) {
                                         modelDefInfos[[modelDefName]] <<- modelDefInfoClass(
                                             labelMaker = labelFunctionCreator(paste0(modelDefName, '_'))
                                             )
                                         if(identical(model$name, character(0))) {
                                             model$name <- modelDefInfos[[modelDefName]]$labelMaker()
                                         } else {
                                             if(!is.null(models[[ model$name ]])) stop('Model provided to project has same name as another one in the same project')
                                         }
                                         if(identical(model$Cname, character(0))) {
                                             model$Cname <- Rname2CppName(model$name)
                                         } 
                                         model$nimbleProject <- .self
                                         models[[ model$name ]] <<- model
                                         modelCppInterfaces[ model$name ] <<- list(NULL)
                                     }
                                 },
                                 addNimbleFunction = function(fun, fromModel = FALSE) {
                                     if(!is.nf(fun)) stop('nimbleFunction provided to project is not a nimbleFunction.', call. = FALSE)
                                     inProjectAlready <- nf_getRefClassObject(fun)[['nimbleProject']]
                                     if(!is.null(inProjectAlready)) {
                                         if(!identical(inProjectAlready, .self)) stop('Trying to add a specialized nimbleFunction to a project but it is already part of another project. If you are recompiling, try redefining models and specialized nimbleFunctions. (The reset option works now for nimbleFunctions but not models.)', call. = FALSE)
                                         else warning('Adding a specialized nimbleFunction to a project to which it already belongs', call. = FALSE)
                                     }
                                     generatorName <- nfGetDefVar(fun, 'name')
                                     if(is.null(nfCompInfos[[generatorName]])) {
                                         ## nfProc could have been created already during makeTypeObject for another nimbleFunction so it knows the types of this one.
                                         nfCompInfos[[generatorName]] <<- nfCompilationInfoClass(nfGenerator = nf_getGeneratorFunction(fun),
                                                                                                 Rcompiled = FALSE, written = FALSE, cppCompiled = FALSE, loaded = FALSE,
                                                                                                 RinitTypesProcessed = FALSE, virtual = FALSE,
                                                                                                 fromModel = fromModel)
                                         nfCompInfos[[generatorName]]$labelMaker <<- labelFunctionCreator(paste0(generatorName,'_'))
                                     }
                                     if(!exists('name', envir = nf_getRefClassObject(fun), inherits = FALSE)) {
                                         assign('name', nfCompInfos[[generatorName]]$labelMaker(), envir = nf_getRefClassObject(fun))
                                     } else {
                                         if(!is.null(nimbleFunctions[[ nf_getRefClassObject(fun)$name ]])) {
                                             stop('nimbleFunction provided to project has same name as another one in the same project', call. = FALSE)
                                         }
                                     }
                                     nimbleFunctions[[ nf_getRefClassObject(fun)$name ]] <<- fun
                                     nfCompInfos[[generatorName]]$addRinstance(fun)
                                     nimbleFunctionCppInterfaces[ nf_getRefClassObject(fun)$name ] <<- list(NULL)
                                    
                                     if(!exists('Cname', envir = nf_getRefClassObject(fun), inherits = FALSE)) {
                                         assign('Cname', Rname2CppName(nf_getRefClassObject(fun)$name), envir = nf_getRefClassObject(fun))
                                     }

                                     assign('nimbleProject', .self, envir = nf_getRefClassObject(fun))
                                     ## could check for duplicate Cnames here, but if the names are unique the Cnames should be too.
                                 },
                                 addRCfun = function(nfmObj, fromModel = FALSE) {
                                     if(!inherits(nfmObj, 'nfMethodRC')) stop("Can't add this function. nfmObj is not an nfMethodRC", call. = FALSE)
                                     className <- nfmObj$uniqueName
                                     if(is.null(RCfunInfos[[className]])) {
                                         RCfunInfos[[className]] <<- RCfunInfoClass(nfMethodRCobj = nfmObj, RCfunProc = NULL, cppClass = NULL, fromModel = fromModel)
                                     }
                                 },
                                 getRCfunCppDef = function(nfmObj, NULLok = FALSE) {
                                     className <- nfmObj$uniqueName
                                     ans <- RCfunInfos[[className]]
                                     if(is.null(ans)) {
                                         if(!NULLok) stop("Requested to get an RCfunCppDef but it is not in the project and NULLok = FALSE", call. = FALSE)
                                     }
                                     ans
                                 },
                                 needRCfunCppClass = function(nfmObj, genNeededTypes = TRUE, control = list(debug = FALSE, debugCpp = FALSE), fromModel = FALSE) {
                                     if(!inherits(nfmObj, 'nfMethodRC')) stop("Can't compile this function. nfmObj is not an nfMethodRC", call. = FALSE)
                                     className <- nfmObj$uniqueName
                                     Cname <- Rname2CppName(className)
                                     RCfunInfo <- RCfunInfos[[className]]
                                     if(is.null(RCfunInfo)) addRCfun(nfmObj, fromModel = fromModel)
                                     if(is.null(RCfunInfos[[className]]$RCfunProc)) {
                                         RCfunInfos[[className]]$RCfunProc <<- RCfunProcessing$new(nfmObj, Cname)
                                         RCfunInfos[[className]]$RCfunProc$process(debug = control$debug, debugCpp = control$debugCpp)
                                     }
                                     cppClass <- RCfunInfos[[className]]$cppClass
                                     if(is.null(cppClass)) {
                                         cppClass <- RCfunctionDef(project = .self)
                                         cppClass$buildFunction(RCfunInfos[[className]]$RCfunProc)
                                         cppClass$buildSEXPinterfaceFun()
                                         if(genNeededTypes) cppClass$genNeededTypes()
                                         RCfunInfos[[className]]$cppClass <<- cppClass
                                     }
                                     cppClass
                                 },
                                 compileRCfun = function( fun, filename = NULL, control = list(debug = FALSE, debugCpp = FALSE, writeFiles = TRUE, returnAsList = FALSE ) ) {
                                     if(is.rcf(fun)) fun <- environment(fun)$nfMethodRCobject
                                     addRCfun(fun) ## checks if it already exists and if it is an nfMethodRC
                                     cppClass <- needRCfunCppClass(fun, genNeededTypes = TRUE, control = control)
                                     className <- fun$uniqueName
                                     if(control$writeFiles) {
                                         cppProj <- cppProjectClass(dirName = dirName)
                                         cppProjects[[ className ]] <<- cppProj
                                         if(is.null(filename)) filename <- paste0(projectName, '_', className)
                                         cppProj$addClass( cppClass, className, filename )
                                         cppProj$writeFiles(filename)
                                     }
                                     if(control$compileCpp) {
                                         cppProj$compileFile(filename)
                                     }
                                     if(control$loadSO) {
                                         cppProj$loadSO(filename)
                                     }
                                     RCfunCppInterfaces[[className]] <<- cppClass$buildRwrapperFunCode(includeLHS = FALSE, eval = TRUE, returnArgsAsList = control$returnAsList, dll = cppProj$dll)
                                     RCfunCppInterfaces[[className]]
                                 },
                                 needModelValuesCppClass = function(mvSpec, fromModel = FALSE) {
                                     if(!isModelValuesSpec(mvSpec)) stop("Can't compileModelValues: mvSpec is not a modelValuesSpec", call. = FALSE)
                                     mvClassName <- environment(mvSpec)$className
                                     mvInfo <- mvInfos[[mvClassName]]
                                     if(is.null(mvInfo)) addModelValuesClass(mvSpec, fromModel)
                                     cppClass <- mvInfos[[mvClassName]]$cppClass
                                     if(is.null(cppClass)) {
                                         cppClass <- cppModelValuesClass(name = mvClassName,
                                                                            vars = environment(mvSpec)$symTab,
                                                                            project = .self)
                                         ##vars = nimbleProject$modelValuesLibrary[[className]]$symTab                               
                                         cppClass$buildAll()
                                          mvInfos[[mvClassName]]$cppClass <<- cppClass
                                         ##nimbleProject$modelValuesLibrary[[className]]$cppClass <<- newCppClass
                                         ##neededTypeDefs[[className]] <<- newCppClass
                                     }
                                     cppClass
                                 },
                                 instantiateCmodelValues = function(mv, dll) {
                                     mvClassName <- class(mv)
                                     cppDef <- mvInfos[[mvClassName]]$cppClass
                                     if(is.null(cppDef)) stop('Trying to instantiate a modelValues type that the project has no record of.')
                                     generatorName <- cppDef$SEXPgeneratorFun$name
                                     sym = if(!is.null(dll))
                                         getNativeSymbolInfo(generatorName, dll)
                                     else {
                                         warning('a nimbleFunctionInterface is about to build a CmodelValues without dll info, based on generatorFun name only.', call. = FALSE)
                                         generatorName
                                     }
                                     ans <- CmodelValues(sym)
                                     mvInfos[[mvClassName]]$addRmv(mv) ## simply a list for later clearing
                                     mv$CobjectInterface <- ans
                                     ans
                                 },
                                 compileModel = function(model, filename = NULL, control = list(debug = FALSE, debugCpp = FALSE), where = globalenv()) {
                                     if(is.character(model)) {
                                         tmp <- models[[model]]
                                         if(is.null(tmp)) stop(paste0("Model provided as name: ", model, " but it is not in this project."), call. = FALSE)
                                         model <- tmp
                                     } else addModel(model)
                                                                          
                                     modelDef <- model$getModelDef()
                                     modelDefName <- modelDef$name
                                     Cname <- Rname2CppName(modelDefName)
                                     if(is.null(filename)) filename <- paste0(projectName, '_', Rname2CppName(modelDefName)) 
                                     
                                     modelCpp <- cppBUGSmodelClass(modelDef = modelDef, model = model,
                                                                   name = Cname, project = .self)
                                     ## buildAll will call back to the project to add its nimbleFunctions 
                                     modelCpp$buildAll(buildNodeDefs = TRUE, where = where, debugCpp = control$debugCpp)
                                     
                                     cppProj <- cppProjectClass(dirName = dirName)
                                     cppProjects[[ modelDefName ]] <<- cppProj
                                     ## genModelValuesCppClass will back to the project to add its mv class
                                     mvc <- modelCpp$genModelValuesCppClass()
                                     if(is.null(filename)) filename <- paste0(projectName, '_', modelDefName)
                                     cppProj$addClass(mvc, filename = filename)
                                     cppProj$addClass(modelCpp, modelDefName, filename)
                                     ##if compileNodes
                                     nfFileName <- paste0(projectName, '_', Rname2CppName(modelDefName),'_nfCode')
                                     for(i in names(modelCpp$nodeFuns)) {
                                         cppProj$addClass(modelCpp$nodeFuns[[i]], filename = nfFileName)
                                     }
                                     ##if writeFiles
                                     cppProj$writeFiles(filename)
                                     cppProj$writeFiles(nfFileName) ## if compileNodes
                                     ##if compileCpp
                                     compileList <- filename
                                     compileList <- c(compileList, nfFileName) ## if compileNodes
                                     cppProj$compileFile(compileList)
                                     ## if loadSO
                                     cppProj$loadSO(filename)
                                     ## if buildInterface
                                     interfaceName <- paste0('C', modelDefName)
                                     
                                     compiledModel <- modelCpp ## cppProj$cppDefs[[2]]
                                     newCall <- paste0('new_',Rname2CppName(modelDefName))
                                     ans <- buildModelInterface(interfaceName, compiledModel, newCall, where = where, project = .self, dll = cppProj$dll)
                                     createModel <- TRUE
                                     if(!createModel) return(ans) else return(ans(model, where, dll = cppProj$dll))
                                 },
                                 getNimbleFunctionCppDef = function(generatorName, nfProc) {
                                     if(missing(generatorName)) {
                                         if(missing(nfProc)) stop('No good information provided to getNimbleFunctionCppDef', call. = FALSE)
                                         generatorName <- environment(nfProc$nfGenerator)[['name']]
                                         if(is.null(generatorName)) stop('Invalid generatorName', call. = FALSE)
                                     }
                                     if(is.null(nfCompInfos[[generatorName]])) return(NULL)
                                     ans <- nfCompInfos[[generatorName]]$cppDef
                                     if(inherits(ans, 'uninitializedField')) return(NULL)
                                     ans
                                 },
                                 getNimbleFunctionNFproc = function(fun) {
                                     generatorName <- nfGetDefVar(fun, 'name')
                                     if(is.null(nfCompInfos[[generatorName]])) return(NULL)
                                     ans <- nfCompInfos[[generatorName]]$nfProc
                                     if(inherits(ans, 'uninitializedField')) return(NULL)
                                     ans
                                 },
                                 buildVirtualNimbleFunctionCompilationInfo = function(vfun, initialTypeInferenceOnly = FALSE, control = list(debug = FALSE, debugCpp = FALSE)) {
                                     if(!is.character(vfun)) {
                                         if(!is.nfGenerator(vfun)) stop("Something provided as a nimbleFunctionVirtual does not appear to be correct.", call. = FALSE)
                                         if(!(environment(vfun)$virtual)) stop("Something provided as a nimbleFunctionVirtual is an nfGenerator but not a virtual one.", call. = FALSE)
                                         if(initialTypeInferenceOnly) stop("Can't do initialTypeInferenceOnly on a virtualNimbleFunction", call. = FALSE)
                                         generatorName <- environment(vfun)$name
                                     } else {
                                         generatorName <- vfun
                                     }
                                     if(is.null(nfCompInfos[[generatorName]])) stop("It doesn't look like nfCompInfos was set up for this generator.  Call setupVirtualNimbleFunction first.", call. = FALSE) 
                                     if(inherits(nfCompInfos[[generatorName]]$nfProc, 'uninitializedField')) {## might always be FALSE by this point in processing
                                         if(is.character(vfun)) stop("vfun given as character but nfProc doesn't exist yet", call. = FALSE)
                                         nfCompInfos[[generatorName]]$nfProc <<- virtualNFprocessing$new(vfun, generatorName)
                                     }
                                     if(!nfCompInfos[[generatorName]]$Rcompiled) {
                                         nfCompInfos[[generatorName]]$nfProc$process(control = control)
                                         nfCompInfos[[generatorName]]$Rcompiled <<- TRUE
                                     }
                                     if(inherits(nfCompInfos[[generatorName]]$cppDef, 'uninitializedField')) {
                                         newCppClass <- cppVirtualNimbleFunctionClass(name = generatorName,
                                                                                      nfProc = nfCompInfos[[generatorName]]$nfProc,
                                                                                      project = .self)
                                         nfCompInfos[[generatorName]]$cppDef <<- newCppClass
                                         newCppClass ## possible return value
                                     } else {
                                         nfCompInfos[[generatorName]]$cppDef ## return value if already exists
                                     }
                                 },
                                 buildNimbleFunctionCompilationInfo = function(funList = NULL, generatorName, initialTypeInferenceOnly = FALSE,
                                     isNode = FALSE, control = list(debug = FALSE, debugCpp = FALSE), where = globalenv(), fromModel = FALSE) {
                                     ## like old makeCppNIMBLEfunction
                                     ## check of make new nfCompInfos item
                                     ## ensure it is build up to the cppNimbleFunctionClass
                                     if(!is.null(funList)) {
                                         generatorName <- nfGetDefVar(funList[[1]], 'name')
                                         name <- nf_getRefClassObject(funList[[1]])$name
                                         Cname <- nf_getRefClassObject(funList[[1]])$Cname
                                         if(is.null(nfCompInfos[[generatorName]])) stop("Requested buildNimbleFunctionCompilationInfo for a generator for which no specialized NF has been added to the project", call. = FALSE)
                                         if(inherits(nfCompInfos[[generatorName]]$nfProc, 'uninitializedField')) 
                                             nfCompInfos[[generatorName]]$nfProc <<- nfProcessing(funList, generatorName, fromModel = fromModel, project = .self)
                                     } else {
                                         if(missing(generatorName)) stop("If funList is omitted, a generator name must be provided to buildNimbleFunctionCompilationInfo", call. = FALSE)
                                         if(inherits(nfCompInfos[[generatorName]]$nfProc, 'uninitializedField')) stop("buildNimbleFunctionCompilationInfo was called with only a generatorName (probably from genNeededTypes), but the nfProc is missing.", call. = FALSE)
                                     }
                                     if(initialTypeInferenceOnly) {
                                         if(!nfCompInfos[[generatorName]]$RinitTypesProcessed) {
                                             nfCompInfos[[generatorName]]$nfProc$setupTypesForUsingFunction() 
                                             nfCompInfos[[generatorName]]$RinitTypesProcessed <<- TRUE
                                         }
                                         return(nfCompInfos[[generatorName]]$nfProc)
                                     }
                                     if(!nfCompInfos[[generatorName]]$Rcompiled) {
                                         nfCompInfos[[generatorName]]$nfProc$process(control = control)
                                         nfCompInfos[[generatorName]]$Rcompiled <<- TRUE
                                     }
                                     if(inherits(nfCompInfos[[generatorName]]$cppDef, 'uninitializedField')) {
                                         newCppClass <- cppNimbleFunctionClass(name = generatorName,
                                                                               nfProc = nfCompInfos[[generatorName]]$nfProc,
                                                                               isNode = isNode,
                                                                               debugCpp = control$debugCpp,
                                                                               project = .self,
                                                                               fromModel = fromModel
                                                                               )
                                         newCppClass$buildAll(where = where)
                                         nfCompInfos[[generatorName]]$cppDef <<- newCppClass
                                         newCppClass ## possible return value
                                     } else {
                                         nfCompInfos[[generatorName]]$cppDef ## return value if already exists
                                     }
                                 },
                                 instantiateNimbleFunction = function(nf, dll) {
                                     if(!is.nf(nf)) stop("Can't instantiateNimbleFunction, nf is not a nimbleFunction")
                                     generatorName <- nfGetDefVar(nf, 'name')
                                     ans <- getNimbleFunctionCppDef(generatorName = generatorName)$Rgenerator(nf, dll = dll, project = .self)
                                     ## This was wrong (goes in the refClassObject now) and it is set in Rgenerator
##                                     environment(nf)$.CobjectInterface <- ans
                                     ans
                                 },
                                 setupVirtualNimbleFunction = function(vfun, fromModel = FALSE) {
                                     if(!is.nfGenerator(vfun)) stop('nimbleFunctionVirtual provided to project is not a nimbleFunction generator.', call. = FALSE)
                                     if(!(environment(vfun)$virtual)) stop("Something provided as a nimbleFunctionVirtual is an nfGenerator but not a virtual one.", call. = FALSE)
                                     
                                     generatorName <- environment(vfun)$name
                                     if(is.null(nfCompInfos[[generatorName]])) {
                                         nfCompInfos[[generatorName]] <<- nfCompilationInfoClass(nfGenerator = vfun,
                                                                                                 Rcompiled = FALSE, written = FALSE, cppCompiled = FALSE, loaded = FALSE,
                                                                                                 RinitTypesProcessed = FALSE, virtual = TRUE, fromModel = fromModel)
                                         nfCompInfos[[generatorName]]$labelMaker <<- NULL ## not needed?
                                     }
                                     if(inherits(nfCompInfos[[generatorName]]$nfProc, 'uninitializedField'))
                                         nfCompInfos[[generatorName]]$nfProc <<- virtualNFprocessing$new(vfun, generatorName)
                                     nfCompInfos[[generatorName]]$nfProc
                                 },
                                 compileNimbleFunctionMulti = function(funList, isNode = FALSE, filename = NULL, initialTypeInferenceOnly = FALSE,
                                     control = list(debug = FALSE, debugCpp = FALSE, compileR = TRUE, writeFiles = TRUE, compileCpp = TRUE, loadSO = TRUE),
                                     reset = FALSE, returnCppClass = FALSE, where = globalenv(), fromModel = FALSE) {
                                     if(!is.list(funList)) stop('funList in compileNimbleFunctionMulti should be a list', call. = FALSE)
                                     allGeneratorNames <- try(lapply(funList, nfGetDefVar, 'name'))
                                     if(inherits(allGeneratorNames, 'try-error')) stop('Problem getting generatorNames in compileNimbleFunctionMulti.', call. = FALSE)
                                     uniqueGeneratorNames <- unique(allGeneratorNames)
                                     if(initialTypeInferenceOnly || returnCppClass) {
                                         ans <- list()
                                         if(initialTypeInferenceOnly) oldKnownGeneratorNames <- names(nfCompInfos)
                                     } else ans <- vector('list', length(funList))
                                     for(uGN in uniqueGeneratorNames) {
                                         thisBool <- allGeneratorNames == uGN
                                         thisAns <- compileNimbleFunction( funList[ thisBool ], isNode = isNode, filename = filename,
                                                                          initialTypeInferenceOnly = initialTypeInferenceOnly,
                                                                          control = control, reset = reset, returnCppClass = returnCppClass, where = where,
                                                                          fromModel = fromModel)
                                         if(initialTypeInferenceOnly || returnCppClass) {
                                             
                                             if(initialTypeInferenceOnly) { ## return only new NFprocs in this case.
                                                 thisGeneratorName <- environment(thisAns$nfGenerator)$name ## should be same as uGN
                                                 if(!(thisGeneratorName %in% oldKnownGeneratorNames)) ans[[ thisGeneratorName ]] <- thisAns
                                             }
                                             else ## they should all be new in this case anyway
                                                 ans[[ uGN ]] <- thisAns
                                         } else {
                                             ans[thisBool] <- if(is.list(thisAns)) thisAns else list(thisAns)
                                         }
                                     }
                                     ans
                                 },
                                 compileNimbleFunction = function(fun, isNode = FALSE, filename = NULL, initialTypeInferenceOnly = FALSE,
                                     control = list(debug = FALSE, debugCpp = FALSE, compileR = TRUE, writeFiles = TRUE, compileCpp = TRUE, loadSO = TRUE),
                                     reset = FALSE, returnCppClass = FALSE, where = globalenv(), fromModel = FALSE) {
                                     ## fundamental difference: fun should be a specialized nf: nfGenerator will not longer track instances
                                     
                                     if(is.character(fun)) {
                                         tmp <- nimbleFunctions[[fun]]
                                         if(is.null(tmp)) stop(paste0("nimbleFunction name ", fun, " not recognized in this project."), call. = FALSE)
                                         if(reset) {
                                             nf_getRefClassObject(tmp)$.CobjectInterface <- NULL
                                             ##                                             if(exists('.CobjectInterface', envir = environment(tmp), inherits = FALSE))
##                                                 rm('.CobjectInterface', envir = environment(tmp))
                                         }
                                         fun <- tmp
                                         funList <- list(fun)
                                         generatorName <- nfGetDefVar(fun, 'name')
                                         if(reset) nfCompInfos[[generatorName]] <<- NULL
                                     } else {
                                         if(is.list(fun)) {
                                             if(length(fun) == 0) stop('Empty list provided to compileNimbleFunction', call. = FALSE)
                                             if(!all(unlist(lapply(fun, function(x) is.nf(x) ) ) ) ) stop('Not all elements in fun list for compileNimbleFunction are specialized nimbleFunctions', call. = FALSE)
                                             generatorName <- unique(unlist(lapply(fun, nfGetDefVar, 'name')))
                                             if(length(generatorName) != 1) stop(paste0('Not all elements in the fun list for compileNimbleFunction are specialized from the same nimbleFunction.  The generator names include:', paste(generatorName, collapse = ' ')), call. = FALSE)
                                             if(reset) nfCompInfos[[generatorName]] <<- NULL
                                             funList <- fun
                                         } else {
                                             if(!is.nf(fun)) stop(paste0("fun argument provided is not a nimbleFunction."), call. = FALSE)
                                             funList <- list(fun)
                                             generatorName <- nfGetDefVar(fun, 'name')
                                             if(reset) nfCompInfos[[generatorName]] <<- NULL
                                         }
                                         for(i in seq_along(funList)) {
                                             addNF <- TRUE
                                             thisName <- nf_getRefClassObject(funList[[i]])[['name']] 
                                             if(!is.null(thisName)) {
                                                 tmp <- nimbleFunctions[[thisName]]
                                                 if(!is.null(tmp)) {
                                                     if(reset) {
                                                         nimbleFunctions[[thisName]] <<- NULL
                                                     } else {
                                                         if(!identical(funList[[i]], tmp)) stop('Trying to compile simething with same name as previously added nimbleFunction that is not the same thing')
                                                         addNF <- FALSE
                                                     }
                                                 }
                                             }
                                             if(addNF) {
                                                 addNimbleFunction(funList[[i]], fromModel = fromModel)
                                             }
                                         }
                                     }
##                                     Cname <- nf_getRefClassObject(funList[[1]])$Cname
                                     
                                     ## what about virtuals?
                                     if(!exists('name', envir = nf_getRefClassObject(funList[[1]]), inherits = FALSE)) stop('Something is wrong if by this point in compileNimbleFunction there is no name.', call. = FALSE)
                                     cppClass <- buildNimbleFunctionCompilationInfo(funList, isNode = isNode, initialTypeInferenceOnly = initialTypeInferenceOnly, control = control, where = where, fromModel = fromModel)
                                     if(initialTypeInferenceOnly || returnCppClass) return(cppClass) ## cppClass is an nfProc in this case
                                     
                                     if(!nfCompInfos[[generatorName]]$written && control$writeFiles) {
                                         cppProj <- cppProjectClass(dirName = dirName)
                                         cppProjects[[ generatorName ]] <<- cppProj
                                         if(is.null(filename)) filename <- paste0(projectName, '_', Rname2CppName(generatorName))
                                         cppProj$addClass(cppClass, generatorName, filename)
                                     ## need cppProj
                                         cppProj$writeFiles(filename)
                                         nfCompInfos[[generatorName]]$written <<- TRUE
                                     } else {
                                         cppProj <- cppProjects[[ generatorName ]]
                                         writeLines('Using previously generated C++ code.  This will not work if the current nimbleFunction specializations use types of modelValues or other nimbleFunctions that have not already been compiled in this project.  If that is the case, you should include these specializiations in the first compilation of the nimbleFunction.  You can compile all the specializations of this nimbleFunction together with reset = TRUE.')
                                     }
                                     if(!nfCompInfos[[generatorName]]$cppCompiled) {
                                         cppProj$compileFile(filename)
                                         nfCompInfos[[generatorName]]$cppCompiled <<- TRUE
                                     } else writeLines('Using previously compiled C++ code.')
                                     if(!nfCompInfos[[generatorName]]$loaded) {
                                         cppProj$loadSO(filename)
                                         nfCompInfos[[generatorName]]$loaded <<- TRUE
                                     } else writeLines('Using previously loaded compilation unit.')
                                     
                                     ans <- vector('list', length(funList))

                                     for(i in seq_along(funList)) {
                                         ans[[i]] <- nfCompInfos[[generatorName]]$cppDef$buildCallable(funList[[i]], cppProj$dll)
                                         nimbleFunctionCppInterfaces[[ nf_getRefClassObject(funList[[i]])$name ]] <<- ans[[i]]
                                     }
                                     if(length(ans) == 1) ans[[1]] else ans
                                 }
                                 )
                             )

#' compile NIMBLE models and nimbleFunctions
#'
#' compile a collection of models and nimbleFunctions: generate C++, compile the C++, load the result, and return an interface object
#'
#' @param ...  An arbitrary set of NIMBLE models and nimbleFunctions, or lists of them.  If given as named parameters, those names may be used in the return list.
#' @param project  Optional NIMBLE model or nimbleFunction already associated with a project, which the current units for compilation should join. If not provided, a new project will be created and the current compilation units will be associated with it.
#' @param dirName  Optional directory name in which to generate the C++ code.  If not provided, a temporary directory will be generated using R's \code{tempdir} function.
#' @param projectName Optional character name for labeling the project if it is new
#' @param control  A list mostly for internal use. See details.
#' @param resetFunctions Logical value stating whether nimbleFunctions associated with an existing project should all be reset for compilation purposes.  See details.
#' @author Perry de Valpine
#' @export
#' @details
#' This is the main function for calling the NIMBLE compiler.  A set of compiler calls and output will be seen.  Compiling in NIMBLE does 4 things:
#' 1. It generates C++ code files for all the model and nimbleFunction components.  2. It calls the system's C++ compiler.  3. It loads the compiled object(s) into R using \code{dyn.load}. And 4. it generates R objects for using the compiled model and nimbleFunctions.
#'
#' When the units for compilation provided in \code{...} include multiple models and/or nimbleFunctions, models are compiled first, in the order in which they are provided.  Groups of nimbleFunctions that were specialized from the same nimbleFunction generator (the result of a call to \code{nimbleFunction}, which then takes setup arguments and returns a specialized nimbleFunction) are then compiled as a group, in the order of first appearance.
#'
#' The behavior of adding new compilation units to an existing project is limited.  For example, one can compile a model in one call to \code{compileNimble} and then compile a nimbleFunction that uses the model (i.e. was given the model as a setup argument) in a second call to \code{compileNimble}, with the model provided as the \code{project} argument.  Either the uncompiled or compiled model can be provided.  However, compiling a second nimbleFunction and adding it to the same project will only work in limited circumstances.  Basically, the limitations occur because it attempts to re-use already compiled pieces, but if these do not have all the necessary information for the new compilation, it gives up.  An attempt has been made to give up in a controlled manner and provide somewhat informative messages.
#'
#' When compilation is not allowed or doesn't work, try using \code{resetFunctions = TRUE}, which will force recompilation of all nimbleFunctions in the new call.  Previously compiled nimbleFunctions will be unaffected, and their R interface objects should continue to work.  The only cost is additional compilation time for the current compilation call.  If that doesn't work, try re-creating the model and/or the nimbleFunctions from their generators.  An alternative possible fix is to compile multiple units in one call, rather than sequentially in multiple calls.
#'
#'#' The control list can contain the following named elements, each with \code{TRUE} or \code{FALSE}: debug, which sets a debug mode for the compiler for development purposes; debugCpp, which inserts an output message before every line of C++ code for debugging purposes; compileR, which determines whether the R-only steps of compilation should be executed; writeCpp, which determines whether the C++ files should be generated; compileCpp, which determines whether the C++ should be compiled;  loadSO, which determines whether the DLL or shared object should be loaded and interfaced; and returnAsList, which determines whether calls to the compiled nimbleFunction should return only the returned value of the call (\code{returnAsList = FALSE}) or whether a list including the input arguments, possibly modified, should be returned in a list with the returned value of the call at the end (\code{returnAsList = TRUE}).  The control list is mostly for developer use, although \code{returnAsArgs} may be useful to a user.  An example of developer use is that one can have the compiler write the C++ files but not compile them, then modify them by hand, then have the C++ compiler do the subsequent steps without over-writing the files.
#'
#' See NIMBLE User Manual for examples
#' 
#' @return If there is only one compilation unit (one model or nimbleFunction), an R interface object is returned.  This object can be used like the uncompiled model or nimbleFunction, but execution will call the corresponding compiled objects or functions.  If there are multiple compilation units, they will be returned as a list of interface objects, in the order provided.  If names were included in the arguments, or in a list if any elements of \code{...} are lists, those names will be used for the corresponding element of the returned list.  Otherwise an attempt will be made to generate names from the argument code.  For example \code{compileNimble(A = fun1, B = fun2, project = myModel)} will return a list with named elements A and B, while \code{compileNimble(fun1, fun2, project = myModel)} will return a list with named elements fun1 and fun2.
#'
#' 
#' 
compileNimble <- function(..., project, dirName = NULL, projectName = '',
                          control = list(),
                          resetFunctions = FALSE) {
## 1. Extract compilation items
    reset <- FALSE
    ## This pulls out ... arguments, makes names from their expressions if names weren't provided, and combines them with any ... arguments that are lists.

    controlDefaults = list(debug = FALSE, debugCpp = FALSE, compileR = TRUE, writeFiles = TRUE, compileCpp = TRUE, loadSO = TRUE, returnAsList = FALSE)
    
    dotsDeparses <- unlist(lapply( substitute(list(...))[-1], deparse ))
    origList <- list(...)
    if(is.null(names(origList))) names(origList) <- rep('', length(origList))
    boolNoName <- names(origList)==''
    origIsList <- unlist(lapply(origList, is.list))
    dotsDeparses[origIsList] <- ''
    names(origList)[boolNoName] <- dotsDeparses[boolNoName]
    units <- do.call('c', origList)
##    units <- do.call('c', list(...)) ## combines any elements that are themselves lists
    unitTypes <- getNimbleTypes(units)
    if(length(grep('unknown', unitTypes)) > 0) stop(paste0('Some items provided for compilation do not have types that can be compiled.  The types provided were: ', paste0(unitTypes, collapse = ' '), '. Be sure only specialized nimbleFunctions are provided, not nimbleFunction generators.'), call. = FALSE)
    if(is.null(names(units))) names(units) <- rep('', length(units))
    if(length(units) == 0) stop('No objects for compilation provided')
    
    ## 2. Get project or make new project
    if(missing(project)) {
        if(reset) warning('reset = TRUE but no project was provided.  If you are trying to re-compiled something into the same project, give it as the project argument as well as a compilation item.  e.g. compileNimble(myFunction, project = myFunction, reset = TRUE)')
        project <- nimbleProjectClass(dirName, name = projectName)
    } else {
        project <- getNimbleProject(project, TRUE)
        if(is.null(project)) stop("Invalid project argument", call. = FALSE)
    }
    if(resetFunctions) project$resetFunctions()

    for(i in names(controlDefaults)) {
        if(!i %in% names(control)) control[[i]] <- controlDefaults[[i]]
    }
    
    ## Units should be either Rmodel, nimbleFunction, or RCfunction (now coming from nimbleFunction with no setup)
    ## Compile models first
    ans <- list()
    rcfUnits <- unitTypes == 'rcf'
    if(sum(rcfUnits) > 0) {
        whichUnits <- which(rcfUnits)
        for(i in whichUnits) {
            ans[[i]] <- project$compileRCfun(units[[i]], control = control)
            if(names(units)[i] != '') names(ans)[i] <- names(units)[i]
        }
    }
    
    modelUnits <- unitTypes == 'model'
    if(sum(modelUnits) > 0) {
        whichUnits <- which(modelUnits)
        for(i in whichUnits) {
            ans[[ i ]] <- project$compileModel(units[[i]], control = control)
            if(names(units)[i] != '') names(ans)[i] <- names(units)[i]
        }
    }
    nfUnits <- unitTypes == 'nf'
    if(sum(nfUnits) > 0) {
        whichUnits <- which(nfUnits)
        nfAns <- project$compileNimbleFunctionMulti(units[whichUnits], control = control, reset = reset)
        ans[whichUnits] <- nfAns
        for(i in whichUnits) if(names(units)[i] != '') names(ans)[i] <- names(units)[i]
    }
    
    
    if(length(ans) == 1) ans[[1]] else ans
}

getNimbleTypes <- function(units) {
    ans <- character(length(units))
    for(i in seq_along(units)) {
        if(inherits(units[[i]], 'RModelBaseClass')) ans[i] <- 'model'
        else if(is.nf(units[[i]])) ans[i] <- 'nf'   ## a nimbleFunction
        else if(is.rcf(units[[i]])) ans[i] <- 'rcf' ## an RCfunction = a nimbleFunction with no setup
        else if(is.nfGenerator(units[[i]])) ans[i] <- 'unknown(nf generator)'
        else ans[i] <- 'unknown'
    }
    ans
}

#' return the nimble project, if any, associated with a model or nimbleFunction object.
getNimbleProject <- function(project, stopOnNull = FALSE) {
    if(inherits(project, 'nimbleProjectClass')) return(project)
    if(is.nf(project)) return(nfVar(project, 'nimbleProject'))
    if(inherits(project, 'RmodelBaseClass')) stop('need to handle RmodelBaseClass in getNimbleProject')
    ans <- try(project$nimbleProject)
    if(inherits(ans, 'try-error') | is.null(ans)) {
        if(stopOnNull) stop(paste0('cannot determine nimbleProject from provided project argument'))
        return(NULL)
    }
    ans
}


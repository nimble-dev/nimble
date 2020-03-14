##.nimbleProjectClassMasterList <- new.env(, emptyenv())

projectNameCreator <- labelFunctionCreator('P')

nfCompilationInfoClass <- setRefClass('nfCompilationInfoClass',
    fields = list(
        nfProc = 		'ANY',      ## an nfProcessing object 
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
        addRinstance = function(nfi) {Rinstances[[ length(Rinstances)+1 ]] <<- nfi},
        addRinstanceList = function(nfList) {Rinstances[length(Rinstances) + seq_along(nfList)] <<- nfList}
    )
)

nlCompilationInfoClass <- setRefClass('nlCompilationInfoClass',
    fields = list(
        nlProc = 'ANY',
        cppDef = 'ANY',       ## a cppNimbleFunctionClass object
        written =  'ANY',		#'logical'
        loaded = 'ANY',
        cppCompiled =  'ANY',		#'logical'
        labelMaker = 'ANY', ## a label maker function
        RinitTypesProcessed = 'ANY',		# 'logical', ## setupTypesForUsingFunction() 
        Rcompiled = 'ANY'   # 'logical'
    ),
    methods = list(
        initialize = function(...){callSuper(...)}
    )
)

mvInfoClass <- setRefClass('mvInfoClass',
    fields = list(
        mvConf = 'ANY', ## a custom modelValues class
        cppClassName =  'ANY',		#'character',
        cppClass = 'ANY', ## a cppModelValuesClass object,
        fromModel =  'ANY',		#'logical',
        RmvObjs =  'ANY'		#'list'
        ),
    methods = list(
        initialize = function(...) {
            RmvObjs <<- list()
            callSuper(...)
        },
        addRmv = function(Rmv) RmvObjs[[length(RmvObjs)+1]] <<- Rmv
    )
)

RCfunInfoClass <- setRefClass('RCfunInfoClass',
    fields = list(
        nfMethodRCobj = 'ANY', ## an mfMethodRC
        RCfunProc     = 'ANY', ## an RCfunProcessing or NULL
        cppClass      = 'ANY',  ## an RCfunctionDef or NULL
        RinitTypesProcessed = 'ANY',
        Rcompiled           = 'ANY',
        fromModel     =  'ANY'		#'logical'
    )
)

modelDefInfoClass <- setRefClass('modelDefInfoClass',
    fields = list(
        labelMaker = 'ANY'
    )
)

nimbleProjectClass <- setRefClass('nimbleProjectClass',
    fields = list(
        RCfunInfos         =  'ANY',		#'list', ## a list of RCfunInfoClass objects
        RCfunCppInterfaces =  'ANY',		#'list', 
        mvInfos            =  'ANY',		#'list', ## a list of mvInfoClass objects
        modelDefInfos      =  'ANY',		#'list',
        models             =  'ANY',		#'list',
        nimbleFunctions    =  'ANY',		#'list',
        nimbleLists        =  'ANY',   #'list',
        nfCompInfos        =  'ANY',		#'list', ## list of nfCompilationInfoClass objects
        nlCompInfos        =  'ANY',   #'list', ## list of nfCompilationInfoClass objects
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
            RCfunInfos <<- new.env()						# list()
            RCfunCppInterfaces <<- new.env()				# list()
            mvInfos <<- new.env()							# list()
            modelDefInfos <<- new.env()						# list()
            models <<- new.env()							# list()
            nimbleFunctions <<- new.env()					# list
            nimbleLists <<- new.env()
            nfCompInfos <<- list()							# list()
            nlCompInfos <<- list()
            cppProjects <<- list()							#new.env()						#list()
            refClassDefsEnv <<- new.env()
            dirName <<- if(is.null(dir)) makeDefaultDirName() else dir
            if(name == '') projectName <<- projectNameCreator() else projectName <<- name
           },
        clearCompiled = function(functions = TRUE, models = TRUE, DLLs = TRUE) {
            if(functions) resetFunctions(finalize = TRUE)
            if(models) resetModels(finalize = TRUE)
            if(DLLs) unloadDLLs()
        },
        unloadDLLs = function() {
            for(i in seq_along(cppProjects)) {
                cppProjects[[i]]$unloadSO(check = TRUE, force = FALSE)
            }
        },
        resetModels = function(finalize = TRUE) {
            for(i in ls(models)) {
                models[[i]]$CobjectInterface$finalizeInternal()
                models[[i]]$CobjectInterface <<- NULL
                models[[i]]$nimbleProject <<- NULL
                models[[i]]$Cname <<- character(0)
                models[[i]] <<- NULL
            }
            for(i in ls(RCfunInfos)) {
                if(length(RCfunInfos[[i]]$fromModel) > 0) {
                    if(RCfunInfos[[i]]$fromModel) {
                        if(!is.null(RCfunCppInterfaces[[i]])) { ## could be null if it was a neededType and an interface was never built 
                            environment(RCfunCppInterfaces[[i]])$CnativeSymbolInfo_ <<- NULL
                            ##RCfunCppInterfaces[[i]] <<- NULL
                            rm(list = i, envir = RCfunCppInterfaces)
                        }
                        assign('nimbleProject', NULL, envir = RCfunInfos[[i]]$nfMethodRCobj)
                        thisName <- RCfunInfos[[i]]$nfMethodRCobj$uniqueName
                        ##if(!is.null(cppProjects[[thisName]])) cppProjects[[thisName]] <<- NULL
                        ##RCfunInfos[[i]] <<- NULL
                        rm(list = i, envir = RCfunInfos)
                    }
                }
            }
            for(i in ls(nfCompInfos)) {
                if(length(nfCompInfos[[i]]$fromModel) > 0) {
                    if(nfCompInfos[[i]]$fromModel) {
                        for(j in seq_along(nfCompInfos[[i]]$Rinstances)) {
                            thisEnv <- environment(nfCompInfos[[i]]$Rinstances[[j]])
                            thisRCO <- nf_getRefClassObject(nfCompInfos[[i]]$Rinstances[[j]])
                            if(exists('name', envir = thisRCO, inherits = FALSE)) {
                                thisname <- thisRCO$name
                                rm(list = thisname, envir = nimbleFunctions)
                                rm('name', envir = thisRCO)
                            }
                            thisRCO[['nimbleProject']] <- NULL
                            ## finalization done internally to the model interface
                        }
                        nfCompInfos[[i]] <<- NULL
                        ## cppProjects[[i]] <<- NULL
                    }
                }
            }
        },
        resetFunctions = function(finalize = FALSE) {
            ## clear everything except models and nimbleFunctions from models
            for(i in ls(mvInfos)) {
                clearThisMV <- TRUE ## It looked like in some situations we'd not want to clear here, but apparently we always should
                ## OK, what happens is we used to check if a cppClass was from a model and then not clear, but that isn't right.
                ## It could be defined from a model but then have objects from nimbleFunctions that need to be cleared.
                if(clearThisMV) {
                    mvInfos[[i]]$cppClass <<- NULL
                    for(j in seq_along(mvInfos[[i]]$RmvObjs)) {
                        if(!is.null(mvInfos[[i]]$RmvObjs[[j]]$CobjectInterface)) {
                            if(finalize) {
                                mvInfos[[i]]$RmvObjs[[j]]$CobjectInterface$finalizeInternal() 
                            }
                            mvInfos[[i]]$RmvObjs[[j]]$CobjectInterface <<- NULL
                        }
                    }
                    rm(list = i, envir = mvInfos)
                }
            }
            
            for(i in ls(RCfunInfos)) {
                if(length(RCfunInfos[[i]]$fromModel) > 0) {
                    if(!RCfunInfos[[i]]$fromModel) {
                        if(!is.null(RCfunCppInterfaces[[i]])) { ## could be null if it was a neededType and an interface was never built 
                            environment(RCfunCppInterfaces[[i]])$CnativeSymbolInfo_ <<- NULL
                            ##RCfunCppInterfaces[[i]] <<- NULL
                            rm(list = i, envir = RCfunCppInterfaces)
                        }
                        assign('nimbleProject', NULL, envir = RCfunInfos[[i]]$nfMethodRCobj)
                        thisName <- RCfunInfos[[i]]$nfMethodRCobj$uniqueName
                        ##if(!is.null(cppProjects[[thisName]])) cppProjects[[thisName]] <<- NULL
                        ##RCfunInfos[[i]] <<- NULL
                        rm(list = i, envir = RCfunInfos)
                    }
                }
            }
            for(i in ls(nfCompInfos)) {
                if(length(nfCompInfos[[i]]$fromModel) > 0) {
                    if(!nfCompInfos[[i]]$fromModel) {
                        for(j in seq_along(nfCompInfos[[i]]$Rinstances)) {
                            thisEnv <- environment(nfCompInfos[[i]]$Rinstances[[j]])
                            thisRCO <- nf_getRefClassObject(nfCompInfos[[i]]$Rinstances[[j]])
                            if(exists('name', envir = thisRCO, inherits = FALSE)) {
                                thisname <- thisRCO$name
                                rm(list = thisname, envir = nimbleFunctions)
                                ##nimbleFunctions[[ thisname ]] <<- NULL
                                rm('name', envir = thisRCO)
                            }
                            thisRCO[['nimbleProject']] <- NULL
                            if(finalize) {
                                if(!is.null(thisRCO$.CobjectInterface)) {
                                    if(is.list( thisRCO$.CobjectInterface)) { ## CmultiNimbleFunctionInterface
                                        thisRCO$.CobjectInterface[[1]]$finalizeInstance(thisRCO$.CobjectInterface[[2]])
                                    } else {
                                        thisRCO$.CobjectInterface$finalizeInternal()
                                    }
                                    thisRCO$.CobjectInterface <- NULL
                                }
                            }
                        }
                        nfCompInfos[[i]] <<- NULL
                        ## cppProjects[[i]] <<- NULL
                    }
                }
            }
        },
            
        addModelValuesClass = function(mvConf, fromModel = FALSE) {
            mvClassName <- environment(mvConf)$className
            if(!is.null(mvInfos[[mvClassName]])) stop('Trying to add a modelValues class with the same name as one already in this project', call. = FALSE)
            mvInfos[[mvClassName]] <<-  mvInfoClass(cppClassName = mvClassName, cppClass = NULL, mvConf = mvConf, fromModel = fromModel)
        },
        getModelValuesCppDef = function(mvConf, NULLok = FALSE) {
            mvClassName <- environment(mvConf)$className
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
            if(!inherits(model, 'RmodelBaseClass')) stop('model provided to project is not an RmodelBaseClass', call. = FALSE)
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
                ##modelCppInterfaces[[ model$name ]] <<- new.env()	#list(NULL)
            }
        },
        addNimbleFunctionMulti = function(funList, fromModel = FALSE, generatorFunNames = NULL) {
            if(length(funList) == 0) return(invisible(NULL))
            if(!is.nf(funList[[1]])) stop('first nimbleFunction provided to project is not a nimbleFunction.', call. = FALSE)
            inProjectAlready <- unlist(lapply(funList, function(x) identical(x[['nimbleProject']], .self)))
            if(any(inProjectAlready)) {
                stop('Trying to add list of nimbleFunctions but some are already part of another project. If you are recompiling, try redefining models and specialized nimbleFunctions. (The reset option works now for nimbleFunctions but not models.)', call. = FALSE)
            }
            allGeneratorNames <- if(is.null(generatorFunNames))
                                     unlist(lapply(funList, function(x) environment(x$.generatorFunction)$name), use.names = FALSE)
                                 else
                                     generatorFunNames
            generatorName2Indices <- split(seq_along(funList), allGeneratorNames) ##uniqueGeneratorNamesIndices <- which(!duplicated(allGeneratorNames))
            for(i in seq_along(generatorName2Indices)) { ##genID in uniqueGeneratorNameIndices) {
                genID <- generatorName2Indices[[i]][1]
                generatorName <- names(generatorName2Indices)[i] ##allGeneratorNames[genID]
                if(is.null(nfCompInfos[[generatorName]])) {
                    ## nfProc could have been created already during makeTypeObject for another nimbleFunction so it knows the types of this one.
                    nfCompInfos[[generatorName]] <<- nfCompilationInfoClass(nfGenerator = nf_getGeneratorFunction(funList[[genID]]),
                                                                            Rcompiled = FALSE, written = FALSE, cppCompiled = FALSE, loaded = FALSE,
                                                                            RinitTypesProcessed = FALSE, virtual = FALSE,
                                                                            fromModel = fromModel)
                    nfCompInfos[[generatorName]]$labelMaker <<- labelFunctionCreator(paste0(generatorName,'_'))
                }
                genIDs <- generatorName2Indices[[i]]
                newLabels <- nfCompInfos[[generatorName]]$labelMaker(count = length(genIDs))
                
                namedFunList <- funList[genIDs]
                names(namedFunList) <- newLabels
                list2env(namedFunList, envir = nimbleFunctions)
                nfCompInfos[[generatorName]]$addRinstanceList(namedFunList)
                
                newEnvs <- lapply(seq_along(newLabels), new.env)
                names(newEnvs) <- newLabels
                CppNewLabels <- Rname2CppName(newLabels)
                for(j in seq_along(newLabels)) {
                    funList[[genIDs[j]]][['Cname']] <- CppNewLabels[j]
                    funList[[genIDs[j]]][['name']] <- newLabels[j] ## skip the checking done in addNimbleFunction
                    funList[[genIDs[j]]][['nimbleProject']] <- .self
                }
            }
        },
        addNimbleFunction = function(fun, fromModel = FALSE) {
            if(!is.nf(fun)) stop('nimbleFunction provided to project is not a nimbleFunction.', call. = FALSE)
            inProjectAlready <- nf_getRefClassObject(fun)[['nimbleProject']]
            if(!is.null(inProjectAlready)) {
                if(!identical(inProjectAlready, .self)) stop('Trying to add a specialized nimbleFunction to a project, but it is already part of another project. \nIf you did not specify a project, this error can occur in trying to create a new project -- you likely need to specify the relevant model as the project.\nIf you are recompiling, try redefining models and specialized nimbleFunctions. (The reset option works now for nimbleFunctions but not models.)', call. = FALSE)
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
           
            if(!exists('Cname', envir = nf_getRefClassObject(fun), inherits = FALSE)) {
                assign('Cname', Rname2CppName(nf_getRefClassObject(fun)$name), envir = nf_getRefClassObject(fun))
            }

            assign('nimbleProject', .self, envir = nf_getRefClassObject(fun))
            ## could check for duplicate Cnames here, but if the names are unique the Cnames should be too.
        },
        addNimbleListGen = function(nlGen) {
            if(!is.nlGenerator(nlGen)) stop('invalid nimbleListGen provided to addNimbleListGen.', call. = FALSE)
            ## get className
            className <- nl.getListDef(nlGen)$className
            
            ## check if there is a nlCompInfos,
            ## add if needed
            if(is.null(nlCompInfos[[className]])) {
                ## nfProc could have been created already during makeTypeObject for another nimbleFunction so it knows the types of this one.
                nlCompInfos[[className]] <<- nlCompilationInfoClass(written = FALSE, cppCompiled = FALSE, Rcompiled = FALSE,
                                                                    RinitTypesProcessed = FALSE, loaded = FALSE)
                nlCompInfos[[className]]$labelMaker <<- labelFunctionCreator(paste0(className,'_'))
            }
            nestedListGens <- nl.getNestedGens(nlGen)
            for(i in seq_along(nestedListGens))
                addNimbleListGen(nestedListGens[[i]])
        },
        addNimbleList = function(nl, fromModel = FALSE, nestedList = FALSE) { ##fromModel never used: clean up. 
          if(!is.nl(nl)) stop('nimbleList provided to project is not a nimbleList.', call. = FALSE)
          inProjectAlready <- nl[['nimbleProject']]
          if(!is.null(inProjectAlready)) {
            if(!identical(inProjectAlready, .self)) stop('Trying to add a specialized nimbleList to a project, but it is already part of another project. \nIf you did not specify a project, this error can occur in trying to create a new project -- you likely need to specify the relevant model as the project.\nIf you are recompiling, try redefining models and specialized nimbleFunctions and nimbleLists.', call. = FALSE)
            else warning('Adding a specialized nimbleList to a project to which it already belongs', call. = FALSE)
          }

          ## Next two lines can become addNimbleListGen, but would have to migrate recursion out of compileNimbleList
          className <- nl$nimbleListDef$className
          if(is.null(nlCompInfos[[className]])) {
            ## nfProc could have been created already during makeTypeObject for another nimbleFunction so it knows the types of this one.
            nlCompInfos[[className]] <<- nlCompilationInfoClass(written = FALSE, cppCompiled = FALSE, Rcompiled = FALSE,
                                                                RinitTypesProcessed = FALSE, loaded = FALSE)
            nlCompInfos[[className]]$labelMaker <<- labelFunctionCreator(paste0(className,'_'))
          }
         
          if(!exists('name', envir = nl, inherits = FALSE)) {
            assign('name', nlCompInfos[[className]]$labelMaker(), envir = nl)
          } else {
            if(!is.null(nimbleLists[[ nl$name ]])) {
              stop('nimbleList provided to project has same name as another one in the same project', call. = FALSE)
            }
          }
          if(!nestedList)   nimbleLists[[ nl$name ]] <<- nl
          # nlCompInfos[[generatorName]]$addRinstance(nl)
          
          if(!exists('Cname', envir = nl, inherits = FALSE)) {
            assign('Cname', Rname2CppName(nl$name), envir = nl)
          }
          
          assign('nimbleProject', .self, envir = nl)
          ## could check for duplicate Cnames here, but if the names are unique the Cnames should be too.
        },
        addRCfun = function(nfmObj, fromModel = FALSE) {
            if(!inherits(nfmObj, 'nfMethodRC')) stop("Can't add this function. nfmObj is not an nfMethodRC", call. = FALSE)
            className <- nfmObj$uniqueName
            if(is.null(RCfunInfos[[className]])) {
                RCfunInfos[[className]] <<- RCfunInfoClass(nfMethodRCobj = nfmObj, RCfunProc = NULL, cppClass = NULL, fromModel = fromModel, RinitTypesProcessed = FALSE, Rcompiled = FALSE)
            }
            assign('nimbleProject', .self, envir = nfmObj) ## needed for clearCompiled(), i.e. safe dyn.unload()
        },
        getRCfunCppDef = function(nfmObj, NULLok = FALSE) {
            className <- nfmObj$uniqueName
            ans <- RCfunInfos[[className]]
            if(is.null(ans)) {
                if(!NULLok) stop("Requested to get an RCfunCppDef but it is not in the project and NULLok = FALSE", call. = FALSE)
                return(NULL)
            }
            ans <- ans$cppClass
            if(inherits(ans, 'uninitializedField') )  return(NULL)                                     	 
            ans
        },
        needRCfunCppClass = function(nfmObj, genNeededTypes = TRUE, initialTypeInference = FALSE, control = list(debug = FALSE, debugCpp = FALSE), fromModel = FALSE) {
            if(!inherits(nfmObj, 'nfMethodRC')) stop("Can't compile this function. nfmObj is not an nfMethodRC", call. = FALSE)
            className <- nfmObj$uniqueName
            Cname <- Rname2CppName(className)
            RCfunInfo <- RCfunInfos[[className]]
            if(is.null(RCfunInfo)) addRCfun(nfmObj, fromModel = fromModel)
            if(is.null(RCfunInfos[[className]]$RCfunProc)) {
                RCfunInfos[[className]]$RCfunProc <<- RCfunProcessing$new(nfmObj, Cname)
            }
            if(!RCfunInfos[[className]]$RinitTypesProcessed) {
                RCfunInfos[[className]]$RCfunProc$process(debug = control$debug, debugCpp = control$debugCpp, initialTypeInferenceOnly = TRUE, nimbleProject = .self)
                RCfunInfos[[className]]$RinitTypesProcessed <<- TRUE
            }
            if(initialTypeInference) return(RCfunInfos[[className]]$RCfunProc)
            if(!RCfunInfos[[className]]$Rcompiled) {
                RCfunInfos[[className]]$RCfunProc$process(debug = control$debug, debugCpp = control$debugCpp, initialTypeInferenceOnly = FALSE, nimbleProject = .self)
                RCfunInfos[[className]]$Rcompiled <<- TRUE
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
        compileRCfun = function( fun, filename = NULL, initialTypeInference = FALSE, control = list(debug = FALSE, debugCpp = FALSE, writeFiles = TRUE, returnAsList = FALSE), showCompilerOutput = nimbleOptions('showCompilerOutput')) {
            disableWrite <- FALSE
            if(nimbleOptions('enableSpecialHandling')) {
                SH <- filenameFromSpecialHandling(fun)
                if(!is.null(filename)) {
                    filename <- SH
                    disableWrite <- TRUE
                }
            }
            if(is.rcf(fun)) fun <- environment(fun)$nfMethodRCobject
            addRCfun(fun) ## checks if it already exists and if it is an nfMethodRC ## redundant? done also in next step.
            cppClass <- needRCfunCppClass(fun, genNeededTypes = TRUE, initialTypeInference = initialTypeInference, control = control)
            if(initialTypeInference) return(cppClass) ## in this case cppClass with be an RCfunProc
            className <- fun$uniqueName
            if(control$writeFiles) {
                cppProj <- cppProjectClass(dirName = dirName)
                cppProjects[[ className ]] <<- cppProj
                if(is.null(filename)) filename <- paste0(projectName, '_', className)
                cppProj$addClass( cppClass, className, filename )
                if(!disableWrite) cppProj$writeFiles(filename)
            }
            if(control$compileCpp) {
                cppProj$compileFile(filename, showCompilerOutput)
            }
            if(control$loadSO) {
                cppProj$loadSO(filename)
            }
            RCfunCppInterfaces[[className]] <<- cppClass$buildRwrapperFunCode(includeLHS = FALSE, eval = TRUE, returnArgsAsList = control$returnAsList, dll = cppProj$dll)
            RCfunCppInterfaces[[className]]
        },
        needModelValuesCppClass = function(mvConf, fromModel = FALSE) {
            if(!isModelValuesConf(mvConf)) stop("Can't compileModelValues: mvConf is not a modelValuesConf", call. = FALSE)
            mvClassName <- environment(mvConf)$className
            mvInfo <- mvInfos[[mvClassName]]
            if(is.null(mvInfo)) addModelValuesClass(mvConf, fromModel)
            cppClass <- mvInfos[[mvClassName]]$cppClass
            if(is.null(cppClass)) {
                cppClass <- cppModelValuesClass(name = mvClassName,
                                                   vars = environment(mvConf)$symTab,
                                                   project = .self)
                cppClass$buildAll()
                 mvInfos[[mvClassName]]$cppClass <<- cppClass
            }
            cppClass
        },
        instantiateCmodelValues = function(mv, dll) {
            mvClassName <- class(mv)
            cppDef <- mvInfos[[mvClassName]]$cppClass
            if(is.null(cppDef)) stop('Trying to instantiate a modelValues type that the project has no record of. Try setting option resetFunctions = TRUE in compileNimble')
            generatorName <- cppDef$SEXPgeneratorFun$name
            sym = if(!is.null(dll))
                getNativeSymbolInfo(generatorName, dll)
            else {
                warning('a nimbleFunctionInterface is about to build a CmodelValues without dll info, based on generatorFun name only.', call. = FALSE)
                generatorName
            }
            ans <- CmodelValues(sym, dll = dll)
            mvInfos[[mvClassName]]$addRmv(mv) ## simply a list for later clearing
            mv$CobjectInterface <- ans
            ans
        },
        compileModel = function(model, filename = NULL,
                                control = list(debug = FALSE, debugCpp = FALSE, writeFiles = TRUE, compileCpp = TRUE, loadSO = TRUE), showCompilerOutput = nimbleOptions('showCompilerOutput'), where = globalenv()) {
            disableWrite <- FALSE
            if(nimbleOptions('enableSpecialHandling')) {
                filenames <- filenameFromSpecialHandling(model)
                if(!is.null(filenames)) {
                    filename <- filenames$filename
                    nfFileName <- filenames$nfFileName
                    disableWrite <- TRUE
                }
            }
            if(is.character(model)) {
                tmp <- models[[model]]
                if(is.null(tmp)) stop(paste0("Model provided as name: ", model, " but it is not in this project."), call. = FALSE)
                model <- tmp
            } else addModel(model)
                                                 
            modelDef <- model$getModelDef()
            modelDefName <- modelDef$name
            Cname <- Rname2CppName(modelDefName)
            if(!disableWrite) {
                if(is.null(filename)) {
                    filename <- paste0(projectName, '_', Rname2CppName(modelDefName)) 
                }
                nfFileName <- paste0(projectName, '_', Rname2CppName(modelDefName),'_nfCode')
            }
            modelCpp <- cppBUGSmodelClass(modelDef = modelDef, model = model,
                                          name = Cname, project = .self)
            ## buildAll will call back to the project to add its nimbleFunctions 
            modelCpp$buildAll(buildNodeDefs = TRUE, where = where, control = control)
            
            cppProj <- cppProjectClass(dirName = dirName)
            cppProjects[[ modelDefName ]] <<- cppProj
            ## genModelValuesCppClass will back to the project to add its mv class
            mvc <- modelCpp$genModelValuesCppClass()
            ##if(is.null(filename)) filename <- paste0(projectName, '_', modelDefName)
            cppProj$addClass(mvc, filename = filename)
            cppProj$addClass(modelCpp, modelDefName, filename)
            ##if compileNodes
            ##nfFileName <- paste0(projectName, '_', Rname2CppName(modelDefName),'_nfCode')
            for(i in names(modelCpp$nodeFuns)) {
                cppProj$addClass(modelCpp$nodeFuns[[i]], filename = nfFileName)
            }
            if(control$writeFiles) {
                if(!disableWrite) {
                    cppProj$writeFiles(filename)
                    cppProj$writeFiles(nfFileName) ## if compileNodes
                }
            } else return(cppProj)
            if(control$compileCpp) {
                compileList <- filename
                compileList <- c(compileList, nfFileName) ## if compileNodes
                cppProj$compileFile(compileList, showCompilerOutput)
            } else return(cppProj)
            if(control$loadSO) {
                ## if loadSO
                cppProj$loadSO(filename)
            } else return(cppProj)
            ## if buildInterface
            interfaceName <- paste0('C', modelDefName)
            
            compiledModel <- modelCpp ## cppProj$cppDefs[[2]]
            newCall <- paste0('new_',Rname2CppName(modelDefName))
            ans <- buildModelInterface(interfaceName, compiledModel, newCall, where = where, project = .self, dll = cppProj$dll)
            createModel <- TRUE
            if(!createModel) return(ans) else return(ans(model, where, dll = cppProj$dll))
            ## creating the model populates model$CobjectInterface
        },
        ## nimbleList functions
        addNestedNls = function(nl){
          for(iNl in names(nl$nestedListGenList)){
            addNimbleList(nl[[iNl]], nestedList = TRUE)
            if(length(nl[[iNl]]$nestedListGenList) > 0){
              addNestedNls(nl[[iNl]])
            }
          }
        },
        compileNimbleList = function(nl, filename = NULL, initialTypeInferenceOnly = FALSE,
            control = list(debug = FALSE, debugCpp = FALSE, compileR = TRUE, writeFiles = TRUE, compileCpp = TRUE, loadSO = TRUE),
            reset = FALSE, returnCppClass = FALSE, className = NULL, alreadyAdded = FALSE) { ## className? alreadyAdded?

            ## add possibility that nl is a generator
            generatorOnly <- FALSE
            if(is.nlGenerator(nl)) {
                ## determine className
                className <- nl.getListDef(nl)$className
                generatorOnly <- TRUE
                nlGen <- nl
                nlList <- nl
            } else {
            
                if(is.list(nl)) {
                    if(is.null(className)) className <- unique(unlist(lapply(nl, function(x) x$nimbleListDef$className)))
                    if(length(className) != 1) stop(paste0('Not all elements in the nimbleList list for compileNimbleList are from the same nimbleFunctionDef.  The class names include:', paste(className, collapse = ' ')), call. = FALSE)
                    nlList <- nl
                    ## set generator
                } else {
                    if(!is.nl(nl)) stop(paste0("nl argument provided is not a nimbleList."), call. = FALSE)
                    nlList <- list(nl)
                    className <- nl$nimbleListDef$className
                    ## set generator
                }
                nlGen <- nl.getGenerator(nlList[[1]])
            }
            if(reset) nlCompInfos[[className]] <<- NULL
            if(!alreadyAdded) {
                if(generatorOnly) {
                    ## check if generator exists and do addNimbleListGen
                    ## and recurse into nestedListGenList
                    addNimbleListGen(nlGen)
                } else {
                    for(i in seq_along(nlList)) {
                        addNL <- TRUE
                        thisName <- nlList[[i]][['name']]
                        if(!is.null(thisName)) {
                            tmp <- nimbleLists[[thisName]]
                            if(!is.null(tmp)) {
                                if(reset) {
                                    nimbleLists[[thisName]] <<- NULL
                                } else {
                                    if(!identical(nlList[[i]], tmp)) stop('Trying to compile something with same name as previously added nimbleList that is not the same thing')
                                    addNL <- FALSE
                                }
                            }
                        }
                        if(addNL){
                            addNimbleList(nlList[[i]])
                            ## if any nested lists, add them too (recursively)
                            if(length(nlList[[i]]$nestedListGenList) > 0){
                                addNestedNls(nlList[[i]])
                            }
                        }
                    }
                }
            }

            ## modify to pull nestedListGenList from generator
            nestedListGens <- nl.getNestedGens(nlGen)
            for(iNestedNL in seq_along(nestedListGens)) {
                compileNimbleList(nestedListGens[[iNestedNL]], initialTypeInferenceOnly = TRUE, alreadyAdded = TRUE)
            }
            ## for(iNestedNl in seq_along(nlList[[1]]$nestedListGenList)){
            ##   ## create cppInfo for any nested list classes 
            ##   compileNimbleList(nlList[[1]][[names(nlList[[1]]$nestedListGenList)[iNestedNl]]], initialTypeInferenceOnly = TRUE, alreadyAdded = TRUE)
            ## }
            cppClass <- buildNimbleListCompilationInfo(nlList, initialTypeInferenceOnly = initialTypeInferenceOnly)
            

            
            if(initialTypeInferenceOnly || returnCppClass) return(cppClass)
            message('Remaining compileNimbleList is not yet adapted')
            if(!nlCompInfos[[className]]$written && control$writeFiles) {
                cppProj <- cppProjectClass(dirName = dirName)
                cppProjects[[ className ]] <<- cppProj
                if(is.null(filename)) filename <- paste0(projectName, '_', Rname2CppName(className))
                cppProj$addClass(cppClass, className, filename)
                cppProj$writeFiles(filename)
                nlCompInfos[[className]]$written <<- TRUE
            } else {
                if(!control$writeFiles) return(cppProj)
                cppProj <- cppProjects[[ className ]]
            }
            if(!nlCompInfos[[className]]$cppCompiled && control$compileCpp) {
                if(control$compileCpp) {
                    cppProj$compileFile(filename)
                    nlCompInfos[[className]]$cppCompiled <<- TRUE
                } else writeLines('Skipping compilation because control$compileCpp is FALSE')
            } else {if(!control$compileCpp) return(cppProj)}#writeLines('Using previously compiled C++ code.')
            if(!nlCompInfos[[className]]$loaded && control$loadSO) {
                cppProj$loadSO(filename)
                nlCompInfos[[className]]$loaded <<- TRUE
            } else{if(!control$loadSO) return(cppProj)}# writeLines('Using previously loaded compilation unit.')
            
            ans <- vector('list', length(nlList))

            for(i in seq_along(nlList)) {
                ans[[i]] <- nlCompInfos[[className]]$cppDef$buildCallable(nlList[[i]], cppProj$dll, asTopLevel = TRUE)
            }
            if(length(ans) == 1) ans[[1]] else ans
        },
        
        ## nimbleFunction functions
        getNimbleFunctionCppDef = function(generatorName, nfProc) {
            if(missing(generatorName)) {
                if(missing(nfProc)) stop('No good information provided to getNimbleFunctionCppDef', call. = FALSE)
                generatorName <- environment(nfProc$nfGenerator)[['name']]
                if(is.null(generatorName)) stop('Invalid generatorName', call. = FALSE)
            }
            if(is.null(nfCompInfos[[generatorName]])){
                 return(NULL)
           }
            ans <- nfCompInfos[[generatorName]]$cppDef
            if(inherits(ans, 'uninitializedField') )  return(NULL)                                     	 
            ans
        },
        getNimbleFunctionNFproc = function(fun) {
            generatorName <- nfGetDefVar(fun, 'name')
            if(is.null(nfCompInfos[[generatorName]])) return(NULL)
            ans <- nfCompInfos[[generatorName]]$nfProc
            if(inherits(ans, 'uninitializedField')) return(NULL)
            ans
        },
        getNimbleListCppDef = function(generatorName, nlProc) {
          if(missing(generatorName)) {
            if(missing(nlProc)) stop('No good information provided to getNimbleListCppDef', call. = FALSE)
            generatorName <- nlProc$nimbleListObj$className
            if(is.null(generatorName)) stop('Invalid generatorName', call. = FALSE)
          }
          if(is.null(nlCompInfos[[generatorName]])){
            return(NULL)
          }
          ans <- nlCompInfos[[generatorName]]$cppDef
          if(inherits(ans, 'uninitializedField') )  return(NULL)                                     	 
          ans
        },
        getNimbleListNLproc = function(fun) {
          generatorName <- fun$name
          if(is.null(nlCompInfos[[generatorName]])) return(NULL)
          ans <- nlCompInfos[[generatorName]]$nlProc
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
                nfCompInfos[[generatorName]]$nfProc <<- virtualNFprocessing$new(vfun, generatorName, project = .self)
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
                    nfCompInfos[[generatorName]]$nfProc <<- nfProcessing(funList, generatorName, fromModel = fromModel, project = .self, isNode = isNode)
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
        buildNimbleListCompilationInfo = function(listList = NULL, className, initialTypeInferenceOnly = FALSE, 
                                                    control = list(debug = FALSE, debugCpp = FALSE), where = globalenv(), fromModel = FALSE
                                                  ) {
            if(!is.null(listList)) {
                ## check for nimbleListGen and get className
                if(is.nlGenerator(listList)) {
                    className <- nl.getListDef(listList)$className
                } else {
                    className <- listList[[1]]$nimbleListDef$className
                    name <- listList[[1]]$name
                    Cname <- listList[[1]]$Cname
                }
                if(is.null(nlCompInfos[[className]])) stop("Requested buildNimbleListCompilationInfo for a generator that has not been added to the project", call. = FALSE)
                if(inherits(nlCompInfos[[className]]$nlProc, 'uninitializedField')) 
                    nlCompInfos[[className]]$nlProc <<- nlProcessing(listList, className, project = .self)
            } else {
                if(missing(className)) stop("If listList is omitted, a class name must be provided to buildNimbleListCompilationInfo", call. = FALSE)
                if(inherits(nlCompInfos[[className]]$nlProc, 'uninitializedField')) stop("buildNimbleListCompilationInfo was called with only a className (probably from genNeededTypes), but the nfProc is missing.", call. = FALSE)
            }
            if(initialTypeInferenceOnly) {
                if(!nlCompInfos[[className]]$RinitTypesProcessed) {
              nlCompInfos[[className]]$nlProc$setupTypesForUsingFunction() 
              nlCompInfos[[className]]$RinitTypesProcessed <<- TRUE
            }
            return(nlCompInfos[[className]]$nlProc)
          }
          if(!nlCompInfos[[className]]$Rcompiled) {
            nlCompInfos[[className]]$nlProc$process(control = control)
            nlCompInfos[[className]]$Rcompiled <<- TRUE
          }
          if(inherits(nlCompInfos[[className]]$cppDef, 'uninitializedField')) {
            newCppClass <- cppNimbleListClass(name = className,
                                              nimCompProc = nlCompInfos[[className]]$nlProc,
                                              debugCpp = control$debugCpp,
                                              project = .self
            )
            newCppClass$buildAll(where = where)
            nlCompInfos[[className]]$cppDef <<- newCppClass
            newCppClass ## possible return value
          } else {
            nfCompInfos[[className]]$cppDef ## return value if already exists
          }
        },
        instantiateNimbleList = function(nl, dll, asTopLevel = TRUE) { ## called by cppInterfaces_models and cppInterfaces_nimbleFunctions
          ## to instantiate neededObjects
          for(nestedNL in names(nl$nestedListGenList)) {
            nestedAns <- instantiateNimbleList(nl[[nestedNL]], dll, asTopLevel)
          }

          if(!is.nl(nl)) stop("Can't instantiateNimbleList, nl is not a nimbleList")
          className <- nl$nimbleListDef$className
          nlCppDef <- getNimbleListCppDef(generatorName = className)
            ok <- TRUE
            dllToUse <- if(isTRUE(nl.getDefinitionContent(nl.getGenerator(nl), 'predefined')))
                            nimbleUserNamespace$sessionSpecificDll
                        else dll
          if(asTopLevel) {
            if(is.null(nlCppDef$Rgenerator)) ok <- FALSE
            else ans <- nlCppDef$Rgenerator(nl, dll = dllToUse, project = .self)
          } else {
            if(is.null(nlCppDef$CmultiInterface)) ok <- FALSE
            else ans <- nlCppDef$CmultiInterface$addInstance(nl, dll = dllToUse)
          }
          if(!ok) stop("Oops, there is something in this compilation job that doesn\'t fit together.  This can happen in some cases if you are trying to compile new pieces into an exising project.  If that is the situation, please try including \"resetFunctions = TRUE\" as an argument to compileNimble.  Alternatively please try rebuilding the project from the beginning with more pieces in the same call to compileNimble.  For example, if you are compiling multiple algorithms for the same model in multiple calls to compileNimble, try compiling them all with one call.", call. = FALSE) 
          ans
        },
        instantiateNimbleFunction = function(nf, dll, asTopLevel = TRUE) { ## called by cppInterfaces_models and cppInterfaces_nimbleFunctions
            ## to instantiate neededObjects
            if(!is.nf(nf)) stop("Can't instantiateNimbleFunction, nf is not a nimbleFunction")
            generatorName <- nfGetDefVar(nf, 'name')

            nfCppDef <- getNimbleFunctionCppDef(generatorName = generatorName)

            ok <- !is.null(nfCppDef)
            if(ok) {
                ans <- nfCppDef$buildCallable(nf, dll = dll, asTopLevel = asTopLevel)
                ok <- !is.null(ans)
            }
            ## ok <- TRUE
            ## if(asTopLevel) {
            ##     if(is.null(nfCppDef$Rgenerator)) ok <- FALSE
            ##     else ans <- nfCppDef$Rgenerator(nf, dll = dll, project = .self)
            ## } else {
            ##     if(is.null(nfCppDef$CmultiInterface)) ok <- FALSE
            ##     else ans <- nfCppDef$CmultiInterface$addInstance(nf, dll = dll)
            ## }
            if(!ok) stop("Oops, there is something in this compilation job that doesn\'t fit together.  This can happen in some cases if you are trying to compile new pieces into an exising project.  If that is the situation, please try including \"resetFunctions = TRUE\" as an argument to compileNimble.  Alternatively please try rebuilding the project from the beginning with more pieces in the same call to compileNimble.  For example, if you are compiling multiple algorithms for the same model in multiple calls to compileNimble, try compiling them all with one call.", call. = FALSE) 

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
                nfCompInfos[[generatorName]]$nfProc <<- virtualNFprocessing$new(vfun, generatorName, project = .self)
            if(!nfCompInfos[[generatorName]]$Rcompiled) { ## to support nimbleLists, this step goes here now, so by size processing the method symbol tables will be set up
                nfCompInfos[[generatorName]]$nfProc$process(control = control) ## there is no need for an initialTypeInference flag because that is *all* that a virtual NF really does anyway
                nfCompInfos[[generatorName]]$Rcompiled <<- TRUE
            }
            nfCompInfos[[generatorName]]$nfProc
        },
        compileNimbleFunctionMulti = function(funList, isNode = FALSE, filename = NULL, initialTypeInferenceOnly = FALSE,
            control = list(debug = FALSE, debugCpp = FALSE, compileR = TRUE, writeFiles = TRUE, compileCpp = TRUE, loadSO = TRUE),
            reset = FALSE, returnCppClass = FALSE, where = globalenv(), fromModel = FALSE, generatorFunNames = NULL, alreadyAdded = FALSE, showCompilerOutput = nimbleOptions('showCompilerOutput')) {
            if(!is.list(funList)) stop('funList in compileNimbleFunctionMulti should be a list', call. = FALSE)
            allGeneratorNames <- if(is.null(generatorFunNames)) lapply(funList, nfGetDefVar, 'name') else generatorFunNames
            uniqueGeneratorNames <- unique(allGeneratorNames)
            if(initialTypeInferenceOnly || returnCppClass) {
                ans <- list()
                if(initialTypeInferenceOnly) oldKnownGeneratorNames <- ls(nfCompInfos)	#names(nfCompInfos)
            } else ans <- vector('list', length(funList))
            for(uGN in uniqueGeneratorNames) {
                thisBool <- allGeneratorNames == uGN
                thisAns <- compileNimbleFunction( funList[ thisBool ], isNode = isNode, filename = filename,
                                                 initialTypeInferenceOnly = initialTypeInferenceOnly,
                                                 control = control, reset = reset, returnCppClass = returnCppClass, where = where,
                                                 fromModel = fromModel, generatorName = uGN, alreadyAdded = alreadyAdded, showCompilerOutput = showCompilerOutput)
                if(initialTypeInferenceOnly || returnCppClass) {
                    
                    if(initialTypeInferenceOnly) { ## return only new NFprocs in this case.
                        thisGeneratorName <- environment(thisAns$nfGenerator)$name ## should be same as uGN
                        if(!(thisGeneratorName %in% oldKnownGeneratorNames)) ans[[ thisGeneratorName ]] <- thisAns
                    }
                    else ## they should all be new in this case anyway
                        ans[[ uGN ]] <- thisAns
                } else {
                    ans[thisBool] <- thisAns ##if(is.list(thisAns)) thisAns else list(thisAns)
                }
            }
            ans
        },
        compileNimbleFunction = function(fun, isNode = FALSE, filename = NULL, initialTypeInferenceOnly = FALSE,
            control = list(debug = FALSE, debugCpp = FALSE, compileR = TRUE, writeFiles = TRUE, compileCpp = TRUE, loadSO = TRUE),
            reset = FALSE, returnCppClass = FALSE, where = globalenv(), fromModel = FALSE, generatorName = NULL, alreadyAdded = FALSE, showCompilerOutput = nimbleOptions('showCompilerOutput')) {
          if(is.character(fun)) {
                tmp <- nimbleFunctions[[fun]]
                if(is.null(tmp)) stop(paste0("nimbleFunction name ", fun, " not recognized in this project."), call. = FALSE)
                if(reset) {
                    nf_getRefClassObject(tmp)$.CobjectInterface <- NULL
                }
                fun <- tmp
                funList <- list(fun)
                generatorName <- nfGetDefVar(fun, 'name')
                if(reset) nfCompInfos[[generatorName]] <<- NULL
            } else {
                if(is.list(fun)) {
                    if(length(fun) == 0) stop('Empty list provided to compileNimbleFunction', call. = FALSE)
                    if(is.null(generatorName)) generatorName <- unique(unlist(lapply(fun, nfGetDefVar, 'name')))
                    if(length(generatorName) != 1) stop(paste0('Not all elements in the fun list for compileNimbleFunction are specialized from the same nimbleFunction.  The generator names include:', paste(generatorName, collapse = ' ')), call. = FALSE)
                    if(reset) nfCompInfos[[generatorName]] <<- NULL
                    funList <- fun
                } else {
                    if(!is.nf(fun)) stop(paste0("fun argument provided is not a nimbleFunction."), call. = FALSE)
                    funList <- list(fun)
                    generatorName <- nfGetDefVar(fun, 'name')
                    if(reset) nfCompInfos[[generatorName]] <<- NULL
                }
                if(!alreadyAdded) {
                    for(i in seq_along(funList)) {
                        addNF <- TRUE
                        thisName <- nf_getRefClassObject(funList[[i]])[['name']] 
                        if(!is.null(thisName)) {
                            tmp <- nimbleFunctions[[thisName]]
                            if(!is.null(tmp)) {
                                if(reset) {
                                    nimbleFunctions[[thisName]] <<- NULL
                                } else {
                                    if(!identical(funList[[i]], tmp)) stop('Trying to compile something with same name as previously added nimbleFunction that is not the same thing')
                                    addNF <- FALSE
                                }
                            }
                        }
                        if(addNF) {
                            addNimbleFunction(funList[[i]], fromModel = fromModel)
                        }
                    }
                }
            }
              Cname <- nf_getRefClassObject(funList[[1]])$Cname

            if(!exists('name', envir = nf_getRefClassObject(funList[[1]]), inherits = FALSE)) stop('Something is wrong if by this point in compileNimbleFunction there is no name.', call. = FALSE)
            cppClass <- buildNimbleFunctionCompilationInfo(funList, isNode = isNode, initialTypeInferenceOnly = initialTypeInferenceOnly, control = control, where = where, fromModel = fromModel)
            if(initialTypeInferenceOnly || returnCppClass) return(cppClass) ## cppClass is an nfProc in this case

            ## At this point we are ready to write, compile, load and instantiate.
            ## However the system for tracking these steps is not perfect.
            ## Specifically, if another nimbleFunction containing an object of the current nimbleFunction
            ## has already been compiled, then the files for the current one will have been written and
            ## compiled, but that status is not recorded in its nfCompInfos object.
            ## Also there could be other objects in the current funList that have not been instantiated.
            ## As a kludgy fix, we determine whether any of the objects already have a .CobjectInterface
            ## If they do, we skip over writing, compiling and loading steps.
            hasCobjectInterface <- unlist(lapply(funList, function(x) {
                RCO <- nf_getRefClassObject(x)
                !inherits(RCO$.CobjectInterface, 'uninitializedField') || is.null(RCO$.CobjectInterface))
            }))

            if(!any(hasCobjectInterface)) {
                if(!nfCompInfos[[generatorName]]$written && control$writeFiles) {
                    cppProj <- cppProjectClass(dirName = dirName)
                    cppProjects[[ generatorName ]] <<- cppProj
                    if(is.null(filename)) filename <- paste0(projectName, '_', Rname2CppName(generatorName))
                    cppProj$addClass(cppClass, generatorName, filename)
                    cppProj$writeFiles(filename)
                    nfCompInfos[[generatorName]]$written <<- TRUE
                } else {
                    if(!control$writeFiles) return(cppProj)
                    cppProj <- cppProjects[[ generatorName ]]
                }
                if(!nfCompInfos[[generatorName]]$cppCompiled && control$compileCpp) {
                    if(control$compileCpp) {
                        cppProj$compileFile(filename, showCompilerOutput)
                        nfCompInfos[[generatorName]]$cppCompiled <<- TRUE
                    } else writeLines('Skipping compilation because control$compileCpp is FALSE')
                } else {if(!control$compileCpp) return(cppProj)}#writeLines('Using previously compiled C++ code.')
                if(!nfCompInfos[[generatorName]]$loaded && control$loadSO) {
                    cppProj$loadSO(filename)
                    nfCompInfos[[generatorName]]$loaded <<- TRUE
                } else{if(!control$loadSO) return(cppProj)}# writeLines('Using previously loaded compilation unit.')
            }
            ans <- vector('list', length(funList))
            
            for(i in seq_along(funList)) {
                if(!(hasCobjectInterface[i]))
                    ans[[i]] <- nfCompInfos[[generatorName]]$cppDef$buildCallable(funList[[i]], cppProj$dll, asTopLevel = TRUE)
                else {
                    ## A curious possibility: If a nf was built first nested in another nf,
                    ## its interface may be a Cmulti interface, which is not directly user friendly
                    ## But if we are here via compileNimbleFunction, the user has included it in the compile request and wants
                    ## a full interface
                    ## Hence we call promote interface which checks and builds a new one if needed
                    ans[[i]] <- nfCompInfos[[generatorName]]$cppDef$promoteCallable(funList[[i]]) ##nf_getRefClassObject(funList[[i]])$.CobjectInterface
                }
            }
            ##if(length(ans) == 1) ans[[1]] else ans
            ans
        }
    )
)

clearCompiled <- function(obj) { # for now just take one obj as input
    np <- getNimbleProject(obj)
    np$clearCompiled()
}

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
#' @param showCompilerOutput Logical value indicating whether details of C++ compilation should be printed. 
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
#' The control list can contain the following named elements, each with \code{TRUE} or \code{FALSE}: debug, which sets a debug mode for the compiler for development purposes; debugCpp, which inserts an output message before every line of C++ code for debugging purposes; compileR, which determines whether the R-only steps of compilation should be executed; writeCpp, which determines whether the C++ files should be generated; compileCpp, which determines whether the C++ should be compiled;  loadSO, which determines whether the DLL or shared object should be loaded and interfaced; and returnAsList, which determines whether calls to the compiled nimbleFunction should return only the returned value of the call (\code{returnAsList = FALSE}) or whether a list including the input arguments, possibly modified, should be returned in a list with the returned value of the call at the end (\code{returnAsList = TRUE}).  The control list is mostly for developer use, although \code{returnAsArgs} may be useful to a user.  An example of developer use is that one can have the compiler write the C++ files but not compile them, then modify them by hand, then have the C++ compiler do the subsequent steps without over-writing the files.
#'
#' See NIMBLE User Manual for examples
#' 
#' @return If there is only one compilation unit (one model or nimbleFunction), an R interface object is returned.  This object can be used like the uncompiled model or nimbleFunction, but execution will call the corresponding compiled objects or functions.  If there are multiple compilation units, they will be returned as a list of interface objects, in the order provided.  If names were included in the arguments, or in a list if any elements of \code{...} are lists, those names will be used for the corresponding element of the returned list.  Otherwise an attempt will be made to generate names from the argument code.  For example \code{compileNimble(A = fun1, B = fun2, project = myModel)} will return a list with named elements A and B, while \code{compileNimble(fun1, fun2, project = myModel)} will return a list with named elements fun1 and fun2.
#'
#' 
#' 
compileNimble <- function(..., project, dirName = NULL, projectName = '',
                          control = list(),
                          resetFunctions = FALSE, 
			  showCompilerOutput = nimbleOptions('showCompilerOutput')) {
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
    if(any(sapply(units, is, "MCMCconf")))
       stop("You have provided an MCMC configuration object, which cannot be compiled. Instead, use run 'buildMCMC' on the configuration object and compile the resulting MCMC object.")
    unitTypes <- getNimbleTypes(units)
    if(length(grep('unknown', unitTypes)) > 0) stop(paste0('Some items provided for compilation do not have types that can be compiled: ', paste0(names(units), collapse = ' '), '.  The types provided were: ', paste0(unitTypes, collapse = ' '), '. Be sure only specialized nimbleFunctions are provided, not nimbleFunction generators.'), call. = FALSE)
    if(is.null(names(units))) names(units) <- rep('', length(units))
    if(length(units) == 0) stop('No objects for compilation provided')

    ## 2. Get project or make new project
    if(missing(project)) {
        if(reset) warning("reset = TRUE but no project was provided.  If you are trying to re-compiled something into the same project, give it as the project argument as well as a compilation item. For example, 'compileNimble(myFunction, project = myFunction, reset = TRUE)'")
        if(!is.null(nimbleOptions()$nimbleProject)) project <- nimbleOptions()$nimbleProject
        else project <- nimbleProjectClass(dirName, name = projectName)

        ## Check for uncompiled models.
        if(!any(sapply(units, is, 'RmodelBaseClass'))) {
            mcmcUnits <- which(sapply(units, class) == "MCMC")
            if(any(sapply(mcmcUnits, function(idx) {
                class(units[[idx]]$model$CobjectInterface) == "uninitializedField"
            })))
                stop("compileNimble: The model associated with an MCMC is not compiled. Please compile the model first.")
        }
    } else {
        project <- getNimbleProject(project, TRUE)
        if(!inherits(project, 'nimbleProjectClass'))
            stop("Invalid project argument; note that models and nimbleFunctions need to be compiled before they can be used to specify a project. Once compiled you can use an R model or nimbleFunction to specify the project.", call. = FALSE)
    }
    if(resetFunctions) project$resetFunctions()

    for(i in names(controlDefaults)) {
        if(!i %in% names(control)) control[[i]] <- controlDefaults[[i]]
    }
    

    ## Units should be either Rmodel, nimbleFunction, or RCfunction (now coming from nimbleFunction with no setup)
    if(nimbleOptions('verbose') && !showCompilerOutput) message("compiling... this may take a minute. Use 'showCompilerOutput = TRUE' to see C++ compilation details.")
    if(nimbleOptions('verbose') && showCompilerOutput) message("compiling... this may take a minute. On some systems there may be some compiler warnings that can be safely ignored.")

    ## Compile models first
    ans <- list()
    rcfUnits <- unitTypes == 'rcf'
    if(sum(rcfUnits) > 0) {
        whichUnits <- which(rcfUnits)
        for(i in whichUnits) {
            ans[[i]] <- project$compileRCfun(units[[i]], control = control, showCompilerOutput = showCompilerOutput)
            if(names(units)[i] != '') names(ans)[i] <- names(units)[i]
        }
    }
    
    modelUnits <- unitTypes == 'model'
    if(sum(modelUnits) > 0) {
        whichUnits <- which(modelUnits)
        for(i in whichUnits) {
            ans[[ i ]] <- project$compileModel(units[[i]], control = control, showCompilerOutput = showCompilerOutput)
            if(names(units)[i] != '') names(ans)[i] <- names(units)[i]
        }
    }
    nfUnits <- unitTypes == 'nf'
    if(sum(nfUnits) > 0) {
        whichUnits <- which(nfUnits)
        nfAns <- project$compileNimbleFunctionMulti(units[whichUnits], control = control,
                                                    reset = reset, showCompilerOutput = showCompilerOutput)
        ans[whichUnits] <- nfAns
        for(i in whichUnits) if(names(units)[i] != '') names(ans)[i] <- names(units)[i]
    }
    nlUnits <- unitTypes == 'nl'
    if(sum(nlUnits) > 0) {
      whichUnits <- which(nlUnits)
      nlAns <- project$compileNimbleList(units[whichUnits], control = control, reset = reset)
      ans[[whichUnits]] <- nlAns
      for(i in whichUnits) if(names(units)[i] != '') names(ans)[i] <- names(units)[i]
    }
    
    if(nimbleOptions('verbose')) message("compilation finished.")

    if(length(ans) == 1) ans[[1]] else ans
}

getNimbleTypes <- function(units) {
    ans <- character(length(units))
    for(i in seq_along(units)) {
        if(inherits(units[[i]], 'RmodelBaseClass')) ans[i] <- 'model'
        else if(is.nf(units[[i]])) ans[i] <- 'nf'   ## a nimbleFunction
        else if(is.rcf(units[[i]])) ans[i] <- 'rcf' ## an RCfunction = a nimbleFunction with no setup
        else if(is.nfGenerator(units[[i]])) ans[i] <- 'unknown(nf generator)'
        else if(is.nl(units[[i]])) ans[i] <- 'nl'  ## a nimbleList
        else ans[i] <- 'unknown'
    }
    ans
}

# return the nimble project, if any, associated with a model or nimbleFunction object.
getNimbleProject <- function(project, stopOnNull = FALSE) {
    if(inherits(project, 'nimbleProjectClass')) return(project)
    if(is.nf(project)) return(nfVar(project, 'nimbleProject'))
    if(is.rcf(project)) return(environment(project)$nfMethodRCobject$nimbleProject)
    ans <- try(project$nimbleProject)
    if(inherits(ans, 'try-error') | is.null(ans)) {
        if(stopOnNull) stop(paste0('cannot determine nimbleProject from provided project argument'))
        return(NULL)
    }
    ans
}

countDllObjects <- function() {
    eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$CountDllObjects))
}

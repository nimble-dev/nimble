## Actual symbolTable class is at the end

buildSymbolTable <- function(vars, type, size){
    symTab <- symbolTable()
    for(i in 1:length(vars) ) {
        if(is.list(size) ) 
            symTab$addSymbol( symbolBasic(name = vars[i], type = type[i], nDim = length(size[[i]]) , size = size[[i]] ) )
        else
            symTab$addSymbol( symbolBasic(name = vars[i], type = type[i], nDim = length(size) , size = size ) )
    }
    return(symTab)
}


argTypeList2symbolTable <- function(ATL, neededTypes, origNames) {
    ## ATL is the type of argument type list from run-time args to a nimble function
    symTab <- symbolTable()
    for(i in seq_along(ATL)) {
        symTab$addSymbol(argType2symbol(ATL[[i]], neededTypes, names(ATL)[i], origNames[i]))
    }
    symTab
}

## This takes an unevaluated expression as input (AT), such as
## double(), double(2), or double(2, c(5,5))
## The first argument is the number of dimensions, defaulting to 0 (indicated scalar)
## The second argument is a vector of the sizes, defaulting to rep(NA, nDim)
## If only some sizes are known, something like double(2, c(5, NA)) should be valid, but we'll have to check later handling to be sure.

argType2symbol <- function(AT, neededTypes, name = character(), origName = "") {
    ans <- try(argType2symbolInternal(AT, neededTypes, name))
    if(inherits(ans, 'try-error')) stop(paste0("Invalid type type declaration for ",origName,"."), call.=FALSE)
    ans
}

argType2symbolInternal <- function(AT, neededTypes, name = character()) {
    if(!is.null(AT$default))    AT$default <- NULL     ## remove the 'default=' argument, if it's present
    type <- as.character(AT[[1]])
    if(type == "internalType") {
      return(symbolInternalType(name = name, type = "internal", argList = as.list(AT[-1]))) ## save all other contents for any custom needs later
    }
    if(type %in% c('double', 'integer', 'character', 'logical', 'void')){
      nDim <- if(length(AT)==1) 0 else AT[[2]]
      if(!is.numeric(nDim) || nDim %% 1 != 0)
          stop("argType2symbol: unexpected dimension, '", AT[[2]], "', found in argument '", deparse(AT), "'. Dimension should be integer-valued.")
      size <- if(nDim == 0) 1 else {
          if(length(AT) < 3)
              as.numeric(rep(NA, nDim))
          else
              eval(AT[[3]])
      }
      if(type == "character") {
          if(nDim > 1) {
              warning(paste("character argument",name," with nDim > 1 will be treated as a vector"))
              nDim <- 1
              size <- if(any(is.na(size))) as.numeric(NA) else prod(size)
          }
          return(symbolString(name = name, type = "character", nDim = nDim, size = size))
      }  else {
          return(symbolBasic(name = name, type = type, nDim = nDim, size = size))
      }
    }
    if(is.list(neededTypes)){
        ##  isANeededType <- unlist(lapply(neededTypes, function(x) return(type == x$name)))
        isANeededType <- unlist(lapply(neededTypes, `[[`, 'name')) == type
      if(any(isANeededType == 1)){
        listST <- neededTypes[[which(isANeededType == 1)[1]]]$copy(shallow  = TRUE)
        listST$name <- name
        return(listST)
      }
    }
    if(name == "return"){
        message('found the return case of interest')
        browser()
      if(exists(as.character(AT), envir = globalenv()) &&
         is.nlGenerator(eval(parse(text = as.character(AT), keep.source = FALSE)))){
        nlList <- eval(parse(text = paste0(as.character(AT), "$new"), keep.source = FALSE))()
      }
      else if(type %in% names(nlEigenReferenceList)){
        nlList <- nlEigenReferenceList[[type]]$eigenNimbleListDef$new() 
      }
      if(exists('nlList')){
        className <- nlList$nimbleListDef$className
        isANeededType <- (className == names(neededTypes))
        if(any(isANeededType == 1)){
          listST <- neededTypes[[which(isANeededType == 1)[1]]]$copy(shallow = TRUE)
        }
        else{
          listST <- recurseGetListST(className, neededTypes)
        }
        listST$name <- name
        return(listST)
      }
    }   
}

recurseGetListST <- function(className, neededTypes){
  listST <- NULL
  for(NT in neededTypes){
    if(NT$type %in% c('symbolNimbleList', 'symbolNimbleListGenerator')){
      if(!inherits(NT$nlProc$neededTypes, 'uninitializedField')){
         if(className %in% names(NT$nlProc$neededTypes)){
          isANeededType <- (className == names(NT$nlProc$neededTypes))
          listST <- NT$nlProc$neededTypes[[which(isANeededType == 1)[1]]]$copy(shallow = TRUE)
          return(listST)
         }
        else listST <- recurseGetListST(className, NT$nlProc$neededTypes)
      }
    }
  }
  return(listST)
}

symbolTable2cppVars <- function(symTab, arguments = character(), include, parentST = NULL) {
    newSymTab <- symbolTable(parentST = parentST)
    if(missing(include)) include <- symTab$getSymbolNames()
    for(s in include) {
        inputArg <- s %in% arguments
        sObj <- symTab$getSymbolObject(s)
        if(inherits(sObj$type, 'uninitializedField')) stop(paste('Error in symbolTable2cppVars for ', symTab, '. type field is not set.'))
        if(length(sObj$type == 'Ronly') == 0) browser()
        if(sObj$type == 'Ronly') next
        newSymOrList <- symTab$getSymbolObject(s)$genCppVar(inputArg)
        if(is.list(newSymOrList)) {
            for(i in seq_along(newSymOrList)) {
                newSymTab$addSymbol(newSymOrList[[i]])
            }
        } else {
            newSymTab$addSymbol(newSymOrList)
        }
    }
    newSymTab
}

symbolBase <- 
    setRefClass(Class  = 'symbolBase',
                fields = list(name = 'ANY', 		#'character',
                              type = 'ANY' 	),	#'character'),
                methods = list(
                    generateUse = function(...) name
                    )
                )

## nDim and size are redundant for convenience with one exception:
## nDim = 0 must have size = 1 and means it is a true scalar -- NOT sure this is correct anymore...
## nDim = 1 with size = 1 means it is a 1D vector that happens to be length 1
symbolBasic <-
    setRefClass(Class    = 'symbolBasic',
                contains = 'symbolBase',
                fields   = list(nDim  = 'ANY', 		#'numeric',
                                size = 'ANY'), 		#'numeric'),
                methods = list(
                    show = function() {
                        if(inherits(size, 'uninitializedField')) {
                            writeLines(paste0(name,': ', type, ' sizes = (uninitializedField), nDim = ', nDim))
                        } else {
                            writeLines(paste0(name,': ', type, ' sizes = (', paste(size, collapse = ', '), '), nDim = ', nDim))
                        }
                    },
                    genCppVar = function(functionArg = FALSE) {
                        if(type == 'void') return(cppVoid())
                        else if(type == 'integer') cType <- 'int'
                        else if(type == 'double') cType <- 'double'
                        else if(type == 'logical') cType <- 'bool'
                        else warning(paste("in genCppVar method for",name,"in symbolBasic class, type", type,"unrecognized\n"), FALSE)
                        
                        if(nDim == 0) {
                            return(cppVar(baseType = cType,
                                          name = name,
                                          ptr = 0,
                                          ref = FALSE))
                        }
                        if(functionArg) {
                            return(cppNimArr(name = name,
                                             nDim = nDim,
                                             type = cType,
                                             ref = TRUE))
                        } else {
                            return(cppNimArr(name = name,
                                             nDim = nDim,
                                             type = cType))
                        }
                        })
    )

symbolPtr <- setRefClass(
    Class = "symbolPtr",
    contains = "symbolBase",
    methods = list(
        show = function() writeLines(paste('symbolPtr', name)),
        genCppVar = function(...) {
            if(type == 'integer') cType <- 'int'
            if(type == 'double') cType <- 'double'
            if(type == 'logical') cType <- 'bool'
            cppVar(name = name,
                   ptr = 1,
                   baseType = cType)
        })
    )

symbolSEXP <- setRefClass(
  Class = "symbolSEXP",
  contains = "symbolBasic", ## inheriting from symbolBasic instead of symbolBase make initSizes work smoothly
  methods = list(
    show = function() writeLines(paste('symbolSEXP', name)),
    genCppVar = function(...) {
        cppVar(name = name,
               ptr = 0,
               baseType = "SEXP")
    })
)


symbolString <- setRefClass(
    Class = "symbolString",
    contains = "symbolBasic", ## inheriting from symbolBasic instead of symbolBase make initSizes work smoothly
    ## fields   = list(
    ##     nDim  = 'ANY',
    ##     size = 'ANY'), 
    methods = list(
        show = function() writeLines(paste('symbolString', name)),
        genCppVar = function(...) {
            if(nDim == 0) {
                cppVar(name = name, baseType = "std::string")
            } else {
                cppVarFull(name = name, baseType = "vector", templateArgs = list(cppVar(baseType = "std::string")))
            }
        })
    )

symbolNimbleTimer <- setRefClass(
    Class = "symbolNimbleTimer",
    contains = "symbolBase",
    methods = list(
        show = function() writeLines(paste('symbolNimbleTimer', name)),
        genCppVar = function(...) {
            cppVar(name = name, baseType = "nimbleTimerClass_")
        }))

symbolNimArrDoublePtr <- 
    setRefClass(Class    = 'symbolNimArrDoublePtr',
                contains = 'symbolBasic',
                methods = list(
                    show = function() writeLines(paste('symbolNimArrDoublePtr', name)),
                    genCppVar = function(...){
                        if(type == 'integer')      cType <- 'int'
                        else if(type == 'double')  cType <- 'double'
                        else if(type == 'logical') cType <- 'bool'
                        else cat(paste0("Warning: in genCppVar method in symbolBasic class, type unrecognized: ", type, '\n'))
                        return(cppNimArrPtr(name = name, ## this is self-dereferencing
                                            nDim = nDim,
                                            ptr = 2,
                                            type = cType))}
                    )
    )

symbolVecNimArrPtr <- 
    setRefClass(Class    = 'symbolVecNimArrPtr',
                contains = 'symbolBase', ## Important that this is not symbolBase or it will be thought to be directly numeric
                fields = list(
                    nDim = 'ANY', 		#'numeric',
                    size = 'ANY'), 		#'numeric'),
                methods = list(
                    show = function() writeLines(paste('symbolVecNimArrPtr', name)),
                    genCppVar = function(...){
                        if(type == 'integer')      cType <- 'int'
                        else if(type == 'double')  cType <- 'double'
                        else if(type == 'logical') cType <- 'bool'
                        else cat(paste0("Warning: in genCppVar method in symbolBasic class, type unrecognized: ", type, '\n'))
                        return(cppVecNimArrPtr(name = name, selfDereference = TRUE,
                                               nDim = nDim,
                                               ptr = 1,
                                               type = cType))}
                    )
                )

symbolNodeFunctionVector <- 
    setRefClass(Class = 'symbolNodeFunctionVector',
                contains = 'symbolBase',
                methods = list(
                    initialize = function(...) {
                        callSuper(...)
                        type <<- 'symbolNodeFunctionVector'  
                    },
                    show = function() writeLines(paste('symbolNodeFunctionVector', name)),
                    genCppVar = function(...) {
                        return(cppNodeFunctionVector(name = name)) 
                    }
                    )
                )

symbolModel <- 
    setRefClass(Class = 'symbolModel',
                contains = 'symbolBase',
                fields = list(className = 'ANY'), 
                methods = list(
                    initialize = function(...) {
                        callSuper(...)
                        ## type == 'local' means it is defined in setupCode and so will need to have an object and be built
                        ## type == 'Ronly' means it is a setupArg and may be a different type for different nimbleFunction specializations
                        ##                 and it will be like a model in C++ code: not there except by extracted pointers inside of it
                    },
                    show = function() writeLines(paste('symbolModel', name)),
                    genCppVar = function(...) {
                        return(cppVar(name = name, baseType = "Ronly", ptr = 1))
                    }
                    )
                )


symbolModelValues <- 
    setRefClass(Class = 'symbolModelValues',
                contains = 'symbolBase',
                fields = list(mvConf = 'ANY'), 
                methods = list(
                    initialize = function(...) {
                        callSuper(...)
                        ## type == 'local' means it is defined in setupCode and so will need to have an object and be built
                        ## type == 'Ronly' means it is a setupArg and may be a different type for different nimbleFunction specializations
                        ##                 and it will be like a model in C++ code: not there except by extracted pointers inside of it
                    },
                    show = function() writeLines(paste('symbolModelValues', name)),
                    genCppVar = function(...) {
                        return(cppVar(name = name, baseType = "Values", ptr = 1))
                    }
                    )
                )

symbolMemberFunction <-
    setRefClass(Class = 'symbolMemberFunction',
                contains = 'symbolBase',
                fields = list(nfMethodRCobj = 'ANY'),
                methods = list(
                    initialize = function(...) {callSuper(...); type <<- 'Ronly'},
                    show = function() writeLines(paste('symbolMemberFunction', name)),
                    genCppVar = function(...) {
                        stop(paste('Error, you should not be generating a cppVar for symbolMemberFunction', name))
                    } ))

symbolNimbleListGenerator <-
  setRefClass(Class = 'symbolNimbleListGenerator',
              contains = 'symbolBase',
              fields = list(nlProc = 'ANY'),
              methods = list(
                initialize = function(...){callSuper(...); type <<- 'Ronly'},
                show = function() writeLines(paste('symbolNimbleListGenerator', name)),
                genCppVar = function(...) {
                  return(  cppVarFull(name = name,
                                      baseType = 'Ronly',
                                      templateArgs = nlProc$name) )
                }
              ))

symbolNimbleList <-
    setRefClass(Class = 'symbolNimbleList',
                contains = 'symbolBase',
                fields = list(nlProc = 'ANY'),
                methods = list(
                    initialize = function(...){callSuper(...); type <<- 'symbolNimbleList'},
                    show = function() writeLines(paste('symbolNimbleList', name)),
                    genCppVar = function(...) {
                        return(  cppVarFull(name = name,
                                            baseType = 'nimSmartPtr',
                                            templateArgs = nlProc$name) )
                    }
                    ))

symbolNimbleFunction <-
    setRefClass(Class = 'symbolNimbleFunction',
                contains = 'symbolBase',
                fields = list(nfProc = 'ANY'), 
                methods = list(
                    initialize = function(...) {callSuper(...)},
                    show = function() writeLines(paste('symbolNimbleFunction', name)),
                    genCppVar = function(...) {
                        return(cppVarFull(name = name, baseType = environment(nfProc$nfGenerator)$name, ptr = 1)) ## selfDerefence idea is not general for A$B$C, selfDereference = TRUE))
                    }
                    ))

symbolOptimReadyFunction <- 
	setRefClass(Class = 'symbolOptimReadyFunction',
				contains = 'symbolBase',
				fields = list(nfName = 'ANY', nfProc = 'ANY', genName = 'ANY'),
				methods = list(
					initialize = function(...){callSuper(...)},
					show = function() writeLines(paste('symbolOptimObject', name)),
					genCppVar = function(...){
						cat('note: in genCppVar for symbolOptimReadyFunction, need to figure out how to get unique names\n')
						return(cppVarFull(name = name, baseType = genName, ptr = 0))
					}
					)
				)


symbolVoidPtr <- setRefClass(Class = 'symbolVoidPtr',
                             contains = 'symbolBase',
                             methods = list(
                                 initialize = function(...) callSuper(...),
                                 show = function() writeLines(paste('symbolVoidPtr', name)),
                                 genCppVar = function(...) {
                                     return(cppVarFull(name = name, baseType = 'void', ptr = 1))
                                 }))


symbolModelVariableAccessorVector <- 
    setRefClass(Class = 'symbolModelVariableAccessorVector',
                contains = 'symbolBase',
                fields = list(lengthName = 'ANY'), 		#'character'),
                methods = list(
                    initialize = function(...) {
                        callSuper(...)
                        type <<- 'symbolModelVariableAccessorVector'  
                    },
                    show = function() writeLines(paste('symbolModelVariableAccessorVector', name)),
                    genCppVar = function(...) {
                        return(cppModelVariableMapAccessorVector(name = name))
                    }
                    )
                )

symbolModelValuesAccessorVector <-  
    setRefClass(Class = 'symbolModelValuesAccessorVector',
                contains = 'symbolBase',
                methods = list(
                    initialize = function(...) {
                        callSuper(...)
                        type <<- 'symbolModelValuesAccessorVector'  
                    },
                    show = function() writeLines(paste('symbolModelValuesAccessorVector', name)),
                    genCppVar = function(...) {
                        return(cppModelValuesMapAccessorVector(name = name)) 
                    }
                    )
                )

symbolGetParamInfo <-
    setRefClass(Class = 'symbolGetParamInfo',
                contains = 'symbolBase',
                fields = list(paramInfo = 'ANY'), ## getParam_info, i.e. simple list
                methods = list(
                    initialize = function(paramInfo, ...) {
                        callSuper(...)
                        paramInfo <<- paramInfo
                        type <<- 'Ronly'
                    },
                    show = function() writeLines(paste('symbolGetParamInfo', name)),
                    genCppVar = function(...) {
                        stop(paste('Error, you should not be generating a cppVar for symbolGetParamInfo', name))
                    } ))

symbolGetBoundInfo <-
    setRefClass(Class = 'symbolGetBoundInfo',
                contains = 'symbolBase',
                fields = list(boundInfo = 'ANY'), ## getBound_info, i.e. simple list
                methods = list(
                    initialize = function(boundInfo, ...) {
                        callSuper(...)
                        boundInfo <<- boundInfo
                        type <<- 'Ronly'
                    },
                    show = function() writeLines(paste('symbolGetBoundInfo', name)),
                    genCppVar = function(...) {
                        stop(paste('Error, you should not be generating a cppVar for symbolGetBoundInfo', name))
                    } ))

symbolNumericList <- 
    setRefClass(Class = 'symbolNumericList',
                contains = 'symbolBase', 
                fields = list(		className = 'ANY', 		#'character',
                					nDim ='ANY'), 		# 'numeric'), 
                methods = list(
                    initialize = function(...) {
                        callSuper(...)
                        type <<- 'symbolNumericList'
                    },
                    show = function() writeLines(paste('symbolNumericList', name)),
                    genCppVar = function(...){
                        if(type == 'integer')      cType <- 'int'
                        else if(type == 'double')  cType <- 'double'
                        else if(type == 'logical') cType <- 'bool'
                        else cat(paste0("Warning: in genCppVar method in symbolBasic class, type unrecognized: ", type, '\n'))
                        return(cppVecNimArrPtr(name = name,
                                         nDim = nDim,
                                         ptr = 0,
                                         type = cType))}
                    )
                )

symbolNimPtrList <-
    setRefClass(Class = 'symbolNimPtrList',
                contains = 'symbolBase',
                fields = list(baseClass = 'ANY'),
                methods = list(
                    initialize = function(...) callSuper(...),
                    show = function() writeLines(paste('symbolNimPtrList', name)),
                    genCppVar = function(...) {
                        componentCclassName <- environment(baseClass)$name
                        return(list(
                            cppVarFull(name = name,
                                       baseType = 'vector',
                                       templateArgs = list(cppVar(ptr = 1, baseType = componentCclassName) )
                                       ),
                            cppVarFull(name = paste0(name,'_setter'),
                                       baseType = 'vectorOfPtrsAccess',
                                       templateArgs = list(cppVar(baseType = componentCclassName) )
                                       )
                            ))
                    }
                    ))


symbolCopierVector <-
    setRefClass(Class = 'symbolCopierVector',
                contains = 'symbolBase',
                methods = list(
                    initialize = function(...) callSuper(...),
                    show = function() writeLines(paste('symbolCopierVector', name)),
                    genCppVar = function(...) {
                        cppVar(name = name, baseType = 'copierVectorClass')
                    }
                ))

symbolNimbleFunctionList <-
    setRefClass(Class = 'symbolNimbleFunctionList',
                contains = 'symbolNimPtrList')

symbolEigenMap <- setRefClass(Class = 'symbolEigenMap',
                              contains = 'symbolBase',
                              fields = list(
                                  eigMatrix = 'ANY', 		#'logical', ## or array
                                  strides ='ANY' 		# 'numeric' ## NA for Dynamic.  length(0) means no strides
                                  ),
                              methods = list(
                                  show = function() {
                                      writeLines(paste0(name,': Eigen ',if(eigMatrix) 'Matrix' else 'Array', ' Map of ', type, if(length(strides) > 0) paste0(' with strides ', paste(strides, collapse = ', ')) else character()))
                                  },
                                  genCppVar = function(functionArg = FALSE) {
                                      if(functionArg) stop('Error: cannot take Eigen Map as a function argument (without more work).')
                                      if(length(strides)==2 & eigMatrix) {
                                          if(all(is.na(strides))) {
                                              baseType <- paste0('EigenMapStr', if(type == 'double') 'd' else if(type == 'integer') 'i' else 'b' )
                                              return(cppVarFull(name = name,
                                                                baseType = baseType,
                                                                constructor = '(0,0,0, EigStrDyn(0, 0))',
                                                                ptr = 0,
                                                                static = FALSE))
                                          }
                                      }
                                      cppEigenMap(name = name,
                                                  type = type,
                                                  strides = strides,
                                                  eigMatrix = eigMatrix)
                                  }
                              )
                              )

symbolIndexedNodeInfoTable <-
    setRefClass(Class = "symbolIndexedNodeInfoTable",
                contains = "symbolBase",
                methods = list(
                    initialize = function(...) callSuper(...),
                    show = function() writeLines(paste('symbolIndexedNodeInfoTable', name)),
                    ## We need this to be copied, but it must be copied to a variable already declared in the nodeFun base class,
                    ## so we don't want any genCppVar.
                    genCppVar = function(...)  {
                        cppVarFull(name = name, silent = TRUE) ## this symbol exists to get a base class member data copied, so it shouldn't be declared
                        ##stop(paste('Error, you should not be generating a cppVar for symbolIndexedNodeInfoTable', name))
                        ## looks like if a copy type is created in makeNimbleFxnCppCopyTypes (based on the symbolXXXclass then it will be copied
                        ## and if the type is Ronly then it will not be turned into a cppVar.  So that bit of design worked out well
                       ## it's in the nodeFun base class as vector<indexedNodeInfo>
                    }))

symbolInternalType <-
    setRefClass(Class = "symbolInternalType",
                contains = "symbolBase",
                fields = list(argList = 'ANY'),
                methods = list(
                    initialize = function(...) callSuper(...),
                    show = function() writeLines(paste('symbolInternalType', name)),
                    genCppVar = function(functionArg = FALSE) {
                        if(length(argList) == 0) stop(paste('No information for outputting C++ type of', name))
                        if(argList[[1]] == 'indexedNodeInfoClass'){
                            if(functionArg) return(cppVarFull(name = name, baseType = "indexedNodeInfo", const = TRUE, ref = TRUE))
                            return(cppVar(name = name, baseType = "indexedNodeInfo"))
                        }
                    })
                )

## nDim is set to length(size) unless provided, which is how scalar (nDim = 0) must be set
symbolDouble <- function(name, size = numeric(), nDim = length(size)) {
    if(is.logical(size)) size <- as.numeric(size)
    if(nDim != length(size)) {
        if(nDim != 0 | !identical(size, 1)) stop('Error in symbolDouble, nDim must be length(size) unless nDim == 0 and size == 1')
    }
    symbolBasic(name = name, type = 'double', nDim = nDim, size = size)
}

symbolInt <- function(name, size = numeric(), nDim = length(size)) {
    if(is.logical(size)) size <- as.numeric(size)
    if(nDim != length(size)) {
        if(nDim != 0 | !identical(size, 1)) stop('Error in symbolInt, nDim must be length(size) unless nDim == 0 and size == 1')
    }
    symbolBasic(name = name, type = 'int', nDim = nDim, size = size)
}

symbolTable <- 
    setRefClass(Class   = 'symbolTable',
                fields  = list(symbols  = 'ANY', 		#'list',
                    parentST = 'ANY',
                    dimAndSizeList = 'ANY',
                    dimAndSizeListMade = 'ANY'),
                methods = list(
                    initialize = function(parentST = NULL, ...) {
                        symbols  <<- list()
                        dots <- list(...)
                        if('symbols' %in% names(dots)) {
                            if(is.list(dots$symbols)) {
                                for(s in dots$symbols) {
                                    if(inherits(s$name, 'uninitializedField')) stop(paste0('Error: all symbols in list must have meaningful name fields'))
                                    if(identical(s$name, character())) stop(paste0('Error: all symbols in list must have meaningful name fields'))
                                    symbols[[s$name]] <<- s
                                }
                            } else stop('Error: symbols provided must be a list')
                        }
                        parentST <<- parentST
                        dimAndSizeListMade <<- FALSE
                    },
                    
                    ## add a symbol RC object to this symbolTable; checks for valid symbolRC object, and duplicate symbol names
                    addSymbol  = function(symbolRCobject) {
                      ##  if(!is(symbolRCobject, 'symbolBase'))   stop('adding non-symbol object to symbolTable')
                        name <- symbolRCobject$name
                        if(name %in% getSymbolNames())            warning(paste0('duplicate symbol name: ', name))
                        symbols[[name]] <<- symbolRCobject
                        if(dimAndSizeListMade) {
                            dimAndSizeList[[name]] <<- {ans <- try(list(symbolRCobject$size, symbolRCobject$nDim)); if(inherits(ans, 'try-error')) NULL else ans}
                        }
                    },
                    ## remove a symbol RC object from this symbolTable; gives warning if symbol isn't in table
                    removeSymbol = function(name) {
                        if(!(name %in% getSymbolNames()))         warning(paste0('removing non-existant symbol name: ', name))
                        symbols[[name]] <<- NULL },
                    
                    ## symbol accessor functions
                    getLength = function() return(length(symbols)),
                    getSymbolObjects = function()      return(symbols),
                    getSymbolNames   = function()      if(is.null(names(symbols))) return(character(0)) else return(names(symbols)),
                    getSymbolObject  = function(name, inherits = FALSE) {
                        ans <- symbols[[name]]
                        if(is.null(ans)) if(inherits) if(!is.null(parentST)) ans <- parentST$getSymbolObject(name, TRUE)
                        return(ans)
                    },
                    symbolExists = function(name, inherits = FALSE) {
                        return(!is.null(getSymbolObject(name, inherits)))
                    },
                    initDimAndSizeList = function() {
                        dimAndSizeList <<- lapply(symbols, function(x) {
                            ans <- try(list(x$size, x$nDim))
                            if(inherits(ans, 'try-error')) NULL else ans
                        })
                        dimAndSizeListMade <<- TRUE
                    },
                    makeDimAndSizeList = function(names) {
                        if(!dimAndSizeListMade) initDimAndSizeList()
                        dimAndSizeList[names]
                    },
                    getSymbolType    = function(name)               return(symbols[[name]]$type),
                    getSymbolField   = function(name, field)        return(symbols[[name]][[field]]),
                    setSymbolField   = function(name, field, value) symbols[[name]][[field]] <<- value,
                    
                    ## parentST accessor functions
                    getParentST = function()   return(parentST),
                    setParentST = function(ST) parentST <<- ST,
                    show = function() {
                        writeLines('symbol table:')
                        for(i in seq_along(symbols)) symbols[[i]]$show()
                        if(!is.null(parentST)) {
                            writeLines('parent symbol table:')
                            parentST$show()
                        }
                    }
                )
    )




areMVSymTabsEqual <- function(symTab1, symTab2) {
  vN1 = symTab1$getSymbolNames()
  vN2 = symTab2$getSymbolNames()
  if(length(vN1) != length(vN2) )
    return(FALSE)
  for(i in seq_along(vN1) ) {
    if(vN1[i] != vN2[i])
      return(FALSE)
    if(symTab1$symbols[[i]]$type != symTab2$symbols[[i]]$type)
      return(FALSE)
    if( !all( is(symTab1$symbols[[i]]) == is(symTab2$symbols[[i]]) ) )
      return(FALSE)
    if( inherits(symTab1$symbols[[i]], "symbolBasic")) {
      nDim = symTab1$symbols[[i]]$nDim
      if(nDim != symTab2$symbols[[i]]$nDim )
        return(FALSE)
      size1 = symTab1$symbols[[i]]$size
      size2 = symTab2$symbols[[i]]$size
      if(any(size1 != size2) ) 
        return(FALSE)
    }
  }
  return(TRUE)
}

nl_refClassLabelMaker <- labelFunctionCreator('nimListClass')

nimbleListDefClass <- setRefClass(
    ## This class holds a list of type information such as
    ## A = double(1), B = integer(2)
    ## The types need not be numeric.
    ## In general, ideally, they could be another nimbleList or a nimbleFunction
    Class = "nimbleListDefClass",
    fields = list(types = 'ANY',
                  className = 'ANY')
)

nimbleListBase <- setRefClass(Class = 'nimbleListBase', 
                                  fields = list(
                                      .CobjectInterface = 'ANY',
                                      .generatorFunction = 'ANY',
                                      nimbleListDef = 'ANY',
                                      nestedListGenList = 'ANY'
                                  ),
                                  methods = list(
                                    initialize = function(...)
                                      callSuper(...)
                                  ))

#' create a nimbleType object
#'
#' create a nimbleType object, with information on the name, type, and dimension of an object to be placed in a \link{nimbleList} 
#'
#' @param name The name of the object
#' @param type The type of the object
#' @param dim  The dimension of the object.  This can be left blank if the object is a nimbleList.
#'
#' @author NIMBLE development team
#'
#' @export
#'
#' @details
#' 
#' This function creates \code{nimbleType} objects, which can be used to define the elements of a \link{nimbleList}.  
#' 
#' The \code{type} argument can be chosen from among \code{character}, \code{double}, \code{integer}, and \code{logical},
#' or can be the name of a previously created \link{nimbleList} definition.
#' 
#' See the NIMBLE User Manual for examples.
#'


nimbleType <- setRefClass(
  Class = 'nimbleType',
  fields = c('name', 'type', 'dim'),
  methods = list(
    initialize = function(name, type, dim = NA){
      name <<- name
      type <<- type
      dim <<- dim
    },
    show = function(){
      cat("nimbleType object with name ", name, ", type ", type, ", dim ",
          dim,"\n", sep = "")
    }
  )
)


#' create a nimbleList
#'
#' create a nimbleList from a nimbleList definition 
#'
#' @param types objects defining the names, types, and dimensions of the nimbleList elements.  
#' @param name An optional name used internally, for example in generated C++ code.  Usually this is left blank and NIMBLE provides a name.
#' @param where An optional \code{where} argument passed to \code{setRefClass} for where the reference class definition generated for this nimbleFunction will be stored.  This is needed due to R package namespace issues but should never need to be provided by a user.
#'
#' @author NIMBLE development team
#'
#' @export
#'
#' @details
#' This function creates a definition for a nimbleList.  The \code{types} argument defines the names, types, and dimensions of the elements of the nimbleList.  Elements of nimbleLists can be either basic types (e.g. integer, double) or other nimbleList definitions.   
#' The \code{types} argument can be either a series of expressions of the form \code{name = type(dim)}, or a list of \link{nimbleType} objects.
#' 
#' \code{nimbleList} returns a definition, which can be used to create instances of this type of nimbleList via the \code{new()} member function. 
#' 
#' Definitions can be created in \code{R}'s general environment or in \code{nimbleFunction} setup code.  Instances can be created using the \code{new()} function in \code{R}'s global environment, in \code{nimbleFunction} setup code, or in \code{nimbleFunction} run code.  
#' 
#' Instances of \code{nimbleList} definitions can be used as arguments to run code of \code{nimbleFunction}s, and as the return type of \code{nimbleFunction}s.
#' 
#' See the NIMBLE User Manual for examples.
#'
nimbleList <- function(...,
                       name = NA,
                       where =  getNimbleFunctionEnvironment()) {
    ## This has a role like nimbleFunction but a much simpler implementation
    ## It returns a function that simply makes a regular R list and
    ## attaches two attributes, one to mark it as a nimbleList (for efficienct checking
    ## compatible with checking of other objects that have a class) and
    ## one that has the nimbleListDefClass object
  
  ## 3 possibilities: arguments as expressions, arguments as list created within call,
  ## arguments as list created outside of call
  
    Call <-  match.call(expand.dots = TRUE)
    if(any(names(Call) == 'name')){
      Call <- Call[-which(names(Call) == 'name')]
    }
    if(length(Call) < 2)
      stop("No arguments specified for nimbleList")
    argList <- list()
    
    ## left side of || statement catches list(...) arguments
    ## right side captures name of a previously defined list
    if((is.call(Call[[2]]) && deparse(Call[[2]][[1]]) == 'list') || 
       (!is.call(Call[[2]]) && is.list(eval(Call[[2]], envir = parent.frame())))){ 
      callList <- eval(Call[[2]], envir = parent.frame())
      for(iArg in seq_along(callList)){
        argList[[iArg]] <- list(name = callList[[iArg]]$name,
                                        type = callList[[iArg]]$type,
                                        dim = callList[[iArg]]$dim)
      }
    }
    else{  ## if arguments are expressions e.g. nimListDouble = double(2)
      for(iArg in 2:length(Call)){
        argList[[iArg-1]] <- list(name = names(Call)[iArg],
                                            type = deparse(Call[[iArg]][[1]]))
        argList[[iArg-1]]$dim  <- if(length(Call[[iArg]])>1) deparse(Call[[iArg]][[2]])
                                            else 0
      }
    }
    
    types <- list(vars = sapply(argList, function(x){return(x$name)}),
                  types =  sapply(argList, function(x){return(x$type)}),
                  dims =  sapply(argList, function(x){return(x$dim)}))
    
    if(is.na(name)) name <- nf_refClassLabelMaker()
    nlDefClassObject <- nimbleListDefClass(types = types, className = name) 
    basicTypes <- c("double", "integer", "character", "logical")
    nestedListGens <- list()
    for(i in seq_along(types$types)){
      if(!(types$types[i] %in% basicTypes)){
        if(types$types[i] %in% c('eigen', 'nimEigen')){  ## eigen() will not have been converted to nimEigen() yet, so
                                                         ## need to check for both
          nestedListGens[[types$vars[i]]] <- nlEigenReferenceList[['nimEigen']]$createListDef()
        }
        if(types$types[i] %in% c('svd', 'nimSvd')){
          nestedListGens[[types$vars[i]]] <- nlEigenReferenceList[['nimSvd']]$createListDef()
        }
        for(searchEnvironment in c(parent.frame(), globalenv())){
          if(is.nlGenerator(get(types$types[i], envir = searchEnvironment))){
            nestedListGens[[types$vars[i]]] <- get(types$types[i], envir = searchEnvironment)
            break
          }
        }
      }
    }
    
    classFields <- as.list(rep('ANY', length(types$vars)))
    names(classFields) <- types$vars
    ## classFields[[length(classFields)+1]] <- "ANY"
    ## names(classFields)[length(classFields)] <- "nimbleListDef" ## initial nl definition stored here
    ## classFields[[length(classFields)+1]] <- "ANY"
    ## names(classFields)[length(classFields)] <- "nestedListGenList" ## nl generators for any nested lists stored here
    
    nlRefClass <- setRefClass(
      Class = name,
      fields = classFields,
      contains = 'nimbleListBase',
      methods = list(
        initialize = function(NLDEFCLASSOBJECT, NESTEDGENLIST, ...){
          nimbleListDef <<- NLDEFCLASSOBJECT
          nestedListGenList <<- NESTEDGENLIST
          for(i in seq_along(nestedListGenList)){
            .self[[names(nestedListGenList)[i]]] <- nestedListGenList[[i]]$new()
          }
          nimListFields <- nimbleListDef$types$vars
          initializeFields <- list(...)
          nonInitializeFields <- which(!(nimListFields %in% c(names(nestedListGenList), names(initializeFields))))
          ## initialize uninitialized fields
          for(i in nonInitializeFields){
            thisType <- nimbleListDef$types$types[i]
            thisDim <-  nimbleListDef$types$dims[i]
            if(thisType == 'character'){
              .self[[nimListFields[i]]] <- ""
            }
            else if(thisType %in% c('integer', 'double')){
              if(thisDim == 0)
                .self[[nimListFields[i]]] <- 0
              if(thisDim == 1)
                .self[[nimListFields[i]]] <- integer(0)
              if(thisDim == 2)
                .self[[nimListFields[i]]] <- matrix(0, 0, 0)
              if(thisDim > 2)
                .self[[nimListFields[i]]] <- array(0, dim = rep(0, thisDim))
            }
            else if(thisType == 'logical'){
              .self[[nimListFields[i]]] <- FALSE
            }
          }     
          callSuper(...)
        },
        show = function(){
          cat("nimbleList object of type ", .self$nimbleListDef$className, 
              "\n", sep = "")
          nimListPrintFields <- nimbleListDef$types$vars
          for(fieldName in nimListPrintFields){
            cat("Field \"", fieldName, "\":\n", sep = "")
            cat(methods::show(field(fieldName)))
          }
        }
      ),
      where = where
    )

    nlGeneratorFunction <-   eval(  substitute(
      function(...){
      return(nlRefClass(NLDEFCLASSOBJECT, NESTEDGENLIST, ..., .generatorFunction = nlGeneratorFunction))},
      list(NLDEFCLASSOBJECT = nlDefClassObject,
           NESTEDGENLIST = nestedListGens)))
    return(list(new = nlGeneratorFunction))
}

## nimbleList processing class
## analogous to but must simpler than NFprocessing
nlProcessing <- setRefClass('nlProcessing',
                            fields = list(
                                cppDef = 'ANY',
                                nimbleListObj = 'ANY',
                                symTab = 'ANY',
                                neededTypes = 'ANY',
                                nimbleProject = 'ANY',
                                name = 'ANY',
                                ##   instances = 'ANY',
                                nlGenerator = 'ANY',
                                nestedListGens = 'ANY',
                                neededObjectNames =  'ANY'		#'character', ## a character vector of the names of objects such as models or modelValues that need to exist external to the nimbleFunction object so their contents can be pointed to 
                            ),
                            methods = list(
                                show = function() {
                                    writeLines(paste0('nlProcessing object ', nimbleListObj$className))
                                },
                                initialize = function(nimLists = NULL, className, project, ...) {
                                    ## modifying this so nimLists is allowed to be a nlGenerator.  That way we don't need to create objects just to access their definition information.
                                    ## nimLists can also be a nimbleList object or list of them (all from same generator)
                                  neededTypes <<- list()
                                  callSuper(...)
                                  if(!is.null(nimLists)) {
                                      nimbleProject <<- project
                                      if(is.nlGenerator(nimLists)) {
                                          sl <- nl.getDefinitionContent(nimLists, 'nlDefClassObject')
                                          nlGenerator <<- nimLists
                                      } else {
                                          if(is.list(nimLists)) {
                                              sl <- nimLists[[1]]$nimbleListDef
                                              nlGenerator <<- nimLists[[1]]$.generatorFunction
                                          } else {
                                              sl <- nimLists$nimbleListDef
                                              nlGenerator <<- nimLists$.generatorFunction
                                          }
                                      }
                                      
                                    nimbleListObj <<- sl
                                    if(missing(className)) {
                                      name <<- sl$className
                                    } else {
                                      name <<- className
                                    }
                                   ## instances <<- if(inherits(nimLists, 'list')) nimLists else list(nimLists)
                                      nestedListGens <<- environment(nlGenerator)$nestedListGenList
                                      ##nestedListGens <<- if(inherits(nimLists, 'list')) nimLists[[1]]$nestedListGenList else nimLists$nestedListGenList
                                  }
                                },
                                setupTypesForUsingFunction= function() buildSymbolTable(), ## required name
                                process = function(control = list(debug = FALSE, debugCpp = FALSE)) {
                                    debug <- control$debug
                                    debugCpp <- control$debugCpp
                                    if(!is.null(nimbleOptions()$debugNFProcessing)) {
                                        if(nimbleOptions()$debugNFProcessing) {
                                            debug <- TRUE
                                            control$debug <- TRUE
                                            writeLines('Debugging nfProcessing (nimbleOptions()$debugRCfunProcessing is set to TRUE)') 
                                        }
                                    }
                                    if(debug) browser()
                                    if(inherits(symTab, "uninitializedField")) buildSymbolTable()
                                },
                                buildSymbolTable = function() {
                                  if(!inherits(symTab, "uninitializedField")) return(warning("Symbol Table for nimbleList already built."))
                                  if(length(nimbleListObj$types$types) != length(nimbleListObj$types$vars))
                                    stop("Number of nimbleList vars provided is not equal to number of nimbleList types provided")
                                  symTab <<- symbolTable()
                                  for(i in seq_along(nimbleListObj$types$vars)){
                                      nimbleListObjVar <- nimbleListObj$types$vars[i]
                                    if(nimbleListObjVar %in% names(nestedListGens)){
                                      ##nlList <- nestedListGens[[nimbleListObj$types$vars[i]]]$new()
                                        thisNestedListGen <- nestedListGens[[nimbleListObjVar]]
                                        thisNestedListDef <- nl.getDefinitionContent(thisNestedListGen, 'nlDefClassObject')
                                        ##className <- nlList$nimbleListDef$className
                                        className <- thisNestedListDef$className
                                      nlp <- nimbleProject$nlCompInfos[[className]]$nlProc
                                      newSym <- symbolNimbleList(name = nimbleListObjVar, type = 'symbolNimbleList', nlProc = nlp)
                                      neededTypes[[className]] <<- newSym  ## if returnType is a NLG, this will ensure that it can be found in argType2symbol()
                                      symTab$addSymbol(newSym)
                                    }
                                    else{
                                     nimbleListObjType <- nimbleListObj$types$types[i]
                                     nimbleListObjDim <-  as.numeric(nimbleListObj$types$dims[i])
                                     symTab$addSymbol(argType2symbol(call(nimbleListObjType, nimbleListObjDim),
                                                                     neededTypes, nimbleListObjVar))
                                    }
                                  }
                                },
                                getSymbolTable = function() symTab
                            ))


nlEigenClass <- setRefClass('nlEigenClass',
                            fields = list(
                              className = 'ANY',
                              funcName = 'ANY',
                              nimFuncName = 'ANY',
                              listElements = 'ANY',
                              eigenNimbleListDef = 'ANY'),
                            methods = list(
                              initialize = function(...){
                                callSuper(...)
                                createListDef()
                              },
                              addEigenListInfo = function(nfProc){
                                thisProj <- nfProc$nimbleProject
                                eigenNimbleList <- eigenNimbleListDef$new() 
                                nlp <- thisProj$compileNimbleList(eigenNimbleList, initialTypeInferenceOnly = TRUE)
                                eigenListSym <- symbolNimbleList(name = className, nlProc = nlp)
                                nfProc$neededTypes[[className]] <- eigenListSym 
                                nfProc$setupSymTab$addSymbol(symbolNimbleListGenerator(name = nimFuncName, nlProc = nlp))
                              },
                              createListDef = function(){
                                eigenNimbleListDef <<- nimbleList(listElements, name = className)
                              }
                            ))

nlEigenEigenInfo <- nlEigenClass(funcName = 'nimEigen',
                                 className = 'EIGEN_EIGENCLASS',
                                 nimFuncName = 'EIGEN_EIGEN',
                                 listElements = list(nimbleType('values', 'double', 1),
                                                     nimbleType('vectors', 'double', 2)))
nlEigenSvdInfo    <- nlEigenClass(funcName = 'nimSvd',
                                  className = 'EIGEN_SVDCLASS',
                                  nimFuncName = 'EIGEN_SVD',
                                  listElements = list(nimbleType('d', 'double', 1),
                                                      nimbleType('u', 'double', 2),
                                                      nimbleType('v', 'double', 2)))

nlEigenReferenceList <- list(nimEigen = nlEigenEigenInfo,
                             nimSvd = nlEigenSvdInfo)



is.nl <- function(f){
  if(inherits(f, 'nimbleListBase')) return(TRUE)
  return(FALSE)
}

is.nlGenerator <- function(x, inputIsName = FALSE) {
    if(inputIsName) x <- get(x)
    if(is.list(x) && is.function(x$new)) {
        if(is.null(environment(x$new))) return(FALSE)
        if(exists('nlDefClassObject', envir = environment(x$new), inherits = FALSE)) return(TRUE)
    }
    FALSE
}

nl.getDefinitionContent <- function(nlGen, name) {
    environment(testNL$new)[[name]]
}

nl.getNestedGens <- function(nlGen) {
    environment(testNL$new)$nestedListGens
}

nl_refClassLabelMaker <- labelFunctionCreator('nimListClass')

nimbleListDefClass <- setRefClass(
    ## This class holds a list of type information such as
    ## A = double(1), B = integer(2)
    ## The types need not be numeric.
    ## In general, ideally, they could be another nimbleList or a nimbleFunction
    Class = "nimbleListDefClass",
    fields = list(types = 'ANY',
                  className = 'ANY',
                  predefined = 'ANY')
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
#' Create a nimbleType object, with information on the name, type, and dimension of an object to be placed in a \code{\link{nimbleList}}.
#'
#' @param name The name of the object, given as a character string.
#' @param type The type of the object, given as a character string.
#' @param dim  The dimension of the object, given as an integer.  This can be left blank if the object is a nimbleList.
#'
#' @author NIMBLE development team
#'
#' @export
#'
#' @details
#' 
#' This function creates \code{nimbleType} objects, which can be used to define the elements of a \code{\link{nimbleList}}.  
#' 
#' The \code{type} argument can be chosen from among \code{character}, \code{double}, \code{integer}, and \code{logical},
#' or can be the name of a previously created \code{\link{nimbleList} definition}.
#' 
#' See the NIMBLE \href{https://r-nimble.org/html_manual/cha-welcome-nimble.html}{User Manual} for additional examples.
#' 
#' @examples 
#' nimbleTypeList <- list()
#' nimbleTypeList[[1]] <- nimbleType(name = 'x', type = 'integer', dim = 0)
#' nimbleTypeList[[2]] <- nimbleType(name = 'Y', type = 'double', dim = 2)
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
#' @param ... arbitrary set of names and types for the elements of the list or a single R list of type \code{nimbleType}.
#' @param name optional character providing a name used internally, for example in generated C++ code.  Usually this is left blank and NIMBLE provides a name.
#' @param predefined logical for internal use only.
#' @param where optional argument passed to \code{setRefClass} for where the reference class definition generated for this nimbleFunction will be stored.  This is needed due to R package namespace issues but should never need to be provided by a user.
#'
#' @author NIMBLE development team
#'
#' @export
#'
#' @details
#' This function creates a definition for a nimbleList.  The \code{types} argument defines the names, types, and dimensions of the elements of the nimbleList.  Elements of nimbleLists can be either basic types (e.g., \code{integer}, \code{double}) or other nimbleList definitions.   
#' The \code{types} argument can be either a series of expressions of the form \code{name = type(dim)}, or a list of \code{\link{nimbleType}} objects.
#' 
#' \code{nimbleList} returns a definition, which can be used to create instances of this type of nimbleList via the \code{new()} member function. 
#' 
#' Definitions can be created in R's general environment or in nimbleFunction setup code.  Instances can be created using the \code{new()} function in R's global environment, in nimbleFunction setup code, or in nimbleFunction run code.  
#' 
#' Instances of \code{nimbleList} definitions can be used as arguments to run code of nimbleFunctions, and as the return type of nimbleFunctions.
#' @examples 
#'  exampleNimListDef <- nimbleList(x = integer(0), Y = double(2))
#'  
#'  nimbleListTypes <- list(nimbleType(name = 'x', type = 'integer', dim = 0),
#'                          nimbleType(name = 'Y', type = 'double', dim = 2))
#'  
#'  ## this nimbleList definition is identical to the one created above
#'  exampleNimListDef <- nimbleList(nimbleListTypes)
nimbleList <- function(...,
                       name = NA,
                       predefined = FALSE,
                       where =  getNimbleFunctionEnvironment()) {
    ## This has a role like nimbleFunction but a much simpler implementation
    ## It returns a function that simply makes a regular R list and
    ## attaches two attributes, one to mark it as a nimbleList (for efficienct checking
    ## compatible with checking of other objects that have a class) and
    ## one that has the nimbleListDefClass object

    ## This manual override allows us to generate static code by temporarily setting
    ## predefined = FALSE for all predefined nimbleLists.
    GENERATE_STATIC_CODE <- FALSE  ## Enable this before using generateStaticCode.R.
    if(GENERATE_STATIC_CODE) predefined <- FALSE
    
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
    nlDefClassObject <- nimbleListDefClass(types = types, className = name, predefined = predefined) 
    basicTypes <- c("double", "integer", "character", "logical")
    nestedListGens <- list()
    for(i in seq_along(types$types)){
        if(!(types$types[i] %in% basicTypes)){
        for(searchEnvironment in c(parent.frame(), globalenv())){
          if(try(is.nlGenerator(get(types$types[i], envir = searchEnvironment)), silent = TRUE)){
            nestedListGens[[types$vars[i]]] <- get(types$types[i], envir = searchEnvironment)
            break
          }
        }
      }
    }
    
    classFields <- as.list(rep('ANY', length(types$vars)))
    names(classFields) <- types$vars
    nlRefClass <- setRefClass(
      Class = name,
      fields = classFields,
      contains = 'nimbleListBase',
      methods = list(
        initialize = function(...){
          callSuper(...)
          nimListFields <- nimbleListDef$types$vars
          initializeFields <- list(...)
          nonInitializeFields <- which(!(nimListFields %in% names(initializeFields)))
          ## initialize uninitialized fields
          for(i in nonInitializeFields){
            thisType <- nimbleListDef$types$types[i]
            thisDim <-  nimbleListDef$types$dims[i]
            if(thisType == 'character'){
              initValue  <- ""
            }
            else if(thisType %in% c('integer', 'double')){
              if(thisDim == 0)
                initValue <- 0
              if(thisDim == 1)
                initValue <- integer(0)
              if(thisDim == 2)
                initValue <- matrix(0, 0, 0)
              if(thisDim > 2)
                initValue <- array(0, dim = rep(0, thisDim))
            }
            else if(thisType == 'logical'){
              initValue <- FALSE
            }
            else if(nimListFields[i] %in% names(nestedListGenList)){
              initValue <- nestedListGenList[[nimListFields[i]]]$new()
            }
            else(stop(paste("unrecognized type given for nimbleList element", nimListFields[i])))
            eval(substitute(.self[[nimListFields[i]]] <<-initValue))
          }   
        },
        show = function(){
          cat("nimbleList object of type ", nimbleListDef$className, 
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
    
    nlGeneratorFunction <- function(...){
        return(nlRefClass(nimbleListDef = nlDefClassObject, nestedListGenList = nestedListGens, .generatorFunction = nlGeneratorFunction,
                          ...))}
    nlGenerator <- list(new = nlGeneratorFunction)
    return(nlGenerator)
}

makeNimbleListTemplateWithBlankFirstArg <- function(nlDef) {
    vars <- c('.LEFTSIDE', nlDef$types$vars)
    functionAsList <- list(as.name('function'))
    functionAsList[2] <- list(NULL)
    if(length(vars) > 0) {
        argsList <- nf_createAList(vars)
        functionAsList[[2]] <- as.pairlist(argsList)
    }
    functionAsList[[3]] <- quote({})
    eval(as.call(functionAsList))
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
                                              nlGenerator <<- nl.getGenerator(nimLists[[1]])
                                          } else {
                                              sl <- nimLists$nimbleListDef
                                              nlGenerator <<- nl.getGenerator(nimLists)
                                          }
                                      }
                                      
                                    nimbleListObj <<- sl
                                    if(missing(className)) {
                                      name <<- sl$className
                                    } else {
                                      name <<- className
                                    }
                                      nestedListGens <<- nl.getNestedGens(nlGenerator)
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
                                        thisNestedListGen <- nestedListGens[[nimbleListObjVar]]
                                        thisNestedListDef <- nl.getDefinitionContent(thisNestedListGen, 'nlDefClassObject')
                                        className <- thisNestedListDef$className
                                      nlp <- nimbleProject$nlCompInfos[[className]]$nlProc
                                      newSym <- symbolNimbleList(name = nimbleListObjVar, nlProc = nlp)
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




## Below are nimbleList definitions for predefined nimbleLists in nimble
## Note that currently, any nimbleList definition that has "predefined = TRUE" must have existing c++ code that defines
## the c++ class.


#' eigenNimbleList definition
#' 
#' \code{nimbleList} definition for the type of \code{nimbleList} returned by \code{\link{nimEigen}}.
#' 
#' @author NIMBLE development team
#'
#' @export
#'
#' @seealso  \code{\link{nimEigen}} 
eigenNimbleList <- nimbleList(list(nimbleType('values', 'double', 1),
                                   nimbleType('vectors', 'double', 2)), name = "EIGEN_EIGENCLASS", predefined = TRUE)


#' svdNimbleList definition
#' 
#' \code{nimbleList} definition for the type of \code{nimbleList} returned by \code{\link{nimSvd}}.
#' 
#' @author NIMBLE development team
#'
#' @export
#' 
#' @seealso  \code{\link{nimSvd}} 
svdNimbleList <-  nimbleList(list(nimbleType('d', 'double', 1),
                                  nimbleType('u', 'double', 2),
                                  nimbleType('v', 'double', 2)), name = "EIGEN_SVDCLASS", predefined = TRUE)


#' waicList definition
#' 
#' \code{waicList} definition for the \code{nimbleList} type returned by WAIC
#' computation.
#'
#' @details
#'
#' See \code{help(waic)} for details on the elements of the list.
#' 
#' @author NIMBLE development team
#'
#' @export
#' 
waicList <- nimbleList(
    list(
        nimbleType('WAIC', 'double', 0),
        nimbleType('lppd', 'double', 0),
        nimbleType('pWAIC', 'double', 0)
    ), name = 'waicList',
    predefined = TRUE
)
    
#' waicDetailsList definition
#' 
#' \code{waicDetailsList} definition for the \code{nimbleList} type returned by WAIC
#' computation.
#'
#' @details
#'
#' See \code{help(waic)} for details on the elements of the list.
#' 
#' @author NIMBLE development team
#'
#' @export
#' 
waicDetailsList <- nimbleList(
    list(
        nimbleType('marginal', 'logical', 0),
        nimbleType('niterMarginal', 'double', 0),
        nimbleType('thin', 'logical', 0),
        nimbleType('online', 'logical', 0),

        ## values for shorter MC runs to assess convergence for marginal calculation
        nimbleType('WAIC_partialMC', 'double', 1),
        nimbleType('lppd_partialMC', 'double', 1),
        nimbleType('pWAIC_partialMC', 'double', 1),
        nimbleType('niterMarginal_partialMC', 'double' , 1),  # checkIts

        ## per data group values potentially useful for SE for contrasting WAIC of two models
        nimbleType('WAIC_elements', 'double', 1),
        nimbleType('lppd_elements', 'double', 1),
        nimbleType('pWAIC_elements', 'double', 1)

    ), name = 'waicDetailsList',
    predefined = TRUE
)


#' EXPERIMENTAL Data type for the return value of \code{\link{nimDerivs}}
#'
#' \code{\link{nimbleList}} definition for the type of \code{\link{nimbleList}} returned by \code{\link{nimDerivs}}.
#'
#' @field value The value of the function evaluated at the given input arguments. 
#' @field gradient	The gradient of the function evaluated at the given input arguments. 
#' @field hessian The Hessian of the function evaluated at the given input arguments. 
#' @field thirdDerivs Currently unused.
#'
#' @export
#' @seealso \code{\link{nimDerivs}}
ADNimbleList <-  nimbleList(list(nimbleType('value', 'double', 1),
                                 nimbleType('gradient', 'double', 2),
                                 nimbleType('hessian', 'double', 3),
                                 nimbleType('thirdDerivs', 'double', 4)),
                            name = "NIMBLE_ADCLASS", predefined = TRUE)

#' EXPERIMENTAL Data type for the return value of \code{\link{nimOptim}}
#'
#' \code{\link{nimbleList}} definition for the type of \code{\link{nimbleList}} returned by \code{\link{nimOptim}}.
#'
#' @field par The best set of parameters found.
#' @field value	The value of fn corresponding to par.
#' @field counts A two-element integer vector giving the number of calls to fn and gr respectively.
#' @field convergence An integer code. 0 indicates successful completion. Possible error codes are
#'        1 indicates that the iteration limit maxit had been reached.
#'        10 indicates degeneracy of the Nelder-Mead simplex.
#'        51 indicates a warning from the "L-BFGS-B" method; see component message for further details.
#'        52 indicates an error from the "L-BFGS-B" method; see component message for further details.
#' @field message A character string giving any additional information returned by the optimizer, or NULL.
#' @field hessian Only if argument hessian is true. A symmetric matrix giving an estimate of the Hessian at the solution found.
#'
#' @export
#' @seealso \code{\link{optim}}, \code{\link{nimOptim}}
optimResultNimbleList <- nimbleList(
    list(
        nimbleType('par', 'double', 1),
        nimbleType('value', 'double', 0),
        nimbleType('counts', 'integer', 1),
        nimbleType('convergence', 'integer', 0),
        nimbleType('message', 'character', 0),
        nimbleType('hessian', 'double', 2)
    ),
    name = "OptimResultNimbleList",
    predefined = TRUE
)

#' EXPERIMENTAL Data type for the \code{control} parameter of \code{\link{nimOptim}}
#'
#' \code{\link{nimbleList}} definition for the type of \code{\link{nimbleList}} input as the \code{control} parameter
#' to \code{\link{nimOptim}}. See \code{\link{optim}} for details.
#' 
#' @export
#' @seealso \code{\link{optim}}, \code{\link{nimOptim}}
optimControlNimbleList <- nimbleList(
    list(
        nimbleType('trace', 'integer', 0),
        nimbleType('fnscale', 'double', 0),
        nimbleType('parscale', 'double', 1),
        nimbleType('ndeps', 'double', 1),
        nimbleType('maxit', 'integer', 0),
        nimbleType('abstol', 'double', 0),
        nimbleType('reltol', 'double', 0),
        nimbleType('alpha', 'double', 0),
        nimbleType('beta', 'double', 0),
        nimbleType('gamma', 'double', 0),
        nimbleType('REPORT', 'integer', 0),
        nimbleType('type', 'integer', 0),
        nimbleType('lmm', 'integer', 0),
        nimbleType('factr', 'double', 0),
        nimbleType('pgtol', 'double', 0),
        nimbleType('temp', 'double', 0),
        nimbleType('tmax', 'integer', 0)
    ),
    name = "OptimControlNimbleList",
    predefined = TRUE
)

## any DSL functions that return nimbleLists should be added to the list below, in the form:
## functionName = list(nlGen = nimbleList definition, cppName = name of cpp function corresponding to dsl function)
nimbleListReturningFunctionList <- list(nimEigen = list(nlGen = eigenNimbleList, cppName = 'EIGEN_EIGEN'),
                                        nimSvd = list(nlGen = svdNimbleList, cppName = "EIGEN_SVD"),
                                        nimDerivs = list(nlGen = ADNimbleList, cppName = "NIM_DERIVS"),
                                        getDerivs = list(nlGen = ADNimbleList, cppName = 'getDerivs'),
                                        nimOptim = list(nlGen = optimResultNimbleList, cppName = "OptimResultNimbleList"),
                                        nimOptimDefaultControl = list(nlGen = optimControlNimbleList, cppName = "OptimControlNimbleList"))


## TODO Add nimbleList definitions for nimOptimResult and nimOptimControl.


#' check if a nimbleList
#'
#' Checks an object to determine if it is a nimbleList (i.e., a list created by \code{nlDef$new()}).
#'
#' @param l object to be tested
#'
#' @seealso \code{\link{nimbleList}} for how to create a nimbleList
#' @export
is.nl <- function(l){
  if(inherits(l, 'nimbleListBase')) return(TRUE)
  return(FALSE)
}

is.nlGenerator <- function(x, inputIsName = FALSE, where = -1) {
    if(inputIsName) x <- get(x, pos = where)
    if(is.list(x) && is.function(x$new)) {
        if(is.null(environment(x$new))) return(FALSE)
        if(exists('nlDefClassObject', envir = environment(x$new), inherits = FALSE)) return(TRUE)
    }
    FALSE
}

nl.getGenerator <- function(nl) {
    environment(nl$.generatorFunction)$nlGenerator
}

nl.getDefinitionContent <- function(nlGen, name) {
    environment(nlGen$new)[[name]]
}

nl.getNestedGens <- function(nlGen) {
    environment(nlGen$new)$nestedListGens
}

nl.getListDef <- function(nlGen) {
    environment(nlGen$new)$nlDefClassObject
}

## makeNewNimListSEXPRESSIONFromC is added to nimbleInternalFunctions and called from c++ function makeNewNimbleList
## nimbleUerNamespace$nimListGens is populated in buildRwrapperFunCode method of RCfunctionDef ref class
makeNewNimListSEXPRESSIONFromC <- function(name){
  returnList <- nimbleUserNamespace$nimListGens[[name]]$new()
  return(returnList)
}



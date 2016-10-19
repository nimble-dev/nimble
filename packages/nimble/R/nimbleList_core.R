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
                                    .CobjectInterface = 'ANY'
                                  ),
                                  methods = list(
                                    initialize = function(...)
                                      callSuper(...)
                                  ))



nimbleList <- function(types,
                       name = NA,
                       where =  getNimbleFunctionEnvironment()) {
    ## This has a role like nimbleFunction but a much simpler implementation
    ## It returns a function that simply makes a regular R list and
    ## attaches two attributes, one to mark it as a nimbleList (for efficienct checking
    ## compatible with checking of other objects that have a class) and
    ## one that has the nimbleListDefClass object
    if(is.na(name)) name <- nf_refClassLabelMaker()
    nlDefClassObject <- nimbleListDefClass(types = types, className = name) 

    fields <- as.list(rep('ANY', length(types$vars)))
    names(fields) <- types$vars
    fields[[length(fields)+1]] <- "ANY"
    names(fields)[length(fields)] <- "nimbleListDef"



    nlGeneratorFunction <-   eval(  substitute(
      function(...){
      nlDefClassObject <- NLDEFCLASSOBJECT
      
      nlRefClassObject <- setRefClass(
          Class = NLREFCLASS_CLASSNAME,
          fields = NLREFCLASS_FIELDS,
          contains = 'nimbleListBase',
          methods = list(
            initialize = function(nlDefClassObject, ...){
              nimbleListDef <<- nlDefClassObject
              callSuper(...)
            }
          ),
          where = where
          )
      return(nlRefClassObject(nlDefClassObject, ...))},
      list(NLREFCLASS_CLASSNAME = name,
           NLREFCLASS_FIELDS = fields,
           NLDEFCLASSOBJECT = nlDefClassObject,
           where =where)))
    return(nlGeneratorFunction)
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
                                instances = 'ANY',
                                neededObjectNames =  'ANY'		#'character', ## a character vector of the names of objects such as models or modelValues that need to exist external to the nimbleFunction object so their contents can be pointed to 
                            ),
                            methods = list(
                                show = function() {
                                    writeLines(paste0('nlProcessing object ', nimbleListObj$className))
                                },
                                initialize = function(nimLists = NULL, className, project, ...) {
                                  neededTypes <<- list()
                                  callSuper(...)
                                  if(!is.null(nimLists)) {
                                    ## in new system, f must be a specialized nf, or a list of them
                                    nimbleProject <<- project
                                    sl <- if(is.list(nimLists)) nimLists[[1]]$nimbleListDef else nimLists$nimbleListDef
                                    nimbleListObj <<- sl
                                    if(missing(className)) {
                                      name <<- sl$className
                                    } else {
                                      name <<- className
                                    }
                                    instances <<- if(inherits(nimLists, 'list')) nimLists else list(nimLists)
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
                                    buildSymbolTable()
                                },
                                buildSymbolTable = function() {
                                  if(length(nimbleListObj$types$types) != length(nimbleListObj$types$vars))
                                    stop("Number of nimbleList vars provided is not equal to number of nimbleList types provided")
                                  if(is.null(nimbleListObj$types$sizes)){  ## if sizes haven't been provided, construct them from types
                                    nimbleListObj$types$sizes <<- list()
                                    varDims <- as.numeric(strsplit(gsub("[^0-9]", "", nimbleListObj$types$types), ""))
                                    for(i in seq_along(nimbleListObj$types$vars)){
                                      if(is.na(varDims[i])) nimbleListObj$types$sizes[[i]] <<-  NA
                                      else nimbleListObj$types$sizes[[i]] <<- rep(NA, varDims[i])
                                    }
                                  }
                                  symTab <<- symbolTable()
                                  for(i in seq_along(nimbleListObj$types$vars)){ ## ensure that types are of form "double", not "double(1)"
                                    nimbleListObj$types$types[i] <<- unlist(strsplit(nimbleListObj$types$types[i], "\\("))[1]
                                    ## first make work for nlDefs defined outside of function
                                    ## also consider what happens if nimbleProject doesn't exist!
                                    if(exists(nimbleListObj$types$types[i], envir = globalenv()) && 
                                       is.nlGenerator(eval(parse(text = nimbleListObj$types$types[i], keep.source = FALSE)))){
                                      nlList <- eval(parse(text = nimbleListObj$types$types[i], keep.source = FALSE))()
                                      nlp <- nimbleProject$compileNimbleList(nlList, initialTypeInferenceOnly = TRUE)
                                      className <- nlList$nimbleListDef$className
                                      newSym <- symbolNimbleList(name = nimbleListObj$types$vars[i], type = 'nimbleList', nlProc = nlp)
                                      symTab$addSymbol(newSym)
                                      if(!(className %in% names(neededTypes))) neededTypes[[className]] <<- newSym
                                    }
                                    else{
                                      symTab$addSymbol(symbolBasic(name = nimbleListObj$types$vars[i], 
                                                                   type = nimbleListObj$types$type[i], 
                                                                   nDim = length(nimbleListObj$types$size[[i]]), 
                                                                   size =nimbleListObj$types$ size[[i]]))
                                    }
                                  }
                                },
                                getSymbolTable = function() symTab
                            ))
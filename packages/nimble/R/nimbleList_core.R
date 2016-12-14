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
    basicTypes <- c("double", "integer", "character", "logical")
    nestedListGens <- list()
    elementTypes <- strsplit(types$types, '\\(')
    for(i in seq_along(types$types)){
      if(!(elementTypes[[i]][1] %in% basicTypes)){
        for(searchEnvironment in c(parent.frame(), globalenv())){
          if(is.nlGenerator(get(elementTypes[[i]][1], envir = searchEnvironment))){
            # nestedListDefs[[types$vars[i]]] <- get(elementTypes[[i]][1], envir = searchEnvironment)$new()$nimbleListDef
            nestedListGens[[types$vars[i]]] <- get(elementTypes[[i]][1], envir = searchEnvironment)
            break
          }
        }
      }
    }
    
    classFields <- as.list(rep('ANY', length(types$vars)))
    names(classFields) <- types$vars
    classFields[[length(classFields)+1]] <- "ANY"
    names(classFields)[length(classFields)] <- "nimbleListDef" ## initial nl definition stored here
    classFields[[length(classFields)+1]] <- "ANY"
    names(classFields)[length(classFields)] <- "nestedListGenList" ## nl generators for any nested lists stored here

    nlRefClassObject <- setRefClass(
      Class = name,
      fields = classFields,
      contains = 'nimbleListBase',
      methods = list(
        initialize = function(NLDEFCLASSOBJECT, NESTEDGENLIST, NLPRINTFIELDS, ...){
          nimbleListDef <<- NLDEFCLASSOBJECT
          nestedListGenList <<- NESTEDGENLIST
          for(i in seq_along(NESTEDGENLIST)){
            .self[[names(NESTEDGENLIST)[i]]] <- NESTEDGENLIST[[i]]$new()
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
      return(nlRefClassObject(NLDEFCLASSOBJECT, NESTEDGENLIST, ...))},
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
                                instances = 'ANY',
                                nestedListGens = 'ANY',
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
                                    
                                    nestedListGens <<- if(inherits(nimLists, 'list')) nimLists[[1]]$nestedListGenList else nimLists$nestedListGenList
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
                                    if(nimbleListObj$types$vars[i] %in% names(nestedListGens)){
                                      nlList <- nestedListGens[[nimbleListObj$types$vars[i]]]$new()
                                      nlp <- nimbleProject$compileNimbleList(nlList, initialTypeInferenceOnly = TRUE)
                                      className <- nlList$nimbleListDef$className
                                      newSym <- symbolNimbleList(name = nimbleListObj$types$vars[i], type = 'nimbleList', nlProc = nlp)
                                      # if(!(className %in% names(neededTypes))) 
                                      neededTypes[[className]] <<- newSym  ## if returnType is a NLG, this will ensure that it can be found in argType2symbol()
                                      symTab$addSymbol(newSym)
                                    }
                                    else{
                                     nimbleListObjType <- strsplit(nimbleListObj$types$types[i], "\\(")[[1]][1]
                                     nimbleListObjDim <-  as.numeric(strsplit(strsplit(nimbleListObj$types$types[i], "\\(")[[1]][2], "\\)")[[1]][1])
                                     symTab$addSymbol(argType2symbol(call(nimbleListObjType, nimbleListObjDim),
                                                                     neededTypes, nimbleListObj$types$vars[i]))
                                    }
                                  }
                                },
                                getSymbolTable = function() symTab
                            ))
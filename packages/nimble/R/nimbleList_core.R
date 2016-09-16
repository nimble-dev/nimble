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




#do we stil want to return an R list?
nimbleList <- function(types,
                       name = NA,
                       where =  getNimbleFunctionEnvironment()) {
    ## This has a role like nimbleFunction but a much simpler implementation
    ## It returns a function that simply makes a regular R list and
    ## attaches two attributes, one to mark it as a nimbleList (for efficienct checking
    ## compatible with checking of other objects that have a class) and
    ## one that has the nimbleListDefClass object
    if(is.na(name)) name <- nimble:::nf_refClassLabelMaker()
    nlDefClassObject <- nimbleListDefClass(types = types, className = name) 

    fields <- as.list(rep('ANY', length(types$vars)))
    names(fields) <- types$vars
    

    nlGeneratorFunction <-   eval(  substitute(
      nlRefClassObject <- setRefClass(
          Class = NLREFCLASS_CLASSNAME,
          fields = NLREFCLASS_FIELDS,
          where = where
      ),       
      list(NLREFCLASS_CLASSNAME = name,
           NLREFCLASS_FIELDS = fields,
           where =where)))
    return(nlGeneratorFunction)
}

## nimbleList processing class
## analogous to but must simpler than NFprocessing
nlProcessing <- setRefClass('nlProcessing',
                            fields = list(
                                cppDef = 'ANY',
                                nimbleListDef = 'ANY',
                                symTab = 'ANY',
                                neededTypes = 'ANY',
                                nimbleProject = 'ANY'
                            ),
                            methods = list(
                                show = function() {
                                    writeLines(paste0('nlProcessing object ', nimbleListDef$className))
                                },
                                initialize = function(...){
                                    callSuper(...)
                                    neededTypes <<- list()
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
                                  symTab <<- nimble:::buildSymbolTable(nimbleListDef$types$vars, nimbleListDef$types$types, nimbleListDef$types$sizes)
                                },
                                getSymbolTable = function() symTab
                            ))
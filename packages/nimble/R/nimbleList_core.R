nl_refClassLabelMaker <- labelFunctionCreator('nimListClass')

nimbleListDefClass <- setRefClass(
    ## This class holds a list of type information such as
    ## A = double(1), B = integer(2)
    ## The types need not be numeric.
    ## In general, ideally, they could be another nimbleList or a nimbleFunction
    Class = "nimbleListDefClass",
    fields = list(vars = 'ANY',
                  types = 'ANY',
                  className = 'ANY')
)

nimbleList <- function(types,
                       name = NA) {
    ## This has a role like nimbleFunction but a much simpler implementation
    ## It returns a function that simply makes a regular R list and
    ## attaches two attributes, one to mark it as a nimbleList (for efficienct checking
    ## compatible with checking of other objects that have a class) and
    ## one that has the nimbleListDefClass object
    if(is.na(name)) name <- nl_refClassLabelMaker()
    
    thisNimbleListDef <- nimbleListDefClass(types = types, className = name)
    ans <- function(...) {
        structure(list(...), class = "nimbleList", nimbleListDef = thisNimbleListDef)
    }
    ans
}

## nimbleList processing class
## analogous to but must simpler than NFprocessing
nlProcessing <- setRefClass('nlProcessing',
                            fields = list(
                                nimbleListDef = 'ANY',
                                symTab = 'ANY',
                                neededTypes = 'ANY',
                                nimbleProject = 'ANY'
                            ),
                            methods = list(
                                show = function() {
                                    writeLines(paste0('nlProcessing object ', name))
                                },
                                initialize = function(...){
                                    callSuper(...)
                                },
                                getSymbolTable = function() symTab, ## required name
                                buildSymbolTable = function() {},
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
                                })
                            )

nlProcessing$methods(buildSymbolTable = function() {
    ##message("nlProcessing::buildSymbolTable not written yet")
    symTab <<- nimble:::buildSymbolTable(nimbleListDef$vars, nimbleListDef$types, nimbleListDef$sizes)
})

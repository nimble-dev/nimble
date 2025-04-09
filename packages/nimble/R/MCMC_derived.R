

####################################################################
### virtual nimbleFunction template for all derived functions ######
####################################################################

#' @rdname derived
#' @export
derived_BASE <- nimbleFunctionVirtual(
    name = 'derived_BASE',
    methods = list(
        getResults = function() { returnType(double(2)) }
    )
)



####################################################################
### derived: logProb ###############################################
####################################################################

#' @rdname derived
#' @export
derived_logProb <- nimbleFunction(
    name = 'derived_logProb',
    contains = derived_BASE,
    setup = function(model, mvSaved, mvSamples, control) {
        ## control list extraction
        nodes <- extractControlElement(control, 'nodes', defaultValue = '.all')
        ## node list generation
        ## MORE CHECKS HERE WOULD BE GOOD
        if(nodes == '.all') {
            lpNodes <- model$getNodeNames(stochOnly = TRUE)
        } else {
            lpNodes <- nodes
        }
        ## create results array
        results <- array(0, c(0,1))    ## 0x1 array
    },
    run = function() {
        ### THIS RESIZING WILL NEED TO CHANGE
        nr <- dim(results)[2]
        resultsNew <- matrix(0, nrow = nr+1, ncol = 1)
        resultsNew[1:nr, 1] <- results[1:nr, 1]
        resultsNew[nr+1, 1] <- model$getLogProb(lpNodes)
        ##
        results <<- resultsNew
    },
    methods = list(
        getResults = function() {
            returnType(double(2))
            return(results)
        }
    )
)



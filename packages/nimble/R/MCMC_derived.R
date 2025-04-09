

####################################################################
### virtual nimbleFunction template for all derived functions ######
####################################################################

#' @rdname derived
#' @export
derived_BASE <- nimbleFunctionVirtual(
    name = 'derived_BASE',
    run = function(iter = double()) { },
    methods = list(
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), nchains = double()) { },
        after_chain  = function() { },
        getResults   = function() { returnType(double(2)) },
        reset        = function() { }
    ),
    methodControl = list(
        before_chain = list(required = FALSE),
        after_chain  = list(required = FALSE)
    )
)



####################################################################
### derived: test ##################################################
####################################################################

#' @rdname derived
#' @export
derived_test <- nimbleFunction(
    name = 'derived_test',
    contains = derived_BASE,
    setup = function(model, mvSaved, mvSamples, mvSamples2, interval, control) {
        ## control list extraction
        ## node list generation
        ## numeric value generation
        results <- matrix(0, nrow = 1, ncol = 1)
        count <- 1
        ## checks
    },
    run = function(iter = double()) {
        results[count, 1] <<- iter
        count <<- count + 1
    },
    methods = list(
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), nchains = double()) {
            nKeep <- floor(niter/interval)
            setSize(results, nKeep, 1)
        },
        after_chain = function() { },
        getResults = function() {
            returnType(double(2))
            return(results)
        },
        reset = function() {
            results <<- results * 0
            count <<- 1
        }
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
    setup = function(model, mvSaved, mvSamples, mvSamples2, interval, control) {
        ## control list extraction
        nodes <- extractControlElement(control, 'nodes', defaultValue = '.all')
        ## node list generation
        ## MORE CHECKS HERE WOULD BE GOOD  XXXXXXX
        if(nodes == '.all') {
            lpNodes <- model$getNodeNames(stochOnly = TRUE)
        } else {
            lpNodes <- nodes
        }
        ## numeric value generation
        results <- array(0, c(0,1))    ## 0x1 array
        ## checks
    },
    run = function(iter = double()) { },
    methods = list(
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), nchains = double()) {
            nKeep <- floor((niter-nburnin) / thin)
            setSize(results, nKeep, 1)
        },
        after_chain = function() { },
        getResults = function() {
            returnType(double(2))
            return(results)
        },
        reset = function() { }
    )
)



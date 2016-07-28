#' Makes the Metropolis-Hastings acceptance decision, based upon the input (log) Metropolis-Hastings ratio
#' 
#' This function returns a logical TRUE/FALSE value, indicating whether the proposed transition should be accepted (TRUE) or rejected (FALSE).
#' 
#' 
#' @param logMetropolisRatio The log of the Metropolis-Hastings ratio, which is calculated from model probabilities and forward/reverse transition probabilities. Calculated as the ratio of the model probability under the proposal to that under the current values multiplied by the ratio of the reverse transition probability to the forward transition probability.
#'
#' @details The Metropolis-Hastings accept/reject decisions is made as follows.  If \code{logMetropolisRatio} is greater than 0, accept (return \code{TRUE}).  Otherwise draw a uniform random number between 0 and 1 and accept if it is less that \code{exp(logMetropolisRatio}.  The proposed transition will be rejected (return \code{FALSE}). If \code{logMetropolisRatio} is NA, NaN, or -Inf, a reject (\code{FALSE}) decision will be returned.
#' 
#' @author Daniel Turek
#' @export
decide <- function(logMetropolisRatio) {
  if(is.na(logMetropolisRatio))	return(FALSE)
  if(logMetropolisRatio > 0) return(TRUE)
  if(runif(1,0,1) < exp(logMetropolisRatio)) return(TRUE)
  return(FALSE)
}

#NOTE: DETAILS(WAS BLANK) REMOVED



#' Creates a nimbleFunction for executing the Metropolis-Hastings jumping decision,
#' and updating values in the model, or in a carbon copy modelValues object, accordingly.
#' 
#' This nimbleFunction generator must be specialized to three required arguments: a model, a modelValues, and a character vector of node names.
#' 
#' @param model An uncompiled or compiled NIMBLE model object.  
#' @param mvSaved A modelValues object containing identical variables and logProb variables as the model. Can be created by \code{modelValues(model)}.
#' @param calcNodes A character vector representing a set of nodes in the model (and hence also the modelValues) object.  
#' @author Daniel Turek
#' @export
#' @details
#' Calling decideAndJump(model, mvSaved, calcNodes) will generate a specialized nimbleFunction with four required numeric arguments:
#' 
#' modelLP1: The model log-probability associated with the newly proposed value(s)
#' 
#' modelLP0: The model log-probability associated with the original value(s)
#' 
#' propLP1: The log-probability associated with the proposal forward-transition
#' 
#' propLP0: The log-probability associated with the proposal reverse-tranisiton
#' 
#' Executing this function has the following effects:
#' -- Calculate the (log) Metropolis-Hastings ratio, as logMHR = modelLP1 - modelLP0 - propLP1 + propLP0
#' -- Make the proposal acceptance decision based upon the (log) Metropolis-Hastings ratio
#' -- If the proposal is accepted, the values and associated logProbs of all calcNodes are copied from the model object into the mvSaved object
#' -- If the proposal is rejected, the values and associated logProbs of all calcNodes are copied from the mvSaved object into the model object
#' -- Return a logical value, indicating whether the proposal was accepted
decideAndJump <- nimbleFunction(
    setup = function(model, mvSaved, calcNodes) { },
    run = function(modelLP1 = double(), modelLP0 = double(), propLP1 = double(), propLP0 = double()) {
        logMHR <- modelLP1 - modelLP0 - propLP1 + propLP0
        jump <- decide(logMHR)
        if(jump) { nimCopy(from = model,   to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else   { nimCopy(from = mvSaved, to = model,   row = 1, nodes = calcNodes, logProb = TRUE) }
        returnType(logical())
        return(jump)
    }, where = getLoadingNamespace()
)





#' Creates a nimbleFunction for setting the value of a scalar model node,
#' calculating the associated deterministic dependents and logProb values,
#' and returning the total sum log-probability.
#' 
#' This nimbleFunction generator must be specialized to any model object and any scalar model node.
#' A specialized instance of this nimbleFunction will set the value of the target node in the specified model,
#' calculate the associated logProb, calculate the values of any deterministic dependents,
#' calculate the logProbs of any stochastic dependents,
#' and return the sum log-probability associated with the target node and all stochastic dependent nodes.
#' 
#' @param model An uncompiled or compiled NIMBLE model.  This argument is required.
#' @param targetNode The character name of any scalar node in the model object.  This argument is required.
#' @author Daniel Turek
#' @export
#' @details
#' Calling setAndCalculateOne(model, targetNode) will return a function with a single, required argument:
#' 
#' targetValue: The numeric value which will be put into the target node, in the specified model object.
#'
#' @examples
#' code <- nimbleCode({ for(i in 1:3) x[i] ~ dnorm(0, 1) })
#' Rmodel <- nimbleModel(code)
#' my_setAndCalc <- setAndCalculateOne(Rmodel, 'x[1]')
#' lp <- my_setAndCalc$run(2)
setAndCalculateOne <- nimbleFunction(
    setup = function(model, targetNode) {
        targetNodeAsScalar <- model$expandNodeNames(targetNode, returnScalarComponents = TRUE)
        if(length(targetNodeAsScalar) > 1)     stop('more than one targetNode; cannot use setAndCalculateOne()')
        calcNodes <- model$getDependencies(targetNode)
    },
    run = function(targetValue = double()) {
        model[[targetNode]] <<- targetValue
        lp <- calculate(model, calcNodes)
        returnType(double())
        return(lp)
    },  where = getLoadingNamespace()
)




#' Creates a nimbleFunction for setting the values of one or more model nodes,
#' calculating the associated deterministic dependents and logProb values,
#' and returning the total sum log-probability.
#' 
#' This nimbleFunction generator must be specialized to any model object and one or more model nodes.
#' A specialized instance of this nimbleFunction will set the values of the target nodes in the specified model,
#' calculate the associated logProbs, calculate the values of any deterministic dependents,
#' calculate the logProbs of any stochastic dependents,
#' and return the sum log-probability associated with the target nodes and all stochastic dependent nodes.
#' 
#' @param model An uncompiled or compiled NIMBLE model.  This argument is required.
#' @param targetNodes A character vector containing the names of one or more nodes or variables in the model.  This argument is required.
#' @author Daniel Turek
#' @aliases setAndCalculateDiff
#' @export
#' @details Calling \code{setAndCalculate(model, targetNodes)} or \code{setAndCalculate(model, targetNodes)} will return a nimbleFunction object whose \code{run} function takes a single, required argument:
#' 
#' targetValues: A vector of numeric values which will be put into the target nodes in the specified model object.  The length of this numeric vector much exactly match the number of target nodes.
#'
#' The difference between \code{setAndCalculate} and \code{setAndCalculateDiff} is the return value of their \code{run} functions.  In the former, \code{run} returns the sum of the log probabilities of the \code{targetNodes} with the provided \code{targetValues}, while the latter returns the difference between that sum with the new \code{targetValues} and the previous values in the \code{model}.
#' 
#' @examples
#' code <- nimbleCode({ for(i in 1:3) { x[i] ~ dnorm(0,1); y[i] ~ dnorm(0, 1)}})
#' Rmodel <- nimbleModel(code)
#' my_setAndCalc <- setAndCalculate(Rmodel, c('x[1]', 'x[2]', 'y[1]', 'y[2]'))
#' lp <- my_setAndCalc$run(c(1.2, 1.4, 7.6, 8.9))
setAndCalculate <- nimbleFunction(
    setup = function(model, targetNodes) {
        targetNodesAsScalar <- model$expandNodeNames(targetNodes, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(targetNodes)
    },
    run = function(targetValues = double(1)) {
        values(model, targetNodesAsScalar) <<- targetValues
        lp <- calculate(model, calcNodes)
        returnType(double())
        return(lp)
    }, where = getLoadingNamespace()
)

#' @rdname setAndCalculate
#' @export
setAndCalculateDiff <- nimbleFunction(
    setup = function(model, targetNodes) {
        targetNodesAsScalar <- model$expandNodeNames(targetNodes, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(targetNodes)
    },
    run = function(targetValues = double(1)) {
        values(model, targetNodesAsScalar) <<- targetValues
        lpD <- calculateDiff(model, calcNodes)
        returnType(double())
        return(lpD)
    }, where = getLoadingNamespace()
)


calcAdaptationFactor <- nimbleFunction(
    setup = function(paramDimension) {
        ## optimal acceptance rates:  (dim=1) .44,    (dim=2) .35,    (dim=3) .32,    (dim=4) .25,    (dim>=5) .234
        acceptanceRates <- c(0.44, 0.35, 0.32, 0.25, 0.234)
        if(paramDimension > 5)     paramDimension <- 5
        optimalAR <- acceptanceRates[paramDimension]
        timesAdapted <- 0
        gamma1       <- 0
    },
    run = function(acceptanceRate = double()) {
        timesAdapted <<- timesAdapted + 1
        gamma1 <<- 1 / ((timesAdapted + 3) ^ 0.8)   ## need this variable separate; it's extracted for use in 'RW_block' sampler
        gamma2 <- 10 * gamma1
        adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
        returnType(double())
        return(adaptFactor)
    },
    methods = list(
        reset = function() {
            timesAdapted <<- 0
            gamma1       <<- 0
        }
    ), where = getLoadingNamespace()
)


## A modified version of calcAdaptationFactor 
## The idea is to replace timesAdapted in the calculation of gamma1 with something more sensitive to the current state of convergence
## This gives an on-the-fly reset during the hill-climbing phase of convergence, 
## and therefore avoids the adaptative kernel becoming stuck when it needs to change direction.
## So timesAdapted is replaced, in that one line, with an effectiveTimesAdapted (ETA or effTimesAdapted)
## Once converged, downscaling of effTimesAdapted become increasingly negligible and it will grow at the same rate as timesAdapted
calcAdaptationFactor_ETA <- nimbleFunction( 
    setup = function(paramDimension, readaptability) {
        ## If user does not specify the targetAcRate then use
        ## the theoretically optimal (under ideal conditions) acceptance rates...
        ## (dim=1) .44, (dim=2) .35, (dim=3) .32, (dim=4) .25, (dim>=5) .234
        acceptanceRates <- c(0.44, 0.35, 0.32, 0.25, 0.234)
        if(paramDimension > 5) paramDimension <- 5
        optimalAR          <- acceptanceRates[paramDimension]
        Scale              <- 1        
        gamma1             <- 0
        timesAdapted       <- 0
        effTimesAdapted    <- 1 ## Effective Times Adapted
        maxEWMA_LogProbOld <- -.Machine$double.xmax
        ApproxNegInf       <- -.Machine$double.xmax
    },
    run = function(acceptanceRate = double(), maxEWMA_LogProb = double()) { ## reAdaptability = double()
        ## maxEWMA_LogProb is the maximum value taken so far in the Exponential Weighted Moving Average of LogProb
        ## It increases monotonically and the size of those increments goes to zero
        ## browser()
        nimPrint("Scale = ", Scale)
        Scale           <<- exp(readaptability*(maxEWMA_LogProbOld - maxEWMA_LogProb)) ## (1-Scale) -> 0 as timesAdapted -> Inf
        nimPrint("readaptability =", readaptability, 
                 " ETA = ", effTimesAdapted,
                 " exp(old-new)", exp(maxEWMA_LogProbOld - maxEWMA_LogProb),
                 " maxEWMA_LogProbOld = ", maxEWMA_LogProbOld,
                 " maxEWMA_LogProb = ", maxEWMA_LogProb,
                 " Scale = ", Scale)
        effTimesAdapted <<- 1 + effTimesAdapted * Scale
        nimPrint("ETA = ", effTimesAdapted)
        timesAdapted    <<- 1 + timesAdapted
        gamma1          <<- 1 / ((effTimesAdapted + 3)^0.8) ## Here timesAdapted was swapped for the more flexible effTimesAdapted 
        gamma2           <- 10 * gamma1                                ## Unchanged
        adaptFactor      <- exp(gamma2 * (acceptanceRate - optimalAR)) ## Unchanged
        ## nimPrint("timesAdapted = ", timesAdapted,
        ##          "readaptability =", readaptability, 
        ##          " ETA = ", effTimesAdapted,
        ##          " exp(old-new)", exp(maxEWMA_LogProbOld - maxEWMA_LogProb),
        ##          " maxEWMA_LogProb = ", maxEWMA_LogProb, " adaptFactor = ", adaptFactor)
        maxEWMA_LogProbOld <<- maxEWMA_LogProb
        returnType(double())
        return(adaptFactor)
    },
    methods = list(
        reset = function() {
            Scale              <<- 1   
            gamma1             <<- 0
            timesAdapted       <<- 0
            effTimesAdapted    <<- 1
            maxEWMA_LogProbOld <<- ApproxNegInf
        }
    ), where = getLoadingNamespace()
)







# for now export this as R<3.1.2 give warnings if don't

#' Class \code{codeBlockClass}
#' @aliases codeBlockClass
#' @export
#' @description
#' Classes used internally in NIMBLE and not expected to be called directly by users.
codeBlockClass <- setRefClass(
    Class   = 'codeBlockClass',
    fields  = list(codeBlock = 'ANY'),
    methods = list(
        
        initialize = function()     codeBlock <<- quote({}),
        
        addCode = function(expr, substituteList = NULL, quote = TRUE) {
            if(quote) expr <- substitute(expr)
            if(!is.null(substituteList)) {
                expr <- eval(substitute(substitute(EXPR, substituteList), list(EXPR=expr)))
            }
            if(is.call(expr) && expr[[1]] == '{') {  ## expr is a block of code
                for(i in seq_along(expr)[-1]) {
                    codeBlock[[length(codeBlock) + 1]] <<- expr[[i]]
                }
            } else {  ## expr is a single expression
                codeBlock[[length(codeBlock) + 1]] <<- expr
            }
        },
        
        getCode = function()     return(codeBlock)
    )
    
)


mcmc_listContentsToStr <- function(ls) {
    if(any(unlist(lapply(ls, is.function)))) warning('probably provided wrong type of function argument')
    ls <- lapply(ls, function(el) if(is.nf(el)) 'function' else el)   ## functions -> 'function'
    ls2 <- list()
    defaultOptions <- getNimbleOption('MCMCcontrolDefaultList')
    for(i in seq_along(ls)) {
        controlName <- names(ls)[i]
        controlValue <- ls[[i]]
        if(length(controlValue) == 0) next   ## remove length 0
        if(controlName %in% names(defaultOptions))   ## skip default control values
            if(identical(controlValue, defaultOptions[[controlName]])) next
        deparsedItem <- deparse(controlValue)
        if(length(deparsedItem) > 1) deparsedItem <- paste0(deparsedItem, collapse='')
        ls2[[i]] <- paste0(controlName, ': ', deparsedItem)
    }
    ls2 <- ls2[unlist(lapply(ls2, function(i) !is.null(i)))]
    str <- paste0(ls2, collapse = ',  ')
    ##if(length(ls2) == 1)
    ##    str <- paste0(str, ', default')
    str <- gsub('\"', '', str)
    str <- gsub('c\\((.*?)\\)', '\\1', str)
    return(str)
}


mcmc_findControlListNamesInCode <- function(code) {
    if(is.function(code))     return(mcmc_findControlListNamesInCode(body(code)))
    if(is.name(code) || is.numeric(code) || is.logical(code) || is.character(code) || is.pairlist(code))     return(character())
    
    if(is.call(code)) {
        if(code[[1]] == '$' && code[[2]] == 'control') {
            if(is.name(code[[3]])) { return(as.character(code[[3]]))
            } else                 { warning(paste0('having trouble processing control list elements: ', deparse(code)))
                                     return(character()) }
        }
        if(code[[1]] == '[[' && code[[2]] == 'control') {
            if(is.character(code[[3]])) { return(as.character(code[[3]]))
            } else                 { warning(paste0('having trouble processing control list elements: ', deparse(code)))
                                     return(character()) }
        }
        ## code is some call, other than $ or [[
        return(unique(unlist(lapply(code, mcmc_findControlListNamesInCode))))
    }
    browser()
    stop('not sure how to handle this code expression')
}





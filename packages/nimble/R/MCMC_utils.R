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
    if(is.na(logMetropolisRatio)) return(FALSE)
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
#' @param target A character vector providing the target node.
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
    name = 'decideAndJump',
    setup = function(model, mvSaved, target, calcNodes) {
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)   ## should be made faster
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
    },
    run = function(modelLP1 = double(), modelLP0 = double(), propLP1 = double(), propLP0 = double()) {
        logMHR <- modelLP1 - modelLP0 - propLP1 + propLP0
        jump <- decide(logMHR)
        if(jump) {
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
        returnType(logical())
        return(jump)
    }
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
    name = 'setAndCalculateOne',
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
    }
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
    name = 'setAndCalculate',
    setup = function(model, targetNodes) {
        targetNodesAsScalar <- model$expandNodeNames(targetNodes, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(targetNodes)
    },
    run = function(targetValues = double(1)) {
        values(model, targetNodesAsScalar) <<- targetValues
        lp <- calculate(model, calcNodes)
        returnType(double())
        return(lp)
    }
)

#' @rdname setAndCalculate
#' @export
setAndCalculateDiff <- nimbleFunction(
    name = 'setAndCalculateDiff',
    setup = function(model, targetNodes) {
        targetNodesAsScalar <- model$expandNodeNames(targetNodes, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(targetNodes)
    },
    run = function(targetValues = double(1)) {
        values(model, targetNodesAsScalar) <<- targetValues
        lpD <- calculateDiff(model, calcNodes)
        returnType(double())
        return(lpD)
    }
)


calcAdaptationFactor <- nimbleFunction(
    name = 'calcAdaptationFactor',
    setup = function(paramDimension, adaptFactorExponent) {
        ## optimal acceptance rates:  (dim=1) .44,    (dim=2) .35,    (dim=3) .32,    (dim=4) .25,    (dim>=5) .234
        acceptanceRates <- c(0.44, 0.35, 0.32, 0.25, 0.234)
        if(paramDimension > 5)     paramDimension <- 5
        optimalAR <- acceptanceRates[paramDimension]
        timesAdapted <- 0
        gamma1       <- 0
    },
    run = function(acceptanceRate = double()) {
        timesAdapted <<- timesAdapted + 1
        gamma1 <<- 1 / ((timesAdapted + 3) ^ adaptFactorExponent)   ## need this variable separate; it's extracted for use in 'RW_block' sampler
        gamma2 <- 10 * gamma1
        adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
        returnType(double())
        return(adaptFactor)
    },
    methods = list(
        getGamma1 = function() {
            returnType(double())
            return(gamma1)
        },
        reset = function() {
            timesAdapted <<- 0
            gamma1       <<- 0
        }
    )
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

mcmc_generateControlListArgument <- function(control, controlDefaults) {
    if(missing(control))           control         <- list()
    if(missing(controlDefaults))   controlDefaults <- list()
    thisControlList <- controlDefaults           ## start with all the defaults
    thisControlList[names(control)] <- control   ## add in any controls provided as an argument
    return(thisControlList)
}



mcmc_listContentsToStr <- function(ls, displayControlDefaults=FALSE, displayNonScalars=FALSE, displayConjugateDependencies=FALSE) {
    ##if(any(unlist(lapply(ls, is.function)))) warning('probably provided wrong type of function argument')
    if(!displayConjugateDependencies) {
        if(grepl('^conjugate_d', names(ls)[1])) ls <- ls[1]    ## for conjugate samplers, remove all 'dep_dnorm', etc, control elements (don't print them!)
        if(grepl('^CRP_cluster_wrapper', names(ls)[1]) && 'wrapped_type' %in% names(ls) &&
           grepl('conjugate_d', ls$wrapped_type)[1]) ls <- ls[!names(ls) %in% c('wrapped_conf')]
    }
    if((names(ls)[1] == 'CRP sampler') && ('clusterVarInfo' %in% names(ls)))   ls[which(names(ls) == 'clusterVarInfo')] <- NULL   ## remove 'clusterVarInfo' from CRP sampler printing
    ls <- lapply(ls, function(el) if(is.nf(el) || is.function(el)) 'function' else el)   ## functions -> 'function'
    ls2 <- list()
    ## to make displayControlDefaults argument work again, would need to code process
    ## setup function of each sampler, find & extract the default control values, and pass
    ## them into this function.  Then set:
    ## defaultOptions <- [the default list for this sampler]
    ## (the commented function below, mcmc_findControlListNamesInCode, is a helpful start).
    ## old roxygen documentation for displayControlDefaults (from conf$print() method):
    ## displayControlDefaults: A logical argument, specifying whether to display default values of control list elements (default FALSE).
    ## -DT July 2017
    for(i in seq_along(ls)) {
        controlName <- names(ls)[i]
        controlValue <- ls[[i]]
        if(length(controlValue) == 0) next   ## remove length 0
        ##if(!displayControlDefaults)
        ##    if(controlName %in% names(defaultOptions))   ## skip default control values
        ##        if(identical(controlValue, defaultOptions[[controlName]])) next
        if(!displayNonScalars)
            if(is.numeric(controlValue) || is.logical(controlValue))
                if(length(controlValue) > 1)
                    controlValue <- ifelse(is.null(dim(controlValue)), 'custom vector', 'custom array')
        if(inherits(controlValue, 'samplerConf'))   controlValue <- controlValue$name
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


#' Extract named elements from MCMC sampler control list
#'
#' @param controlList control list object, which is passed as an argument to all MCMC sampler setup functions.
#' @param elementName character string, giving the name of the element to be extracted from the control list.
#' @param defaultValue default value of the control list element, giving the value to be used when the elementName does not exactly match the name of an element in the controlList.
#' @param error character string, giving the error message to be printed if no defaultValue is provided and elementName does not match the name of an element in the controlList.
#' @return The element of controlList with name matching elementName; or, if no controlList name matches elementName, rather the defaultValue is returned.
#' @author Daniel Turek
#' @export
extractControlElement <- function(controlList, elementName, defaultValue, error) {
    if(missing(controlList) | !is.list(controlList))      stop('extractControlElement: controlList argument must be a list')
    if(missing(elementName) | !is.character(elementName)) stop('extractControlElement: elementName argument must be a character string variable name')
    if(missing(defaultValue) & missing(error))            stop('extractControlElement: must provide either defaultValue or error argument')
    if(elementName %in% names(controlList)) {
        return(controlList[[elementName]])
    } else {
        if(!missing(defaultValue)) return(defaultValue) else stop(error)
    }
}


## obselete, since control defaults were moved to sampler function setup code
## -DT July 2017
##mcmc_findControlListNamesInCode <- function(code) {
##    if(is.function(code))     return(mcmc_findControlListNamesInCode(body(code)))
##    if(is.name(code) || is.numeric(code) || is.logical(code) || is.character(code) || is.pairlist(code))     return(character())
##    if(is.call(code)) {
##        if(code[[1]] == '$' && code[[2]] == 'control') {
##            if(is.name(code[[3]])) { return(as.character(code[[3]]))
##            } else                 { warning(paste0('having trouble processing control list elements: ', deparse(code)))
##                                     return(character()) }
##        }
##        if(code[[1]] == '[[' && code[[2]] == 'control') {
##            if(is.character(code[[3]])) { return(as.character(code[[3]]))
##            } else                 { warning(paste0('having trouble processing control list elements: ', deparse(code)))
##                                     return(character()) }
##        }
##        ## code is some call, other than $ or [[
##        return(unique(unlist(lapply(code, mcmc_findControlListNamesInCode))))
##    }
##    browser()
##    stop('not sure how to handle this code expression')
##}



#' @export
samplesSummary <- function(samples, round) {
    summary <- try(cbind(
        `Mean`      = apply(samples, 2, mean),
        `Median`    = apply(samples, 2, median),
        `St.Dev.`   = apply(samples, 2, sd),
        `95%CI_low` = apply(samples, 2, function(x) quantile(x, 0.025)),
        `95%CI_upp` = apply(samples, 2, function(x) quantile(x, 0.975))),
                   silent = TRUE)
    if(inherits(summary, 'try-error')) {
        warning('Could not calculate the full summary of posterior samples, possibly due to NA or NaN values present in the samples array', call. = FALSE)
        summary <- array(as.numeric(NA), dim = c(ncol(samples), 5))
        rownames(summary) <- colnames(samples)
        colnames(summary) <- c('Mean','Median','St.Dev.','95%CI_low','95%CI_upp')
        if(ncol(samples) > 0) for(i in 1:ncol(samples)) {
            theseSamples <- samples[,i]
            if(isValid(theseSamples)) {
                summary[i, 1] <- mean(theseSamples)
                summary[i, 2] <- median(theseSamples)
                summary[i, 3] <- sd(theseSamples)
                summary[i, 4] <- quantile(theseSamples, 0.025)
                summary[i, 5] <- quantile(theseSamples, 0.975)
            }
        }
    }
    if(!missing(round)) summary <- round(summary, digits = round)
    return(summary)
}




## weed out missing indices from the monitors
mcmc_processMonitorNames <- function(model, nodes) {
    isLogProbName <- grepl('logProb', nodes)
    expandedNodeNames <- model$expandNodeNames(nodes[!isLogProbName])
    origLogProbNames <- nodes[isLogProbName]
    expandedLogProbNames <- character()
    if(length(origLogProbNames) > 0) {
        nodeName_fromLogProbName <- gsub('logProb_', '', origLogProbNames)
        expandedLogProbNames <- model$modelDef$nodeName2LogProbName(nodeName_fromLogProbName)
    }
    return(c(expandedNodeNames, expandedLogProbNames))
}

mcmc_checkWAICmonitors_conditional <- function(model, monitors, dataNodes) {
    parentNodes <- getParentNodes(dataNodes, model, stochOnly = TRUE)
    parentVars <- model$getVarNames(nodes = parentNodes)
    wh <- which(!parentVars %in% monitors)
    if(length(wh)) {
        if(length(wh) > 10)
            badVars <- c(parentVars[wh[1:10]], "...") else badVars <- parentVars[wh]
        stop(paste0("In calculate WAIC in NIMBLE, all parameters of",
                    " data nodes in the model must be monitored.", "\n", 
                    "  Currently, the following parameters are not monitored: ",
                    paste0(badVars, collapse = ", ")))
    }
    message('Monitored nodes are valid for WAIC.')
}

## Used through version 0.10.0 and likely to be used in some form once we re-introduce mWAIC
mcmc_checkWAICmonitors <- function(model, monitors, dataNodes) {
    monitoredDetermNodes <- model$expandNodeNames(monitors)[model$isDeterm(model$expandNodeNames(monitors))]
    if(length(monitoredDetermNodes) > 0) {
        monitors <- monitors[- which(monitors %in% model$getVarNames(nodes = monitoredDetermNodes))]
    }
    thisNodes <- model$getNodeNames(stochOnly = TRUE, topOnly = TRUE)
    thisVars <- model$getVarNames(nodes = thisNodes)
    thisVars <- thisVars[!(thisVars %in% monitors)]
    while(length(thisVars) > 0) {
        nextNodes <- model$getDependencies(thisVars, stochOnly = TRUE, omit = monitoredDetermNodes, self = FALSE, includeData = TRUE)
        if(any(nextNodes %in% dataNodes)) {
            badDataNodes <- dataNodes[dataNodes %in% nextNodes]
            if(length(badDataNodes) > 10) {
                badDataNodes <- c(badDataNodes[1:10], "...")
            }
            stop(paste0("In order for a valid WAIC calculation, all parameters of",
                        " data nodes in the model must be monitored, or be", 
                        " downstream from monitored nodes.", 
                        " See help(buildMCMC) for more information on valid sets of",
                        " monitored nodes for WAIC calculations.", "\n",
                        " Currently, the following data nodes have un-monitored",
                        " upstream parameters:", "\n ", 
                        paste0(badDataNodes, collapse = ", ")))
        }
        thisVars <- model$getVarNames(nodes = nextNodes)
        thisVars <- thisVars[!(thisVars %in% monitors)]
    }
    message('Monitored nodes are valid for WAIC')
}


## formerly used in conf$printSamplers(byType = TRUE)
## to compress scalar index ranges of variables.
## deprecated Nov 2019.
##mcmc_compressIndexRanges <- function(nums) {
##    nums <- sort(unique(nums))
##    rangeStartInd <- c(1, which(diff(nums) != 1) + 1)
##    ranges <- vector('list', length(rangeStartInd))
##    for(i in seq_along(rangeStartInd)) {
##        startInd <- rangeStartInd[i]
##        if(i == length(rangeStartInd)) {
##            if(rangeStartInd[i] == length(nums)) {
##                ranges[[i]] <- nums[startInd]
##            } else {
##                ranges[[i]] <- substitute(START:END, list(START = as.numeric(nums[startInd]),
##                                                          END = as.numeric(nums[length(nums)])))
##            }
##        } else {
##            if(startInd+1 < rangeStartInd[i+1]) {
##                ranges[[i]] <- substitute(START:END, list(START = as.numeric(nums[startInd]),
##                                                          END = as.numeric(nums[rangeStartInd[i+1]-1])))
##            } else ranges[[i]] <- nums[startInd]
##        }
##    }
##    return(ranges)
##}
##
##
##mcmc_getIndexNumberFromNodeNames <- function(nodeNames, indNumber) {
##    (nodeInsideBrackets <- gsub('^[[:alpha:]]+\\[(.+)\\]$', '\\1', nodeNames))
##    splitList <- strsplit(nodeInsideBrackets, ',')
##    if(length(splitList[[1]]) < indNumber) stop('in mcmc_getIndexNumberFromNodeNames: indNumber exceeds node dimension', call. = FALSE)
##    indices <- as.numeric(sapply(splitList, function(el) el[indNumber]))
##    return(indices)
##}









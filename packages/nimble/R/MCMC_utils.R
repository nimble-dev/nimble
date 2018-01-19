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
    name = 'decideAndJump',
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
    }, where = getLoadingNamespace()
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
    }, where = getLoadingNamespace()
)


calcAdaptationFactor <- nimbleFunction(
    name = 'calcAdaptationFactor',
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
    if(!displayConjugateDependencies)
        if(grepl('^conjugate_d', names(ls)[1])) ls <- ls[1]    ## for conjugate samplers, remove all 'dep_dnorm', etc, control elements (don't print them!)
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
samplesSummary <- function(samples) {
    cbind(
        `Mean`      = apply(samples, 2, mean),
        `Median`    = apply(samples, 2, median),
        `St.Dev.`   = apply(samples, 2, sd),
        `95%CI_low` = apply(samples, 2, function(x) quantile(x, 0.025)),
        `95%CI_upp` = apply(samples, 2, function(x) quantile(x, 0.975)))
}



#' @export
samplesPlot <- function(samples, var=colnames(samples), ind=NULL, burnin=NULL, width=7, height=4, legend=TRUE, legend.location='topright', traceplot=TRUE, densityplot=TRUE, file=NULL) {
    if(!is.null(file)) pdf(file, width=width, height=height) else
    if(inherits(try(eval(parse(text=paste0('knitr','::','opts_chunk$get(\'dev\')'))[[1]]), silent=TRUE), 'try-error') || is.null(eval(parse(text=paste0('knitr','::','opts_chunk$get(\'dev\')'))[[1]])))   ## if called from Rmarkdown/knitr
        dev.new(height=height, width=width)
    par.save <- par(no.readonly = TRUE)
    par(mfrow=c(1,traceplot+densityplot), cex=0.7, cex.main=1.5, cex.axis=0.9, lab=c(3,3,7), mgp=c(0,0.4,0), mar=c(1.6,1.6,2,0.6), oma=c(0,0,0,0), tcl=-0.3, bty='l')
    ## process samples
    var <- gsub('\\[', '\\\\\\[', gsub('\\]', '\\\\\\]', var))   ## add \\ before any '[' or ']' appearing in var
    var <- unlist(lapply(var, function(n) grep(paste0('^', n,'(\\[.+\\])?$'), colnames(samples), value=TRUE)))  ## expanded any indexing
    samples <- samples[, var, drop=FALSE]
    if(!is.null(ind) && !is.null(burnin)) stop('only specify either ind or burnin')
    if(!is.null(ind))     samples <- samples[ind, , drop=FALSE]
    if(!is.null(burnin))  samples <- samples[(burnin+1):dim(samples)[1], , drop=FALSE]
    nparam <- ncol(samples)
    rng <- range(samples)
    if(!traceplot & !densityplot) stop('both traceplot and densityplot are false')
    if(traceplot) {  ## traceplot
        plot(1:nrow(samples), ylim=rng, type='n', main='Traceplots', xlab='', ylab='')
        for(i in 1:nparam)
            lines(samples[,i], col=rainbow(nparam, alpha=0.75)[i])
        if(legend & !densityplot & !is.null(dimnames(samples)) & is.character(dimnames(samples)[[2]]))
            legend(legend=dimnames(samples)[[2]], fill=rainbow(nparam, alpha=0.5), bty='n', x=legend.location)
    }  ## finish traceplot
    if(densityplot) {  ## denstyplot
        xMin <- xMax <- yMax <- NULL
        for(i in 1:nparam) {
            d <- density(samples[,i])
            xMin <- min(xMin,d$x); xMax <- max(xMax,d$x); yMax <- max(yMax, d$y) }
        plot(1, xlim=c(xMin,xMax), ylim=c(0,yMax), type='n', main='Posterior Densities', xlab='', ylab='', yaxt='n')
        for(i in 1:nparam)
            polygon(density(samples[,i]), col=rainbow(nparam, alpha=0.2)[i], border=rainbow(nparam, alpha=0.2)[i])
        if(legend & !is.null(dimnames(samples)) & is.character(dimnames(samples)[[2]]))
            legend(legend=dimnames(samples)[[2]], fill=rainbow(nparam, alpha=0.5), bty='n', x=legend.location)
    }  ## finish densityplot
    if(!is.null(file)) dev.off()
    invisible(par(par.save))
}



#' @export
chainsPlot <- function(samplesList, var=NULL, nrows=3, width=7, height=min(1+3*nrows,7), legend=!is.null(names(samplesList)), legend.location='topright', jitter=1, buffer.right=0, buffer.left=0, cex=1, file=NULL) {
    if(!is.null(file)) pdf(file, width=width, height=height) else
    if(inherits(try(eval(parse(text=paste0('knitr','::','opts_chunk$get(\'dev\')'))[[1]]), silent=TRUE), 'try-error') || is.null(eval(parse(text=paste0('knitr','::','opts_chunk$get(\'dev\')'))[[1]])))   ## if called from Rmarkdown/knitr
        dev.new(height=height, width=width)
    par.save <- par(no.readonly = TRUE)
    par(mfrow=c(nrows,1), oma=c(3,1,1,1), mar=c(4,1,0,1), mgp=c(3,0.5,0))
    if(class(samplesList) != 'list') samplesList <- list(samplesList)
    if(!is.null(var)) samplesList <- lapply(samplesList, function(samples) {
        var <- gsub('\\[', '\\\\\\[', gsub('\\]', '\\\\\\]', var))   ## add \\ before any '[' or ']' appearing in var
        theseVar <- unlist(lapply(var, function(n) grep(paste0('^', n,'(\\[.+\\])?$'), colnames(samples), value=TRUE)))  ## expanded any indexing
        samples[, theseVar, drop=FALSE]
    })
    chainParamNamesList <- lapply(samplesList, function(s) colnames(s))
    nChains <- length(samplesList)
    paramNamesAll <- unique(unlist(lapply(samplesList, function(s) colnames(s))))
    nParamsAll <- length(paramNamesAll)
    cols <- rainbow(nChains)
    ## construct 3D summary array:
    summary <- array(as.numeric(NA), dim = c(nChains, 3, nParamsAll))
    if(!is.null(names(samplesList))) dimnames(summary)[[1]] <- names(samplesList)
    dimnames(summary)[[2]] <- c('mean','low','upp')
    dimnames(summary)[[3]] <- paramNamesAll
    for(iChain in 1:nChains) {
        theseSamples <- samplesList[[iChain]]
        thisSummary <- rbind(mean = apply(theseSamples, 2, mean),
                             low  = apply(theseSamples, 2, function(x) quantile(x, 0.025)),
                             upp  = apply(theseSamples, 2, function(x) quantile(x, 0.975)))
        summary[iChain,c('mean','low','upp'),colnames(thisSummary)] <- thisSummary
    }
    nParamsPerRow <- ceiling(nParamsAll/nrows)
    sq <- if(nChains==1) 0 else seq(-1,1,length=nChains)
    scale <- width/nParamsPerRow * jitter * 0.1  ## adjust jitter scale factor at end
    for(iRow in 1:nrows) {
        rowParamInd <- (1+(iRow-1)*nParamsPerRow) : ifelse(iRow==nrows,nParamsAll,iRow*nParamsPerRow)
        nRowParams <- length(rowParamInd)
        rowParamNames <- paramNamesAll[rowParamInd]
        xs <- 1:nRowParams
        names(xs) <- rowParamNames
        ylim <- range(summary[,c('low','upp'),rowParamNames], na.rm=TRUE)
        plot(x=-100, y=0, xlim=c(1-buffer.left,nParamsPerRow+buffer.right), ylim=ylim, xaxt='n', ylab='', xlab='', tcl=-0.3, cex.axis=cex)
        axis(1, at=1:nRowParams, labels=FALSE, tcl=-0.3)
        text(x=1:nRowParams, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), labels=rowParamNames, srt=45, adj=1, xpd=TRUE, cex=0.9*cex)
        for(iChain in 1:nChains) {
            ps <- intersect(rowParamNames, chainParamNamesList[[iChain]])
            xsJittered <- xs + sq[iChain]*scale
            points(x=xsJittered[ps], y=summary[iChain,'mean',ps], pch=16, col=cols[iChain])
            segments(x0=xsJittered[ps], y0=summary[iChain,'low',ps], y1=summary[iChain,'upp',ps], lwd=1, col=cols[iChain])
        }
        if(legend) legend(legend.location, legend=names(samplesList), pch=16, col=cols, cex=cex)
    }
    if(!is.null(file)) dev.off()
    invisible(par(par.save))
}







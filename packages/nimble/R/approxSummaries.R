marginalSplineR <- nimbleRcall(function(theta = double(1), logdens = double(1)){},
  Rfun = 'marginalSpline',
  returnType = double(2)
)

marginalSpline = function(theta, logdens){
    n <- length(theta)
    rn <- range(theta)
    rnl <- diff(rn)
    thetarange <- c(min(rn) - rnl/2,max(rn) + rnl/2)
    finegrid <- c(seq(thetarange[1],thetarange[2],length.out=1000)) ## This is based off of Stringer for a fine grid.
    if (n <= 3) {
        log_pdf <- as.function(polynom::poly.calc(x = theta, y = logdens))
        logPDF <- log_pdf(finegrid)             
    }else{
        ss <- splines::interpSpline(theta, logdens, bSpline = TRUE, sparse = FALSE)
        if(isS4(co <- ss[["coefficients"]])) ss[["coefficients"]] <- as.vector(co)  ## Not sure why this might be necessary.
        logPDF <- as.numeric(stats::predict(ss, finegrid)$y)
    }
    ## Normalize the PDF:
    pdf <- exp(logPDF)
    trapezoids <- diff(finegrid)*(pdf[-n]+pdf[-1])/2  # trapezoidal rule (could use Simpson as (2M+T)/3
    norm <- sum(trapezoids)
    pdf <- pdf/norm
    cdf <- c(0, cumsum(trapezoids)/norm)
    ## Return the gridded distribution information.
    return(cbind(finegrid, pdf, cdf))
}

## Trapezoidal rule for estimating expectation.
estimateExpectation <- function(grid, pdf, functional, n = length(grid)) {
    return(sum(diff(grid)*(pdf[-n]*functional[-n]+pdf[-1]*functional[-1])/2))
}

## Works only for transformations that are 1:1 between theta and p.
## For non-1:1 report on theta scale probably and probably provide a simulation method like inla.hyperpar.sample.
## See also what Stringer does. 
## Need to wire in the reTrans transformation and logDetJacobian.

summarizeMarginal <- function(marginalApprox, transform = inverseTransform, logDetJacobian = logDetJacobian,
                              quantiles = c(.025,.25,.5,.75,.975), functionals = NULL, ...) {
    n <- nrow(marginalApprox)

    finegridTrans <- marginalApprox[,'finegrid']
    pdfTrans <- marginalApprox[,'pdf']
    cdf <- marginalApprox[,'cdf']  # same for both scales
    
    finegrid <- transform(finegridTrans)
    pdf <- pdfTrans[,'pdf']*logDetJac(finegrid)
    # sum(diff(finegrid)*(pdf[-n]+pdf[-1])/2) # trapezoidal rule showing normalization is preserved

    ## For now use smoothing spline on quantile function on transformed (theta) scale.
    ## How do INLA and Stringer/aghq get their quantiles? I think Paul and Chris discussed this.
    used <- cdf > .001 & cdf < .999
    ss <- splines::interpSpline(cdf[used], finegridTrans[used],
                                bSpline = TRUE, sparse = FALSE)
    quantsTrans <- stats::predict(ss, quantiles)$y
    quants <- transform(quantsTrans)

    ## Posterior expectations: defaults
    ## Use pdf on transformed (theta) scale.
    postMean <- estimateExpectation(finegridTrans, pdfTrans, finegrid, n)
    functional <- (finegrid - postMean)^2
    postVar <- estimateExpectation(finegridTrans, pdfTrans, functional, n)
    postSD <- sqrt(postVar)

    ## transformed scale
    ## postMean <- sum(diff(finegrid)*(pdf[-n]*finegrid[-n]+pdf[-1]*finegrid[-1])/2)
    ## functional <- (finegrid - postMean)^2
    ## postVar <- sum(diff(finegrid)*(pdf[-n]*functional[-n]+pdf[-1]*functional[-1])/2)

    ## Posterior expectations: user-defined
    ## Not sure about use of `...`.
    if(!is.null(functionals)) {
        userExpectations <- sapply(functionals, function(fun) {
            functional <- fun(finegrid, ...)
            return(estimateExpectation(finegridTrans, pdfTrans, functional, n))
        })
    } else userExpectations <- NULL
    ## Also decide what to return in terms of pdf information. (See INLA/Stringer.)
    return(list(quantiles = quants, postMean = postMean, postVar = postVar, postSD = postSD, userExpectations = userExpectations))
}

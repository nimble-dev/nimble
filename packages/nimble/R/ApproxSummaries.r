
## Following aghq interpolation methods for summaries:
marginalSpline = function(theta, logdens){
    n <- length(theta)
    rn <- range(theta)
    rnl <- diff(rn)
    thetarange <- c(min(rn) - rnl/2,max(rn) + rnl/2)
    finegrid <- c(seq(thetarange[1],thetarange[2],length.out=1000))
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
  d <- (finegrid[2]-finegrid[1])
  pdf <- pdf/(sum(pdf)*d)
  return(cbind(finegrid, pdf))
}

marginalSplineR <- nimbleRcall(function(theta = double(1), logdens = double(1)){},
  Rfun = 'marginalSpline',
  returnType = double(2)
)

library(rjags, quietly=TRUE)

read.jagsdata <- function(file)
{
  e <- new.env()
  eval(parse(file), e)
  return(as.list(e))
}

check.fun <- function() {
   checkstats <- summary(x)$statistics
   if (is.matrix(benchstats)) {
     rnb <- order(rownames(benchstats))
     rnc <- order(rownames(checkstats))
     z <- (benchstats[rnb,1] - checkstats[rnc,1])
     z <- ifelse(z==0, 0, z/benchstats[rnb,2])
   } else {
     z <- (benchstats[1] - checkstats[1])
     z <- ifelse(z==0, 0, z/benchstats[2])
   }
   print(z)
   if (any(abs(z) > 0.15)) {
     stop("FAIL")
     quit(save="no", status=1)
   } else {
     cat("OK\n")
   }
}

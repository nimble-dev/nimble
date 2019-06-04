### INSTRUCTIONS:
## enter each test as a list, with an informative name, NIMBLE expression to evaluate, vector of input dimensions, value of output dimension, and (if NIMBLE expression cannot be directly evaluated in R) the equivalent pure R expression whose result should match the NIMBLE result

testsVaried = list(
  list(name = "matrix direct product", expr = quote(out <- arg1 * arg2), inputDim = c(2,2), outputDim = 2),
  list(name = "matrix direct product with scalar addition", expr = quote(out <- (arg1+1) * (arg2+1)), inputDim = c(2,2), outputDim = 2),
  list(name = "matrix absolute value", expr = quote(out <- abs(arg1)), inputDim = c(2), outputDim = 2),
  list(name = "matrix absolute value with scalar addition", expr = quote(out <- abs(arg1 - 2)), inputDim = c(2), outputDim = 2),
  list(name = "vector pmin", expr = quote(out <- pmin(arg1, arg2)), inputDim = c(1,1), outputDim = 1),
  list(name = "vector pmax", expr = quote(out <- pmax(arg1, arg2)), inputDim = c(1,1), outputDim = 1),
  list(name = "sd with addition", expr = quote(out <- sd(arg1) + 3), inputDim = c(1), outputDim = 0),
  list(name = "sd of vector with addition", expr = quote(out <- sd(arg1 + 3)), inputDim = c(1), outputDim = 0),
  list(name = "sd of matrix-vector multiply", expr = quote(out <- sd(arg1 %*% arg2)), inputDim = c(2,1), outputDim = 0),
  list(name = "var of vector", expr = quote(out <- var(arg1)), inputDim = c(1), outputDim = 0),
  list(name = "log determinant", expr = quote(out <- logdet(arg1)), inputDim = c(2), outputDim = 0)
  )

testsBasicMath = list(
  list(name = 'exp of scalar', expr = quote(out <- exp(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'log of scalar', expr = quote(out <- log(abs(arg1))), inputDim = 0, outputDim = 0),
  list(name = 'sqrt of scalar', expr = quote(out <- sqrt(abs(arg1))), inputDim = 0, outputDim = 0),
  list(name = 'abs of scalar', expr = quote(out <- abs(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'step of scalar', expr = quote(out <- step(arg1)), inputDim = 0, outputDim = 0, Rcode = quote( out <- as.numeric(arg1 > 0))),
  list(name = 'cube of scalar', expr = quote(out <- cube(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'cos of scalar', expr = quote(out <- cos(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'acos of cos of scalar', expr = quote(out <- acos(cos(arg1))), inputDim = 0, outputDim = 0),
  list(name = 'sin of scalar', expr = quote(out <- sin(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'asin of sin of scalar', expr = quote(out <- asin(sin(arg1))), inputDim = 0, outputDim = 0),
  list(name = 'tan of scalar', expr = quote(out <- tan(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'atan of tan of scalar', expr = quote(out <- atan(tan(arg1))), inputDim = 0, outputDim = 0),
  list(name = 'cosh of scalar', expr = quote(out <- cosh(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'sinh of scalar', expr = quote(out <- sinh(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'tanh of scalar', expr = quote(out <- tanh(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'acosh of scalar', expr = quote(out <- acosh(1 + abs(arg1))), inputDim = 0, outputDim = 0),
  list(name = 'asinh of scalar', expr = quote(out <- asinh(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'atanh of scalar', expr = quote(out <- atanh(abs(arg1)%%1)), inputDim = 0, outputDim = 0),
  ###
  list(name = 'exp of vector', expr = quote(out <- exp(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'log of vector', expr = quote(out <- log(abs(arg1))), inputDim = 1, outputDim = 1),
  list(name = 'sqrt of vector', expr = quote(out <- sqrt(abs(arg1))), inputDim = 1, outputDim = 1),
  list(name = 'abs of vector', expr = quote(out <- abs(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'step of vector', expr = quote(out <- step(arg1)), inputDim = 1, outputDim = 1, Rcode = quote(out <- as.numeric(arg1 > 0)), returnType = "integer"), # failing as of 10/2/17 with C result for input in (-1,0) being 1; formerly: # knownFailure = '.*compiles'), ## FAILS on compileNimble(nfR) with Eigen error.
  list(name = 'cube of vector', expr = quote(out <- cube(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'cos of vector', expr = quote(out <- cos(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'acos of cos of vector', expr = quote(out <- acos(cos(arg1))), inputDim = 1, outputDim = 1),
  list(name = 'sin of vector', expr = quote(out <- sin(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'asin of sin of vector', expr = quote(out <- asin(sin(arg1))), inputDim = 1, outputDim = 1),
  list(name = 'tan of vector', expr = quote(out <- tan(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'atan of tan of vector', expr = quote(out <- atan(tan(arg1))), inputDim = 1, outputDim = 1),
  list(name = 'cosh of vector', expr = quote(out <- cosh(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'sinh of vector', expr = quote(out <- sinh(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'tanh of vector', expr = quote(out <- tanh(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'acosh of vector', expr = quote(out <- acosh(1 + abs(arg1))), inputDim = 1, outputDim = 1),
  list(name = 'asinh of vector', expr = quote(out <- asinh(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'atanh of vector', expr = quote(out <- atanh(abs(arg1)-floor(abs(arg1)))), inputDim = 1, outputDim = 1), ## formerly failed as using modulo
  ###
  list(name = 'scalar + scalar', expr = quote(out <- arg1 + arg2), inputDim = c(0,0), outputDim = 0),
  list(name = 'diff of scalars', expr = quote(out <- arg1 - arg2), inputDim = c(0,0), outputDim = 0),
  list(name = 'product of scalars', expr = quote(out <- arg1 * arg2), inputDim = c(0,0), outputDim = 0),
  list(name = 'ratio of scalars', expr = quote(out <- arg1 / arg2), inputDim = c(0,0), outputDim = 0),
  list(name = 'power of scalars via ^', expr = quote(out <- arg1 ^ arg2), inputDim = c(0,0), outputDim = 0),
  list(name = 'power of scalars via pow', expr = quote(out <- pow(arg1, arg2)), inputDim = c(0,0), outputDim = 0),
  list(name = 'power of scalars via ^ with positive first arg', expr = quote(out <- exp(arg1) ^ arg2), inputDim = c(0,0), outputDim = 0),
  list(name = 'power of scalars via pow with positive first arg', expr = quote(out <- pow(exp(arg1), arg2)), inputDim = c(0,0), outputDim = 0),
  list(name = 'modulo of positive scalars', expr = quote(out <- abs(arg1) %% abs(arg2)), inputDim = c(0,0), outputDim = 0), 
  list(name = 'modulo of positive and negative scalars', expr = quote(out <- -abs(arg1) %% abs(arg2)), inputDim = c(0,0), outputDim = 0, knownFailure = ".*Cpp"),  ## C and R differ in handling when one input is positive and one negative
  list(name = 'min of scalars', expr = quote(out <- min(arg1, arg2)), inputDim = c(0,0), outputDim = 0),
  list(name = 'max of scalars', expr = quote(out <- max(arg1, arg2)), inputDim = c(0,0), outputDim = 0),
  ###
  list(name = 'vector + vector', expr = quote(out <- arg1 + arg2), inputDim = c(1,1), outputDim = 1),
  list(name = 'diff of vectors', expr = quote(out <- arg1 - arg2), inputDim = c(1,1), outputDim = 1),
  list(name = 'product of vectors', expr = quote(out <- arg1 * arg2), inputDim = c(1,1), outputDim = 1),
  list(name = 'ratio of vectors', expr = quote(out <- arg1 / arg2), inputDim = c(1,1), outputDim = 1),
  list(name = 'power of vectors via ^', expr = quote(out <- arg1 ^ arg2), inputDim = c(1,1), outputDim = 1, knownFailure = 'math.*compiles', knownFailureReport = TRUE), ## fails with Eigen casting: second argument to pow must be a scalar
  list(name = 'power of vectors via pow', expr = quote(out <- pow(arg1, arg2)), inputDim = c(1,1), outputDim = 1, knownFailure = 'math.*compiles', knownFailureReport = TRUE), ## fails with Eigen casting: second argument to pow must be a scalar## FAILS with Eigen casting
  list(name = 'modulo of vectors', expr = quote(out <- abs(arg1) %% abs(arg2)), inputDim = c(1,1), outputDim = 1, knownFailure = '.*compiles'), ## not set up to work for vectors
  list(name = 'pmin of vectors', expr = quote(out <- pmin(arg1, arg2)), inputDim = c(1,1), outputDim = 1),
  list(name = 'pmax of vectors', expr = quote(out <- pmax(arg1, arg2)), inputDim = c(1,1), outputDim = 1),
  ###
  list(name = 'vector + scalar', expr = quote(out <- arg1 + arg2), inputDim = c(1,0), outputDim = 1),
  list(name = 'diff of vector and scalar', expr = quote(out <- arg1 + arg2), inputDim = c(1,0), outputDim = 1),
  list(name = 'product of vector and scalar', expr = quote(out <- arg1 + arg2), inputDim = c(1,0), outputDim = 1),
  list(name = 'ratio of vector and scalar', expr = quote(out <- arg1 + arg2), inputDim = c(1,0), outputDim = 1),
  list(name = 'power of vector and scalar via ^', expr = quote(out <- arg1 ^ arg2), inputDim = c(1,0), outputDim = 1),
  list(name = 'power of vector and scalar via pow', expr = quote(out <- pow(arg1, arg2)), inputDim = c(1,0), outputDim = 1),
  list(name = 'power of vector and constant via ^', expr = quote(out <- arg1 ^ 2), inputDim = c(1,0), outputDim = 1),
  list(name = 'power of vector and constant via pow', expr = quote(out <- pow(arg1, 2)), inputDim = c(1,0), outputDim = 1),
  list(name = 'power of vector and scalar via ^ with positive first arg', expr = quote(out <- exp(arg1) ^ arg2), inputDim = c(1,0), outputDim = 1),
  list(name = 'power of vector and scalar via pow with positive first arg', expr = quote(out <- pow(exp(arg1), arg2)), inputDim = c(1,0), outputDim = 1),
  list(name = 'modulo of vector and scalar', expr = quote(out <- abs(arg1) %% abs(arg2)), inputDim = c(1,0), outputDim = 1, knownFailure = '.*compiles'), ## not set up to work for vectors
  list(name = 'besselK of scalar', expr = quote(out <- besselK(exp(arg1), arg2)), inputDim = c(0, 0), outputDim = 0),
  list(name = 'besselK of scalar, expon.scaled', expr = quote(out <- besselK(exp(arg1), arg2, 1)), inputDim = c(0, 0), outputDim = 0),
  list(name = 'besselK of vector with vector arg', expr = quote(out <- besselK(exp(arg1), arg2)), inputDim = c(1, 1), outputDim = 1)
  )

testsMoreMath = list(
  list(name = 'inverse cloglog of scalar', expr = quote(out <- icloglog(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'cloglog/inverse cloglog of scalar', expr = quote(out <- cloglog(icloglog(arg1))), inputDim = 0, outputDim = 0),
  list(name = 'inverse logit of scalar', expr = quote(out <- ilogit(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'expit of scalar', expr = quote(out <- expit(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'logit/expit of scalar', expr = quote(out <- logit(expit(arg1))), inputDim = 0, outputDim = 0),
  list(name = 'inverse probit of scalar', expr = quote(out <- iprobit(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'inverse probit of scalar via phi', expr = quote(out <- phi(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'probit/iprobit of scalar', expr = quote(out <- probit(iprobit(arg1))), inputDim = 0, outputDim = 0),
  ###
  list(name = 'ceiling of scalar', expr = quote(out <- ceiling(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'floor of scalar', expr = quote(out <- floor(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'round of scalar', expr = quote(out <- round(arg1)), inputDim = 0, outputDim = 0),
  list(name = 'trunc of scalar', expr = quote(out <- trunc(arg1)), inputDim = 0, outputDim = 0),
  ###
  list(name = 'gamma of scalar', expr = quote(out <- gamma(abs(arg1))), inputDim = 0, outputDim = 0),
  list(name = 'lgamma of scalar', expr = quote(out <- lgamma(abs(arg1))), inputDim = 0, outputDim = 0),
  list(name = 'loggam of scalar', expr = quote(out <- loggam(abs(arg1))), inputDim = 0, outputDim = 0),
  list(name = 'log1p of scalar', expr = quote(out <- log1p(abs(arg1))), inputDim = 0, outputDim = 0),
  list(name = 'factorial of scalar', expr = quote(out <- factorial(ceiling(abs(arg1)))), inputDim = 0, outputDim = 0),
  list(name = 'lfactorial of scalar', expr = quote(out <- lfactorial(ceiling(abs(arg1)))), inputDim = 0, outputDim = 0),
  ###
  list(name = 'inverse cloglog of vector', expr = quote(out <- icloglog(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'cloglog/inverse cloglog of vector', expr = quote(out <- cloglog(icloglog(arg1))), inputDim = 1, outputDim = 1),
  list(name = 'inverse logit of vector', expr = quote(out <- ilogit(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'expit of vector', expr = quote(out <- expit(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'logit/expit of vector', expr = quote(out <- logit(expit(arg1))), inputDim = 1, outputDim = 1),
  list(name = 'inverse probit of vector', expr = quote(out <- iprobit(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'inverse probit of vector via phi', expr = quote(out <- phi(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'probit/iprobit of vector', expr = quote(out <- probit(iprobit(arg1))), inputDim = 1, outputDim = 1),
  ###
  list(name = 'ceiling of vector', expr = quote(out <- ceiling(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'floor of vector', expr = quote(out <- floor(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'round of vector', expr = quote(out <- round(arg1)), inputDim = 1, outputDim = 1),
  list(name = 'trunc of vector', expr = quote(out <- trunc(arg1)), inputDim = 1, outputDim = 1),
  ###
  list(name = 'gamma of vector', expr = quote(out <- gamma(abs(arg1))), inputDim = 1, outputDim = 1),
  list(name = 'lgamma of vector', expr = quote(out <- lgamma(abs(arg1))), inputDim = 1, outputDim = 1),
  list(name = 'loggam of vector', expr = quote(out <- loggam(abs(arg1))), inputDim = 1, outputDim = 1),
  list(name = 'log1p of vector', expr = quote(out <- log1p(abs(arg1))), inputDim = 1, outputDim = 1),
  list(name = 'factorial of vector', expr = quote(out <- factorial(ceiling(abs(arg1)))), inputDim = 1, outputDim = 1),
  list(name = 'lfactorial of vector', expr = quote(out <- lfactorial(ceiling(abs(arg1)))), inputDim = 1, outputDim = 1)
  )

testsReduction = list(
  ### vector
  list(name = 'min of vector', expr = quote(out <- min(arg1)), inputDim = 1, outputDim = 0),
  list(name = 'max of vector', expr = quote(out <- max(arg1)), inputDim = 1, outputDim = 0),
  list(name = 'sum of vector', expr = quote(out <- sum(arg1)), inputDim = 1, outputDim = 0),
  list(name = 'mean of vector', expr = quote(out <- mean(arg1)), inputDim = 1, outputDim = 0),
  list(name = 'sd of vector', expr = quote(out <- sd(arg1)), inputDim = 1, outputDim = 0),
  list(name = 'var of vector', expr = quote(out <- var(arg1)), inputDim = 1, outputDim = 0),
  list(name = 'prod of vector', expr = quote(out <- prod(arg1)), inputDim = 1, outputDim = 0),
  list(name = 'norm of vector', expr = quote(out <- norm(arg1)), inputDim = 1, outputDim = 0, knownFailure = '(.*compiles|.*runs)', expectWarnings = 'builds'),  ## norm doesn't work on vector in R and, in addition, is disabled because of C-R inconsistency so compilation fails as well
  ### matrix
  list(name = 'min of matrix', expr = quote(out <- min(arg1)), inputDim = 2, outputDim = 0),
  list(name = 'max of matrix', expr = quote(out <- max(arg1)), inputDim = 2, outputDim = 0),
  list(name = 'sum of matrix', expr = quote(out <- sum(arg1)), inputDim = 2, outputDim = 0),
  list(name = 'mean of matrix', expr = quote(out <- mean(arg1)), inputDim = 2, outputDim = 0),
  list(name = 'sd of matrix', expr = quote(out <- sd(arg1)), inputDim = 2, outputDim = 0),
  list(name = 'var of matrix', expr = quote(out <- var(arg1)), inputDim = 2, outputDim = 0, knownFailure = '.*compiles'),  # Not supported
  list(name = 'prod of matrix', expr = quote(out <- prod(arg1)), inputDim = 2, outputDim = 0),
  list(name = 'norm of matrix', expr = quote(out <- norm(arg1)), inputDim = 2, outputDim = 0, Rcode = quote(out <- norm(arg1, "F")), knownFailure = '(.*compiles|.*runs)', expectWarnings = 'builds') ## NIMBLE's C norm is apparently Frobenius, so R and C nimble functions, and, in addition, is disabled because of C-R inconsistency so compilation fails as well
  )

testsComparison = list(
  ## scalar
  list(name = 'greater than, scalar', expr = quote(out <- arg1 > arg2), inputDim = c(0,0), outputDim = 0, returnType = "logical"),
  list(name = 'equals, scalar', expr = quote(out <- arg1 == arg2), inputDim = c(0,0), outputDim = 0, returnType = "logical"),
  list(name = 'not equals, scalar', expr = quote(out <- arg1 != arg2), inputDim = c(0,0), outputDim = 0, returnType = "logical"),
  ## vector
  list(name = 'greater than, vector', expr = quote(out <- arg1 > arg2), inputDim = c(1,1), outputDim = 1, returnType = "logical"), 
  list(name = 'equals, vector', expr = quote(out <- arg1 == arg2), inputDim = c(1,1), outputDim = 1, returnType = "logical"),  
  list(name = 'not equals, vector', expr = quote(out <- arg1 != arg2), inputDim = c(1,1), outputDim = 1, returnType = "logical"),  
  ## logical scalar
  list(name = 'and operator, scalar', expr = quote(out <- arg1 & arg2), inputDim = c(0,0), outputDim = 0, logicalArgs = c(TRUE, TRUE), returnType = "logical"),
  list(name = 'or operator, scalar', expr = quote(out <- arg1 | arg2), inputDim = c(0,0), outputDim = 0, logicalArgs = c(TRUE, TRUE), returnType = "logical"),
  list(name = 'not operator, scalar', expr = quote(out <- !arg1), inputDim = c(0), outputDim = 0, logicalArgs = c(TRUE), returnType = "logical"),
  ## logical vector
  list(name = 'and operator, vector', expr = quote(out <- arg1 & arg2), inputDim = c(1,1), outputDim = 1, logicalArgs = c(TRUE, TRUE), returnType = "logical"),
  list(name = 'or operator, vector', expr = quote(out <- arg1 | arg2), inputDim = c(1,1), outputDim = 1, logicalArgs = c(TRUE, TRUE), returnType = "logical"),
  list(name = 'not operator, vector', expr = quote(out <- !arg1), inputDim = c(1), outputDim = 1, logicalArgs = c(TRUE), returnType = "logical")
)


testsMatrix = list(
    list(name = 'forwardsolve matrix-vector', expr = quote(out <- forwardsolve(arg1, arg2)), inputDim = c(2, 1), outputDim = 1),
    list(name = 'forwardsolve matrix-matrix', expr = quote(out <- forwardsolve(arg1, arg2)), inputDim = c(2, 2), outputDim = 2),
    list(name = 'backsolve matrix-vector', expr = quote(out <- backsolve(arg1, arg2)), inputDim = c(2, 1), outputDim = 1),
    list(name = 'backsolve matrix-matrix', expr = quote(out <- backsolve(arg1, arg2)), inputDim = c(2, 2), outputDim = 2),

    list(name = 'forwardsolve matrix-vector with indices', expr = quote(out <- forwardsolve(arg1[1:2,1:2], arg2[1:2])), inputDim = c(2, 1), outputDim = 1),
    list(name = 'forwardsolve matrix-matrix with indices', expr = quote(out <- forwardsolve(arg1[1:2,1:2], arg2[1:2,1:2])), inputDim = c(2, 2), outputDim = 2),
    list(name = 'backsolve matrix-vector with indices', expr = quote(out <- backsolve(arg1[1:2,1:2], arg2[1:2])), inputDim = c(2, 1), outputDim = 1),
    list(name = 'backsolve matrix-matrix with indices', expr = quote(out <- backsolve(arg1[1:2,1:2], arg2[1:2,1:2])), inputDim = c(2, 2), outputDim = 2),

    list(name = 'forwardsolve matrix-vector amid expr', expr = quote(out <- arg2[1:2] + forwardsolve(arg1[1:2,1:2], arg2[1:2] + arg2[1:2])), inputDim = c(2, 1), outputDim = 1),
    list(name = 'forwardsolve matrix-matrix amid expr', expr = quote(out <- arg2[1:2,1:2] + forwardsolve(arg1[1:2,1:2], arg2[1:2,1:2] + arg2[1:2,1:2])), inputDim = c(2, 2), outputDim = 2),
    list(name = 'backsolve matrix-vector amid expr', expr = quote(out <- arg2[1:2] + backsolve(arg1[1:2,1:2], arg2[1:2] + arg2[1:2])), inputDim = c(2, 1), outputDim = 1),
    list(name = 'backsolve matrix-matrix amid expr', expr = quote(out <- arg2[1:2,1:2] + backsolve(arg1[1:2,1:2], arg2[1:2,1:2] + arg2[1:2,1:2])), inputDim = c(2, 2), outputDim = 2),

    list(name = 'chol', expr = quote({ A <- arg1; for(i in 1:dim(A)[1]) A[i,i] <- A[i,i] + 10; out <- chol(A) }), inputDim = c(2), outputDim = 2),
    list(name = 'matrix-vector multiply', expr = quote(out <- arg1 %*% arg2), inputDim = c(2, 1), outputDim = 2),
    list(name = 'vector-matrix multiply', expr = quote(out <- t(arg1) %*% arg2), inputDim = c(1, 2), outputDim = 2),
    list(name = 'matrix-matrix multiply', expr = quote(out <- arg1 %*% arg2), inputDim = c(2, 2), outputDim = 2)
)





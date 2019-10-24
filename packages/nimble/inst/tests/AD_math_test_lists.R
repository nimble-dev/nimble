###################
## create R aliases
###################

## could instead use an inverted version of nimble:::specificCallReplacements in
## test_AD but for now we need to ensure there is an R function with the same
## name as every operator tested

gammafn <- gamma
lgammafn <- lgamma
ceil <- ceiling
ftrunc <- trunc
squaredNorm <- function(x) sum(x^2)

##################
## unary cwise ops
##################

unaryArgs <- c('double(0)', 'double(1, 4)', 'double(2, c(3, 4))')
unaryOps <- c(
  '-', nimble:::unaryDoubleOperators,
  nimble:::unaryPromoteNoLogicalOperators
)

unaryOpTests <- make_AD_test_batch(
  unaryOps, unaryArgs
)

## tranform args
modify_on_match(unaryOpTests, '(log|sqrt) .+', 'input_gen_funs', function(x) abs(rnorm(x)))
modify_on_match(unaryOpTests, 'log1p .+', 'input_gen_funs', function(x) abs(rnorm(x)) - 1)
modify_on_match(unaryOpTests, '(logit|probit|cloglog) .+', 'input_gen_funs', runif)
modify_on_match(unaryOpTests, '(acos|asin|atanh) .+', 'input_gen_funs', function(x) runif(x, -1, 1))
modify_on_match(unaryOpTests, 'acosh .+', 'input_gen_funs', function(x) abs(rnorm(x)) + 1)

## see AD_knownFailures entry for factorial
modify_on_match(unaryOpTests, '^(gammafn|factorial) .+', 'input_gen_funs', function(x) abs(rnorm(x)))

## abs fails on negative inputs
modify_on_match(unaryOpTests, 'abs .+', 'input_gen_funs', function(x) abs(rnorm(x)))

## e.g. of how to set tolerances
## defaults are:
##   tol1 = 1e-8 (values)
##   tol2 = 1e-7 (jacobians)
##   tol3 = 1e-6 (hessians)
modify_on_match(unaryOpTests, 'ilogit .+', 'tol2', 1e-7)

## lgammafn gets large Hessians and needs looser tolerance
modify_on_match(unaryOpTests, 'lgammafn .+', 'tol3', 0.0003)

######################
## unary reduction ops
######################

unaryReductionArgs <- c('double(1, 4)', 'double(2, c(3, 4))')

unaryReductionOps <- c(
  'sum', nimble:::reductionUnaryDoubleOperatorsEither,
  nimble:::reductionUnaryOperatorsArray
)

unaryReductionOpTests <- make_AD_test_batch(
  unaryReductionOps, unaryReductionArgs
)

#############
## binary ops
#############

## does not include combinations of vector and matrix
binaryArgs <- as.list(
  cbind(
    data.frame(t(expand.grid(unaryArgs[1], unaryArgs)), stringsAsFactors=FALSE),
    data.frame(t(expand.grid(unaryArgs[2:3], unaryArgs[1])), stringsAsFactors=FALSE)
  )
)
names(binaryArgs) <- NULL
binaryArgs[[length(binaryArgs) + 1]] <- rep(unaryArgs[2], 2)
binaryArgs[[length(binaryArgs) + 1]] <- rep(unaryArgs[3], 2)

binaryOps <- c(
  nimble:::binaryOrUnaryOperators,
  '/', '*', '%%'
)

binaryOpTests <- make_AD_test_batch(
  binaryOps, binaryArgs
)

## set tolerances (default is tol1 = 0.00001 and tol2 = 0.0001)
##modify_on_match(binaryOpTests, '(\\+|-) double\\(0\\) double\\(1, 4\\)', 'tol2', 0.001)
##modify_on_match(binaryOpTests, '(\\+|-) double\\(1, 4\\) double\\(0\\)', 'tol2', 0.001)
##modify_on_match(binaryOpTests, '(\\+|-) double\\(1, 4\\) double\\(1, 4\\)', 'tol2', 0.001)

## runtime failures

## example of specifying when a particular method fails:
## modify_on_match(
##   binaryOpTests,
##   '\\+ double\\(0\\) double\\(0\\)',
##   'knownFailures',
##   list(
##     method3 = list( ## arg1, arg2
##       jacobian = expect_failure,
##       hessian = expect_failure
##     ),
##     method4 = list( ## no wrt
##       jacobian = expect_failure,
##       hessian = expect_failure
##     )
##   )
## )

###############
## pow-like ops
###############

powArgs <- list(
  c('double(0)', 'double(0)'),
  c('double(1, 4)', 'double(0)')
)

powOps <- c(
  'pow', '^'
)

powOpTests <- make_AD_test_batch(
  powOps, powArgs
)

#######################
## binary reduction ops
#######################

binaryReductionArgs <- list(
  c('double(1, 4)', 'double(1, 4)')
)

binaryReductionOps <- nimble:::reductionBinaryOperatorsEither

binaryReductionOpTests <- make_AD_test_batch(
  binaryReductionOps, binaryReductionArgs
)

##########################
## unary square matrix ops
##########################

squareMatrixArgs <- list('double(2, c(2, 2))', 'double(2, c(5, 5))')

squareMatrixOps <- c(nimble:::matrixSquareOperators,
                     nimble:::matrixSquareReductionOperators)

squareMatrixOpTests <- make_AD_test_batch(
  squareMatrixOps, squareMatrixArgs
)

modify_on_match(
  squareMatrixOpTests, 'chol .+', 'input_gen_funs',
  gen_pos_def_matrix ## see AD_test_utils.R
)

## ## compilation failures
## modify_on_match(
##   squareMatrixOpTests, '(logdet|det) .+',
##   'knownFailures', list(compilation = TRUE)
## )

####################
## binary matrix ops
####################

binaryMatrixOps <- nimble:::matrixMultOperators

binaryMatrixArgs <- as.list(
  data.frame(t(expand.grid(
    c('double(2, c(3, 4))', 'double(2, c(4, 4))', 'double(2, c(1, 4))'),
    c('double(2, c(4, 3))', 'double(2, c(4, 4))', 'double(2, c(4, 1))'))),
    stringsAsFactors=FALSE)
)

binaryMatrixOpTests <- make_AD_test_batch(
  binaryMatrixOps, binaryMatrixArgs
)

## set tolerances (default is tol1 = 0.00001 and tol2 = 0.0001)
## modify_on_match(binaryMatrixOpTests, '%*% .+', 'tol2', 0.001)

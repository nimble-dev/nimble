## To-do
## Check atomics with some care.
## replace CppADround with nimDerivs_round

## This file had a major update from an earlier "version 1" to "version 2" in spring 2022.
## Version 1 portions are commented out but retained for easy reference.

###################
## create R aliases
###################

## could instead use an inverted version of nimble:::specificCallReplacements in
## test_AD but for now we need to ensure there is an R function with the same
## name as every operator tested.
## In addition there are now three variants of all tests:
## 1. Just the function itself, f(x) (potentially with variations in argument types)
## 2. The function with a nested argument, f(g(x))
## 3. The function as a nested argument, g(f(x))
## The nested cases are important because the invoke less trivial use of CppAD
## forward and backward mode calculations.  The choice of c() in both cases
## is to have non-trivial (and non-constant) derivatives itself to second order
## while not throwing result values into difficult regions.
## Testing relies on comparisons among R-vs-R (to a minimal extent in these tests),
##   R-vs-C, and C-vs-C (in various modes of obtaining derivatives).
## The R-vs-C in particular makes it tricky to choose tolerances for verification
##   that are not too loose or too tight.

gammafn <- gamma
lgammafn <- lgamma
ceil <- ceiling
ftrunc <- trunc
squaredNorm <- function(x) sum(x^2)

##################
## unary cwise ops
##################

# unaryArgs <- c('double(0)', 'double(1, 4)', 'double(2, c(3, 4))')
unaryOps <- c(
  '-', nimble:::unaryDoubleOperators,
  nimble:::unaryPromoteNoLogicalOperators
)

## # old versions
## unaryOpTests <- make_AD_test_batch(
##   unaryOps, unaryArgs
## )

## ## tranform args
## modify_on_match(unaryOpTests, '(log|sqrt) .+', 'input_gen_funs', function(x) abs(rnorm(x)))
## modify_on_match(unaryOpTests, 'log1p .+', 'input_gen_funs', function(x) abs(rnorm(x)) - 1)
## modify_on_match(unaryOpTests, '(logit|probit|cloglog) .+', 'input_gen_funs', runif)
## modify_on_match(unaryOpTests, '(acos|asin|atanh) .+', 'input_gen_funs', function(x) runif(x, -1, 1))
## modify_on_match(unaryOpTests, 'acosh .+', 'input_gen_funs', function(x) abs(rnorm(x)) + 1)


unaryArgs2 <- c('double(0)', 'double(1)', 'double(2)')
set.seed(123) # seed for randomly generating seeds for each test
unaryOpTests2 <- make_AD_test_batch(
  unaryOps, unaryArgs2, maker = make_AD_test2
)


modify_on_match(unaryOpTests2, '(log|sqrt) .+', 'input_gen_funs', function(x) abs(rnorm(x)))
modify_on_match(unaryOpTests2, 'log1p .+', 'input_gen_funs', function(x) abs(rnorm(x)) - 1)
modify_on_match(unaryOpTests2, '(logit|probit|cloglog) .+', 'input_gen_funs', function(x) runif(x, 0.05, .95))
modify_on_match(unaryOpTests2, '(acos|asin|atanh) .+', 'input_gen_funs', function(x) runif(x, -0.95, 0.95))
modify_on_match(unaryOpTests2, 'acosh .+', 'input_gen_funs', function(x) abs(rnorm(x)) + 1)
modify_on_match(unaryOpTests2, '(factorial|factorial) .+', 'input_gen_funs', function(x) sample(3:10, size = x, replace = TRUE))
modify_on_match(unaryOpTests2, '^cos .+', 'input_gen_funs', function(x) runif(x, 1, 3)) # Avoid high 2nd derivs b/c numerical gradient is inaccurate three
modify_on_match(unaryOpTests2, '^tan .+', 'input_gen_funs', function(x) runif(x, -1.5, 1.5)) #Avoid pi/2=1.57
modify_on_match(unaryOpTests2, '^cosh .+', 'input_gen_funs', function(x) sample(c(-1,1), size = x, replace = TRUE) * runif(x, .2, 1.2)) # Avoid high 2nd derivs b/c numerical gradient is inaccurate there
modify_on_match(unaryOpTests2, '^tanh .+', 'input_gen_funs', function(x) runif(x, -0.8, 0.8)) # Avoid flat rmodify_on_match(unaryOpTests2, '^sin .+', 'input_gen_funs', function(x) runif(x, -1, 1))
gamma_x_vals <- c(-1.8, -1.4, -0.9, -0.2, 0.2, 0.8, 1, 1.2, 2, 2.5)
modify_on_match(unaryOpTests2, '^gammafn .+', 'input_gen_funs', function(x) sample(gamma_x_vals, size = x, replace = TRUE))
modify_on_match(unaryOpTests2, '^cube .+', 'input_gen_funs', function(x) sample(c(seq(-1.5, -.2, length = 20), seq(.2, 1.5, length = 20)), size = x, replace = TRUE))

###########################################################
# f(g(x)) where f is function being tested and g(x) = x^3 

set.seed(123) # Ok to use same as above
unaryOpTests2_inner <- make_AD_test_batch(
  unaryOps, unaryArgs2, maker = make_AD_test2, inner_codes = list(quote(X*X*X))
)

# ADtestEnv$RCrelTol <- c(1e-15, 1e-4, 1e-2) ## Loosen tols because there are more operations

sCbRt <- function(x) sign(x) * abs(x)^(1/3) #signed cube root. This preserves sign and scale of input to f since g(x) = x^3 
modify_on_match(unaryOpTests2_inner, '(log|sqrt) .+', 'input_gen_funs', function(x) abs(rnorm(x)))
modify_on_match(unaryOpTests2_inner, 'log1p .+', 'input_gen_funs', function(x) sCbRt(abs(rnorm(x)) - 1))
modify_on_match(unaryOpTests2_inner, '(logit|probit|cloglog) .+', 'input_gen_funs', function(x) sCbRt(runif(x, 0.05, 0.95)))
modify_on_match(unaryOpTests2_inner, '(acos|asin|atanh) .+', 'input_gen_funs', function(x) sCbRt(runif(x, -0.95, 0.95)))
modify_on_match(unaryOpTests2_inner, 'acosh .+', 'input_gen_funs', function(x) sCbRt(abs(rnorm(x)) + 1.05))
modify_on_match(unaryOpTests2_inner, '(factorial|factorial) .+', 'input_gen_funs', function(x) sample(1:3, size = x, replace = TRUE))
modify_on_match(unaryOpTests2_inner, '^cos .+', 'input_gen_funs', function(x) sCbRt(runif(x, 1, 3))) # Avoid high 2nd derivs b/c numerical gradient is inaccurate there
modify_on_match(unaryOpTests2_inner, '^cosh .+', 'input_gen_funs', function(x) sCbRt(sample(c(-1,1), size = x, replace = TRUE) * runif(x, .2, 1.2))) # Avoid high 2nd derivs b/c numerical gradient is inaccurate there
modify_on_match(unaryOpTests2_inner, '^tanh .+', 'input_gen_funs', function(x) sCbRt(runif(x, -0.8, 0.8))) # Avoid flat regions
modify_on_match(unaryOpTests2_inner, '^sin .+', 'input_gen_funs', function(x) sCbRt(runif(x, -1, 1)))
modify_on_match(unaryOpTests2_inner, '^tan .+', 'input_gen_funs', function(x) sCbRt(runif(x, -1.5, 1.5))) #Avoid pi/2=1.57
modify_on_match(unaryOpTests2_inner, '^gammafn .+', 'input_gen_funs', function(x) sCbRt(sample(gamma_x_vals, size = x, replace = TRUE)))
modify_on_match(unaryOpTests2_inner, '^cube .+', 'input_gen_funs', function(x) sCbRt(sample(c(seq(-1.5, -.2, length = 20), seq(.2, 1.5, length = 20)), size = x, replace = TRUE)))



###########################################################
# g(f(x)) where f is function being tested and g(y) is exp(.5 * y)

set.seed(123) # Ok to use same as above
unaryOpTests2_outer <- make_AD_test_batch(
  unaryOps, unaryArgs2, maker = make_AD_test2, outer_code = quote(exp(0.5 * Y))
)

#ADtestEnv$RCrelTol <- c(1e-15, 1e-6, 1e-2)

modify_on_match(unaryOpTests2_outer, '(log|sqrt) .+', 'input_gen_funs', function(x) abs(rnorm(x)))
modify_on_match(unaryOpTests2_outer, 'log1p .+', 'input_gen_funs', function(x) abs(rnorm(x)) - 1)
modify_on_match(unaryOpTests2_outer, '(logit|probit|cloglog) .+', 'input_gen_funs', function(x) runif(x, 0.05, .95))
modify_on_match(unaryOpTests2_outer, '(acos|asin|atanh) .+', 'input_gen_funs', function(x) runif(x, -0.95, 0.95))
modify_on_match(unaryOpTests2_outer, 'acosh .+', 'input_gen_funs', function(x) abs(rnorm(x)) + 1)
modify_on_match(unaryOpTests2_outer, '(factorial|factorial) .+', 'input_gen_funs', function(x) sample(1:3, size = x, replace = TRUE))
modify_on_match(unaryOpTests2_outer, '^cos .+', 'input_gen_funs', function(x) runif(x, 1, 3)) # Avoid high 2nd derivs b/c numerical gradient is inaccurate three
modify_on_match(unaryOpTests2_outer, '^cosh .+', 'input_gen_funs', function(x) sample(c(-1,1), size = x, replace = TRUE) * runif(x, .2, 1.2)) # Avoid high 2nd derivs b/c numerical gradient is inaccurate there
modify_on_match(unaryOpTests2_outer, '^tanh .+', 'input_gen_funs', function(x) runif(x, -0.8, 0.8)) # Avoid flat regions
modify_on_match(unaryOpTests2_outer, '^sin .+', 'input_gen_funs', function(x) runif(x, -1, 1))
modify_on_match(unaryOpTests2_outer, '^tan .+', 'input_gen_funs', function(x) runif(x, -1.5, 1.5)) #Avoid pi/2=1.57
modify_on_match(unaryOpTests2_outer, '^gammafn .+', 'input_gen_funs', function(x) sample(gamma_x_vals, size = x, replace = TRUE))
modify_on_match(unaryOpTests2_outer, '^cube .+', 'input_gen_funs', function(x) sample(c(seq(-1.5, -.2, length = 20), seq(.2, 1.5, length = 20)), size = x, replace = TRUE))


######################
## unary reduction ops
######################


# nimble:::reductionUnaryOperatorsArray contains only sd and var
# but var(matrix) is not supported because, per note in size processing code,
# in R var(matrix) is interpeted as cov(data.frame) and so it really different.

unaryReductionOps <- c(
  'sum', nimble:::reductionUnaryDoubleOperatorsEither,
  'sd' # Add 'var' below without matrix input
)

## unaryReductionArgs <- c('double(1, 4)', 'double(2, c(3, 4))')
## unaryReductionOpTests <- make_AD_test_batch(
##   unaryReductionOps, unaryReductionArgs
## )
## unaryReductionOpTests <- c(unaryReductionOpTests,
##                            make_AD_test_batch('var', unaryReductionArgs[1]))

unaryReductionArgs2 <- c('double(1)', 'double(2)')
set.seed(123) # Ok to use same as above
unaryReductionOpTests2 <- make_AD_test_batch(
  unaryReductionOps, unaryReductionArgs2, maker = make_AD_test2
)
unaryReductionOpTests2 <- c(unaryReductionOpTests2,
                            make_AD_test_batch(
                              'var', unaryReductionArgs2[1], maker = make_AD_test2))
modify_on_match(unaryReductionOpTests2, 'arg1 = double\\(2\\)',
                'RCrelTol', c(1e-13, 1e-5, 1e-3) ) ## lower tol for matrix (larger) inputs

#res <- lapply(unaryReductionOpTests2, test_AD2)


## f(g(x))
set.seed(123) # Ok to use same as above
unaryReductionOpTests2_inner <- make_AD_test_batch(
    unaryReductionOps, unaryReductionArgs2, maker = make_AD_test2, inner_codes = list(quote(X*X*X))
)
unaryReductionOpTests2_inner <- c(unaryReductionOpTests2_inner,
                            make_AD_test_batch(
                              'var', unaryReductionArgs2[1], maker = make_AD_test2, inner_codes = list(quote(X*X*X))))
modify_on_match(unaryReductionOpTests2_inner, 'sum .+', 'input_gen_funs', function(x) runif(x, .5, 1.5))
modify_on_match(unaryReductionOpTests2_inner, 'mean .+', 'input_gen_funs', function(x) runif(x, .5, 1.5))
modify_on_match(unaryReductionOpTests2_inner, 'squaredNorm .+', 'input_gen_funs', function(x) runif(x, .5, 1.5))
modify_on_match(unaryReductionOpTests2_inner, 'prod .+', 'input_gen_funs', function(x) runif(x, .5, 1.5))
modify_on_match(unaryReductionOpTests2_inner, 'sd .+', 'input_gen_funs', function(x) runif(x, .5, 1.5))
modify_on_match(unaryReductionOpTests2_inner, 'var .+', 'input_gen_funs', function(x) runif(x, .5, 1.5))

## g(f(x)
set.seed(123) # Ok to use same as above
unaryReductionOpTests2_outer <- make_AD_test_batch(
    unaryReductionOps, unaryReductionArgs2, maker = make_AD_test2, outer_code = quote(exp(0.5 * Y))
)
unaryReductionOpTests2_outer <- c(unaryReductionOpTests2_outer,
                            make_AD_test_batch(
                              'var', unaryReductionArgs2[1], maker = make_AD_test2,  outer_code = quote(exp(0.5 * Y))))
modify_on_match(unaryReductionOpTests2_outer, 'sum .+', 'input_gen_funs', function(x) runif(x, .5, 1.5))
modify_on_match(unaryReductionOpTests2_outer, 'mean .+', 'input_gen_funs', function(x) runif(x, .5, 1.5))
modify_on_match(unaryReductionOpTests2_outer, 'squaredNorm .+', 'input_gen_funs', function(x) runif(x, .5, 1.5))
modify_on_match(unaryReductionOpTests2_outer, 'prod .+', 'input_gen_funs', function(x) runif(x, .5, 1.5))
modify_on_match(unaryReductionOpTests2_outer, 'sd .+', 'input_gen_funs', function(x) runif(x, .5, 1.5))
modify_on_match(unaryReductionOpTests2_outer, 'var .+', 'input_gen_funs', function(x) runif(x, .5, 1.5))


#############
## binary ops
#############

## does not include combinations of vector and matrix
## binaryArgs <- as.list(
##   cbind(
##     data.frame(t(expand.grid(unaryArgs[1], unaryArgs)), stringsAsFactors=FALSE),
##     data.frame(t(expand.grid(unaryArgs[2:3], unaryArgs[1])), stringsAsFactors=FALSE)
##   )
## )
## names(binaryArgs) <- NULL
## binaryArgs[[length(binaryArgs) + 1]] <- rep(unaryArgs[2], 2)
## binaryArgs[[length(binaryArgs) + 1]] <- rep(unaryArgs[3], 2)

binaryArgs2 <- as.list(
  cbind(
    data.frame(t(expand.grid(unaryArgs2[1], unaryArgs2)), stringsAsFactors=FALSE),
    data.frame(t(expand.grid(unaryArgs2[2:3], unaryArgs2[1])), stringsAsFactors=FALSE)
  )
)
names(binaryArgs2) <- NULL
binaryArgs2[[length(binaryArgs2) + 1]] <- rep(unaryArgs2[2], 2)
binaryArgs2[[length(binaryArgs2) + 1]] <- rep(unaryArgs2[3], 2)

binaryOps <- c(
  nimble:::binaryOrUnaryOperators,
  '/', '*', '%%' # %% is not supported and is on knownFailures list
)

## binaryOpTests <- make_AD_test_batch(
##   binaryOps, binaryArgs
## )
## modify_on_match(binaryOpTests, "/ arg1 = double\\(0\\) arg2 = double\\(2, c\\(3, 4\\)\\)", 'input_gen_funs', function(x) {res <- rnorm(x);res <- res + sign(res)*0.1; res}) ## This generator avoids small numbers that, in the denominator, give huge hes

set.seed(123)
binaryOpTests2 <- make_AD_test_batch(
  binaryOps, binaryArgs2, maker = make_AD_test2
)
modify_on_match(binaryOpTests2, "/ arg1 = double\\(0\\) arg2 = double\\(2\\)", 'input_gen_funs',
                function(x) {res <- rnorm(x);res <- res + sign(res)*0.1; res}) ## This generator avoids small numbers that, in the denominator, give huge hessians that are hard to numerically match from finite elements

## f(g(x1), g(x2))
set.seed(123)
binaryOpTests2_inner <- make_AD_test_batch(
  binaryOps, binaryArgs2, maker = make_AD_test2, inner_codes = list(quote(X*X*X), quote(X*X*X))
)
modify_on_match(binaryOpTests2_inner, "", 'input_gen_funs', function(x) {res <- rnorm(x);res <- res + sign(res)*0.1; res}) ## This generator avoids small numbers that, in the denominator, give huge

## g(f(x1, x2))
set.seed(123)
binaryOpTests2_outer <- make_AD_test_batch(
  binaryOps, binaryArgs2, maker = make_AD_test2, outer_code = quote(exp(0.5 * Y))
)
modify_on_match(binaryOpTests2_outer, "", 'input_gen_funs', function(x) {res <- rnorm(x);res <- res + sign(res)*0.1; res}) ## This generator avoids small numbers that, in the denominator, give huge

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

## powArgs <- list(
##   c('double(0)', 'double(0)'),
##   c('double(1, 4)', 'double(0)')
## )
## powOpTests <- make_AD_test_batch(
##   powOps, powArgs
## )
## pow_int_OpTests <- list(
##   make_AD_test(c('double(0)', 'integer(0)'), op = 'pow_int', wrt_args = 'arg1'),
##   make_AD_test(c('double(1, 4)', 'integer(0)'), op = 'pow_int', wrt_args = 'arg1')
## )

powOps <- c(
  'pow', '^'
)
powArgs2 <- list(
  c('double(0)', 'double(0)'),
  c('double(1)', 'double(0)'),
  c('double(2)', 'double(0)')
)
powOpTests2 <- make_AD_test_batch(
  powOps, powArgs2, maker = make_AD_test2
)
modify_on_match(powOpTests2, '', 'input_gen_funs', function(x) runif(x, .5, 1.5)) # a^b now defined only valid for both > 0

## STOPPED HERE


# TO-DO: TEST NON-INTEGER INPUT FOR b
# TO-DO: Add matrix^scalar case above

# This is a little delicated and should be explained in the manual.
# pow_int should still be given a double as a second argument it will
# not end up flowing through the derivatives system correctly.
# There could be a crash, or, worse, incorrect results when the value is baked in.
pow_int_OpTests2 <- list(
  make_AD_test2(c('double(0)', 'double(0)'), op = 'pow_int', wrt_args = c('arg1', 'arg2'),
                input_gen_funs = list(arg1 = function(x) rnorm(x),
                                      arg2 = function(x) sample(-3:3, size = x, replace = TRUE))),
  make_AD_test2(c('double(1)', 'double(0)'), op = 'pow_int', wrt_args = c('arg1', 'arg2'),
                input_gen_funs = list(arg1 = function(x) rnorm(x),
                                      arg2 = function(x) sample(-3:3, size = x, replace = TRUE))),
  make_AD_test2(c('double(2)', 'double(0)'), op = 'pow_int', wrt_args = c('arg1', 'arg2'),
                input_gen_funs = list(arg1 = function(x) rnorm(x),
                                      arg2 = function(x) sample(-3:3, size = x, replace = TRUE)))
)

resetTols()

res <- lapply(powOpTests2, test_AD2)
res <- lapply(pow_int_OpTests2, test_AD2)

# f(g(x)).  Use a different g than above because that is basically a power function

powOpTests2 <- make_AD_test_batch(
  powOps, powArgs2, maker = make_AD_test2, inner_codes = list(quote(0.5*exp(X)), quote(X*X*X))
)

modify_on_match(powOpTests2, '', 'input_gen_funs', function(x) runif(x, .5, 1.5)) # a^b now defined only valid for both > 0
pow_int_OpTests2 <- list(
  make_AD_test2(c('double(0)', 'double(0)'), op = 'pow_int', wrt_args = c('arg1', 'arg2'),
                input_gen_funs = list(arg1 = function(x) rnorm(x),
                                      arg2 = function(x) sample(-3:3, size = x, replace = TRUE)),
                inner_codes = list(quote(0.5*(exp(X/2)-1)), quote(1.1*X))),
  make_AD_test2(c('double(1)', 'double(0)'), op = 'pow_int', wrt_args = c('arg1', 'arg2'),
                input_gen_funs = list(arg1 = function(x) rnorm(x),
                                      arg2 = function(x) sample(-3:3, size = x, replace = TRUE)),
                inner_codes = list(quote(0.5*(exp(X/2)-1)), quote(1.1*X))),
  make_AD_test2(c('double(2)', 'double(0)'), op = 'pow_int', wrt_args = c('arg1', 'arg2'),
                input_gen_funs = list(arg1 = function(x) rnorm(x),
                                      arg2 = function(x) sample(-3:3, size = x, replace = TRUE)),
                inner_codes = list(quote(0.5*(exp(X/2)-1)), quote(1.1*X)))
)
# The 0.5*exp(X/2-1) is made up to give non-trivial derivs and some <0 and >0 values.
# The 1.1*X doesn't do much but checks that round(b) is used and at least isn't nothing.

resetTols()

res <- lapply(powOpTests2, test_AD2)
res <- lapply(pow_int_OpTests2, test_AD2)

# g(f(x))

powOpTests2 <- make_AD_test_batch(
  powOps, powArgs2, maker = make_AD_test2, outer_code = quote(exp(0.5*Y))
)

modify_on_match(powOpTests2, '', 'input_gen_funs', function(x) runif(x, .5, 1.5)) 
pow_int_OpTests2 <- list(
  make_AD_test2(c('double(0)', 'double(0)'), op = 'pow_int', wrt_args = c('arg1', 'arg2'),
                input_gen_funs = list(arg1 = function(x) rnorm(x),
                                      arg2 = function(x) sample(-3:3, size = x, replace = TRUE)),
                outer_code = quote(exp(0.5*Y))),
  make_AD_test2(c('double(1)', 'double(0)'), op = 'pow_int', wrt_args = c('arg1', 'arg2'),
                input_gen_funs = list(arg1 = function(x) rnorm(x),
                                      arg2 = function(x) sample(-3:3, size = x, replace = TRUE)),
                outer_code = quote(exp(0.5*Y))),
  make_AD_test2(c('double(2)', 'double(0)'), op = 'pow_int', wrt_args = c('arg1', 'arg2'),
                input_gen_funs = list(arg1 = function(x) rnorm(x),
                                      arg2 = function(x) sample(-3:3, size = x, replace = TRUE)),
                outer_code = quote(exp(0.5*Y)))
)

resetTols()

ADtestEnv$RCrelTol <- c(1e-12, 1e-4, 1e-2) # This seems to need looser tols, based on one case.

res <- lapply(powOpTests2, test_AD2)
res <- lapply(pow_int_OpTests2, test_AD2)

resetTols()


#######################
## binary reduction ops
#######################

binaryReductionArgs <- list(
  c('double(1, 4)', 'double(1, 4)')
)

binaryReductionArgs2 <- list(
  c('double(1)', 'double(1)')
)

binaryReductionOps <- nimble:::reductionBinaryOperatorsEither

binaryReductionOpTests <- make_AD_test_batch(
  binaryReductionOps, binaryReductionArgs
)

binaryReductionOpTests2 <- make_AD_test_batch(
  binaryReductionOps, binaryReductionArgs2, maker = make_AD_test2
)
resetTols()
res <- lapply(binaryReductionOpTests2, test_AD2)

binaryReductionOpTests2 <- make_AD_test_batch(
  binaryReductionOps, binaryReductionArgs2, maker = make_AD_test2,
  inner_codes = list(quote(X*X*X), quote(X*X*X))
)
resetTols()
ADtestEnv$RCrelTol <- c(1e-12, 1e-4, 1e-2)
res <- lapply(binaryReductionOpTests2, test_AD2)

binaryReductionOpTests2 <- make_AD_test_batch(
  binaryReductionOps, binaryReductionArgs2, maker = make_AD_test2,
  outer_code = quote(exp(0.5*Y))
)
resetTols()
ADtestEnv$RCrelTol <- c(1e-12, 1e-4, 1e-2)
res <- lapply(binaryReductionOpTests2, test_AD2)


##########################
## unary square matrix ops
##########################

squareMatrixArgs <- list('double(2, c(2, 2))', 'double(2, c(5, 5))')
squareMatrixArgs2 <- list('double(2)')
# To-Do: run with atomics on and off.

squareMatrixOps <- c(nimble:::matrixSquareOperators,
                     nimble:::matrixSquareReductionOperators)

squareMatrixOpTests <- make_AD_test_batch(
  squareMatrixOps, squareMatrixArgs
)

modify_on_match(
  squareMatrixOpTests, 'chol .+', 'input_gen_funs',
  gen_pos_def_matrix ## see AD_test_utils.R
)

squareMatrixOpTests2a <- make_AD_test_batch(squareMatrixOps, squareMatrixArgs2, maker = make_AD_test2)
modify_on_match(squareMatrixOpTests2a, '', 'size', c(2, 2)) ## see AD_test_utils.R
modify_on_match(squareMatrixOpTests2a, 'chol .+', 'input_gen_funs', function(x) gen_pos_def_matrix(x))

squareMatrixOpTests2b <- make_AD_test_batch(squareMatrixOps, squareMatrixArgs2, maker = make_AD_test2)
modify_on_match(squareMatrixOpTests2b, '', 'size', c(5, 5)) ## see AD_test_utils.R
modify_on_match(squareMatrixOpTests2b, 'chol .+', 'input_gen_funs', function(x) gen_pos_def_matrix(x))

squareMatrixOpTests2 <- c(squareMatrixOpTests2a, squareMatrixOpTests2b)

# modify_on_match after c() would fail because only the first of redundant names would get handled

ADtestEnv$RCrelTol <- c(1e-12, 1e-4, 1e-2)

lapply(squareMatrixOpTests2, test_AD2)
debug(test_AD2)
lapply(squareMatrixOpTests2[6], test_AD2) # Chol 5x5 can have some crazy unstable 2nd derivs.
# Hessian of det is numerically nearly zero in different ways - comparison disaster
# This affects log det too
# We evidently don't really support trace?


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

## To-Do: Add solve, forward-solve and back-solve

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
nimbleOptions(experimentalEnableDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)

context("Testing of derivatives for nimbleFunctions.")

test_that('Derivatives of dnorm function correctly.',
  {
    ADfun1 <- nimbleFunction(
      setup = function(){},
      run = function(y = double(1)) {
        outList <- derivs(testMethod(y), wrt = c('x'))
        returnType(ADNimbleList())
        return(outList)
      },
      methods = list(
        testMethod = function(x = double(1, 2)) {
          out <- dnorm(x[1],0,1)
          returnType(double())
          return(out)
        }
      ), enableDerivs = list('testMethod')
    )
    ADfunInst <- ADfun1()
    x <- matrix(c(2, -2))
    Rderivs <- ADfunInst$run(x)
    temporarilyAssignInGlobalEnv(ADfunInst)  
    cADfunInst <- compileNimble(ADfunInst)
    cderivs <- cADfunInst$run(x)
    expect_equal(cderivs$value, Rderivs$value)
    expect_equal(cderivs$jacobian, Rderivs$jacobian)
    expect_equal(cderivs$hessian, Rderivs$hessian)
  }
)



test_that('Derivatives of x^2 function correctly.',
          {
            listOfADLists <- nimbleList(list1 = ADNimbleList(),
                                        list2 = ADNimbleList())
            temporarilyAssignInGlobalEnv(listOfADLists)  
            ADfun2 <- nimbleFunction(
              setup = function(){},
              run = function(y = double(1, c(2))) {
                ADlistList <- listOfADLists$new()
                ADlistList$list1 <- derivs(testMethod(y, y))
                ADlistList$list2 <- derivs(testMethod2(y, y))
                returnType(listOfADLists())
                return(ADlistList)
              },
              methods = list(
                testMethod = function(x = double(1, c(2)), y = double(1, c(2))) {
                  returnType(double(1, c(2)))
                  return(x[]^2)
                },
                testMethod2 = function(x = double(1, c(2)), y = double(1, c(2))) {
                  returnType(double(0))
                  return(2^x[1])
                }
                
              ), enableDerivs = list('testMethod', 'testMethod2')
            )
            ADfunInst <- ADfun2()
            x <- c(1, 1)
            Rderivs <- ADfunInst$run(x)
            temporarilyAssignInGlobalEnv(ADfunInst)  
            cADfunInst <- compileNimble(ADfunInst)
            cderivs <- cADfunInst$run(x)
            expect_equal(cderivs$list1$value, Rderivs$list1$value)
            expect_equal(cderivs$list1$jacobian, Rderivs$list1$jacobian, tolerance = 0.01)
            expect_equal(cderivs$list1$hessian, Rderivs$list1$hessian, tolerance = 0.01)
            expect_equal(cderivs$list2$value, Rderivs$list2$value)
            expect_equal(cderivs$list2$jacobian, Rderivs$list2$jacobian, tolerance = 0.01)
            expect_equal(cderivs$list2$hessian, Rderivs$list2$hessian, tolerance = 0.01)
          }
)

test_that('Derivatives of sum(log(x)) function correctly.',
          {
            ADfun3 <- nimbleFunction(
              setup = function(){},
              run = function(x = double(1, c(2))) {
                outVals <- 2*derivs(testMethod(x))$jacobian
                returnType(double(2))
                return(outVals)
              },
              methods = list(
                testMethod = function(x = double(1, c(2))) {
                  returnType(double(0))
                  return(sum(log(x[])))
                }
              ), enableDerivs = list('testMethod')
            )
            ADfunInst <- ADfun3()
            x <- c(1, 1)
            Rderivs <- ADfunInst$run(x)
            temporarilyAssignInGlobalEnv(ADfunInst)  
            cADfunInst <- compileNimble(ADfunInst)
            cderivs <- cADfunInst$run(x)
            expect_equal(cderivs, Rderivs, tolerance = 0.01)
          }
)


test_that('Derivatives of model$calculate() work in expressions.',
          {
            testModelCode <- nimbleCode({
              x[1:2] ~ dmnorm(zeroVec[1:2], diagMat[1:2, 1:2])
            })
            testModel <- nimbleModel(code = testModelCode,
                                     inits = list(x = rep(0, 2)),
                                     constants = list(zeroVec = rep(0, 2),
                                                      diagMat = diag(2)))
            ADfun4 <- nimbleFunction(
              setup = function(model){},
              run = function() {
                outVals <- 2*sum(derivs(model$calculate('x'), 
                                        wrt = 'x')$jacobian)
                returnType(double())
                return(outVals)
              }
            )
            ADfunInst <- ADfun4(testModel)
            Rderivs <- ADfunInst$run()
            temporarilyAssignInGlobalEnv(testModel)  
            temporarilyAssignInGlobalEnv(ADfunInst)  
            compileNimble(testModel)
            cADfunInst <- compileNimble(ADfunInst)
            cderivs <- cADfunInst$run()
            print(cderivs)
            expect_equal(cderivs, Rderivs, tolerance = 0.01)
          }
)

test_that('Derivatives of matrix multiplication function correctly.',
          {
            ADfun5 <- nimbleFunction(
              setup = function(){},
              run = function(x = double(2, c(2, 2)),
                             y = double(2, c(2, 4))) {
                outVals <- derivs(testMethod(x, y), wrt = 'x', order = 1)$jacobian
                returnType(double(2))
                return(outVals)
              },
              methods = list(
                testMethod = function(x = double(2, c(2, 2)),
                                      y = double(2, c(2, 4))) {
                  returnType(double(2, c(2, 4)))
                  return(x%*%y)
                }
              ), enableDerivs = list('testMethod')
            )
            ADfunInst <- ADfun5()
            x <- diag(2)
            set.seed(1)
            y <- matrix(rnorm(8), ncol = 4, nrow = 2)
            Rderivs <- ADfunInst$run(x, y)
            temporarilyAssignInGlobalEnv(ADfunInst)  
            cADfunInst <- compileNimble(ADfunInst)
            cderivs <- cADfunInst$run(x, y)
            expect_equal(cderivs, Rderivs, tolerance = 0.01)
          }
)



test_that('Derivatives of matrix exponentiation function correctly.',
          {
            ADfun6 <- nimbleFunction(
              setup = function(){},
              run = function(x = double(2, c(2, 2))) {
                outVals <- derivs(testMethod(x), wrt = 'x', order = 1)$jacobian
                returnType(double(2))
                return(outVals)
              },
              methods = list(
                testMethod = function(x = double(2, c(2, 2))) {
                  returnType(double(2, c(2, 2)))
                  return(exp(x))
                }
              ), enableDerivs = list('testMethod')
            )
            ADfunInst <- ADfun6()
            x <- diag(2)
            Rderivs <- ADfunInst$run(x)
            temporarilyAssignInGlobalEnv(ADfunInst)  
            cADfunInst <- compileNimble(ADfunInst)
            cderivs <- cADfunInst$run(x)
            expect_equal(cderivs, Rderivs, tolerance = 0.01)
          }
)


## one for exp(mat)
## get order to owrk without c()!

test_AD <- function(param, size = 3,
                    dir = file.path(tempdir(), "nimble_generatedCode"),
                    control = list(), verbose = nimbleOptions('verbose')) {
  if (!is.null(param$debug) && param$debug) browser()
  if (verbose) cat(paste0("### Testing ", param$name, "\n"))

  nf <- nimbleFunction(
    setup = function() {},
    run = param$run,
    methods = param$methods,
    enableDerivs = param$enableDerivs
  )
  nfInst <- nf()

  opParam <- param$opParam
  input <- lapply(opParam$args, argType2input, dist = param$input_dist)

  if (verbose) cat("## Calling R versions of nimbleFunction methods\n")
  Rderivs <- try(
    sapply(names(param$methods), function(method) {
      do.call(
        ## can't access nfInst[[method]] until $ has been used :(
        eval(substitute(nfInst$method, list(method = as.name(method)))), input
      )
    }, USE.NAMES = TRUE), silent = TRUE
  )
  if (inherits(Rderivs, 'try-error')) {
    warning(
      paste(
        'Calling R version of test', opParam$name,
        'resulted in an error:', Rderivs[1]
      ),
      immediate. = TRUE
    )
    return(invisible(NULL))
  }

  temporarilyAssignInGlobalEnv(nfInst)

  wrap_if_true(param$knownFailures$compilation, expect_error, {
    if (verbose) cat("## Compiling nimbleFunction \n")
    if (!is.null(param$dir)) dir <- param$dir
    CnfInst <- compileNimble(nfInst, dirName = dir)
    Cderivs <- sapply(names(param$methods), function(method) {
      do.call(
        ## same issue as with Rderivs
        eval(substitute(CnfInst$method, list(method = as.name(method)))),
        input
      )
    }, USE.NAMES = TRUE)
    tol1 <- if (is.null(param$tol1)) 0.00001 else param$tol1
    tol2 <- if (is.null(param$tol2)) 0.0001 else param$tol2
    for (method_name in names(param$methods)) {
      if (verbose) {
        cat(paste0(
          "## Testing ", method_name, ': ', param$wrts[[method_name]], '\n'
        ))
      }
      wrap_if_true(
        !is.null(param$knownFailures[[method_name]]$value),
        param$knownFailures[[method_name]]$value, {
          if (verbose) cat("## Checking values\n")
          expect_equal(
            Cderivs[[method_name]]$value,
            Rderivs[[method_name]]$value
          )
        }
      )
      wrap_if_true(
        !is.null(param$knownFailures[[method_name]]$jacobian),
        param$knownFailures[[method_name]]$jacobian, {
          if (verbose) cat("## Checking jacobians\n")
          expect_equal(
            Cderivs[[method_name]]$jacobian,
            Rderivs[[method_name]]$jacobian,
            tolerance = tol1
          )
        }
      )
      wrap_if_true(
        !is.null(param$knownFailures[[method_name]]$hessian),
        param$knownFailures[[method_name]]$hessian, {
          if (verbose) cat("## Checking hessians\n")
          expect_equal(
            Cderivs[[method_name]]$hessian,
            Rderivs[[method_name]]$hessian,
            tolerance = tol2
          )
        }
      )
    }
    if (verbose) cat("### Test successful \n\n")
  })
  invisible(NULL)
}

## Takes an argType call and returns up to two components
## of the arg with respect to which to take derivatives
## in test_AD. One of the components will always be the
## argument itself and the other will be a random element
## within the arg (if the arg is non-scalar).
make_wrt <- function(argType, argname) {
  argSymbol <- nimble:::argType2symbol(argType)
  wrt <- argname
  if (argSymbol$nDim == 2)
    c(
      wrt, paste0(
             argname, '[',
             sample(1:argSymbol$size[1], size = 1) , ',',
             sample(1:argSymbol$size[2], size = 1),
             ']'
           )
    )
  else ## nDim is 0 or 1
    c(wrt, paste0(argname, '[', sample(1:argSymbol$size, size = 1) , ']'))
}

makeADtest <- function(op, argTypes) {
  opType <- switch(
    length(argTypes),
    'unary',
    'binary'
  )
  if (is.null(opType))
    stop('Only unary and binary opeartors are supported.', call. = FALSE)

  opParam <- makeOperatorParam(op, argTypes)

  run <- gen_runFunCore(opParam)
  method <- function() {
    {}
    returnType(ADNimbleList())
    return(outList)
  }
  if (opType == 'unary')
    formals(method) <- list(
      arg1 = parse(text = argTypes[1])[[1]]
    )
  else if (opType == 'binary')
    formals(method) <- list(
      arg1 = parse(text = argTypes[1])[[1]],
      arg2 = parse(text = argTypes[2])[[1]]
    )

  methods <- list()
  wrts <- list()

  wrt <- unlist(mapply(make_wrt, opParam$args, names(opParam$args)))
  for (i in seq_along(wrt)) {
    method_i <- paste0('method', i)
    methods[[method_i]] <- method
    wrts[[method_i]] <- paste0('wrt ', wrt[[i]])
    body(methods[[method_i]])[[2]] <-
      if (opType == 'unary')
        substitute(
          outList <- derivs(run(arg1), wrt = this_wrt),
          list(this_wrt = wrt[[i]])
        )
      else
        substitute(
          outList <- derivs(run(arg1, arg2), wrt = this_wrt),
          list(this_wrt = wrt[[i]])
        )
  }

  method_no_wrt <- paste0('method', length(wrt) + 1)
  wrts[[method_no_wrt]] <- 'no wrt'

  method_all_wrt <- paste0('method', length(wrt) + 2)
  wrts[[method_all_wrt]] <- paste0('wrt ', paste0(wrt, collapse = ', '))
  
  methods[[method_all_wrt]] <- methods[[method_no_wrt]] <-  method

  if (opType == 'unary') {
    body(methods[[method_no_wrt]])[[2]] <- quote(
      outList <- derivs(run(arg1))
    )
    body(methods[[method_all_wrt]])[[2]] <- substitute(
      outList <- derivs(run(arg1), wrt = this_wrt),
      list(this_wrt = wrt)
    )
  } else {
    body(methods[[method_no_wrt]])[[2]] <- quote(
      outList <- derivs(run(arg1, arg2))
    )
    body(methods[[method_all_wrt]])[[2]] <- substitute(
      outList <- derivs(run(arg1, arg2), wrt = this_wrt),
      list(this_wrt = wrt)
    )
  }

  list(
    name = opParam$name,
    opType = opType,
    opParam = opParam,
    run = run,
    methods = methods,
    enableDerivs = list('run'),
    wrts = wrts
  )
}

############
## unary ops
############

gammafn <- gamma
lgammafn <- lgamma
ceil <- ceiling

unaryArgs <- c('double(0)', 'double(1, 4)', 'double(2, c(3, 4))')
unaryOps <- c('-', 'sum',
              nimble:::unaryDoubleOperators,
              nimble:::unaryPromoteNoLogicalOperators,
              nimble:::reductionUnaryDoubleOperatorsEither,
              nimble:::reductionUnaryOperatorsArray)

unaryOpTests <- unlist(
  recursive = FALSE,
  x = lapply(
    unaryOps,
    function(x) {
      mapply(
        makeADtest,
        argTypes = unaryArgs,
        MoreArgs = list(op = x),
        SIMPLIFY = FALSE
      )
    })
)
names(unaryOpTests) <- sapply(unaryOpTests, `[[`, 'name')

## tranform args
modifyOnMatch(unaryOpTests, '(log|sqrt) .+', 'input_dist', function(x) abs(rnorm(x)))
modifyOnMatch(unaryOpTests, 'log1p .+', 'input_dist', function(x) abs(rnorm(x)) - 1)
modifyOnMatch(unaryOpTests, '(logit|probit|cloglog) .+', 'input_dist', runif)
modifyOnMatch(unaryOpTests, '(acos|asin|atanh) .+', 'input_dist', function(x) runif(x, -1, 1))
modifyOnMatch(unaryOpTests, 'acosh .+', 'input_dist', function(x) abs(rnorm(x)) + 1)

## set tolerances (default is tol1 = 0.00001 and tol2 = 0.0001)
modifyOnMatch(unaryOpTests, 'ilogit .+', 'tol2', 0.1)
modifyOnMatch(unaryOpTests, '(factorial|exp|log1p|tan|acos) .+', 'tol2', 0.001)

## runtime failures
## abs fails on negative inputs
modifyOnMatch(unaryOpTests, 'abs .+', 'input_dist', function(x) -abs(rnorm(x)))
modifyOnMatch(
  unaryOpTests,
  'abs (double\\(1, 4\\)|double\\(2, c\\(3, 4\\)\\))',
  'knownFailure', '.*runs'
)

## wrt failures
modifyOnMatch(
  unaryOpTests, '.+ double\\(0\\)', 'knownFailures',
  list(
    method4 = list(
      jacobian = expect_error,
      hessian = expect_error
    )
  )
)
modifyOnMatch(
  unaryOpTests, '.+ (double\\(1, 4\\)|double\\(2, c\\(3, 4\\)\\))',
  'knownFailures',
  list(
    method2 = list(
      jacobian = expect_error,
      hessian = expect_error
    ),
    method4 = list(
      jacobian = expect_error,
      hessian = expect_error
    )
  )
)

## compilation failures
modifyOnMatch(
  unaryOpTests,
  paste(
    '(logit|log1p|lfactorial|factorial|cloglog|atan|cosh|sinh|tanh|acosh|asinh|atanh)',
    '(double\\(1, 4\\)|double\\(2, c\\(3, 4\\)\\))'
  ),
  'knownFailures',
  list(compilation = TRUE)
)

modifyOnMatch(
  unaryOpTests, '(nimRound|ceil|floor|sd|var|probit|lgammafn|gammafn) .+',
  'knownFailures', list(compilation = TRUE)
)

ops_regex <- paste0(
  c('sum', nimble:::reductionUnaryDoubleOperatorsEither,
    nimble:::reductionUnaryOperatorsArray),
  collapse = '|'
)
modifyOnMatch(
  unaryOpTests,
  paste0('(', ops_regex, ') double\\(0\\)'),
  'knownFailures', list(compilation = TRUE)
)

set.seed(0)
sapply(unaryOpTests, test_AD)

#############
## binary ops
#############

binaryArgs <- as.list(
  cbind(
    data.frame(t(expand.grid(unaryArgs[1], unaryArgs)), stringsAsFactors=FALSE),
    data.frame(t(expand.grid(unaryArgs[2:3], unaryArgs[1])), stringsAsFactors=FALSE)
  )
)
names(binaryArgs) <- NULL
binaryArgs[[length(binaryArgs) + 1]] <- rep(unaryArgs[2], 2)
binaryArgs[[length(binaryArgs) + 1]] <- rep(unaryArgs[3], 2)

binaryOps <- c(nimble:::binaryOrUnaryOperators,
               nimble:::binaryMidOperators,
               nimble:::binaryLeftDoubleOperators,
               nimble:::reductionBinaryOperatorsEither)

binaryOpTests <- unlist(
  recursive = FALSE,
  x = lapply(
    binaryOps,
    function(x) {
      mapply(
        makeADtest,
        argTypes = binaryArgs,
        MoreArgs = list(op = x, opType = 'binary'),
        SIMPLIFY = FALSE
      )
    })
)
names(binaryOpTests) <- sapply(binaryOpTests, `[[`, 'name')

## set tolerances (default is tol1 = 0.00001 and tol2 = 0.0001)
modifyOnMatch(binaryOpTests, '(\\+|-) double\\(0\\) double\\(1, 4\\)', 'tol2', 0.001)
modifyOnMatch(binaryOpTests, '(\\+|-) double\\(1, 4\\) double\\(0\\)', 'tol2', 0.001)
modifyOnMatch(binaryOpTests, '(\\+|-) double\\(1, 4\\) double\\(1, 4\\)', 'tol2', 0.001)

## compilation failures
ops_regex <- paste0(
  c('\\+', '-'),
  collapse = '|'
)

modifyOnMatch(
  binaryOpTests,
  paste0('\\+ (',
         'double\\(2, c\\(3, 4\\)\\) double\\(0\\)|',
         'double\\(0\\) double\\(2, c\\(3, 4\\)\\))'),
  'knownFailure', '.*compiles'
)

modifyOnMatch(
  binaryOpTests,
  '/ double\\(0\\) double\\(2, c\\(3, 4\\)\\)',
  'knownFailure', '.*compiles'
)

sapply(binaryOpTests, test_AD, info = 'binary')

##########################
## unary square matrix ops
##########################

squareMatrixArgs <- c('double(2, c(2, 2))', 'double(2, c(5, 5))')

squareMatrixOps <- c(nimble:::matrixSquareReductionOperators,
                     'inverse')

squareMatrixOpTests <- unlist(
  recursive = FALSE,
  x = lapply(
    squareMatrixOps,
    function(x) {
      mapply(
        makeADtest,
        argTypes = squareMatrixArgs,
        MoreArgs = list(op = x),
        SIMPLIFY = FALSE
      )
    })
)
names(squareMatrixOpTests) <- sapply(squareMatrixOpTests, `[[`, 'name')

## compilation failures
modifyOnMatch(squareMatrixOpTests, '(logdet|det) .+', 'knownFailure', '.*compiles')

## runtime failures
modifyOnMatch(
  squareMatrixOpTests, 'inverse double\\(2, c\\(5, 5\\)\\)',
  'knownFailure', '.*compiles'
)

sapply(squareMatrixOpTests, test_AD, info = 'unary square matrix')

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

binaryMatrixOpTests <- unlist(
  recursive = FALSE,
  x = lapply(
    binaryMatrixOps,
    function(x) {
      mapply(
        makeADtest,
        argTypes = binaryMatrixArgs,
        MoreArgs = list(op = x, opType = 'binary'),
        SIMPLIFY = FALSE
      )
    })
)
names(binaryMatrixOpTests) <- sapply(binaryMatrixOpTests, `[[`, 'name')

## set tolerances (default is tol1 = 0.00001 and tol2 = 0.0001)
modifyOnMatch(binaryMatrixOpTests, '%*% .+', 'tol2', 0.001)

sapply(binaryMatrixOpTests, test_AD, info = 'binary matrix')

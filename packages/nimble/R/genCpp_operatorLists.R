## This takes a character vector as the first argument and length-1
## character vector as the second argument.  It returns a list with
## the first vector as names and the second argument as the value of
## each element.  E.g. makeCallList(c('A','B'), 'foo') is equivalent
## to list(A = 'foo', B = 'foo')
makeCallList <- function(opList, call) {
    ans <- rep(list(call), length(opList))
    names(ans) <- opList
    ans
}

ifOrWhile <- c('if','while')

## Following are a set of operators organized into categories that
## various processing steps use
binaryMidLogicalOperatorsLogical <- c('&','|')
binaryMidLogicalOperatorsComparison <- c('==','!=','<=','>=','<','>')
binaryMidLogicalOperators <- c(binaryMidLogicalOperatorsLogical,
                               binaryMidLogicalOperatorsComparison)

binaryMidDoubleOperators <- c('/', '^')
binaryMidPromoteNoLogicalOperators <- c('*','%%')
binaryMidOperators <- c(binaryMidDoubleOperators,
                        binaryMidPromoteNoLogicalOperators)

binaryLeftDoubleOperators <- c('pow','nimMod', 'pow_int')
binaryLeftPromoteOperators <- c('pmin','pmax','pairmin','pairmax')
binaryLeftLogicalOperators <- c( 'nimEquals')
binaryLeftOperators <- c(binaryLeftDoubleOperators,
                         binaryLeftPromoteOperators,
                         binaryLeftLogicalOperators)

binaryOperators <- c(binaryMidOperators,
                     binaryLeftOperators)

binaryOrUnaryOperators <- c('+','-')
unaryPromoteNoLogicalOperators <- c('abs','cube')
unaryIntegerOperators <- 'nimStep'
unaryLogicalOperators <- '!'
unaryDoubleOperators <- c('exp',
                          'log',
                          'logit',
                          'ilogit',
                          'probit',
                          'iprobit',
                          'sqrt', ## these do not go directly into cppOutputCalls.  They should be direct C++ names or go through eigProxyCalls or eigProxyCallsExternalUnary
                          'gammafn',
                          'lgammafn',                    ## these also do not go direclty into eigenizeCalls but rather should be entered directly there for eigenize_cWiseUnaryEither, eigenize_cWiseUnaryArray or eigenize_cWiseUnaryMatrix
                          ## 'lgamma1p',
                          'log1p',
                          'lfactorial',
                          'factorial',
                          'cloglog',
                          'icloglog',
                          'nimRound',
                          'ftrunc',
                          'ceil',
                          'floor', 
                          'cos',
                          'sin',
                          'tan',
                          'acos',
                          'asin',
                          'atan',
                          'cosh',
                          'sinh',
                          'tanh',
                          'acosh',
                          'asinh',
                          'atanh')
unaryOperators <- c(unaryPromoteNoLogicalOperators,
                    unaryIntegerOperators,
                    unaryDoubleOperators,
                    unaryLogicalOperators)
unaryOrNonaryOperators <- list() 
assignmentOperators <- c('<-','<<-','=')

reductionUnaryDoubleOperatorsEither <- c('mean', 'prod','squaredNorm')
reductionUnaryPromoteOperatorsEither <-  c('min','max', 'sum')
reductionUnaryLogicalOperatorsEither <- c('any','all')

reductionUnaryOperatorsEither <- c(reductionUnaryDoubleOperatorsEither,
                                   reductionUnaryPromoteOperatorsEither,
                                   reductionUnaryLogicalOperatorsEither)  # removed norm as not consistent between R and C

reductionUnaryOperatorsArray <- c('sd','var')
reductionUnaryOperators <- c(reductionUnaryOperatorsEither,
                             reductionUnaryOperatorsArray)
matrixSquareReductionOperators <- c('det','logdet')
reductionBinaryOperatorsEither <- c('inprod')
reductionBinaryOperators <- reductionBinaryOperatorsEither

coreRnonSeqBlockCalls <- c('nimNonseqIndexedd',
                           'nimNonseqIndexedi',
                           'nimNonseqIndexedb')
coreRmanipulationCalls <- c('nimC',
                            'nimRepd',
                            'nimRepi',
                            'nimRepb',
                            'nimSeqByD',
                            'nimSeqLenD',
                            'nimSeqByLenD',
                            'nimSeqByI',
                            'nimSeqLenI',
                            'nimSeqByLenI',
                            'nimDiagonalD',
                            'nimDiagonalI',
                            'nimDiagonalB',
                            'nimNewMatrixD',
                            'nimNewMatrixI',
                            'nimNewMatrixB')
nonNativeEigenCalls <- c('logdet',
                         'sd',
                         'var',
                         'inprod',
                         coreRmanipulationCalls,
                         coreRnonSeqBlockCalls)

matrixMultOperators <- c('%*%')
matrixFlipOperators <- c('t')
matrixSquareOperators <- c('chol','inverse')
nimbleListReturningOperators <- c('nimEigen',
                                  'nimSvd',
                                  'getDerivs_wrapper',
                                  'nimDerivs_dummy')  ## These use sizeNimbleListReturningFunction. Note that nimOptim and nimDerivs_calculate are handled separately.
matrixSolveOperators <- c('solve','forwardsolve','backsolve')
passThroughOperators <- c('return')

returnTypeCodes <- list(
    double = 1L,
    integer = 2L,
    logical = 3L,
    promote = 4L,
    promoteNoLogical = 5L)

returnTypeHandling <- with(
    returnTypeCodes,
    c(list('(' = promote),
      makeCallList(binaryMidLogicalOperators, logical),
      makeCallList(binaryMidDoubleOperators, double),
      makeCallList(binaryMidPromoteNoLogicalOperators, promoteNoLogical),
      makeCallList(binaryLeftDoubleOperators, double),
      makeCallList(binaryLeftPromoteOperators, promoteNoLogical),
      makeCallList(binaryLeftLogicalOperators, logical),
      makeCallList(binaryOrUnaryOperators, promoteNoLogical),
      makeCallList(unaryPromoteNoLogicalOperators, promoteNoLogical),
      makeCallList(unaryLogicalOperators, logical),
      makeCallList(unaryIntegerOperators, integer),
      makeCallList(unaryDoubleOperators, double),
      makeCallList(reductionUnaryDoubleOperatorsEither, double),
      makeCallList(reductionUnaryPromoteOperatorsEither, promoteNoLogical),
      makeCallList(reductionUnaryLogicalOperatorsEither, logical),
      makeCallList(reductionUnaryOperatorsArray, double),
      makeCallList(matrixSquareReductionOperators, double),
      makeCallList(reductionBinaryOperatorsEither, promoteNoLogical),
      makeCallList(c(matrixMultOperators, matrixSquareOperators, matrixSolveOperators), double)))
## deliberately omitted (so they just return same type as input):
## matrixFlipOperators ('t')


midOperators <- as.list(
    paste0(' ',
           c(binaryMidOperators,
             binaryMidLogicalOperators,
             binaryOrUnaryOperators,
             assignmentOperators),
           ' ')
)

names(midOperators) <- c(binaryMidOperators,
                         binaryMidLogicalOperators,
                         binaryOrUnaryOperators,
                         assignmentOperators)
midOperators <- c(midOperators,
                  list('$' = '$', '%*%' = ' %*% ', ':' = ':', '%o%' = '%o%'))

brackOperators <- list('[' = c('[',']'),
                       '[[' = c('[[',']]'))

## see distributions_processInputList for some relevant lists of distributions functions 

callToSkipInEigenization <- c('copy',
                              'setValues',
                              'setValuesIndexRange',
                              'getValues',
                              'getValuesIndexRange',
                              'setSize',
                              'resize',
                              'getsize',
                              'size',
                              'resizeNoPtr',
                              'assert',
                              'return',
                              'blank',
                              'rankSample',
                              'nimArr_dmnorm_chol',
                              'nimArr_dmvt_chol',
                              'nimArr_dlkj_corr_cholesky',
                              'nimArr_dwish_chol',
                              'nimArr_dinvwish_chol',
                              'nimArr_dcar_normal',
                              'nimArr_dcar_proper',
                              'nimArr_dmulti',
                              'nimArr_dcat',
                              'nimArr_dinterval',
                              'nimArr_ddirch',
                              'nimArr_rmnorm_chol',
                              'nimArr_rmvt_chol',
                              'nimArr_rlkj_corr_cholesky',
                              'nimArr_rwish_chol',
                              'nimArr_rinvwish_chol',
                              'nimArr_rcar_normal',
                              'nimArr_rcar_proper',
                              'nimArr_rmulti',
                              'nimArr_rcat',
                              'nimArr_rinterval',
                              'nimArr_rdirch',
                              'calculate',
                              'calculateDiff',
                              'simulate',
                              'getLogProb',
                              'nimEquals',
                              'startNimbleTimer',
                              'endNimbleTimer')

## used for cppOutputs
## things here should have the inverse listing in the eigenizeTranslate list
eigProxyTranslate <- c(eigTranspose = 'transpose',
                       eigCos = 'cos',
                       eigSin = 'sin',
                       eigTan = 'tan',
                       eigAcos = 'acos',
                       eigAsin = 'asin',
                       eigExp = 'exp',
                       eigPow = 'pow',
                       eigLog = 'log',
                       eigCube = 'cube',
                       eiginprod = 'eigenInprod', ## need lowercase i for nonNativeEigCalls files
                       cwiseProduct = 'cwiseProduct',
                       cwiseQuotient = 'cwiseQuotient',
                       eigArray = 'array',
                       eigMatrix = 'matrix',
                       eigInverse = 'inverse',
                       setAll = 'setConstant',
                       eigEval = 'eval',
                       eigDiagonal = 'diagonal',
                       eigenBlock = 'block') ## created in makeEigenBlockExprFromBrackets called from sizeIndexingBracket

newEPT <- reductionUnaryOperators
names(newEPT) <- paste0('eig', reductionUnaryOperators)
eigProxyTranslate <- c(eigProxyTranslate, newEPT)

eigProxyTranslate[['eigmax']] <- 'maxCoeff'
eigProxyTranslate[['eigmin']] <- 'minCoeff'
eigProxyTranslate[['eigpmax']] <- 'max'
eigProxyTranslate[['eigpmin']] <- 'min'


newEPT <- reductionBinaryOperators
names(newEPT) <- paste0('eig', reductionBinaryOperators)
eigProxyTranslate <- c(eigProxyTranslate,
                       newEPT)

newEPT <- c(coreRmanipulationCalls,
            coreRnonSeqBlockCalls)
names(newEPT) <- paste0('eig',
                        c(coreRmanipulationCalls,
                          coreRnonSeqBlockCalls)
                        )
eigProxyTranslate <- c(eigProxyTranslate, newEPT)

newEPT <- matrixSquareReductionOperators
names(newEPT) <- paste0('eig', matrixSquareReductionOperators)
eigProxyTranslate <- c(eigProxyTranslate, newEPT)
eigProxyTranslate[['eigdet']] <- 'determinant'

## nonNativeEigenProxyCalls should each appear in one of the other
## operatorLists feeding into eigProxyTranslate Then they are removed
## from eigProxyCalls so that in cppOutputCalls the
## nonNativeEigenProxyCalls can be generated differently (as a
## function, not a member function)
nonNativeEigenProxyCalls <- paste0('eig', nonNativeEigenCalls)
eigProxyCalls <- setdiff(names(eigProxyTranslate), nonNativeEigenProxyCalls)

## things here should have the inverse entry in the eigenizeTranslate
## list.  Those are created automatically in genCpp_eigenization.R The
## C++ name here is (unfortunately) also used supposed to match the
## DSL keyword to be included in eigenizeCalls correctly
eigProxyTranslateExternalUnary <- list(
    ## (C++ name, arg type, return type, DSL name if different from C++ name) for std::ptr_fun<argtype, returntype>(fun name)
    eigAtan = c('atan', 'double', 'double'), 
    eigCosh = c('cosh', 'double', 'double'),
    eigSinh = c('sinh', 'double', 'double'),
    eigTanh = c('tanh', 'double', 'double'),
    eigAcosh = c('acosh', 'double', 'double'),
    eigAsinh = c('asinh', 'double', 'double'),
    eigAtanh = c('atanh', 'double', 'double'),
    eigLogit = c('logit', 'double', 'double'),
    eigIlogit = c('ilogit', 'double', 'double'),
    eigProbit = c('probit', 'double', 'double'),
    eigIprobit = c('iprobit', 'double', 'double'),
    eigGammafn = c('gammafn', 'double', 'double'),
    eigLgammafn = c('lgammafn', 'double', 'double'),
    eigLog1p = c('log1p', 'double', 'double'),
    eigLfactorial = c('lfactorial', 'double', 'double'),
    eigFactorial = c('factorial', 'double', 'double'),
    eigCloglog = c('cloglog', 'double', 'double'),
    eigIcloglog = c('icloglog', 'double', 'double'),
    eigNimRound = c('nimRound', 'double', 'double'),
    eigFtrunc = c('ftrunc', 'double', 'double'),
    eigCeil = c('ceil', 'double', 'double'),
    eigFloor = c('floor', 'double', 'double'),
    eigNimStep = c('nimStep', 'double', 'int'),
    'eig!' = c('nimNot','bool','bool', '!')
)
eigProxyCallsExternalUnary <- names(eigProxyTranslateExternalUnary)

eigCalls <- c('matrix', 'array')
cppCasts = list(as.numeric = 'double', as.integer = 'int')

##http://en.cppreference.com/w/cpp/language/operator_precedence

## Used to decide when to put parentheses around LHS or RHS based on operator precendence.
operatorRank <- c(
    list('<-' = 100, '^' = 4, '::' = 3),
    makeCallList(c('*','/','%*%', '%%'), 5),
    makeCallList(c('+', '-'), 6),
    makeCallList(c('>','<','<=', '>='), 7),
    makeCallList(c('==','!='), 8),
    list('!' = 10), ## follows R's precedence order, not C's
    list('&' = 13,
         '|' = 14,
         '&&' = 13,
         '||' = 14)
)

nimDerivsPrependTypeOperators <- c("dnorm", "dpois", "dgamma", "dinvgamma", "dsqrtinvgamma",
                                   "dexp_nimble", "dexp", "ddexp", "dlnorm", "dweibull",
                                   "dbinom", "dbeta", "dchisq", "dlogis", "dt",
                                   "dt_nonstandard", "nimArr_dmulti", "nimArr_dcat", "dnbinom", "dunif", "pairmax", "pairmin", 
                                   "nimArr_ddirch", "nimArr_dmvt_chol", "nimArr_dmnorm_chol", 
                                   "nimArr_dwish_chol", "nimArr_dinvwish_chol", "nimArr_dlkj_corr_cholesky",
                                   "nimArr_dcar_normal", "nimArr_dcar_proper",
                                   "nimStep", 'ilogit', 'icloglog', 'iprobit', 'probit', 'cloglog',
                                   "nimEquals", "lgammafn", "gammafn", "lfactorial", "factorial",
                                   "logit", "floor", "ceil", "nimRound", "ftrunc",
                                   "cube", "inprod",
                                   "nimStep", "dflat", "dhalfflat")

## Reflects distribution funs that support recycling rule with AD -- see nimDerivs_TMB.h.
recyclingRuleOperatorsAD <- c(
  'dbinom', 'dexp_nimble', 'dnbinom', 'dpois', 'dchisq', 'dbeta', 'dnorm', 'dgamma',
  'dinvgamma', 'ddexp', 'dlnorm', 'dlogis', 'dunif', 'dweibull',
  ## The following operators should work but are currently broken for various
  ## reasons which are explained in test-ADfunctions.R.
  'dexp', 'dsqrtinvgamma', 'dt_nonstandard', 'dt', 'pow_int'
)
recyclingRuleOperatorsAD <- paste0(
    recyclingRuleOperatorsAD,
    '_RR_impl<MatrixXd>::',
    recyclingRuleOperatorsAD,
    '_RecyclingRule'
)

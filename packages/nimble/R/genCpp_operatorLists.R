ifOrWhile <- c('if','while')

## Following are a set of operators organized into categories that various processing steps use
binaryMidLogicalOperators <- c('==','!=','<=','>=','<','>','&','|')
binaryMidOperators <- c('/','*','%%','^')
binaryLeftOperators <- c('pow','pmin','pmax', 'nimMod', 'nimEquals','pairmin','pairmax')
binaryOperators <- c(binaryMidOperators, binaryLeftOperators)
binaryOrUnaryOperators <- c('+','-')
unaryOperators <- c('exp','log', 'cube', 'logit','ilogit','probit','iprobit', 'sqrt', ## these do not go directly into cppOutputCalls.  They should be direct C++ names or go through eigProxyCalls or eigProxyCallsExternalUnary
                    'gammafn','lgammafn',                    ## these also do not go direclty into eigenizeCalls but rather should be entered directly there for eigenize_cWiseUnaryEither, eigenize_cWiseUnaryArray or eigenize_cWiseUnaryMatrix
                    'lgamma1p', 'log1p', 'lfactorial', 'factorial', 'cloglog', 'icloglog',
                    'abs','nimRound','ftrunc','ceil','floor','nimStep', 
                    'cos', 'sin', 'tan', 'acos', 'asin', 'atan', 'cosh', 'sinh', 'tanh', 'acosh', 'asinh', 'atanh')
unaryOrNonaryOperators <- list() 
assignmentOperators <- c('<-','<<-','=')
reductionUnaryOperatorsEither <- c('min','max','sum','mean','any','all','prod','squaredNorm')  # removed norm as not consistent between R and C
reductionUnaryOperatorsArray <- c('sd','var')
reductionUnaryOperators <- c(reductionUnaryOperatorsEither, reductionUnaryOperatorsArray)
matrixSquareReductionOperators <- c('det','logdet','trace')
reductionBinaryOperatorsEither <- c('inprod')
reductionBinaryOperators <- reductionBinaryOperatorsEither

nonNativeEigenCalls <- c('logdet','sd','var','inprod')

matrixMultOperators <- c('%*%')
matrixFlipOperators <- c('t')
matrixSquareOperators <- c('chol','inverse')
matrixEigenListOperators <- c('EIGEN_EIGEN')
matrixSolveOperators <- c('solve','forwardsolve','backsolve')
matrixVectorOperators <- c('svdd')
matrixMatrixOperators <- c('svdu', 'svdv')
passThroughOperators <- c('return')
##keywordOperators <- c('for','if', 'while')

midOperators <- as.list(paste0(' ',c(binaryMidOperators,  binaryMidLogicalOperators, binaryOrUnaryOperators, assignmentOperators),' '))
names(midOperators) <- c(binaryMidOperators, binaryMidLogicalOperators, binaryOrUnaryOperators, assignmentOperators)
midOperators <- c(midOperators, list('$' = '$', '%*%' = ' %*% ', ':' = ':', '%o%' = '%o%'))

brackOperators <- list('[' = c('[',']'), '[[' = c('[[',']]'))

callToSkipInEigenization <- c('copy','setValues', 'getValues', 'setSize', 'resize', 'getsize', 'size', 'resizeNoPtr','assert', 'return', 'blank', 'rankSample', 'nimArr_dmnorm_chol', 'nimArr_dmvt_chol', 'nimArr_dwish_chol', 'nimArr_dmulti', 'nimArr_dcat', 'nimArr_dinterval', 'nimArr_ddirch', 'nimArr_rmnorm_chol', 'nimArr_rmvt_chol', 'nimArr_rwish_chol', 'nimArr_rmulti', 'nimArr_rcat', 'nimArr_rinterval', 'nimArr_rdirch', 'calculate', 'calculateDiff', 'simulate', 'getLogProb', 'nimEquals', 'startNimbleTimer', 'endNimbleTimer')

## This takes a character vector as the first argument and length-1 character vector as the second argument.
## It returns a list with the first vector as names and the second argument as the value of each element.
## E.g. makeCallList(c('A','B'), 'foo') is equivalent to list(A = 'foo', B = 'foo')
makeCallList <- function(opList, call) {
    ans <- rep(list(call), length(opList))
    names(ans) <- opList
    ans
}

## used for nimDeparse &/or cppOutputs

## used for cppOutputs
## eigProxyCalls <- c('eigTranspose', 'eigCos', 'eigSin', 'eigTan', 'eigAcos', 'eigAsin', 'eigExp', 'eigLog', 'eigCube', 'cwiseProduct', 'cwiseQuotient', 'eigArray', 'eigMatrix', 'eigInverse', 'setAll', 'eigEval')
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
                       eigEval = 'eval')

newEPT <- reductionUnaryOperators
names(newEPT) <- paste0('eig', reductionUnaryOperators)
eigProxyTranslate <- c(eigProxyTranslate, newEPT)

eigProxyTranslate[['eigmax']] <- 'maxCoeff'
eigProxyTranslate[['eigmin']] <- 'minCoeff'
eigProxyTranslate[['eigpmax']] <- 'max'
eigProxyTranslate[['eigpmin']] <- 'min'


newEPT <- reductionBinaryOperators
names(newEPT) <- paste0('eig', reductionBinaryOperators)
eigProxyTranslate <- c(eigProxyTranslate, newEPT)

newEPT <- matrixSquareReductionOperators
names(newEPT) <- paste0('eig', matrixSquareReductionOperators)
eigProxyTranslate <- c(eigProxyTranslate, newEPT)
eigProxyTranslate[['eigdet']] <- 'determinant'

nonNativeEigenProxyCalls <- paste0('eig', nonNativeEigenCalls)
eigProxyCalls <- setdiff(names(eigProxyTranslate), nonNativeEigenProxyCalls)

## things here shuold have the inverse listing in the eigenizeTranslate list
eigProxyTranslateExternalUnary <- list(eigAtan = c('atan', 'double', 'double'), ## (C++ name, arg type, return type) for std::ptr_fun<argtype, returntype>(fun name)
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
                                       eigLgamma1p = c('lgamma1p', 'double', 'double'),
                                       eigLog1p = c('log1p', 'double', 'double'),
                                       eigLfactorial = c('lfactorial', 'double', 'double'),
                                       eigFactorial = c('factorial', 'double', 'double'),
                                       eigCloglog = c('cloglog', 'double', 'double'),
                                       eigIcloglog = c('icloglog', 'double', 'double'),
                                       eigNimRound = c('nimRound', 'double', 'double'),
                                       eigFtrunc = c('ftrunc', 'double', 'double'),
                                       eigCeil = c('ceil', 'double', 'double'),
                                       eigFloor = c('floor', 'double', 'double'),
                                       eigNimStep = c('nimStep', 'double', 'int')
                                       )
eigProxyCallsExternalUnary <- names(eigProxyTranslateExternalUnary)

eigOtherMemberFunctionCalls <- c('cwiseSqrt', 'cwiseAbs')
eigCalls <- c('matrix', 'array')
cppCasts = list(as.numeric = 'double', as.integer = 'int')

##http://en.cppreference.com/w/cpp/language/operator_precedence

## Used to decide when to put parentheses around LHS or RHS based on operator precendence.
operatorRank <- c(
    list('<-' = 100, '^' = 4),
    makeCallList(c('*','/','%*%', '%%'), 5),
    makeCallList(c('+', '-'), 6),
    makeCallList(c('>','<','<=', '>='), 7),
    makeCallList(c('==','!='), 8),
    list('&' = 13,
         '|' = 14,
         '&&' = 13,
         '||' = 14)                  
)



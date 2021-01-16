source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
context("Testing of numeric type handling and casting")

RwarnLevel <- options('warn')$warn
## There are a bunch of NaN warnings we want to ignore.
options(warn = -1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)


inverseCallReplacements <- as.list(names(nimble:::specificCallReplacements))
names(inverseCallReplacements) <- unlist(nimble:::specificCallReplacements)
inverseReplace <- function(x) {
    replacement <- inverseCallReplacements[[x]]
    if(is.null(replacement)) x else replacement
}

makeUnaryCwiseTypeTest <- function(name, funName, type, nDim) {
    outputHandling <- nimble:::returnTypeHandling[[funName]]
    
    if(is.null(outputHandling)) outputType <- type
    else outputType <- nimble:::setReturnType(funName, type) ## see sizeProcessing and returnTypeCodes

    seqFrom <- if(type == 'double') 0.1 else 1L
    seqBy <- if(type == 'double') 0.1 else 1L
    
    list(name = paste(name, type, nDim),
         expr = substitute(out <- FOO(arg1), list(FOO = as.name(inverseReplace(funName)))),
         args = list(arg1 = substitute(TYPE(NDIM), list(TYPE = as.name(type), NDIM = nDim))),
         setArgVals = substitute( {arg1 <- as(seq(from = SEQFROM, by = SEQBY, length.out = NDIM + 1), TYPE);
             if(NDIM == 2) arg1 <- matrix(arg1, nrow = 1)}, list(NDIM = nDim, TYPE = type, SEQFROM = seqFrom, SEQBY = seqBy)),
         outputType = substitute(OUTPUTTYPE(NDIM), list(OUTPUTTYPE = as.name(outputType), NDIM = nDim)))
}


unaryCwiseTypeTests <- unlist(recursive = FALSE,
                              x= lapply(c('-', nimble:::unaryOperators),
                                        function(x) {
                                            mapply(makeUnaryCwiseTypeTest,
                                                   type = rep(c('double','integer','logical'), 3),
                                                   nDim = rep(0:2, each = 3),
                                                   MoreArgs = list(name = x, funName = x),
                                                   SIMPLIFY = FALSE) 
                                        }
                                        )
                              )

unaryCwiseTypeTests <- indexNames(unaryCwiseTypeTests)

makeBinaryCwiseTypeTest <- function(name, funName, LHStype, RHStype, nDim, outputNDim = nDim) {
    outputHandling <- nimble:::returnTypeHandling[[funName]]
    
    if(is.null(outputHandling)) outputType <- LHStype
    else outputType <- nimble:::setReturnType(funName, nimble:::arithmeticOutputType(LHStype, RHStype)) ## see sizeProcessing and returnTypeCodes

    LHSseqFrom <- if(LHStype == 'double') 0.1 else 1L
    LHSseqBy <- if(LHStype == 'double') 0.1 else 1L

    RHSseqFrom <- if(RHStype == 'double') 0.1 else 1L
    RHSseqBy <- if(RHStype == 'double') 0.1 else 1L
    
    lenOut <- if(nDim == 2) 4 else nDim + 1
    
    list(name = paste(name, LHStype, RHStype, nDim),
         expr = substitute(out <- FOO(arg1, arg2), list(FOO = as.name(inverseReplace(funName)))),
         args = list(arg1 = substitute(TYPE(NDIM), list(TYPE = as.name(LHStype), NDIM = nDim)),
                     arg2 = substitute(TYPE(NDIM), list(TYPE = as.name(RHStype), NDIM = nDim))),
         setArgVals = substitute( {
             arg1 <- as(seq(from = LHSSEQFROM, by = LHSSEQBY, length.out = LENOUT), LHSTYPE);
             if(NDIM == 2) arg1 <- matrix(arg1, nrow = 2);
             arg2 <- as(seq(from = RHSSEQFROM, by = RHSSEQBY, length.out = LENOUT), RHSTYPE);
             if(NDIM == 2) arg2 <- matrix(arg2, nrow = 2);
         }, list(LENOUT = lenOut, LHSTYPE = LHStype, RHSTYPE = RHStype, LHSSEQFROM = LHSseqFrom, LHSSEQBY = LHSseqBy, RHSSEQFROM = RHSseqFrom, RHSSEQBY = RHSseqBy, NDIM = nDim)),
         outputType = substitute(OUTPUTTYPE(NDIM), list(OUTPUTTYPE = as.name(outputType), NDIM = outputNDim)))
}


binaryCwiseTypeTests <- unlist(recursive = FALSE,
                              x = lapply(nimble:::binaryMidLogicalOperatorsComparison,
                                         function(x) {
                                             mapply(makeBinaryCwiseTypeTest,
                                                    LHStype = rep(c('double','integer','logical'), 9),
                                                    RHStype = rep(rep(c('double','integer','logical'), each = 3), 3),
                                                    nDim = rep(0:2, each = 9),
                                                    MoreArgs = list(name = x, funName = x),
                                                    SIMPLIFY = FALSE) 
                                         }
                                         )
                              )

binaryCwiseTypeTests <- indexNames(binaryCwiseTypeTests)


binaryCwiseTypeTestsLogicals <- unlist(recursive = FALSE,
                                       x= lapply(nimble:::binaryMidLogicalOperatorsLogical,
                                                 function(x) {
                                                     mapply(makeBinaryCwiseTypeTest,
                                                            LHStype = rep('logical', 3),
                                                            RHStype = rep('logical', 3),
                                                            nDim = 0:2,
                                                            MoreArgs = list(name = x, funName = x),
                                                            SIMPLIFY = FALSE) 
                                                 }
                                                 )
                                       )

binaryCwiseTypeTestsLogicals <- indexNames(binaryCwiseTypeTestsLogicals)

makeReductionTypeTest <- function(name, funName, type, nDim, checkEqual = FALSE) {
    outputHandling <- nimble:::returnTypeHandling[[funName]]
    
    if(is.null(outputHandling)) outputType <- type
    else outputType <- nimble:::setReturnType(funName, type) ## see sizeProcessing and returnTypeCodes

    seqFrom <- if(type == 'double') 0.1 else 1L
    seqBy <- if(type == 'double') 0.1 else 1L

    lenOut <- if(nDim == 2) 4 else nDim + 1

    list(name = paste(name, type, nDim),
         expr = substitute(out <- FOO(arg1), list(FOO = as.name(inverseReplace(funName)))),
         args = list(arg1 = substitute(TYPE(NDIM), list(TYPE = as.name(type), NDIM = nDim))),
         setArgVals = substitute( {arg1 <- as(seq(from = SEQFROM, by = SEQBY, length.out = LENOUT), TYPE);
                                   if(NDIM == 2) arg1 <- matrix(arg1, nrow = 2)}, list(NDIM = nDim, TYPE = type, SEQFROM = seqFrom, SEQBY = seqBy, LENOUT = lenOut)),
         outputType = substitute(OUTPUTTYPE(0), list(OUTPUTTYPE = as.name(outputType))),
         checkEqual = checkEqual)
}

reductionTypeTests <- unlist(recursive = FALSE,
                              x= lapply(nimble:::reductionUnaryOperators[ !(nimble:::reductionUnaryOperators %in% c('any','all', 'squaredNorm')) ],
                                        function(x) {
                                            mapply(makeReductionTypeTest,
                                                   type = rep(c('double','integer','logical'), if(x == 'var') 1 else 2),
                                                   nDim = if(x == 'var') 1 else rep(1:2, each = 3), ## We error-trap var(matrix) because we don't have cov() yet, which is what it should do to match R
                                                   MoreArgs = list(name = x, funName = x, checkEqual = grepl('prod',x)), ## prod does not produce numerically identical results between R and Eigen
                                                   SIMPLIFY = FALSE)
                                        }
                                        )
                              )

reductionTypeTests <- indexNames(reductionTypeTests)



reductionTypeTestsLogical <- unlist(recursive = FALSE,
                                    x= lapply(c('any','all'),
                                              function(x) {
                                                  mapply(makeReductionTypeTest,
                                                         type = rep('logical', 2),
                                                         nDim = rep(1:2),
                                                         MoreArgs = list(name = x, funName = x),
                                                         SIMPLIFY = FALSE) 
                                              }
                                              )
                                    )

reductionTypeTestsLogical <- indexNames(reductionTypeTestsLogical)

reductionTypeTestsMatrixSquare <- unlist(recursive = FALSE,
                                    x= lapply(nimble:::matrixSquareReductionOperators,
                                              function(x) {
                                                  mapply(makeReductionTypeTest,
                                                         ## As of 2015-03-13, Eigen does not allow integer determinants.
                                                         ## https://bitbucket.org/eigen/eigen/commits/678c42a8
                                                         type = 'double',
                                                         nDim = rep(2, 2),
                                                         MoreArgs = list(name = x, funName = x),
                                                         SIMPLIFY = FALSE) 
                                              }
                                              )
                                    )

reductionTypeTestsMatrixSquare <- indexNames(reductionTypeTestsMatrixSquare)


binaryCwiseTypeTestsMidOps <- unlist(recursive = FALSE,
                              x= lapply(c('+','-','/','*'),
                                        function(x) {
                                            mapply(makeBinaryCwiseTypeTest,
                                                   LHStype = rep(c('double','integer','logical'), 9),
                                                   RHStype = rep(rep(c('double','integer','logical'), each = 3), 3),
                                                   nDim = rep(0:2, each = 9),
                                                   MoreArgs = list(name = x, funName = x),
                                                   SIMPLIFY = FALSE) 
                                        }
                                        )
                              )

binaryCwiseTypeTestsMidOps <- indexNames(binaryCwiseTypeTestsMidOps)

binaryCwiseTypeTestsInprod <- unlist(recursive = FALSE,
                              x= lapply('inprod',
                                        function(x) {
                                            mapply(makeBinaryCwiseTypeTest,
                                                   LHStype = rep(c('double','integer','logical'), 3),
                                                   RHStype = rep(rep(c('double','integer','logical'), each = 3), 1),
                                                   nDim = rep(1, 9),
                                                   MoreArgs = list(name = x, funName = x, outputNDim = 0),
                                                   SIMPLIFY = FALSE) 
                                        }
                                        )
                              )

binaryCwiseTypeTestsInprod <- indexNames(binaryCwiseTypeTestsInprod)

binaryCwiseTypeTestsLeftPromotOps <- unlist(recursive = FALSE,
                              x= lapply(c('pmin','pmax'), ## not for scalars
                                        function(x) {
                                            mapply(makeBinaryCwiseTypeTest,
                                                   LHStype = rep(c('double','integer','logical'), 6),
                                                   RHStype = rep(rep(c('double','integer','logical'), each = 3), 2),
                                                   nDim = rep(1:2, each = 9),
                                                   MoreArgs = list(name = x, funName = x),
                                                   SIMPLIFY = FALSE) 
                                        }
                                        )
                              )

binaryCwiseTypeTestsLeftPromotOps <- indexNames(binaryCwiseTypeTestsLeftPromotOps)


makeBinaryMatrixOpTypeTest <- function(name, funName, LHStype, LHSdims, RHStype, RHSdims, outputNDim, checkEqual = FALSE) {
    outputHandling <- nimble:::returnTypeHandling[[funName]]
    
    if(is.null(outputHandling)) outputType <- LHStype
    else outputType <- nimble:::setReturnType(funName, nimble:::arithmeticOutputType(LHStype, RHStype)) ## see sizeProcessing and returnTypeCodes

    LHSseqFrom <- if(LHStype == 'double') 0.1 else 1L
    LHSseqBy <- if(LHStype == 'double') 0.1 else 1L
    LHSnDim <- length(LHSdims)
    LHSlenOut <- prod(LHSdims)    

    RHSseqFrom <- if(RHStype == 'double') 0.1 else 1L
    RHSseqBy <- if(RHStype == 'double') 0.1 else 1L
    RHSnDim <- length(RHSdims)
    RHSlenOut <- prod(RHSdims)
    
    list(name = paste(name, LHStype, RHStype, LHSnDim, RHSnDim),
         expr = substitute(out <- FOO(arg1, arg2), list(FOO = as.name(inverseReplace(funName)))),
         args = list(arg1 = substitute(TYPE(NDIM), list(TYPE = as.name(LHStype), NDIM = LHSnDim)),
                     arg2 = substitute(TYPE(NDIM), list(TYPE = as.name(RHStype), NDIM = RHSnDim))),
         setArgVals = substitute( {
             arg1 <- as(seq(from = LHSSEQFROM, by = LHSSEQBY, length.out = LHSLENOUT), LHSTYPE);
             if(LHSNDIM == 2) dim(arg1) <- LHSDIMS;
             arg2 <- as(seq(from = RHSSEQFROM, by = RHSSEQBY, length.out = RHSLENOUT), RHSTYPE);
             if(RHSNDIM == 2) dim(arg2) <- RHSDIMS
         }, list(LHSDIMS = LHSdims, RHSDIMS = RHSdims,
                 RHSNDIM = RHSnDim, LHSNDIM = LHSnDim,
                 RHSLENOUT = RHSlenOut, LHSLENOUT = LHSlenOut,
                 LHSTYPE = LHStype, RHSTYPE = RHStype,
                 LHSSEQFROM = LHSseqFrom, LHSSEQBY = LHSseqBy,
                 RHSSEQFROM = RHSseqFrom, RHSSEQBY = RHSseqBy)),
         outputType = substitute(OUTPUTTYPE(NDIM), list(OUTPUTTYPE = as.name(outputType), NDIM = outputNDim)),
         checkEqual = checkEqual)
}


## pmin, pmax
binaryMatrixOpTypeTests <- unlist(recursive = FALSE,
                                x= lapply(c('%*%'),
                                          function(x) {
                                              mapply(makeBinaryMatrixOpTypeTest, ## 3 types * 3 types * 4 shapes  
                                                     LHStype = rep(rep(c('double', 'integer', 'logical'), each = 4), 3),
                                                     RHStype = rep(c('double', 'integer', 'logical'), each = 12),
                                                     LHSdims = rep(list(c(3, 2), c(3, 2), 3      , 3), 9),
                                                     RHSdims = rep(list(c(2, 2), 2      , c(3, 1), 3), 9),
                                                     MoreArgs = list(name = x, funName = x, outputNDim = 2),
                                                     SIMPLIFY = FALSE, checkEqual = TRUE) 
                                          }
                                          )
                              )

binaryMatrixOpTypeTests <- indexNames(binaryMatrixOpTypeTests)


unaryCwiseResults <- test_coreRfeature_batch(unaryCwiseTypeTests, 'unaryCwiseTypeTests') ## lapply(unaryCwiseTypeTests, test_coreRfeature)
binaryCwiseResults <- test_coreRfeature_batch(binaryCwiseTypeTests, 'binaryCwiseTypeTests') ## lapply(binaryCwiseTypeTests, test_coreRfeature)
binaryCwiseLogicalResults <- test_coreRfeature_batch(binaryCwiseTypeTestsLogicals, 'binaryCwiseTypeTestsLogicals') ## lapply(binaryCwiseTypeTestsLogicals, test_coreRfeature)
reductionResults <- test_coreRfeature_batch(reductionTypeTests, 'reductionTypeTests') ## lapply(reductionTypeTests, test_coreRfeature)
reductionLogicalResults <- test_coreRfeature_batch(reductionTypeTestsLogical, 'reductionTypeTestsLogical') ## lapply(reductionTypeTestsLogical, test_coreRfeature)
reductionMatrixSquareResults <- test_coreRfeature_batch(reductionTypeTestsMatrixSquare[3:4], 'reductionTypeTestsMatrixSquare[3:4]') ## lapply(reductionTypeTestsMatrixSquare[3:4], test_coreRfeature)
binaryCwiseMidOpsResults <- test_coreRfeature_batch(binaryCwiseTypeTestsMidOps, 'binaryCwiseTypeTestsMidOps') ## lapply(binaryCwiseTypeTestsMidOps, test_coreRfeature)
binaryCwiseInProdResults <- test_coreRfeature_batch(binaryCwiseTypeTestsInprod, 'binaryCwiseTypeTestsInprod') ## lapply(binaryCwiseTypeTestsInprod, test_coreRfeature)
binaryCwiseLeftPromoteResults <- test_coreRfeature_batch(binaryCwiseTypeTestsLeftPromotOps, 'binaryCwiseTypeTestsLeftPromotOps') ## lapply(binaryCwiseTypeTestsLeftPromotOps, test_coreRfeature)
binaryMatrixOpResults <- test_coreRfeature_batch(binaryMatrixOpTypeTests, 'binaryMatrixOpTypeTests') ## lapply(binaryMatrixOpTypeTests, test_coreRfeature)

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)

## tests of numeric, integer, logical, matrix and array

source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context('Testing of numeric, integer, logical, matrix and array allocation')

numericTests <- list(
    list(name = 'numeric: length',
         expr = quote(out <- numeric(5)),
         outputType = quote(double(1) ),
         args = list(),
         setArgVals = quote({})
         ),
    list(name = 'numeric: length 0',
         expr = quote({
             out <- nimNumeric(length = 0)
         }),
         outputType = quote(double(1))
         ),
    list(name = 'numeric: length, scalar value',
         expr = quote(out <- nimNumeric(5, value = 1.5)),
         outputType = quote(double(1))
         ),
    list(name = 'numeric: length, value scalar, recycle = FALSE, fillZeros = TRUE',
         expr = quote(out <- nimNumeric(5, value = 2, recycle = FALSE, fillZeros = TRUE)),
         outputType = quote(double(1) ),
         args = list(),
         setArgVals = quote({})
         ),
    list(name = 'numeric: length, scalar value, init=FALSE',
         expr = quote({out <- nimNumeric(5, value = 1.5, init = FALSE);
             out[1:5] <- .123}),
         outputType = quote(double(1))
         ),
    list(name = 'numeric: length, init = FALSE', ## not really a good test b/c we need to insert values to make them match. So the test ensures it compiles but not its initialization behavior
         expr = quote({out <- nimNumeric(5, init = FALSE);
             out[1:5] <- .123}),
         outputType = quote(double(1))
         ),
    list(name = 'numeric: length expression, scalar value expression, init expression',
         expr = quote({
             a <- 2:4
             b <- rnorm(5)
             c <- 1:10
             out <- nimNumeric(sum(a) + 1, value = mean(b) + 0.5, init = (c > 2)[4]);
         }),
         outputType = quote(double(1))
         ),
    list(name = 'numeric: length, vector value',
         expr = quote({
             v <- rnorm(5)
             out <- nimNumeric(5, value = v);
         }),
         outputType = quote(double(1))
         ),
    list(name = 'numeric: length, vector value expression',
         expr = quote({
             v <- rnorm(5)
             out <- nimNumeric(5, value = exp(v) + 1);
         }),
         outputType = quote(double(1))
         ),
    list(name = 'numeric: length, vector value with too few elements',
         expr = quote({
             v <- rnorm(4)
             out <- nimNumeric(5, value = v);
         }),
         outputType = quote(double(1))
         ),
    list(name = 'numeric: length vector value with different types',
         expr = quote({
             v <- nimInteger(5, value = 1:5)
             out <- nimNumeric(5, value = v);
         }),
         outputType = quote(double(1))
         ),
    list(name = 'numeric: length vector value with too many elements',
         expr = quote({
             v <- rnorm(6)
             out <- nimNumeric(5, value = v);
         }),
         outputType = quote(double(1))
         ),
    list(name = 'numeric: length,  matrix value',
         expr = quote({
             v <- matrix(rnorm(6), nrow = 2)
             out <- nimNumeric(6, value = v);
         }),
         outputType = quote(double(1))
         )
)

integerTests <- list(
    list(name = 'integer: length',
         expr = quote(out <- integer(5)),
         outputType = quote(integer(1) ),
         args = list(),
         setArgVals = quote({})
         ),
    list(name = 'integer: length, scalar value',
         expr = quote(out <- nimInteger(5, value = 1.5)),
         outputType = quote(integer(1))
         ),
    list(name = 'integer: length, value scalar, recycle = FALSE, fillZeros = TRUE',
         expr = quote(out <- nimInteger(5, value = 2, recycle = FALSE, fillZeros = TRUE)),
         outputType = quote(integer(1) ),
         args = list(),
         setArgVals = quote({})
         ),
    list(name = 'integer: length, scalar value, init=FALSE',
         expr = quote({out <- nimInteger(5, value = 1.5, init = FALSE);
             out[1:5] <- 123}),
         outputType = quote(integer(1)),
         checkEqual = TRUE
         ),
    list(name = 'integer: length, init = FALSE', ## not really a good test b/c we need to insert values to make them match. So the test ensures it compiles but not its initialization behavior
         expr = quote({out <- nimInteger(5, init = FALSE);
             out[1:5] <- 123}),
         outputType = quote(integer(1)),
         checkEqual = TRUE
         ),
    list(name = 'integer: length expression, scalar value expression, init expression',
         expr = quote({
             a <- 2:4
             b <- rnorm(5)
             c <- 1:10
             out <- nimInteger(sum(a) + 1, value = round(mean(b)), init = (c > 2)[4]);
         }),
         outputType = quote(integer(1))
         ),
    list(name = 'integer: length, vector value',
         expr = quote({
             v <- rpois(5, lambda = 12.3)
             out <- nimInteger(5, value = v);
         }),
         outputType = quote(integer(1))
         ),
    list(name = 'integer: length, vector value expression',
         expr = quote({
             v <- rpois(5, lambda = 12.3)
             out <- nimInteger(5, value = exp(v) + 1);
         }),
         outputType = quote(integer(1))
         ),
    list(name = 'integer: length, vector value with too few elements',
         expr = quote({
             v <- rpois(4, lambda = 12.3)
             out <- nimInteger(5, value = v);
         }),
         outputType = quote(integer(1))
         ),
    list(name = 'integer: length vector value with different types',
         expr = quote({
             v <- nimInteger(5, value = 1:5)
             out <- nimInteger(5, value = v);
         }),
         outputType = quote(integer(1))
         ),
    list(name = 'integer: length vector value with too many elements',
         expr = quote({
             v <- rpois(6, lambda = 12.3)
             out <- nimInteger(5, value = v);
         }),
         outputType = quote(integer(1))
         ),
    list(name = 'integer: length,  matrix value',
         expr = quote({
             v <- matrix(rpois(6, lambda = 12.3), nrow = 2)
             out <- nimInteger(6, value = v);
         }),
         outputType = quote(integer(1))
         )
)

logicalTests <- list(
    list(name = 'logical: length',
         expr = quote(out <- logical(5)),
         outputType = quote(logical(1) ),
         args = list(),
         setArgVals = quote({})
         ),
    list(name = 'logical: length, scalar value',
         expr = quote(out <- nimLogical(5, value = 1.5)),
         outputType = quote(logical(1))
         ),
    list(name = 'logical: length, value scalar, recycle = FALSE, fillZeros = TRUE',
         expr = quote(out <- nimLogical(5, value = TRUE, recycle = FALSE, fillZeros = TRUE)),
         outputType = quote(logical(1) ),
         args = list(),
         setArgVals = quote({})
         ),
    list(name = 'logical: length, scalar value, init=FALSE',
         expr = quote({out <- nimLogical(5, value = FALSE, init = FALSE);
             out[1:5] <- TRUE}),
         outputType = quote(logical(1))
         ),
    list(name = 'logical: length, init = FALSE', ## not really a good test b/c we need to insert values to make them match. So the test ensures it compiles but not its initialization behavior
         expr = quote({out <- nimLogical(5, init = FALSE);
             out[1:5] <- TRUE}),
         outputType = quote(logical(1))
         ),
    list(name = 'logical: length expression, scalar value expression, init expression',
         expr = quote({
             a <- 2:4
             b <- rnorm(5)
             c <- 1:10
             out <- nimLogical(sum(a) + 1, value = mean(b) > 0, init = (c > 2)[4]);
         }),
         outputType = quote(logical(1))
         ),
    list(name = 'logical: length, vector value',
         expr = quote({
             v <- rpois(5, lambda = 12.3) > 12.3
             out <- nimLogical(5, value = v);
         }),
         outputType = quote(logical(1))
         ),
    list(name = 'logical: length, vector value expression',
         expr = quote({
             v <- rpois(5, lambda = 12.3) > 12.3
             out <- nimLogical(5, value = !v);
         }),
         outputType = quote(logical(1))
         ),
    list(name = 'logical: length, vector value with too few elements',
         expr = quote({
             v <- rpois(4, lambda = 12.3) > 12.3
             out <- nimLogical(5, value = v);
         }),
         outputType = quote(logical(1))
         ),
    list(name = 'logical: length vector value with different types',
         expr = quote({
             v <- nimLogical(5, value = c(TRUE, TRUE, rep(FALSE, 2), TRUE))
             out <- nimLogical(5, value = v);
         }),
         outputType = quote(logical(1))
         ),
    list(name = 'logical: length vector value with too many elements',
         expr = quote({
             v <- rpois(6, lambda = 12.3) > 12.3
             out <- nimLogical(5, value = v);
         }),
         outputType = quote(logical(1))
         ),
    list(name = 'logical: length,  matrix value',
         expr = quote({
             v <- matrix(rpois(6, lambda = 12.3) > 12.3, nrow = 2)
             out <- nimLogical(6, value = v);
         }),
         outputType = quote(logical(1))
         )
)

matrixTests <- list(
    list(name = 'matrix',
         expr = quote({
             out <- nimMatrix(6, nrow = 2);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix',
         expr = quote({
             out <- nimMatrix(6, ncol = 2);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix size 0',
         expr = quote({
             out <- nimMatrix(type = 'double', nrow = 0, ncol = 0)
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix, recycle = FALSE',
         expr = quote({
             out <- nimMatrix(6, ncol = 2, recycle = FALSE);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix, recycle = FALSE',
         expr = quote({
             out <- nimMatrix(6, nrow = 2, recycle = FALSE);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix, recycle = FALSE',
         expr = quote({
             out <- nimMatrix(6, ncol = 3, nrow = 2, recycle = FALSE);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix',
         expr = quote({
             out <- nimMatrix(6, ncol = 2, nrow = 3);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix no init',
         expr = quote({
             test <- nimMatrix(6, ncol = 3, init = FALSE);
             out <- dim(test)
         }),
         outputType = quote(integer(1))
         ),
    list(name = 'matrix, scalar value, in exression',
         expr = quote({
             out <- nimMatrix(6, ncol = 3, nrow = 2)^2;
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix vector init',
         expr = quote({
             out <- nimMatrix(1:6, ncol = 3);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix, vector value, in exression',
         expr = quote({
             out <- nimMatrix(1:6, ncol = 3, nrow = 2)^2;
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix matrix init',
         expr = quote({
             initMat <- nimMatrix(1:6, nrow = 3);
             out <- nimMatrix(initMat, nrow = 2);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix vector no sizes given',
         expr = quote({
             out <- nimMatrix(1:6)
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix vector init too big',
         expr = quote({
             out <- nimMatrix(1:6, nrow = 2)
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix vector init',
         expr = quote({
             out <- nimMatrix(1:6, nrow = 3);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix vector init, extraneous arg',
         expr = quote({
             out <- nimMatrix(1:6, nrow = 2, recycle = FALSE);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix vector init',
         expr = quote({
             out <- nimMatrix(1:6, ncol = 2, nrow = 3);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix vector init too small',
         expr = quote({
             out <- nimMatrix(1:6, nrow = 3, ncol = 3);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix vector init too big',
         expr = quote({
             out <- nimMatrix(1:6, nrow = 2, ncol = 2)
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix vector init incongruous size',
         expr = quote({
             out <- nimMatrix(1:7, nrow = 3);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'matrix vector init too small, recycle=FALSE',
         expr = quote({
             out <- nimMatrix(1:6, nrow = 3, ncol = 3, recycle = FALSE);
         }),
         outputType = quote(double(2))
         ),
    ## for nimNewMatrix case, zeros will be filled in anyway
    list(name = 'matrix vector init too small, recycle=FALSE, fillZeros = FALSE (ignored)',
         expr = quote({
             out <- nimMatrix(1:6, nrow = 3, ncol = 3, recycle = FALSE, fillZeros = FALSE);
         }),
         outputType = quote(double(2))
         )
)

arrayTests1D <- list(
    list(name = 'array 1D',
         expr = quote({
             out <- nimArray(6, dim = 2);
         }),
         outputType = quote(double(1))
         ),
    list(name = 'array 1D, recycle = FALSE',
         expr = quote({
             out <- nimArray(6, dim = 2, recycle = FALSE);
         }),
         outputType = quote(double(1))
         ),
    list(name = 'array 1D, value scalar, fillZeros = TRUE (but recycle=TRUE takes precedence)',
         expr = quote({
             out <- nimArray(1.23, dim = 2, fillZeros = TRUE);
         }),
         outputType = quote(double(1))
         ),
    list(name = 'array 1D, value scalar, fillZeros = TRUE, recycle = FALSE',
         expr = quote({
             out <- nimArray(1.23, dim = 2, fillZeros = TRUE, recycle = FALSE);
         }),
         outputType = quote(double(1))
         ),
    list(name = 'array 1D (c() notation)',
         expr = quote({
             out <- nimArray(6, dim = c(2));
         }),
         outputType = quote(double(1))
         ),
    list(name = 'array 1D (dim vector)',
         expr = quote({
             dim <- c(2)
             out <- nimArray(6, dim = dim, nDim = 1);
         }),
         outputType = quote(double(1))
         ),
    list(name = 'array 1D in expression',
         expr = quote({
             out <- nimArray(6, dim = 2)^2;
         }),
         outputType = quote(double(1))
         ),
    list(name = 'array 1D, value vector',
         expr = quote({
             out <- nimArray(rnorm(2), dim = 2);
         }),
         outputType = quote(double(1))
         ),
    list(name = 'array 1D, value vector, in expression',
         expr = quote({
             out <- nimArray(rnorm(2), dim = 2)^2;
         }),
         outputType = quote(double(1))
         ),
    list(name = 'array 1D, value vector, recycle = TRUE',
         expr = quote({
             out <- nimArray(rnorm(2), dim = 4, recycle = TRUE);
         }),
         outputType = quote(double(1))
         ),
    list(name = 'array 1D, value vector, fillZeros = TRUE (but recycle=TRUE has precedence)',
         expr = quote({
             out <- nimArray(rnorm(2), dim = 4, fillZeros = TRUE);
         }),
         outputType = quote(double(1))
         ),
    list(name = 'array 1D, value vector, fillZeros = TRUE, recycle=FALSE',
         expr = quote({
             out <- nimArray(rnorm(2), dim = 4, fillZeros = TRUE, recycle = FALSE);
         }),
         outputType = quote(double(1))
     ) 
)
arrayTests3D <- list(
    list(name = 'array 3D',
         expr = quote({
             out <- nimArray(6, dim = c(2, 3, 4));
         }),
         outputType = quote(double(3))
         ),
    list(name = 'array 3D size 0',
         expr = quote({
             out <- nimArray(type = 'double', dim = c(0,0,0))
         }),
         outputType = quote(double(3))
         ),
    list(name = 'array 3D, recycle = FALSE',
         expr = quote({
             out <- nimArray(6, dim = c(2, 3, 4), recycle = FALSE);
         }),
         outputType = quote(double(3))
         ),
    list(name = 'array 3D, value scalar, fillZeros = TRUE (but recycle=TRUE takes precedence)',
         expr = quote({
             out <- nimArray(1.23, dim = c(2, 3, 4), fillZeros = TRUE);
         }),
         outputType = quote(double(3))
         ),
    list(name = 'array 3D, value scalar, fillZeros = TRUE, recycle = FALSE',
         expr = quote({
             out <- nimArray(1.23, dim = c(2, 3, 4), fillZeros = TRUE, recycle = FALSE);
         }),
         outputType = quote(double(3))
         ),
    list(name = 'array 3D (dim vector)',
         expr = quote({
             dim <- c(2, 3, 4)
             out <- nimArray(6, dim = dim, nDim = 3);
         }),
         outputType = quote(double(3))
         ),
    list(name = 'array 3D in expression',
         expr = quote({
             out <- nimArray(6, dim = c(2, 3, 4))[,2,]^2;
         }),
         outputType = quote(double(2))
         ),
    list(name = 'array 3D, value vector in bracket expression',
         expr = quote({
             out <- nimArray(1:7, dim = c(3, 4, 2))[,2,]
         }),
         outputType = quote(double(2))
         ),
    list(name = 'array 3D, value vector',
         expr = quote({
             out <- nimArray(rnorm(2*3*4), c(2, 3, 4));
         }),
         outputType = quote(double(3))
         ),
    list(name = 'array 3D, value vector, recycle = TRUE',
         expr = quote({
             out <- nimArray(rnorm(2), dim = c(2, 3, 4), recycle = TRUE);
         }),
         outputType = quote(double(3))
         ),
    list(name = 'array 3D, value vector, fillZeros = TRUE (but recycle=TRUE has precedence)',
         expr = quote({
             out <- nimArray(rnorm(2), dim = c(2, 3, 4), fillZeros = TRUE);
         }),
         outputType = quote(double(3))
         ),
    list(name = 'array 3D, value vector, fillZeros = TRUE, recycle=FALSE',
         expr = quote({
             out <- nimArray(rnorm(2), dim = c(2, 3, 4), fillZeros = TRUE, recycle = FALSE);
         }),
         outputType = quote(double(3))
         ),
    list(name = 'array matrix init',
         expr = quote({
             initMat <- nimMatrix(1:18, nrow = 3);
             out <- nimArray(initMat, dim = c(3,3,2));
         }),
         outputType = quote(double(3))
         ),
    list(name = 'array 3D too few values',
         expr = quote({
             out <- nimArray(1:9, dim = c(3,3,2));
         }),
         outputType = quote(double(3))
         ),
    list(name = 'array 3D, value vector incongruous to dims',
         expr = quote({
             out <- nimArray(1:7, dim = c(3, 4, 2))
         }),
         outputType = quote(double(3))
         )
    )

arrayTests2D <- list(
    list(name = 'array 2D',
         expr = quote({
             out <- nimArray(6, dim = c(2, 3));
         }),
         outputType = quote(double(2))
         ),
    list(name = 'array 2D, recycle = FALSE',
         expr = quote({
             out <- nimArray(6, dim = c(2, 3), recycle = FALSE);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'array 2D, value scalar, fillZeros = TRUE (but recycle=TRUE takes precedence)',
         expr = quote({
             out <- nimArray(1.23, dim = c(2, 3), fillZeros = TRUE);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'array 2D, value scalar, fillZeros = TRUE, recycle = FALSE',
         expr = quote({
             out <- nimArray(1.23, dim = c(2, 3), fillZeros = TRUE, recycle = FALSE);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'array 2D (dim vector)',
         expr = quote({
             dim <- c(2, 3)
             out <- nimArray(6, dim = dim, nDim = 2);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'array 2D in expression',
         expr = quote({
             out <- nimArray(6, dim = c(2, 3))^2;
         }),
         outputType = quote(double(2))
         ),
    list(name = 'array 2D, value vector',
         expr = quote({
             out <- nimArray(rnorm(2), dim = c(2, 3));
         }),
         outputType = quote(double(2))
         ),
    list(name = 'array 2D, value vector, recycle = TRUE',
         expr = quote({
             out <- nimArray(rnorm(2), dim = c(2, 3), recycle = TRUE);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'array 2D, value vector, fillZeros = TRUE (but recycle=TRUE has precedence)',
         expr = quote({
             out <- nimArray(rnorm(2), dim = c(2, 3), fillZeros = TRUE);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'array 2D, value vector, fillZeros = TRUE, recycle=FALSE',
         expr = quote({
             out <- nimArray(rnorm(2), dim = c(2, 3), fillZeros = TRUE, recycle = FALSE);
         }),
         outputType = quote(double(2))
         ),
    list(name = 'array 2D, value vector with indexed block expression',
         expr = quote({
             out <- nimArray(1:7, dim = c(3, 4))[,2]^2
         }),
         outputType = quote(double(1))
         ),
    list(name = 'array matrix init',
         expr = quote({
             initMat <- nimMatrix(1:6, nrow = 3);
             out <- nimArray(initMat, dim = c(3,3));
         }),
         outputType = quote(double(2))
         ),
    list(name = 'array 2D too few values',
         expr = quote({
             out <- nimArray(1:6, dim = c(3,3));
         }),
         outputType = quote(double(2))
         )
    )


setSize1D <- list(
    list(name = 'setSize 1D',
         expr = quote({
             out <- nimNumeric(4, value = 1:2)
             setSize(out, 7)
         }),
         outputType = quote(double(1))
         ),
    list(name = 'setSize 1D dim expresion',
         expr = quote({
             out <- nimNumeric(4, value = 1:2)
             z <- rpois(2, lambda = 1)
             setSize(out, 1 + sum(z))
         }),
         outputType = quote(double(1))
         ),
    list(name = 'setSize 1D to 0',
         expr = quote({
             out <- nimNumeric(4, value = 1:2)
             setSize(out, 0)
         }),
         outputType = quote(double(1))
         ),
    list(name = 'setSize 1D, copy=FALSE (default fillValues = TRUE)', 
         expr = quote({
             out <- nimNumeric(4, value = 1:2)
             setSize(out, 7, copy=FALSE)
         }),
         outputType = quote(double(1))
         ),
    list(name = 'setSize 1D, (default copy = TRUE), fillZeros = FALSE', ## weak test since we don't check on the non-filled zeros since they shouldn't be equal
         expr = quote({
             v <- nimNumeric(4, value = 1:2)
             setSize(v, 7, fillZeros = FALSE)
             out <- v[1:4]
         }),
         outputType = quote(double(1))
         ),
    list(name = 'setSize 1D, copy=FALSE, fillZeros = FALSE', ## weaker test since we can only check length of resized object
         expr = quote({
             v <- nimNumeric(4, value = 1:2)
             setSize(v, 7, copy = FALSE, fillZeros = FALSE)
             out <- length(v)
         }),
         outputType = quote(integer(0))
         )
    )

setSize1Dinteger <- list(
    list(name = 'setSize integer 1D',
         expr = quote({
             out <- nimInteger(4, value = 1:2)
             setSize(out, 7)
         }),
         outputType = quote(integer(1))
         ),
    list(name = 'setSize integer 1D dim expresion',
         expr = quote({
             out <- nimInteger(4, value = 1:2)
             z <- rpois(2, lambda = 1)
             setSize(out, 1 + sum(z))
         }),
         outputType = quote(integer(1))
         ),
    list(name = 'setSize integer 1D to 0',
         expr = quote({
             out <- nimInteger(4, value = 1:2)
             setSize(out, 0)
         }),
         outputType = quote(integer(1))
         ),
    list(name = 'setSize integer 1D, copy=FALSE (default fillValues = TRUE)', 
         expr = quote({
             out <- nimInteger(4, value = 1:2)
             setSize(out, 7, copy=FALSE)
         }),
         outputType = quote(integer(1))
         ),
    list(name = 'setSize integer 1D, (default copy = TRUE), fillZeros = FALSE', ## weak test since we don't check on the non-filled zeros since they shouldn't be equal
         expr = quote({
             v <- nimInteger(4, value = 1:2)
             setSize(v, 7, fillZeros = FALSE)
             out <- v[1:4]
         }),
         outputType = quote(integer(1))
         ),
    list(name = 'setSize integer 1D, copy=FALSE, fillZeros = FALSE', ## weaker test since we can only check length of resized object
         expr = quote({
             v <- nimInteger(4, value = 1:2)
             setSize(v, 7, copy = FALSE, fillZeros = FALSE)
             out <- length(v)
         }),
         outputType = quote(integer(0))
         )
    )

setSize2D <- list(
    list(name = 'setSize 2D',
         expr = quote({
             out <- nimMatrix(1:5, nrow = 2, ncol = 3)
             setSize(out, c(3, 7))
         }),
         outputType = quote(double(2))
         ),
    list(name = 'setSize 2D dim expresion',
         expr = quote({
             out <- nimMatrix(1:5, nrow = 2, ncol = 3)
             z <- rpois(2, lambda = 1)
             setSize(out, c(3, 1 + sum(z)))
         }),
         outputType = quote(double(2))
         ),
    list(name = 'setSize 2D to 0',
         expr = quote({
             out <- nimMatrix(1:5, nrow = 2, ncol = 3)
             setSize(out, c(0, 0))
         }),
         outputType = quote(double(2))
         ),
    list(name = 'setSize 2D, copy=FALSE (default fillValues = TRUE)', 
         expr = quote({
             out <- nimMatrix(1:5, nrow = 2, ncol = 3)
             setSize(out, c(3,7), copy=FALSE)
         }),
         outputType = quote(double(2))
         ),
    list(name = 'setSize 2D, (default copy = TRUE), fillZeros = FALSE', ## weak test since we don't check on the non-filled zeros since they shouldn't be equal
         expr = quote({
             v <- nimMatrix(1:5, nrow = 2, ncol = 3)
             setSize(v, c(3,7), fillZeros = FALSE)
             out <- v[1:3, 1:2]
         }),
         outputType = quote(double(2))
         ),
    list(name = 'setSize 2D, copy=FALSE, fillZeros = FALSE', ## weaker test since we can only check length of resized object
         expr = quote({
             v <- nimMatrix(1:5, nrow = 2, ncol = 3)
             setSize(v, c(3,7), copy = FALSE, fillZeros = FALSE)
             out <- nimDim(v)
         }),
         outputType = quote(integer(1))
         )
    )

setSize3D <- list(
    list(name = 'setSize 3D',
         expr = quote({
             out <- nimArray(1:5, dim = c(2,3,1))
             setSize(out, c(7, 2, 2))
         }),
         outputType = quote(double(3))
         ),
    list(name = 'setSize 3D dim expresion',
         expr = quote({
             out <- nimArray(1:5, dim = c(2,3,1))
             z <- rpois(2, lambda = 1)
             setSize(out, c(3, 1 + sum(z), 2))
         }),
         outputType = quote(double(3))
         ),
    list(name = 'setSize 3D to 0',
         expr = quote({
             out <- nimArray(1:5, dim = c(2,3,1))
             setSize(out, c(0, 0, 0))
         }),
         outputType = quote(double(3))
         ),
    list(name = 'setSize 3D, copy=FALSE (default fillValues = TRUE)', 
         expr = quote({
             out <- nimArray(1:5, dim = c(2,3,1))
             setSize(out, c(7, 2, 2), copy=FALSE)
         }),
         outputType = quote(double(3))
         ),
    list(name = 'setSize 3D, (default copy = TRUE), fillZeros = FALSE', ## weak test since we don't check on the non-filled zeros since they shouldn't be equal
         expr = quote({
             v <- nimArray(1:5, dim = c(2,3,1))
             setSize(v, c(3, 2, 2), fillZeros = FALSE)
             out <- v[1:3, 1:2, 1]
         }),
         outputType = quote(double(2))
         ),
    list(name = 'setSize 3D, copy=FALSE, fillZeros = FALSE', ## weaker test since we can only check length of resized object
         expr = quote({
             v <- nimArray(1:5, dim = c(2,3,1))
             setSize(v, c(3, 7, 2), copy = FALSE, fillZeros = FALSE)
             out <- nimDim(v)
         }),
         outputType = quote(integer(1))
         )
    )

expectedCompilerErrors <- list(
    list(name = 'numeric: fail for length given as vector',
         expr = quote({
             a <- 2:4
             out <- nimNumeric(a, value = 1.5);
         }),
         outputType = quote(double(1)),
         expectedCompilerError = TRUE
         ),
    list(name = 'numeric: fail for init given as vector',
         expr = quote({
             b <- c(TRUE, FALSE, TRUE)
             out <- nimNumeric(3, value = 1.5, init = b);
         }),
         outputType = quote(double(1)),
         expectedCompilerError = TRUE
         ),
    list(name = 'integer: fail for length given as vector',
         expr = quote({
             a <- 2:4
             out <- nimInteger(a, value = 15);
         }),
         outputType = quote(integer(1)),
         expectedCompilerError = TRUE
         ),
    list(name = 'integer: fail for init given as vector',
         expr = quote({
             b <- c(TRUE, FALSE, TRUE)
             out <- nimInteger(3, value = 15, init = b);
         }),
         outputType = quote(integer(1)),
         expectedCompilerError = TRUE
         ),
    list(name = 'logical: fail for length given as vector',
         expr = quote({
             a <- 2:4
             out <- nimLogical(a, value = TRUE);
         }),
         outputType = quote(logical(1)),
         expectedCompilerError = TRUE
         ),
    list(name = 'logical: fail for init given as vector',
         expr = quote({
             b <- c(TRUE, FALSE, TRUE)
             out <- nimLogical(3, value = TRUE, init = b);
         }),
         outputType = quote(logical(1)),
         expectedCompilerError = TRUE
         ),
    list(name = 'array 1D (dim vector without nDim: safe compiler fail)',
         expr = quote({
             dim <- c(2)
             out <- nimArray(6, dim = dim);
         }),
         outputType = quote(double(1)),
         expectedCompilerError = TRUE
         ),
    list(name = 'array 3D (dim vector without nDim: safe compiler fail)',
         expr = quote({
             dim <- c(2, 3, 4)
             out <- nimArray(6, dim = dim);
         }),
         outputType = quote(double(3)),
         expectedCompilerError = TRUE
         ),
    list(name = 'array 2D (dim vector without nDim: safe compiler fail)',
         expr = quote({
             dim <- c(2, 3)
             out <- nimArray(6, dim = dim);
         }),
         outputType = quote(double(2)),
         expectedCompilerError = TRUE
         )
    )

numericTestResults <- test_coreRfeature_batch(numericTests, 'numericTests') ## lapply(numericTests, test_coreRfeature)
numericTestResults <- test_coreRfeature_batch(integerTests, 'integerTests') ## lapply(integerTests, test_coreRfeature)
logicalTestResults <- test_coreRfeature_batch(logicalTests, 'logicalTests') ## lapply(logicalTests, test_coreRfeature)
matrixTestResults <- test_coreRfeature_batch(matrixTests, 'matrixTests') ## lapply(matrixTests, test_coreRfeature)
array1DTestResults <- test_coreRfeature_batch(arrayTests1D, 'arrayTests1D') ## lapply(arrayTests1D, test_coreRfeature)
array2DTestResults <- test_coreRfeature_batch(arrayTests2D, 'arrayTests2D') ## lapply(arrayTests2D, test_coreRfeature)
array3DTestResults <- test_coreRfeature_batch(arrayTests3D, 'arrayTests3D') ## lapply(arrayTests3D, test_coreRfeature)
setSize1DResults <- test_coreRfeature_batch(setSize1D, 'setSize1D') ## lapply(setSize1D, test_coreRfeature)
setSize1DintegerResults <- test_coreRfeature_batch(setSize1Dinteger, 'setSize1Dinteger') ## lapply(setSize1Dinteger, test_coreRfeature)
setSize2DResults <- test_coreRfeature_batch(setSize2D, 'setSize2D') ## lapply(setSize2D, test_coreRfeature)
setSize3DResults <- test_coreRfeature_batch(setSize3D, 'setSize3D') ## lapply(setSize3D, test_coreRfeature)

allocationExpectedCompilerFailures <- lapply(expectedCompilerErrors, test_coreRfeature)

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)

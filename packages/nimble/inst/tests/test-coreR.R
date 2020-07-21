source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context("Testing of core R functions in NIMBLE code")

## fix result_type in nimbleEigen.h

cTests <- list(
    list(name = "c(double, double)", expr = quote(out <- c(arg1, arg2)), args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(4:5)}), outputType = quote(double(1))),
    list(name = "c(double, integer)", expr = quote(out <- c(arg1, arg2)), args = list(arg1 = quote(double(1)), arg2 = quote(integer(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.integer(4:5)}), outputType = quote(double(1))),
    list(name = "c(double, logical)", expr = quote(out <- c(arg1, arg2)), args = list(arg1 = quote(double(1)), arg2 = quote(logical(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- c(TRUE, FALSE, TRUE)}), outputType = quote(double(1))),


    list(name = "c(integer, double)", expr = quote(out <- c(arg1, arg2)), args = list(arg1 = quote(integer(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- as.integer(1:3); arg2 <- as.numeric(4:5)}), outputType = quote(double(1))),
    list(name = "c(integer, integer)", expr = quote(out <- c(arg1, arg2)), args = list(arg1 = quote(integer(1)), arg2 = quote(integer(1))),
         setArgVals = quote({arg1 <- as.integer(1:3); arg2 <- as.integer(4:5)}), outputType = quote(integer(1))),
    list(name = "c(integer, logical)", expr = quote(out <- c(arg1, arg2)), args = list(arg1 = quote(integer(1)), arg2 = quote(logical(1))),
         setArgVals = quote({arg1 <- as.integer(1:3); arg2 <- c(TRUE, FALSE, TRUE)}), outputType = quote(integer(1))),

    list(name = "c(logical, double)", expr = quote(out <- c(arg1, arg2)), args = list(arg1 = quote(logical(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- c(FALSE, TRUE, FALSE); arg2 <- as.numeric(4:5)}), outputType = quote(double(1))),

    
    list(name = "c(double(2), double)", expr = quote(out <- c(arg1, arg2)), args = list(arg1 = quote(double(2)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:4), nrow = 2); arg2 <- as.numeric(10:11)}), outputType = quote(double(1))),

        list(name = "c(double, double, double)", expr = quote(out <- c(arg1, arg2, arg3)),
         args = list(arg1 = quote(double(1)), arg2 = quote(double(1)), arg3 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3);
                             arg2 <- as.numeric(4:5);
                             arg3 <- as.numeric(10:15)}), outputType = quote(double(1))),
    list(name = "c(double, double, double, double)", expr = quote(out <- c(arg1, arg2, arg3, arg4)),
         args = list(arg1 = quote(double(1)), arg2 = quote(double(1)), arg3 = quote(double(1)), arg4 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3);
                             arg2 <- as.numeric(4:5);
                             arg3 <- as.numeric(10:15);
                             arg4 <- as.numeric(100:120)}), outputType = quote(double(1))),
    list(name = "c(double, double, double, double, double)", expr = quote(out <- c(arg1, arg2, arg3, arg4, arg5)),
         args = list(arg1 = quote(double(1)), arg2 = quote(double(1)), arg3 = quote(double(1)), arg4 = quote(double(1)), arg5 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3);
                             arg2 <- as.numeric(4:5);
                             arg3 <- as.numeric(10:15);
                             arg4 <- as.numeric(100:120);
                             arg5 <- as.numeric(200:203)}), outputType = quote(double(1))),
    list(name = "c(double, double, 1, 2, 3, double)", expr = quote(out <- c(arg1, arg2, 1, 2, 3, arg5)),
         args = list(arg1 = quote(double(1)), arg2 = quote(double(1)), arg3 = quote(double(1)), arg4 = quote(double(1)), arg5 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3);
                             arg2 <- as.numeric(4:5);
                             arg3 <- as.numeric(10:15);
                             arg4 <- as.numeric(100:120);
                             arg5 <- as.numeric(200:203)}), outputType = quote(double(1))),
    list(name = "c(double, double, double, 1, 2, 3, double)", expr = quote(out <- c(arg1, arg2, arg3, 1, 2, 3, arg5)),
         args = list(arg1 = quote(double(1)), arg2 = quote(double(1)), arg3 = quote(double(1)), arg4 = quote(double(1)), arg5 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3);
                             arg2 <- as.numeric(4:5);
                             arg3 <- as.numeric(10:15);
                             arg4 <- as.numeric(100:120);
                             arg5 <- as.numeric(200:203)}), outputType = quote(double(1))),
    list(name = "c(double, double, double, 1, 2, 3, double, 5, 6, 7, double)", expr = quote(out <- c(arg1, arg2, arg3, 1, 2, 3, arg4, 5, 6, 7, arg5)),
         args = list(arg1 = quote(double(1)), arg2 = quote(double(1)), arg3 = quote(double(1)), arg4 = quote(double(1)), arg5 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3);
                             arg2 <- as.numeric(4:5);
                             arg3 <- as.numeric(10:15);
                             arg4 <- as.numeric(100:120);
                             arg5 <- as.numeric(200:203)}), outputType = quote(double(1))),
    list(name = "expressions: c(double, double)", expr = quote(out <- log(c(arg1 + 1, arg2 + 2)) + 1), args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(4:5)}), outputType = quote(double(1)))
)

blockTests <- list(
    ##1
    ## basics
    list(name = "3x3 block simple copy", expr = quote(out <- arg1[2:4, 2:4]), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(2))),
    list(name = "3x3 block simple copy non-arg", expr = quote({temp <- arg1; out <- temp[2:4, 2:4]}), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(2))),
    list(name = "3x3 block", expr = quote(out <- arg1[2:4, 2:4] + 2), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(2))),
    list(name = "3xfull block", expr = quote(out <- arg1[2:4, ] + 2), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(2))),
    list(name = "fullx3 block", expr = quote(out <- arg1[, 2:4] + 2), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(2))),
    ##6
    list(name = "fullxfull block", expr = quote(out <- arg1[, ] + 2), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(2))),
    ## expressions in index ranges (should be lifted) 
    list(name = "3x3 block variable index range", expr = quote({i <- 1; j <- 3; out <- arg1[(j-1):(j+1), 2:4] + 2}), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(2))),
    ## dropping a dimension
    list(name = "3x1 block", expr = quote(out <- arg1[2:4, 3] + 2), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(1))),
    list(name = "1x3 block", expr = quote(out <- arg1[3, 2:4] + 2), args = list(arg1 = quote(double(2))), ## OOPS THIS IS NOT OK
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(1))),
    list(name = "3x1 block non-arg", expr = quote({temp <- arg1;out <- temp[2:4, 3] + 2}), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(1))),
    #11
    list(name = "fullx1 block", expr = quote(out <- arg1[, 3] + 2), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(1))),
    list(name = "1xfull block", expr = quote(out <- arg1[3, ] + 2), args = list(arg1 = quote(double(2))), ## OOPS THIS IS NOT OK
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(1))),
    list(name = "fullx1 block non-arg", expr = quote({temp <- arg1;out <- temp[, 3] + 2}), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(1))),
    ## not dropping a dimension (scalar index but with drop  = FALSE)
    list(name = "3x1 block drop = FALSE", expr = quote(out <- arg1[2:4, 3, drop = FALSE] + 2), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(2))),
    list(name = "1x3 block drop = FALSE", expr = quote(out <- arg1[3, 2:4, drop = FALSE] + 2), args = list(arg1 = quote(double(2))), ## OOPS THIS IS NOT OK
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(2))),
    ##16
    list(name = "3x1 block non-arg drop = FALSE", expr = quote({temp <- arg1;out <- temp[2:4, 3, drop = FALSE] + 2}), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(2))),
    list(name = "1xfull block drop = FALSE", expr = quote(out <- arg1[3, , drop = FALSE] + 2), args = list(arg1 = quote(double(2))), ## OOPS THIS IS NOT OK
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(2))),
    list(name = "fullx1 block non-arg drop = FALSE", expr = quote({temp <- arg1;out <- temp[, 3, drop = FALSE] + 2}), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(2))),
    ## scalar indices via 2:2 - THERE IS NOTHING WE CAN DO AT COMPILE TIME TO DISTINGUISH 2:2 and 2:j
    list(name = "3x1 block with 2:2 index", expr = quote(out <- arg1[2:4, 3:3] + 2), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(1))),
    list(name = "1x3 block with 2:2 index", expr = quote(out <- arg1[3:3, 2:4] + 2), args = list(arg1 = quote(double(2))), ## OOPS THIS IS NOT OK
         setArgVals = quote(arg1 <- matrix(as.numeric(1:25), nrow = 5)), outputType = quote(double(1))),
    ##21
    ## higher dimensions --> 2D or 1D
    list(name = "3x3x1 block", expr = quote(out <- arg1[2:4, 3:5, 2] + 2), args = list(arg1 = quote(double(3))),
         setArgVals = quote(arg1 <- array(as.numeric(1:120), dim = c(4, 5, 6))), outputType = quote(double(2))),
    list(name = "3x1x3 block", expr = quote(out <- arg1[2:4, 2, 3:5] + 2), args = list(arg1 = quote(double(3))),
         setArgVals = quote(arg1 <- array(as.numeric(1:120), dim = c(4, 5, 6))), outputType = quote(double(2))),
    list(name = "1x3x3 block", expr = quote(out <- arg1[2, 2:4, 3:5] + 2), args = list(arg1 = quote(double(3))),
         setArgVals = quote(arg1 <- array(as.numeric(1:120), dim = c(4, 5, 6))), outputType = quote(double(2))),
    list(name = "3x1x1 block", expr = quote(out <- arg1[2:4, 4, 2] + 2), args = list(arg1 = quote(double(3))),
         setArgVals = quote(arg1 <- array(as.numeric(1:120), dim = c(4, 5, 6))), outputType = quote(double(1))),
    list(name = "1x1x3 block", expr = quote(out <- arg1[3, 2, 3:5] + 2), args = list(arg1 = quote(double(3))),
         setArgVals = quote(arg1 <- array(as.numeric(1:120), dim = c(4, 5, 6))), outputType = quote(double(1))),
    ##26
    list(name = "1x3x1 block", expr = quote(out <- arg1[3, 2:4, 5] + 2), args = list(arg1 = quote(double(3))),
         setArgVals = quote(arg1 <- array(as.numeric(1:120), dim = c(4, 5, 6))), outputType = quote(double(1))),
    ## chained
    list(name = "3x3 chained to 2x2 block simple copy", expr = quote(out <- arg1[2:5, 3:6][2:3, 3:4]), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:36), nrow = 6)), outputType = quote(double(2))),
    list(name = "3x3 chained to 2x2 block simple copy non-arg", expr = quote({temp <- arg1; out <- temp[2:5, 3:6][2:3, 3:4]}), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:36), nrow = 6)), outputType = quote(double(2))),
    list(name = "3x3 chained to 2x2  block", expr = quote(out <- arg1[2:5, 3:6][2:3, 3:4] + 2), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:36), nrow = 6)), outputType = quote(double(2))),
    list(name = "3xfull chained to 2x2  block", expr = quote(out <- arg1[2:5, ][2:3, 3:4] + 2), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:36), nrow = 6)), outputType = quote(double(2))),
    ##31
    list(name = "fullx3 chained to 2x2  block", expr = quote(out <- arg1[, 2:5][2:3, 3:4] + 2), args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:36), nrow = 6)), outputType = quote(double(2))),
    ## chained from map to block
    list(name = "3x3x1 chained to 2x2 block", expr = quote(out <- arg1[2:4, 2:5, 2][2:3, 3:4] + 2), args = list(arg1 = quote(double(3))),
         setArgVals = quote(arg1 <- array(as.numeric(1:120), dim = c(4, 5, 6))), outputType = quote(double(2))),
    list(name = "3x3x1 chained to 2x2 block non-arg", expr = quote({temp <- arg1; out <- temp[2:4, 2:5, 2][2:3, 3:4] + 2}), args = list(arg1 = quote(double(3))),
         setArgVals = quote(arg1 <- array(as.numeric(1:120), dim = c(4, 5, 6))), outputType = quote(double(2)))

## Following is not currently supported: 
##    list(name = "5d nimArray map copy", expr = quote({temp <- arg1; out <- temp[2:4, 3:6, 2, 4, 1:3]}), args = list(arg1 = quote(double(5))),
##         setArgVals = quote(arg1 <- array(as.numeric(1:(5^5)), dim = c(5, 5, 5, 5, 5))), outputType = quote(double(3)))
    ## tests to add:
    ## integer and logical types
    ## input passed non-trivially as map
    ## all scalar indices
)

repTests <- list(
    ##1
    ## basic cases with x and times
    list(name = "rep(1, 3)", expr = quote(out <- rep(1, 3)), args = list(arg1 = quote(double(0))),
         setArgVals = quote({arg1 <- 3}), outputType = quote(double(1))),
    list(name = "rep(vector double, 3)", expr = quote(out <- rep(arg1, 3)), args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3)}), outputType = quote(double(1))),
    list(name = "rep(vector integer, 3)", expr = quote(out <- rep(arg1, 3)), args = list(arg1 = quote(integer(1))),
         setArgVals = quote({arg1 <- as.integer(1:3)}), outputType = quote(integer(1))),
    list(name = "rep(vector logical, 3)", expr = quote(out <- rep(arg1, 3)), args = list(arg1 = quote(logical(1))),
         setArgVals = quote({arg1 <- c(TRUE, FALSE, FALSE)}), outputType = quote(logical(1))),
    list(name = "rep(vector double, variable)", expr = quote(out <- rep(arg1, arg2)), args = list(arg1 = quote(double(1)), arg2 = quote(integer())),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- 4}), outputType = quote(double(1))),
    ##6
    ## cases with x and each
    list(name = "rep(1, 3)", expr = quote(out <- rep(1, each = 3)), args = list(arg1 = quote(double(0))),
         setArgVals = quote({arg1 <- 3}), outputType = quote(double(1))),
    list(name = "rep(vector double, 3)", expr = quote(out <- rep(arg1, each = 3)), args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3)}), outputType = quote(double(1))),
    list(name = "rep(vector integer, 3)", expr = quote(out <- rep(arg1, each = 3)), args = list(arg1 = quote(integer(1))),
         setArgVals = quote({arg1 <- as.integer(1:3)}), outputType = quote(integer(1))),
    list(name = "rep(vector logical, 3)", expr = quote(out <- rep(arg1, each = 3)), args = list(arg1 = quote(logical(1))),
         setArgVals = quote({arg1 <- c(TRUE, FALSE, FALSE)}), outputType = quote(logical(1))),
    list(name = "rep(vector double, variable)", expr = quote(out <- rep(arg1, each = arg2)), args = list(arg1 = quote(double(1)), arg2 = quote(integer())),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- 4}), outputType = quote(double(1))),
    ##11
    list(name = "rep(vector double, first arg)", expr = quote(out <- rep(arg1, each = arg2)), args = list(arg1 = quote(double(1)), arg2 = quote(integer(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.integer(4:5)}), outputType = quote(double(1)), expectWarnings = list("R eval" = 'Expected warning: vector each', "R run" = "Expected warning: vector each")),

    ## basic cases with x, times and each
    list(name = "rep(vector double, variable, each = 2)", expr = quote(out <- rep(arg1, times = arg2, each = 2)), args = list(arg1 = quote(double(1)), arg2 = quote(double(0))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(c(4))}), outputType = quote(double(1))),
    list(name = "rep(vector double, variable, variable)", expr = quote(out <- rep(arg1, times = arg2, each = arg3)), args = list(arg1 = quote(double(1)), arg2 = quote(double(0)), arg3 = quote(double(0))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(4); arg3 <- 5}), outputType = quote(double(1))),
    list(name = "rep(vector double, variable, first arg)", expr = quote(out <- rep(arg1, times = arg2, each = arg3)), args = list(arg1 = quote(double(1)), arg2 = quote(double(0)), arg3 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- 4; arg3 <- c(5, 7)}), outputType = quote(double(1)), expectWarnings = list("R eval" = 'Expected warning: vector each', "R run" = 'Expected warning: vector each')),
    ## basic cases with x and length.out
    list(name = "rep(1, length.out = 3)", expr = quote(out <- rep(1, length.out = 3)), args = list(arg1 = quote(double(0))),
         setArgVals = quote({arg1 <- 3}), outputType = quote(double(1))),
    ## 16
    list(name = "rep(vector double, length.out = larger than vector)", expr = quote(out <- rep(arg1, length.out = 5)), args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3)}), outputType = quote(double(1))),
    list(name = "rep(vector double, length.out = smaller than vector)", expr = quote(out <- rep(arg1, length.out = 2)), args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3)}), outputType = quote(double(1))),
    list(name = "rep(vector double, length.out = 0)", expr = quote(out <- rep(arg1, length.out = 0)), args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3)}), outputType = quote(double(1))),
    list(name = "rep(vector double, length.out = first arg)", expr = quote(out <- rep(arg1, length.out = arg2)), args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(c(7, 8))}), outputType = quote(double(1)), expectWarnings = list("R eval" = 'Expected warning: vector each', "R run" = "Expected warning: vector each")),
    list(name = "rep(vector double, length.out = scalar from vectors)", expr = quote(out <- rep(arg1, length.out = sum(arg2))), args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(c(7, 8))}), outputType = quote(double(1))),
    #21
    list(name = "rep(vector double, times to ignore, length.out = scalar from vectors)", expr = quote(out <- rep(arg1, times = 5, length.out = sum(arg2))), args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(c(7, 8))}), outputType = quote(double(1))),
    list(name = "rep(vector double, each, length.out = scalar from vectors)", expr = quote(out <- rep(arg1, each = 3, length.out = sum(arg2))), args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(c(7, 8))}), outputType = quote(double(1))),
    list(name = "rep(vector double, times to ignore, each, length.out = scalar from vectors)", expr = quote(out <- rep(arg1, each = 3, times = 10, length.out = sum(arg2))), args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(c(7, 8))}), outputType = quote(double(1))),
    list(name = "rep(vector double, times to ignore, each = first arg, length.out = first arg)", expr = quote(out <- rep(arg1, each = arg3, times = 10, length.out = arg2)), args = list(arg1 = quote(double(1)), arg2 = quote(double(1)), arg3 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(c(7, 8)); arg3 <- as.numeric(4:5)}), outputType = quote(double(1)), expectWarnings = list("R eval" = 'Expected warning: vector each', "R run" = "Expected warning: vector each")),

    ## x, times expressions
    list(name = "rep(vector double expression, expression)", expr = quote(out <- rep(exp(arg1), arg2^2)), args = list(arg1 = quote(double(1)), arg2 = quote(integer())),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- 4}), outputType = quote(double(1))),
    ##26
    list(name = "rep(vector double expression, non-scalar expression)", expr = quote(out <- rep(exp(arg1), sum(arg2^2))), args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(c(2,3))}), outputType = quote(double(1))),
        

    list(name = "rep(matrix, 3)", expr = quote(out <- rep(arg1, 3)), args = list(arg1 = quote(double(2))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:9), nrow = 3)}),outputType = quote(double(1))),

    list(name = "rep(vector, vector)", expr = quote(out <- rep(arg1, arg2)), args = list(arg1 = quote(double(1)), arg2 = quote(double(1))), 
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(2:4)}), outputType = quote(double(1))),

     list(name = "rep(vector double, 3) in expression", expr = quote(out <- log(rep(arg1, 3))^2 + c(arg1, arg1, arg1)), args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3)}), outputType = quote(double(1)))

)

diagTests <- list(
    ## could add some where non-scalar inputs are copied in order to get different map behavior
    ## 
    ##1
    ## diag(scalar)
    list(name = "diag(scalar)", expr = quote(out <- diag(arg1)), args = list(arg1 = quote(double(0))),
         setArgVals = quote({arg1 <- 3}), outputType = quote(double(2))),
    list(name = "diag(scalar expression)", expr = quote(out <- diag(arg1 + arg2)), args = list(arg1 = quote(double(0)), arg2 = quote(double(0))),
         setArgVals = quote({arg1 <- 3; arg2 <- 2}), outputType = quote(double(2))),
    list(name = "diag(scalar-producing vector expression)", expr = quote(out <- diag(sum(arg1))), args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3);}), outputType = quote(double(2))),    
    list(name = "diag(scalar) with expr", expr = quote(out <- exp(diag(arg1)) + arg2), args = list(arg1 = quote(double(0)), arg2 = quote(double(2))),
         setArgVals = quote({arg1 <- 3; arg2 = matrix(1:9, nrow = 3)}), outputType = quote(double(2))),
    list(name = "diag(0)", expr = quote(out <- diag(arg1)), args = list(arg1 = quote(double(0))),
         setArgVals = quote({arg1 <- 0}), outputType = quote(double(2))),
    ## 6
    ## diag(vector)
    list(name = "diag(vector)", expr = quote(out <- diag(arg1)), args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3)}), outputType = quote(double(2))),
    list(name = "diag(vector expression)", expr = quote(out <- diag(arg1 + arg2)), args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(c(10,20, 30))}), outputType = quote(double(2))),
    list(name = "diag(vector) with expression", expr = quote(out <- exp(diag(arg1)) + arg2), args = list(arg1 = quote(double(1)), arg2 = quote(double(2))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- matrix(as.numeric(11:19), nrow = 3)}), outputType = quote(double(2))), 

    ## diag(matrix)
    list(name = "diag(square matrix)", expr = quote(out <- diag(arg1)), args = list(arg1 = quote(double(2))),
         setArgVals = quote({arg1 <- matrix(rnorm(25), nrow = 5)}), outputType = quote(double(1))),
    list(name = "diag(square matrix from expression)", expr = quote(out <- diag(exp(arg1) + arg2)), args = list(arg1 = quote(double(2)), arg2 = quote(double(2))),
         setArgVals = quote({arg1 <- matrix(rnorm(25), nrow = 5); arg2 <- matrix(1:25, nrow = 5)}), outputType = quote(double(1))),
    ## 11
    list(name = "diag(non-square matrix)", expr = quote(out <- diag(arg1)), args = list(arg1 = quote(double(2))),
         setArgVals = quote({arg1 <- matrix(rnorm(12), nrow = 3)}), outputType = quote(double(1))),

    list(name = "diag(square matrix) <-", expr = quote({diag(arg1) <- arg2; out <- arg1}), args = list(arg1 = quote(double(2)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- matrix(rnorm(25), nrow = 5); arg2 <- as.numeric(101:105)}), outputType = quote(double(2))),
    list(name = "copy, then diag(square matrix) <-", expr = quote({A1 <- arg1; diag(A1) <- arg2; out <- A1}), args = list(arg1 = quote(double(2)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- matrix(rnorm(25), nrow = 5); arg2 <- as.numeric(101:105)}), outputType = quote(double(2))),
    list(name = "diag(square matrix)[subset]", expr = quote(out <- diag(arg1)[2:4]), args = list(arg1 = quote(double(2))),
         setArgVals = quote({arg1 <- matrix(rnorm(25), nrow = 5)}), outputType = quote(double(1))),
    list(name = "diag(square matrix)[subset] <-", expr = quote({diag(arg1)[2:4] <- arg2[1:3]; out <- arg1}), args = list(arg1 = quote(double(2)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- matrix(rnorm(25), nrow = 5); arg2 <- rnorm(3)}), outputType = quote(double(2))),
    ## 16
    ## aliasing
    list(name = "diag(matrix)[3:5] <- diag(matrix[1:3])", expr = quote({diag(arg1)[3:5] <- diag(arg1)[1:3]; out <- arg1}), args = list(arg1 = quote(double(2))),
         setArgVals = quote({arg1 <- matrix(rnorm(25), nrow = 5)}), outputType = quote(double(2)))
    
)

recyclingRuleTests <- list(
    list(name = "dnorm all vector", expr = quote(out <- dnorm(arg1, arg2, arg3)), args = list(arg1 = quote(double(1)), arg2 = quote(double(1)), arg3 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(4:2); arg3 <- as.numeric(c(2, 3, 5, 8))}), outputType = quote(double(1))),
    list(name = "dnorm case 1", expr = quote(out <- dnorm(arg1, arg2, arg3)), args = list(arg1 = quote(double(0)), arg2 = quote(double(1)), arg3 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(2); arg2 <- as.numeric(4:2); arg3 <- as.numeric(c(2, 3, 5, 8))}), outputType = quote(double(1))),
    list(name = "dnorm case 2", expr = quote(out <- dnorm(arg1[1], arg2, arg3)), args = list(arg1 = quote(double(1)), arg2 = quote(double(1)), arg3 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:2); arg2 <- as.numeric(4:1); arg3 <- as.numeric(c(2, 3, 5, 8))}), outputType = quote(double(1))),
    list(name = "dnorm case 3", expr = quote(out <- dnorm(arg1, arg2, arg3)), args = list(arg1 = quote(double(1)), arg2 = quote(double(0)), arg3 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:2); arg2 <- as.numeric(3.5); arg3 <- as.numeric(c(2, 3, 5, 8))}), outputType = quote(double(1))),
    list(name = "dnorm case 4", expr = quote(out <- dnorm(arg1, arg2, arg3)), args = list(arg1 = quote(double(1)), arg2 = quote(double(0)), arg3 = quote(double(0))),
         setArgVals = quote({arg1 <- as.numeric(1:2); arg2 <- as.numeric(3.5); arg3 <- as.numeric(4.1)}), outputType = quote(double(1))),
    list(name = "dnorm case 4 [with expressions]", expr = quote(out <- (dnorm(arg1 + 1.5, arg2 + 1.5, arg3) + 1)^2), args = list(arg1 = quote(double(1)), arg2 = quote(double(0)), arg3 = quote(double(0))),
         setArgVals = quote({arg1 <- as.numeric(1:2); arg2 <- as.numeric(3.5); arg3 <- as.numeric(4.1)}), outputType = quote(double(1))),
    list(name = "dlogis all vector",
         expr = quote(out <- dlogis(arg1, arg2, arg3)),
         args = list(arg1 = quote(double(1)), arg2 = quote(double(1)), arg3 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(4:2); arg3 <- as.numeric(c(2, 3, 5, 8))}),
         outputType = quote(double(1)))
)

pRecyclingRuleTests <- list(
    list(name = "plogis case 1", expr = quote(out <- plogis(arg1, arg2, arg3)),
         args = list(arg1 = quote(double(1)), arg2 = quote(double(1)), arg3 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2  <- as.numeric(1:4); arg3 <- as.numeric(1:5);}),
         outputType = quote(double(1)))
)

qRecyclingRuleTests <- list(
    list(name = "qlogis case 1", expr = quote(out <- qlogis(arg1, arg2, arg3)),
         args = list(arg1 = quote(double(1)), arg2 = quote(double(1)), arg3 = quote(double(1))),
         setArgVals = quote({arg1 <- seq(.1, .4, length = 3); arg2  <- as.numeric(1:4); arg3 <- as.numeric(1:5);}),
         outputType = quote(double(1)))
)

rRecyclingRuleTests <- list(
 list(name = "rnorm case 1", expr = quote(out <- rnorm(5, arg1, arg2)),
             args = list(arg1 = quote(double(1)), arg2 = quote(double(0))),
             setArgVals = quote({arg1 <- as.numeric(1:2); arg2 <- as.numeric(3.5);}),
             outputType = quote(double(1))),
 list(name = "rnorm case 2 (with expressions)", expr = quote(out <- (1 + rnorm(5, arg1 + 1, arg2/2))^2),
             args = list(arg1 = quote(double(1)), arg2 = quote(double(0))),
             setArgVals = quote({arg1 <- as.numeric(1:2); arg2 <- as.numeric(3.5);}),
             outputType = quote(double(1))),
 list(name = "rnorm case 3 (with assignment block)", expr = quote({out <- numeric(10); out[2:6] <- (1 + rnorm(5, arg1 + 1, arg2/2))^2}),
             args = list(arg1 = quote(double(1)), arg2 = quote(double(0))),
             setArgVals = quote({arg1 <- as.numeric(1:2); arg2 <- as.numeric(3.5);}),
             outputType = quote(double(1))),
 list(name = "rbinom case 1", expr = quote(out <- rbinom(5, prob = arg1, size = arg2)),
             args = list(arg1 = quote(double(1)), arg2 = quote(integer(1))),
             setArgVals = quote({arg1 <- seq(.1, .4, length = 10); arg2 <- 1:3}),
             outputType = quote(double(1)), checkEqual = TRUE),
 list(name = "exp all vector", expr = quote(out <- rexp(arg1, arg2)),
             args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
             setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(4:2);}),
             outputType = quote(double(1))),
 list(name = "rexp_nimble case 1", expr = quote(out <- rexp_nimble(5, arg1)),
             args = list(arg1 = quote(double(1))),
             setArgVals = quote({arg1 <- seq(.1, .4, length = 10)}),
             outputType = quote(double(1))),

 list(name = "rexp case 1", expr = quote(out <- rexp(5, arg1)),
             args = list(arg1 = quote(double(1))),
             setArgVals = quote({arg1 <- seq(.1, .4, length = 10)}),
             outputType = quote(double(1))),

 list(name = "rnbinom case 1", expr = quote(out <- rnbinom(5, prob = arg1, size = arg2)),
             args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
             setArgVals = quote({arg1 <- seq(.1, .4, length = 10); arg2 <- seq(4, 1, length = 10)}),
             outputType = quote(double(1)), checkEqual = TRUE),

 list(name = "rpois case 1", expr = quote(out <- rpois(5, arg1)),
             args = list(arg1 = quote(double(1))),
             setArgVals = quote({arg1 <- seq(.1, .4, length = 10)}),
             outputType = quote(double(1)), checkEqual = TRUE),

 list(name = "rchisq case 1", expr = quote(out <- rchisq(5, arg1)),
             args = list(arg1 = quote(integer(1))),
             setArgVals = quote({arg1 <- 1:10}),
             outputType = quote(double(1))),

 list(name = "rbeta case 1", expr = quote(out <- rbeta(5, arg1, arg2)),
             args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
             setArgVals = quote({arg1 <- seq(0.1, 0.4, length = 10); arg2 <- seq(0.9, 0.7, length = 3)}),
             outputType = quote(double(1))),

 list(name = "rgamma case 1", expr = quote(out <- rgamma(5, arg1, arg2)),
             args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
             setArgVals = quote({arg1 <- seq(0.1, 0.4, length = 10); arg2 <- seq(0.9, 0.7, length = 3)}),
             outputType = quote(double(1))),

  list(name = "rinvgamma case 1", expr = quote(out <- rinvgamma(5, arg1, arg2)),
             args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
             setArgVals = quote({arg1 <- seq(0.1, 0.4, length = 10); arg2 <- seq(0.9, 0.7, length = 3)}),
             outputType = quote(double(1))),

 list(name = "rlnorm case 1", expr = quote(out <- rlnorm(5, arg1, arg2)),
             args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
             setArgVals = quote({arg1 <- seq(0.1, 0.4, length = 10); arg2 <- seq(0.9, 0.7, length = 3)}),
             outputType = quote(double(1))),

 list(name = "rlogis case 1", expr = quote(out <- rlogis(5, arg1, arg2)),
             args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
             setArgVals = quote({arg1 <- seq(0.1, 0.4, length = 10); arg2 <- seq(0.9, 0.7, length = 3)}),
             outputType = quote(double(1))),

 list(name = "runif case 1", expr = quote(out <- runif(5, arg1, arg2)),
             args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
             setArgVals = quote({arg1 <- seq(0.1, 0.4, length = 10); arg2 <- seq(0.9, 0.7, length = 3)}),
             outputType = quote(double(1))),

 list(name = "rweibull case 1", expr = quote(out <- rweibull(5, arg1, arg2)),
             args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
             setArgVals = quote({arg1 <- seq(0.1, 0.4, length = 10); arg2 <- seq(0.9, 0.7, length = 3)}),
             outputType = quote(double(1))),

 list(name = "rt case 1", expr = quote(out <- rt(5, arg1)),
             args = list(arg1 = quote(integer(1)), arg2 = quote(double(1))),
             setArgVals = quote({arg1 <- 5:14; arg2 <- seq(0.9, 0.7, length = 3)}),
      outputType = quote(double(1)))

)

seqTests <- list(
    ##1
    list(name = "1:5", expr = quote(out <- 1:5), args = list(),
         setArgVals = quote({}), outputType = quote(double(1)), checkEqual = TRUE),
    list(name = "seq(.1, 10, by = .1)", expr = quote(out <- seq(.1, 10, by = .1)), args = list(),
         setArgVals = quote({}), outputType = quote(double(1))),
    list(name = "seq(.1, 10, length.out = 11)", expr = quote(out <- seq(.1, 10, length.out = 11)), args = list(),
         setArgVals = quote({}), outputType = quote(double(1))),
    list(name = "seq(.1, by = 10, length.out = 11)", expr = quote(out <- seq(.1, by = 10, length.out = 11)), args = list(),
         setArgVals = quote({}), outputType = quote(double(1))),
    list(name = "seq(.1, 10, length.out = 11) in expression", expr = quote(out <- log(seq(.1, 10, length.out = 11)) + 2 + rep(1, 11)), args = list(),
         setArgVals = quote({}), outputType = quote(double(1)))
)

nonSeqIndexTests <- list(
    ##1
    list(name = "non-sequential indexing: out <- arg1[arg2, arg3]", expr = quote(out <- arg1[arg2, arg3]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:25), nrow = 5);
                             arg2 <- c(2, 4);
                             arg3 <- c(1, 3, 4)}),
         outputType = quote(double(2))),
    list(name = "non-sequential indexing: out <- arg1[arg2, 2:4]", expr = quote(out <- arg1[arg2, 2:4]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:25), nrow = 5);
                             arg2 <- c(2, 4);
                             arg3 <- c(1, 3, 4)}),
         outputType = quote(double(2))),
    list(name = "non-sequential indexing: out <- arg1[2:4, arg3]", expr = quote(out <- arg1[2:4, arg3]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:25), nrow = 5);
                             arg2 <- c(2, 4);
                             arg3 <- c(1, 3, 4)}),
         outputType = quote(double(2))),
    list(name = "non-sequential indexing: out <- arg1[2, arg3, drop = FALSE]", expr = quote(out <- arg1[2, arg3, drop = FALSE]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:25), nrow = 5);
                             arg2 <- c(2, 4);
                             arg3 <- c(1, 3, 4)}),
         outputType = quote(double(2))),
    list(name = "non-sequential indexing: out <- arg1[2, arg3]", expr = quote(out <- arg1[2, arg3]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:25), nrow = 5);
                             arg2 <- c(2, 4);
                             arg3 <- c(1, 3, 4)}),
         outputType = quote(double(1))),
    ##6
    list(name = "non-sequential indexing: out <- arg1[arg2, arg3] with scalar arg2", expr = quote(out <- arg1[arg2, arg3]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(0)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:25), nrow = 5);
                             arg2 <- 2;
                             arg3 <- c(1, 3, 4)}),
         outputType = quote(double(1))),
    list(name = "non-sequential indexing: out <- arg1[arg2, arg3] with scalar arg3", expr = quote(out <- arg1[arg2, arg3]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(0))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:25), nrow = 5);
                             arg2 <- c(2, 4);
                             arg3 <- 3}),
         outputType = quote(double(1))),
    list(name = "non-sequential indexing: out <- arg1[arg2]", expr = quote(out <- arg1[arg2]),
         args = list(arg1 = quote(double(1)), arg2 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:25), nrow = 5);
                             arg2 <- c(2, 4)}),
         outputType = quote(double(1))),
    list(name = "non-sequential indexing: out[arg2, arg3] <- arg1", expr = quote({out <- matrix(100, nrow = 5, ncol = 5); out[arg2, arg3] <- arg1}),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:6), nrow = 2);
                             arg2 <- c(2, 4);
                             arg3 <- c(1, 3, 4)}),
         outputType = quote(double(2))),
    list(name = "non-sequential indexing: out[2:3, arg3] <- arg1", expr = quote({out <- matrix(100, nrow = 5, ncol = 5); out[2:3, arg3] <- arg1}),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:6), nrow = 2);
                             arg2 <- c(2, 4);
                             arg3 <- c(1, 3, 4)}),
         outputType = quote(double(2))),
##11
    list(name = "non-sequential indexing: out[arg2, 3:5] <- arg1", expr = quote({out <- matrix(100, nrow = 5, ncol = 5); out[arg2, 3:5] <- arg1}),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:6), nrow = 2);
                             arg2 <- c(2, 4);
                             arg3 <- c(1, 3, 4)}),
         outputType = quote(double(2))),
    list(name = "non-sequential indexing: out <- log(arg1[arg2, arg3]) + 1 [in expression]", expr = quote(out <- log(arg1[arg2, arg3]) + 1),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:25), nrow = 5);
                             arg2 <- c(2, 4);
                             arg3 <- c(1, 3, 4)}),
         outputType = quote(double(2))),
    list(name = "non-sequential indexing: out <- arg1[arg2 - 1, arg3 + 1] [indices in expressions]", expr = quote(out <- arg1[arg2 - 1, arg3 + 1]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:25), nrow = 5);
                             arg2 <- c(2, 4);
                             arg3 <- c(1, 3, 4)}),
         outputType = quote(double(2))),
    list(name = "non-sequential indexing: out[arg2] <- arg1", expr = quote({out <- numeric(5); out[1:5] <- 100; out[arg2] <- arg1}),
         args = list(arg1 = quote(double(1)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- as.numeric(1:2);
                             arg2 <- c(2, 4);
                             arg3 <- c(1, 3, 4)}),
         outputType = quote(double(1))),
    list(name = "non-sequential indexing: out[arg2] <- 200", expr = quote({out <- numeric(5); out[1:5] <- 100; out[arg2] <- 200}),
         args = list(arg1 = quote(double(1)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- as.numeric(1:2);
                             arg2 <- c(2, 4);
                             arg3 <- c(1, 3, 4)}),
         outputType = quote(double(1))),
    list(name = "non-sequential indexing: out[2, arg3] <- arg1 (row matrix)", expr = quote({out <- matrix(100, nrow = 5, ncol = 5); out[2, arg3] <- arg1}),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:3), nrow = 1);
                             arg2 <- c(2, 4);
                             arg3 <- c(1, 3, 4)}),
         outputType = quote(double(2))),
    list(name = "non-sequential indexing: out[2, arg3] <- arg1 (col matrix)", expr = quote({out <- matrix(100, nrow = 5, ncol = 5); out[2, arg3] <- arg1}),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:3), ncol = 1);
                             arg2 <- c(2, 4);
                             arg3 <- c(1, 3, 4)}),
         outputType = quote(double(2))),
    list(name = "non-sequential indexing: out <- arg1[, arg3]", expr = quote(out <- arg1[, 3:5]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:25), nrow = 5);
                             arg2 <- c(2, 4);
                             arg3 <- c(1, 3, 4)}),
         outputType = quote(double(2))),
    list(name = "non-sequential indexing: out <- arg1[arg2, ]", expr = quote(out <- arg1[, 3:5]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:25), nrow = 5);
                             arg2 <- c(2, 4);
                             arg3 <- c(1, 3, 4)}),
         outputType = quote(double(2)))
)


indexChainTests <- list(
    list(name = "block chaining 1", expr = quote(out <- arg1[arg2, arg3][2:3, 2:4]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:100), nrow = 10);
                             arg2 <- c(2, 3, 5);
                             arg3 <- c(4, 6, 7, 8)}),
         outputType = quote(double(2))),

    list(name = "block chaining 1b", expr = quote(out <- arg1[, arg3][2:3, 2:4]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:100), nrow = 10);
                             arg2 <- c(2, 3, 5);
                             arg3 <- c(4, 6, 7, 8)}),
         outputType = quote(double(2))),

    list(name = "block chaining 1c", expr = quote(out <- arg1[arg2, ][2:3, 2:4]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:100), nrow = 10);
                             arg2 <- c(2, 3, 5);
                             arg3 <- c(4, 6, 7, 8)}),
         outputType = quote(double(2))),

    list(name = "block chaining 1d", expr = quote(out <- arg1[arg2, arg3][, 2:4]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:100), nrow = 10);
                             arg2 <- c(2, 3, 5);
                             arg3 <- c(4, 6, 7, 8)}),
         outputType = quote(double(2))),

    list(name = "block chaining 1e", expr = quote(out <- arg1[arg2, arg3][2:3, ]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:100), nrow = 10);
                             arg2 <- c(2, 3, 5);
                             arg3 <- c(4, 6, 7, 8)}),
         outputType = quote(double(2))),


    list(name = "block chaining 2", expr = quote(out <- arg1[2:8, 3:6][arg2, arg3]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:100), nrow = 10);
                             arg2 <- c(2, 3, 5);
                             arg3 <- c(2, 4)}),
         outputType = quote(double(2))),

    list(name = "block chaining 3", expr = quote(out <- arg1[arg2, arg3][arg4, arg5]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1)), arg4 = quote(integer(1)), arg5 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:100), nrow = 10);
                             arg2 <- c(2, 3, 5, 7);
                             arg3 <- c(2, 4, 6, 8, 9);
                             arg4 <- c(2, 3);
                             arg5 <- c(1, 4, 5)}),
         outputType = quote(double(2))),

    list(name = "block chaining 4", expr = quote(out <- arg1[1, arg3][arg5]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1)), arg4 = quote(integer(1)), arg5 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:100), nrow = 10);
                             arg2 <- c(2, 3, 5, 7);
                             arg3 <- c(2, 4, 6, 8);
                             arg4 <- c(2, 3);
                             arg5 <- c(1, 4)}),
         outputType = quote(double(1))),

    list(name = "block chaining 4b", expr = quote(out <- arg1[1, arg3, drop = FALSE][1, arg5]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1)), arg4 = quote(integer(1)), arg5 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:100), nrow = 10);
                             arg2 <- c(2, 3, 5, 7);
                             arg3 <- c(2, 4, 6, 8);
                             arg4 <- c(2, 3);
                             arg5 <- c(1, 4)}),
         outputType = quote(double(1))),


    list(name = "block chaining 4c", expr = quote(out <- arg1[arg2, 2][arg5]),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1)), arg4 = quote(integer(1)), arg5 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:100), nrow = 10);
                             arg2 <- c(2, 3, 5, 7);
                             arg3 <- c(2, 4, 6, 8);
                             arg4 <- c(2, 3);
                             arg5 <- c(1, 4)}),
         outputType = quote(double(1))),

    list(name = "block chaining assignment 1", expr = quote({out <- matrix(-1, nrow = 10, ncol = 10);
                                                             out[arg2, arg3][2:3, 2:4] <- arg1[arg2, arg3][2:3, 2:4]}),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:100), nrow = 10);
                             arg2 <- c(2, 3, 5);
                             arg3 <- c(4, 6, 7, 8)}),
         outputType = quote(double(2))),


    list(name = "block chaining assignment 1b", expr = quote({out <- matrix(-1, nrow = 10, ncol = 10);
                                                              out[, arg3][2:3, 2:4] <- arg1[, arg3][2:3, 2:4]}),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:100), nrow = 10);
                             arg2 <- c(2, 3, 5);
                             arg3 <- c(4, 6, 7, 8)}),
         outputType = quote(double(2))),

    list(name = "block chaining assignment 1c", expr = quote({out <- matrix(-1, nrow = 10, ncol = 10);
                                                              out[arg2, ][2:3, 2:4] <- arg1[arg2, ][2:3, 2:4]}),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:100), nrow = 10);
                             arg2 <- c(2, 3, 5);
                             arg3 <- c(4, 6, 7, 8)}),
         outputType = quote(double(2))),

    list(name = "block chaining assignment 1d", expr = quote({out <- matrix(-1, nrow = 10, ncol = 10);
                                                              out[arg2, arg3][, 2:4] <- arg1[arg2, arg3][, 2:4]}),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:100), nrow = 10);
                             arg2 <- c(2, 3, 5);
                             arg3 <- c(4, 6, 7, 8)}),
         outputType = quote(double(2))),

    list(name = "block chaining assignment 1e", expr = quote({out <- matrix(-1, nrow = 10, ncol = 10);
                                                              out[arg2, arg3][2:3, ] <- arg1[arg2, arg3][2:3, ]}),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:100), nrow = 10);
                             arg2 <- c(2, 3, 5);
                             arg3 <- c(4, 6, 7, 8)}),
         outputType = quote(double(2))),

    list(name = "block chaining assignment 2", expr = quote({out <- matrix(-1, nrow = 10, ncol = 10);
                                                             out[2:8, 3:6][arg2, arg3] <- arg1[2:8, 3:6][arg2, arg3]}),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:100), nrow = 10);
                             arg2 <- c(2, 3, 5);
                             arg3 <- c(2, 4)}),
         outputType = quote(double(2))),

    list(name = "block chaining assignment 3", expr = quote({out <- matrix(-1, nrow = 10, ncol = 10);
                                                             out[arg2, arg3][arg4, arg5] <- arg1[arg2, arg3][arg4, arg5]}),
         args = list(arg1 = quote(double(2)), arg2 = quote(integer(1)), arg3 = quote(integer(1)), arg4 = quote(integer(1)), arg5 = quote(integer(1))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:100), nrow = 10);
                             arg2 <- c(2, 3, 5, 7);
                             arg3 <- c(2, 4, 6, 8);
                             arg4 <- c(2, 3);
                             arg5 <- c(1, 4)}),
         outputType = quote(double(2)))
)

logicalTests <- list(
    list(name = "create boolean vector", expr = quote(out <- arg1 > 3 & arg1 < 6),
         args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- seq(1, 8, length = 100)}),
         outputType = quote(logical(1))),

    list(name = "create boolean vector with expressions", expr = quote(out <- arg1 > 3 & arg1 + 1 < 6),
         args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- seq(1, 8, length = 100)}),
         outputType = quote(logical(1))),

    list(name = "use boolean vector with expressions",
         expr = quote({out <- arg1 > 3 & arg1 + 1 < 6; out <- out | arg1 > 7}),
         args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- seq(1, 8, length = 100)}),
         outputType = quote(logical(1))),

    list(name = "index from boolean vector 1",
         expr = quote({out <- arg1[arg2 < 5]}),
         args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- seq(1, 8, length = 100); arg2 <- seq(2, 9, length = 100)}),
         outputType = quote(double(1))),

    list(name = "index from boolean vector 2 (in expression)",
         expr = quote({out <- (arg1[arg2 < 5]^2) + 1}),
         args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- seq(1, 8, length = 100); arg2 <- seq(2, 9, length = 100)}),
         outputType = quote(double(1))),

    list(name = "index from boolean vector 3 (2D)",
         expr = quote({out <- arg1[arg2 < 5, arg2 > 4]}),
         args = list(arg1 = quote(double(2)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- matrix(seq(1, 8, length = 10000), nrow = 100); arg2 <- seq(2, 9, length = 100)}),
         outputType = quote(double(2))),

    list(name = "index from boolean vector 3 (2D with mixed types)",
         expr = quote({out <- arg1[arg2 < 5, 30:50]}),
         args = list(arg1 = quote(double(2)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- matrix(seq(1, 8, length = 10000), nrow = 100); arg2 <- seq(2, 9, length = 100)}),
         outputType = quote(double(2))),


    list(name = "index assignment from boolean vector 1",
         expr = quote({out <- rep(100, length(arg1)); out[arg2 < 5] <- (arg1[arg2 < 5]^2) + 1}),
         args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- seq(1, 8, length = 100); arg2 <- seq(2, 9, length = 100)}),
         outputType = quote(double(1))),

    list(name = "index assignment from boolean vector 2 (2D)",
         expr = quote({out <- matrix(rep(100, length(arg1)), nrow = dim(arg1)[1]); out[arg2 < 5, arg2 > 4] <- (arg1[arg2 < 5, arg2 > 4]^2) + 1}),
         args = list(arg1 = quote(double(2)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- matrix(seq(1, 8, length = 10000), nrow = 100); arg2 <- seq(2, 9, length = 100)}),
         outputType = quote(double(2)), checkEqual = TRUE), ## small numerical differences

    list(name = "index assignment from boolean vector 2 (2D, mixed types)",
         expr = quote({out <- matrix(rep(100, length(arg1)), nrow = dim(arg1)[1]); out[arg2 < 5, 30:40] <- (arg1[arg2 < 5, 30:40]^2) + 1}),
         args = list(arg1 = quote(double(2)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- matrix(seq(1, 8, length = 10000), nrow = 100); arg2 <- seq(2, 9, length = 100)}),
         outputType = quote(double(2)), checkEqual = TRUE)
)

anyNaTests <- list(
    list(name = "use any_na", expr = quote(out <- any_na(arg1)),
         args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- rnorm(5); arg1[2] <- NA}),
         outputType = quote(logical(0)))
)

returnTests <- list(
    list(name = "return(rnorm scalar)",
         expr = quote({}),
         return = quote(return(rnorm(1))),
         args = list(),
         setArgVals = quote({}),
         outputType = quote(double())),
    list(name = "return(rnorm vector)",
         expr = quote({}),
         return = quote(return(rnorm(4))),
         args = list(),
         setArgVals = quote({}),
         outputType = quote(double(1))),
    list(name = "return(rep(...))",
         expr = quote({}),
         return = quote(return(rep(1.23, 4))),
         args = list(),
         setArgVals = quote({}),
         outputType = quote(double(1))),
    list(name = "return(seq(...))",
         expr = quote({}),
         return = quote(return(seq(from = .1, to = .5, by = .15))),
         args = list(),
         setArgVals = quote({}),
         outputType = quote(double(1))),
    list(name = "return(A + B scalar)",
         expr = quote({A <- .1; B <- .2}),
         return = quote(return(A + B)),
         args = list(),
         setArgVals = quote({}),
         outputType = quote(double(0))),
    list(name = "return(A + B vector)",
         expr = quote({A <- rep(.1, 3); B <- rep(.2, 3)}),
         return = quote(return(A + B)),
         args = list(),
         setArgVals = quote({}),
         outputType = quote(double(1))) 
)

## Regression test for Issue #563
test_that('unary_function( inprod(vector1, vector2) ) compiles and works', {
    nfDef <- nimbleFunction(
        setup = function() {},
        run = function() {
            a <- rep(0, 5)
            b <- rep(0, 5)
            c <- step(inprod(a, b))
            return(c)
            returnType(integer())
        }
    )
    Rnf <- nfDef()
    ## safe use of try followed by expectation of type
    Cnf <- try(compileNimble(Rnf))
    expect_false(inherits(Cnf, 'try-error'),
                 info = 'step(inprod(a, b)) does not compile')
    expect_equal(Cnf$run(),
                 Rnf$run(),
                 info = 'step(inprod(a, b)) compiles but gives wrong answer')
}
)

simpleCopyTests <- list(
    list(name = "1d copy",
         expr = quote(out <- arg1),
         args = list(arg1 = quote(double(1))),
         setArgVals = quote(arg1 <- as.numeric(1:5)),
         outputType = quote(double(1))),
    list(name = "2d copy",
         expr = quote(out <- arg1),
         args = list(arg1 = quote(double(2))),
         setArgVals = quote(arg1 <- matrix(as.numeric(1:20), nrow = 5)),
         outputType = quote(double(2))),
    list(name = "3d copy",
         expr = quote(out <- arg1),
         args = list(arg1 = quote(double(3))),
         setArgVals = quote(arg1 <- array(as.numeric(1:(3*5*7)), dim = c(3, 5, 7))),
         outputType = quote(double(3))),
    list(name = "4d copy", ## bug fixed from Issue #834
         expr = quote(out <- arg1),
         args = list(arg1 = quote(double(4))),
         setArgVals = quote(arg1 <- array(as.numeric(1:(3*5*7*9)), dim = c(3, 5, 7, 9))),
         outputType = quote(double(4))),
    list(name = "5d copy",
         expr = quote(out <- arg1),
         args = list(arg1 = quote(double(5))),
         setArgVals = quote(arg1 <- array(as.numeric(1:(3*5*7*9*11)), dim = c(3, 5, 7, 9, 11))),
         outputType = quote(double(5)))
)

higherDimBlockTests <- list(
    list(name = "3x4x1x1 block",
         expr = quote(out <- arg1[2:4, 3:6, 2, 4] + 2),
         args = list(arg1 = quote(double(4))),
         setArgVals = quote(arg1 <- array(as.numeric(1:(5*7*3*6)), dim = c(5, 7, 3, 6))),
         outputType = quote(double(2))),
    list(name = "1x4x1x3 block",
         expr = quote(out <- arg1[2, 3:6, 5, 2:4] + 2),
         args = list(arg1 = quote(double(4))),
         setArgVals = quote(arg1 <- array(as.numeric(1:(5*7*9*6)), dim = c(5, 7, 9, 6))),
         outputType = quote(double(2))),
    list(name = "1x4x1x1 block",
         expr = quote(out <- arg1[1, 3:6, 2, 4] + 2),
         args = list(arg1 = quote(double(4))),
         setArgVals = quote(arg1 <- array(as.numeric(1:(5*7*3*6)), dim = c(5, 7, 3, 6))),
         outputType = quote(double(1))),
    list(name = "3x4x1x1x1 block",
         expr = quote(out <- arg1[2:4, 3:6, 2, 4, 5] + 2),
         args = list(arg1 = quote(double(5))),
         setArgVals = quote(arg1 <- array(as.numeric(1:(5*7*3*6*7)), dim = c(5, 7, 3, 6, 7))),
         outputType = quote(double(2))),
    list(name = "1x4x1x3x1 block",
         expr = quote(out <- arg1[3, 3:6, 2, 4:6, 5] + 2),
         args = list(arg1 = quote(double(5))),
         setArgVals = quote(arg1 <- array(as.numeric(1:(5*7*3*8*7)), dim = c(5, 7, 3, 8, 7))),
         outputType = quote(double(2))),
   list(name = "1x1x1x3x1 block",
         expr = quote(out <- arg1[3, 4, 2, 4:6, 5] + 2),
         args = list(arg1 = quote(double(5))),
         setArgVals = quote(arg1 <- array(as.numeric(1:(5*7*3*8*7)), dim = c(5, 7, 3, 8, 7))),
         outputType = quote(double(1))),
    list(name = "1x4x1x3x1 block from input ranges",
         expr = quote(out <- arg1[3, arg2:arg3, 2, 4:6, 5] + 2),
         args = list(arg1 = quote(double(5)),
                     arg2 = quote(integer()),
                     arg3 = quote(integer())),
         setArgVals = quote(
         {
             arg1 <- array(as.numeric(1:(5*7*3*8*7)), dim = c(5, 7, 3, 8, 7))
             arg2 <- 3
             arg3 <- 6
         }),
         outputType = quote(double(2))),
     list(name = "block chaining from 1x4x1x3x1",
         expr = quote(out <- arg1[3, 3:6, 2, 4:7, 5][2:3, 2:4] + 2),
         args = list(arg1 = quote(double(5))),
         setArgVals = quote(arg1 <- array(as.numeric(1:(5*7*3*8*7)), dim = c(5, 7, 3, 8, 7))),
         outputType = quote(double(2)))
)

aliasTests <- list(
    list(name = "x <- rep(x, 2)",   ## Fails due to lack of eval with eigenBlock map
         expr = quote({
             x <- arg1
             x <- rep(x, 2)
             out <- x}),
         args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3)}),
         outputType = quote(double(1))),
    list(name = "x <- rep(x, length = 6)",   ## Fails due to lack of eval with eigenBlock map
         expr = quote({
             x <- arg1
             x <- rep(x, length = 6)
             out <- x}),
         args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3)}),
         outputType = quote(double(1))),
    list(name = "x <- rep(x, each = 2)",   ## Fails due to lack of eval with eigenBlock map
         expr = quote({
             x <- arg1
             x <- rep(x, each = 2)
             out <- x}),
         args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3)}),
         outputType = quote(double(1))),
    list(name = "x <- c(x, arg2)",   ## Fails due to lack of eval with eigenBlock map
         expr = quote({
             x <- arg1
             x <- c(x, arg2)
             out <- x}),
         args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(4:6)}),
         outputType = quote(double(1))),    
    list(name = "x[2:3] <- x[1:2]",   ## Fails due to lack of eval with eigenBlock map
         expr = quote({
             x <- arg1
             x[2:3] <- x[1:2]
             out <- x}),
         args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3)}),
         outputType = quote(double(1))),
    list(name = "x[c(2, 3)] <- x[1:2]",   ## Fails due to lack of eval with eigenBlock map
         expr = quote({
             x <- arg1
             x[c(2, 3)] <- x[1:2]
             out <- x}),
         args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3)}),
         outputType = quote(double(1))),
    list(name = "x[c(FALSE, TRUE, TRUE)] <- x[1:2]",   ## Fails due to lack of eval with eigenBlock map
         expr = quote({
             x <- arg1
             x[c(FALSE, TRUE, TRUE)] <- x[1:2]
             out <- x}),
         args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3)}),
         outputType = quote(double(1))),
    list(name = "x[2:3, 1:2] <- x[1:2, 1:2]",   ## Fails due to lack of eval with eigenBlock map
         expr = quote({
             x <- arg1
             x[2:3, 1:2] <- x[1:2, 1:2]
             out <- x}),
         args = list(arg1 = quote(double(2))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:9), nrow = 3)}),
         outputType = quote(double(2))),
    ## No blocking in >2D supported
    list(name = "x <- x[1:2]",   ## Fails due to x being resized before the Eigen assignment
         expr = quote({
             x <- arg1
             x <- x[1:2]
             out <- x}),
         args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3)}),
         outputType = quote(double(1))),
    list(name = "x <- t(x)",   ## Fails due to x being resized before the Eigen assignment
         expr = quote({
             x <- arg1
             x <- t(x)
             out <- x}),
         args = list(arg1 = quote(double(2))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:6), nrow = 2)}),
         outputType = quote(double(2))),
    list(name = "x <- t(x) + 1",   ## Fails due to x being resized before the Eigen assignment
         expr = quote({
             x <- arg1
             x <- t(x) + 1
             out <- x}),
         args = list(arg1 = quote(double(2))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:6), nrow = 2)}),
         outputType = quote(double(2))),
    list(name = "x <- dnorm(x, arg2) with length(x) < length(arg2)",   ## Fails due to x being resized before the Eigen assignment
         expr = quote({
             x <- arg1
             x <- dnorm(x, arg2) ## should handle recycling rule on x, but x is resized first
             out <- x}),
         args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 = 0.1; arg2 <- 0.1 * (1:6)}),
         outputType = quote(double(1)))
)

cTestsResults <- test_coreRfeature_batch(cTests, 'cTests') ##lapply(cTests, test_coreRfeature)
blockTestsResults <- test_coreRfeature_batch(blockTests, 'blockTests') ##lapply(blockTests, test_coreRfeature)
repTestsResults <- test_coreRfeature_batch(repTests, 'repTests') ## lapply(repTests, test_coreRfeature)
diagTestsResults <- test_coreRfeature_batch(diagTests, 'diagTests') ## lapply(diagTests, test_coreRfeature)
recyclingRuleTestsResults <- test_coreRfeature_batch(recyclingRuleTests, 'recyclingRuleTests') ## lapply(recyclingRuleTests, test_coreRfeature)
pRecyclingRuleTestsResults <- test_coreRfeature_batch(pRecyclingRuleTests, 'pRecyclingRuleTests')
qRecyclingRuleTestsResults <- test_coreRfeature_batch(qRecyclingRuleTests, 'qRecyclingRuleTests')
rRecyclingRuleTestsResults <- test_coreRfeature_batch(rRecyclingRuleTests, 'rRecyclingRuleTests') ## lapply(rRecyclingRuleTests, test_coreRfeature)
seqTestsResults <- test_coreRfeature_batch(seqTests, 'seqTests') ## lapply(seqTests, test_coreRfeature)
nonSeqIndexTestsResults <- test_coreRfeature_batch(nonSeqIndexTests, 'nonSeqIndexTests') ## lapply(nonSeqIndexTests, test_coreRfeature)
indexChainTestsResults <- test_coreRfeature_batch(indexChainTests, 'indexChainTests') ## lapply(indexChainTests, test_coreRfeature)
logicalTestsResults <- test_coreRfeature_batch(logicalTests, 'logicalTests') ## lapply(logicalTests, test_coreRfeature)
anyNaTestResults <- test_coreRfeature_batch(anyNaTests, 'anyNaTests') 
returnTestResults <- test_coreRfeature_batch(returnTests, 'returnTests') ## lapply(returnTests, test_coreRfeature)
simpleCopyTestResults <- test_coreRfeature_batch(simpleCopyTests, 'simpleCopyTests') 
higherDimBlockTestResults <- test_coreRfeature_batch(higherDimBlockTests, 'higherDimBlockTests') 
aliasTestResults <- test_coreRfeature_batch(aliasTests, 'aliasTests')

## basic seq_along test

test_that('seq_along works in nimbleFunctions', {
    nf <- nimbleFunction(
        run = function(x = double(1)) {
            for(i in seq_along(x))
                x[i] <- x[i]+1
            returnType(double(1))
            return(x)
        })
    cnf <- compileNimble(nf)
    x <- rnorm(5)
    expect_identical(nf(x), x+1, 'problem with uncompiled seq_along')
    expect_identical(nf(x), cnf(x), 'problem with compiled seq_along')
})

## Some tests of using coreR features in BUGS models

test_that('c(a, 1.1) in BUGS works', {
    mc <- nimbleCode({
        a ~ dnorm(0,1)
        b[1:2] <- c(a, 1.1)
    })
    
    m <- nimbleModel(mc, inits = list(a = 2))
    expect_identical(as.numeric(m$b), c(2, 1.1))
    m$b <- as.numeric(rep(NA, 2))
    cm <- compileNimble(m)
    cm$calculate()
    expect_identical(as.numeric(cm$b), c(2, 1.1))
}
)

##

test_that('c(1.2, 1.1) in BUGS works', {
    mc <- nimbleCode({
        b[1:2] <- c(1.2, 1.1)
    })
    m <- nimbleModel(mc)
    expect_identical(as.numeric(m$b), c(1.2, 1.1))
    m$b <- as.numeric(rep(NA, 2))
    cm <- compileNimble(m)
    cm$calculate()
    expect_identical(as.numeric(cm$b), c(1.2, 1.1))
}
)

##

test_that('rep(a, 2) in BUGS works', {
    mc <- nimbleCode({
        a ~ dnorm(0,1)
        b[1:2] <- rep(a, 2)
    })
    
    m <- nimbleModel(mc, inits = list(a = 1.2))
    expect_identical(as.numeric(m$b), c(1.2, 1.2))
    m$b <- as.numeric(rep(NA, 2))
    cm <- compileNimble(m)
    cm$calculate()
    expect_identical(as.numeric(cm$b), c(1.2, 1.2))
}
)

##

test_that('rep(1,2)  in BUGS works', {
    mc <- nimbleCode({
        b[1:2] <- rep(1, 2)
    })
    m <- nimbleModel(mc)
    expect_identical(as.numeric(m$b), rep(1, 2))
    m$b <- as.numeric(rep(NA, 2))
    cm <- compileNimble(m)
    cm$calculate()
    expect_identical(as.numeric(cm$b), rep(1, 2))
}
)


##

test_that('2:3   in BUGS works', {
    mc <- nimbleCode({
        b[1:2] <- 2:3 
    })
    m <- nimbleModel(mc)
    expect_equal(as.numeric(m$b), 2:3 )
    m$b <- as.numeric(rep(NA, 2))
    cm <- compileNimble(m)
    cm$calculate()
    expect_equal(as.numeric(cm$b), 2:3 )
}
)

##

test_that('seq(1.2, 2.3, length = 3) in BUGS works', {
    mc <- nimbleCode({
        b[1:3] <- seq(1.2, 2.3, length = 3)
    })
    m <- nimbleModel(mc)
    expect_identical(as.numeric(m$b), seq(1.2, 2.3, length = 3) )
    m$b <- as.numeric(rep(NA, 3))
    cm <- compileNimble(m)
    cm$calculate()
    expect_identical(as.numeric(cm$b), seq(1.2, 2.3, length = 3) )
}
)

##

test_that('diag(3) in BUGS works', {
    mc <- nimbleCode({
        b[1:3, 1:3] <- diag(3)
    })
    m <- nimbleModel(mc)
    expect_equal(m$b, diag(3))
    m$b <- matrix(100, nrow = 3, ncol = 3)
    cm <- compileNimble(m)
    cm$calculate()
    expect_identical(cm$b, diag(3))
}
)

## Tests of slices of objects with nDim > 2 as arguments to other functions
test_that('slice of 3d passed as 2d arg works', {
    f1 <- nimbleFunction(
        run = function(x = double(2)) {
            ans <- numeric(value = x, length = length(x))
            return(ans)
            returnType(double(1))
        })
    temporarilyAssignInGlobalEnv(f1)
    f2 <- nimbleFunction(
        run = function(x = double(3)) {
            ans <- f1(x[3, 2:4, 3:6])
            return(ans)
            returnType(double(1))
        })
    c12 <- compileNimble(f1, f2)
    x <- array(as.numeric(1:(4*6*8)), dim = c(4, 6, 8))
    expect_equal(f2(x), c12$f2(x))
})

test_that('slice of 3d copy works', {
    f1 <- nimbleFunction(
        run = function(x = double(3)) {
            ans <- x
            return(ans)
            returnType(double(3))
        })
    temporarilyAssignInGlobalEnv(f1)
    f2 <- nimbleFunction(
        run = function(x = double(3)) {
            ans <- f1(x[2:5, 2:4, 3:6])
            return(ans)
            returnType(double(3))
        })
    c12 <- compileNimble(f1, f2)
    x <- array(as.numeric(1:(7*6*8)), dim = c(7, 6, 8))
    expect_equal(f2(x), c12$f2(x))
})


test_that('slice of 3d passed as 2d arg works in %*%', {
    f1 <- nimbleFunction(
        run = function(x = double(2), b = double(1)) {
            ans <- x %*% b
            return(ans)
            returnType(double(2))
        })
    temporarilyAssignInGlobalEnv(f1)
    f2 <- nimbleFunction(
        run = function(x = double(3), b = double(1)) {
            ans <- f1(x[3, 2:4, 3:6], b)
            return(ans)
            returnType(double(2))
        })
    c12 <- compileNimble(f1, f2)
    x <- array(as.numeric(1:(4*6*8)), dim = c(4, 6, 8))
    b <- 11:14
    expect_equal(f2(x, b), c12$f2(x, b))
})

## 4D 4-dimensional slice tests

test_that('slice of 4d passed as 2d arg works', {
    f1 <- nimbleFunction(
        run = function(x = double(2)) {
            ans <- numeric(value = x, length = length(x))
            return(ans)
            returnType(double(1))
        })
    temporarilyAssignInGlobalEnv(f1)
    f2 <- nimbleFunction(
        run = function(x = double(4)) {
            ans <- f1(x[3, 2:4, 5, 3:6])
            return(ans)
            returnType(double(1))
        })
    c12 <- compileNimble(f1, f2)
    x <- array(as.numeric(1:(4*6*7*8)), dim = c(4, 6, 7, 8))
    expect_equal(f2(x), c12$f2(x))
})

test_that('slice of 4d passed as 3d arg works', {
    f1 <- nimbleFunction(
        run = function(x = double(3)) {
            ans <- numeric(value = x, length = length(x))
            return(ans)
            returnType(double(1))
        })
    temporarilyAssignInGlobalEnv(f1)
    f2 <- nimbleFunction(
        run = function(x = double(4)) {
            ans <- f1(x[3, 2:4, 2:5, 3:6])
            return(ans)
            returnType(double(1))
        })
    c12 <- compileNimble(f1, f2)
    x <- array(as.numeric(1:(4*6*7*8)), dim = c(4, 6, 7, 8))
    expect_equal(f2(x), c12$f2(x))
})

test_that('slice of 4dcopy works', {
    f1 <- nimbleFunction(
        run = function(x = double(4)) {
            ans <- x
            return(ans)
            returnType(double(4))
        })
    temporarilyAssignInGlobalEnv(f1)
    f2 <- nimbleFunction(
        run = function(x = double(4)) {
            ans <- f1(x[3:7, 2:4, 1:5, 3:6])
            return(ans)
            returnType(double(4))
        })
    c12 <- compileNimble(f1, f2)
    x <- array(as.numeric(1:(8*6*7*8)), dim = c(8, 6, 7, 8))
    expect_equal(f2(x), c12$f2(x))
})

## 5D 5-dimensional slice tests

test_that('slice of 5d passed as 4d arg works', {
    f1 <- nimbleFunction(
        run = function(x = double(4)) {
            ans <- numeric(value = x, length = length(x))
            return(ans)
            returnType(double(1))
        })
    temporarilyAssignInGlobalEnv(f1)
    f2 <- nimbleFunction(
        run = function(x = double(5)) {
            ans <- f1(x[3, 2:4, 2:5, 2:7, 3:6])
            return(ans)
            returnType(double(1))
        })
    c12 <- compileNimble(f1, f2)
    x <- array(as.numeric(1:(4*6*7*9*8)), dim = c(4, 6, 7, 9, 8))
    expect_equal(f2(x), c12$f2(x))
})

test_that('slice of 5d passed as 3d arg works', {
    f1 <- nimbleFunction(
        run = function(x = double(3)) {
            ans <- numeric(value = x, length = length(x))
            return(ans)
            returnType(double(1))
        })
    temporarilyAssignInGlobalEnv(f1)
    f2 <- nimbleFunction(
        run = function(x = double(5)) {
            ans <- f1(x[3, 2:4, 2:5, 7, 3:6])
            return(ans)
            returnType(double(1))
        })
    c12 <- compileNimble(f1, f2)
    x <- array(as.numeric(1:(4*6*7*9*8)), dim = c(4, 6, 7, 9, 8))
    expect_equal(f2(x), c12$f2(x))
})

test_that('slice of 5d passed as 2d arg works with %*%', {
    f1 <- nimbleFunction(
        run = function(x = double(2), b = double(1)) {
            ans <- x %*% b
            return(ans)
            returnType(double(2))
        })
    temporarilyAssignInGlobalEnv(f1)
    f2 <- nimbleFunction(
        run = function(x = double(5), b = double(1)) {
            ans <- f1(x[3, 2:4, 4, 7, 3:6], b)
            return(ans)
            returnType(double(2))
        })
    c12 <- compileNimble(f1, f2)
    x <- array(as.numeric(1:(4*6*7*9*8)), dim = c(4, 6, 7, 9, 8))
    b <- 11:14
    expect_equal(f2(x, b), c12$f2(x, b))
})

test_that('slice of 5d copy works', {
    f1 <- nimbleFunction(
        run = function(x = double(5)) {
            ans <- x
            return(ans)
            returnType(double(5))
        })
    temporarilyAssignInGlobalEnv(f1)
    f2 <- nimbleFunction(
        run = function(x = double(5)) {
            ans <- f1(x[3:7, 2:4, 1:5, 3:6, 4:8])
            return(ans)
            returnType(double(5))
        })
    c12 <- compileNimble(f1, f2)
    x <- array(as.numeric(1:(8*6*7*8*9)), dim = c(8, 6, 7, 8, 9))
    expect_equal(f2(x), c12$f2(x))
})

## 6D 6-dimensional slice tests
test_that('slice of 6d passed as 5d arg works', {
    f1 <- nimbleFunction(
        run = function(x = double(5)) {
            ans <- numeric(value = x, length = length(x))
            return(ans)
            returnType(double(1))
        })
    temporarilyAssignInGlobalEnv(f1)
    f2 <- nimbleFunction(
        run = function(x = double(6)) {
            ans <- f1(x[3, 2:4, 3:7, 2:5, 2:7, 3:6])
            return(ans)
            returnType(double(1))
        })
    c12 <- compileNimble(f1, f2)
    x <- array(as.numeric(1:(4*6*7*9*8*11)), dim = c(4, 6, 7, 9, 8, 11))
    expect_equal(f2(x), c12$f2(x))
})

test_that('slice of 6d passed as 4d arg works', {
    f1 <- nimbleFunction(
        run = function(x = double(4)) {
            ans <- numeric(value = x, length = length(x))
            return(ans)
            returnType(double(1))
        })
    temporarilyAssignInGlobalEnv(f1)
    f2 <- nimbleFunction(
        run = function(x = double(6)) {
            ans <- f1(x[3, 2:4, 3, 2:5, 2:7, 3:6])
            return(ans)
            returnType(double(1))
        })
    c12 <- compileNimble(f1, f2)
    x <- array(as.numeric(1:(4*6*7*9*8*11)), dim = c(4, 6, 7, 9, 8, 11))
    expect_equal(f2(x), c12$f2(x))
})

test_that('slice of 6d passed as 3d arg works', {
    f1 <- nimbleFunction(
        run = function(x = double(3)) {
            ans <- numeric(value = x, length = length(x))
            return(ans)
            returnType(double(1))
        })
    temporarilyAssignInGlobalEnv(f1)
    f2 <- nimbleFunction(
        run = function(x = double(6)) {
            ans <- f1(x[3, 2:4, 3, 2:5, 7, 3:6])
            return(ans)
            returnType(double(1))
        })
    c12 <- compileNimble(f1, f2)
    x <- array(as.numeric(1:(4*6*7*9*8*11)), dim = c(4, 6, 7, 9, 8, 11))
    expect_equal(f2(x), c12$f2(x))
})

test_that('slice of 6d passed as 2d arg works with %*%', {
    f1 <- nimbleFunction(
        run = function(x = double(2), b = double(1)) {
            ans <- x %*% b
            return(ans)
            returnType(double(2))
        })
    temporarilyAssignInGlobalEnv(f1)
    f2 <- nimbleFunction(
        run = function(x = double(6), b = double(1)) {
            ans <- f1(x[3, 2:4, 3, 4, 7, 3:6], b)
            return(ans)
            returnType(double(2))
        })
    c12 <- compileNimble(f1, f2)
    x <- array(as.numeric(1:(4*6*7*9*8*11)), dim = c(4, 6, 7, 9, 8, 11))
    b <- 11:14
    expect_equal(f2(x, b), c12$f2(x, b))
})

test_that('slice of 6d copy works', {
    f1 <- nimbleFunction(
        run = function(x = double(6)) {
            ans <- x
            return(ans)
            returnType(double(6))
        })
    temporarilyAssignInGlobalEnv(f1)
    f2 <- nimbleFunction(
        run = function(x = double(6)) {
            ans <- f1(x[3:7, 2:4, 3:5, 1:5, 3:6, 4:8])
            return(ans)
            returnType(double(6))
        })
    c12 <- compileNimble(f1, f2)
    x <- array(as.numeric(1:(8*6*7*8*9*11)), dim = c(8, 6, 7, 8, 9, 11))
    expect_equal(f2(x), c12$f2(x))
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)

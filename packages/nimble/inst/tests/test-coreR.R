##source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of core R functions in NIMBLE code")

writeLines(c("need to get Scalar instead of result_type",
             "fix nbinom"))

gen_runFunCore <- function(input) {
    runFun <- function() {}
    formalsList <- input$args
    if(is.null(names(formalsList)))
        if(length(formalsList) > 0)
            names(formalsList) <- paste0('arg', seq_along(input$args))
    formals(runFun) <- formalsList
    tmp <- quote({})
    tmp[[2]] <- input$expr
    tmp[[3]] <- quote(return(out))
    tmp[[4]] <- substitute(returnType(OUT), list(OUT = input$outputType))
    body(runFun) <- tmp
    return(runFun)
}

test_coreRfeature <- function(input, verbose = TRUE, dirName = NULL) { ## a lot like test_math but a bit more flexible
  if(verbose) cat("### Testing", input$name, "###\n")
  runFun <- gen_runFunCore(input)
  nfR <- nimbleFunction(run = runFun)
  nfC <- compileNimble(nfR, dirName = dirName)
  nArgs <- length(input$args)
  evalEnv <- new.env()
  eval(input$setArgVals, envir = evalEnv)
  savedArgs <- as.list(evalEnv)
  seedToUse <- if(is.null(input[['seed']])) 31415927 else input[['seed']]
  set.seed(seedToUse)
  eval(input$expr, envir = evalEnv)
  savedOutputs <- as.list(evalEnv)
  list2env(savedArgs, envir = evalEnv)
  if(nArgs == 5) {
      set.seed(seedToUse)
      out_nfR = nfR(evalEnv$arg1, evalEnv$arg2, evalEnv$arg3, evalEnv$arg4, evalEnv$arg5)
      list2env(savedArgs, envir = evalEnv)
      set.seed(seedToUse)
      out_nfC = nfC(evalEnv$arg1, evalEnv$arg2, evalEnv$arg3, evalEnv$arg4, evalEnv$arg5)
  }  
  if(nArgs == 4) {
      set.seed(seedToUse)
      out_nfR = nfR(evalEnv$arg1, evalEnv$arg2, evalEnv$arg3, evalEnv$arg4)
      list2env(savedArgs, envir = evalEnv)
      set.seed(seedToUse)
      out_nfC = nfC(evalEnv$arg1, evalEnv$arg2, evalEnv$arg3, evalEnv$arg4)
  }  

  if(nArgs == 3) {
      set.seed(seedToUse)
      out_nfR = nfR(evalEnv$arg1, evalEnv$arg2, evalEnv$arg3)
      list2env(savedArgs, envir = evalEnv)
      set.seed(seedToUse)
      out_nfC = nfC(evalEnv$arg1, evalEnv$arg2, evalEnv$arg3)
  }  
  if(nArgs == 2) {
      set.seed(seedToUse)
      out_nfR = nfR(evalEnv$arg1, evalEnv$arg2)
      list2env(savedArgs, envir = evalEnv)
      set.seed(seedToUse)
      out_nfC = nfC(evalEnv$arg1, evalEnv$arg2)
  }
  if(nArgs == 1) {
      set.seed(seedToUse)
      out_nfR = nfR(evalEnv$arg1)
      list2env(savedArgs, envir = evalEnv)
      set.seed(seedToUse)
      out_nfC = nfC(evalEnv$arg1)
  }
  if(nArgs == 0) {
      set.seed(seedToUse)
      out_nfR = nfR()
      list2env(savedArgs, envir = evalEnv)
      set.seed(seedToUse)
      out_nfC = nfC()
  }
  out <- savedOutputs$out
  attributes(out) <- attributes(out_nfR) <- attributes(out_nfC) <- NULL
  checkEqual <- input[['checkEqual']]
  if(is.null(checkEqual)) checkEqual <- FALSE
  if(!checkEqual) {
      try(test_that(paste0("Identical test of coreRfeature (direct R vs. R nimbleFunction): ", input$name),
                    expect_identical(out, out_nfR)))
      try(test_that(paste0("Identical test of math (direct R vs. C++ nimbleFunction): ", input$name),
                    expect_identical(out, out_nfC)))
  } else {
      try(test_that(paste0("Equal test of coreRfeature (direct R vs. R nimbleFunction): ", input$name),
                    expect_equal(out, out_nfR)))
      try(test_that(paste0("Equal test of math (direct R vs. C++ nimbleFunction): ", input$name),
                    expect_equal(out, out_nfC)))
  }
  # unload DLL as R doesn't like to have too many loaded
  if(.Platform$OS.type != 'windows') nimble:::clearCompiled(nfR) ##dyn.unload(project$cppProjects[[1]]$getSOName())
  invisible(NULL)

}

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
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.integer(4:5)}), outputType = quote(double(1))),

    ## basic cases with x, times and each
    list(name = "rep(vector double, variable, each = 2)", expr = quote(out <- rep(arg1, times = arg2, each = 2)), args = list(arg1 = quote(double(1)), arg2 = quote(double(0))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(c(4))}), outputType = quote(double(1))),
    list(name = "rep(vector double, variable, variable)", expr = quote(out <- rep(arg1, times = arg2, each = arg3)), args = list(arg1 = quote(double(1)), arg2 = quote(double(0)), arg3 = quote(double(0))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(4); arg3 <- 5}), outputType = quote(double(1))),
    list(name = "rep(vector double, variable, firstarg)", expr = quote(out <- rep(arg1, times = arg2, each = arg3)), args = list(arg1 = quote(double(1)), arg2 = quote(double(0)), arg3 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- 4; arg3 <- c(5, 7)}), outputType = quote(double(1))),
    ## basic cases with x and length.out
    list(name = "rep(1, lenght.out = 3)", expr = quote(out <- rep(1, length.out = 3)), args = list(arg1 = quote(double(0))),
         setArgVals = quote({arg1 <- 3}), outputType = quote(double(1))),
    ## 16
    list(name = "rep(vector double, length.out = larger than vector)", expr = quote(out <- rep(arg1, length.out = 5)), args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3)}), outputType = quote(double(1))),
    list(name = "rep(vector double, length.out = smaller than vector)", expr = quote(out <- rep(arg1, length.out = 2)), args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3)}), outputType = quote(double(1))),
    list(name = "rep(vector double, length.out = 0)", expr = quote(out <- rep(arg1, length.out = 0)), args = list(arg1 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3)}), outputType = quote(double(1))),
    list(name = "rep(vector double, length.out = first arg)", expr = quote(out <- rep(arg1, length.out = arg2)), args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(c(7, 8))}), outputType = quote(double(1))),
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
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(c(7, 8)); arg3 <- as.numeric(4:5)}), outputType = quote(double(1))),
    

    ## x, times expressions
    list(name = "rep(vector double expression, expression)", expr = quote(out <- rep(exp(arg1), arg2^2)), args = list(arg1 = quote(double(1)), arg2 = quote(integer())),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- 4}), outputType = quote(double(1))),
    ##26
    list(name = "rep(vector double expression, non-scalar expression)", expr = quote(out <- rep(exp(arg1), sum(arg2^2))), args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(c(2,3))}), outputType = quote(double(1))),
        

    list(name = "rep(matrix, 3)", expr = quote(out <- rep(arg1, 3)), args = list(arg1 = quote(double(2))),
         setArgVals = quote({arg1 <- matrix(as.numeric(1:9), nrow = 3)}),outputType = quote(double(1))),

    list(name = "rep(vector, vector)", expr = quote(out <- rep(arg1, arg2)), args = list(arg1 = quote(double(1)), arg2 = quote(double(1))), ## not built yet
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
    list(name = "diag(vector) with epxression", expr = quote(out <- exp(diag(arg1)) + arg2), args = list(arg1 = quote(double(1)), arg2 = quote(double(1))),
         setArgVals = quote({arg1 <- as.numeric(1:3); arg2 <- as.numeric(11:13)}), outputType = quote(double(2))),

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
         setArgVals = quote({arg1 <- as.numeric(1:2); arg2 <- as.numeric(3.5); arg3 <- as.numeric(4.1)}), outputType = quote(double(1)))
)

rRecyclingRuleTests <- list(

)

seqTests <- list(
    ##1
    list(name = "1:5", expr = quote(out <- 1:5), args = list(),
         setArgVals = quote({}), outputType = quote(double(1)), checkEqual = TRUE),
    list(name = "seq(.1, 10, by = .1)", expr = quote(out <- seq(.1, 10, by = .1)), args = list(),
         setArgVals = quote({}), outputType = quote(double(1))),
    list(name = "seq(.1, 10, length.out = 11)", expr = quote(out <- seq(.1, 10, length.out = 11)), args = list(),
         setArgVals = quote({}), outputType = quote(double(1))),
    list(name = "seq(.1, 10, length.out = 11) in expression", expr = quote(out <- log(seq(.1, 10, length.out = 11)) + 2 + rep(1, 11)), args = list(),
         setArgVals = quote({}), outputType = quote(double(1)))

    ## need to handle this case
##    list(name = "seq(.1, 10, by = .1)", expr = quote(out <- seq(.1, by = 0.1, length.out = 11)), args = list(),
##         setArgVals = quote({}), outputType = quote(double(1)))
    ## need to handle decreasing colon sequenences
)
## STATUS: need to handle by and length.out case.
## need to handle decreasing colon sequences
## need to cast from integer to double.  


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
         outputType = quote(double(2)))
)




## lapply(cTests, test_coreRfeature)
## lapply(blockTests, test_coreRfeature)
## lapply(repTests, test_coreRfeature)
## lapply(diagTests, test_coreRfeature)
## lapply(recyclingRuleTests, test_coreRfeature)
## lapply(seqTests, test_coreRfeature)
## lapply(nonSeqIndexTests, test_coreRfeature)

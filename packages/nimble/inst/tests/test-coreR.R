##source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of core R functions in NIMBLE code")

writeLines(c("need to get Scalar instead of result_type",
             "rep needs to lift each and length expressions and handle vectors",
             "c needs to handle multiple arguments",
             "see notes in sizeConcatenate"))

gen_runFunCore <- function(input) {
    runFun <- function() {}
    formalsList <- input$args
    if(is.null(names(formalsList))) names(formalsList) <- paste0('arg', seq_along(input$args))
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
  eval(input$expr, envir = evalEnv)
  savedOutputs <- as.list(evalEnv)
  list2env(savedArgs, envir = evalEnv)
  if(nArgs == 3) {
      out_nfR = nfR(evalEnv$arg1, evalEnv$arg2, evalEnv$arg3)
      list2env(savedArgs, envir = evalEnv)
      out_nfC = nfC(evalEnv$arg1, evalEnv$arg2, evalEnv$arg3)
  }  
  if(nArgs == 2) {
      out_nfR = nfR(evalEnv$arg1, evalEnv$arg2)
      list2env(savedArgs, envir = evalEnv)
      out_nfC = nfC(evalEnv$arg1, evalEnv$arg2)
  }
  if(nArgs == 1) {
      out_nfR = nfR(evalEnv$arg1)
      list2env(savedArgs, envir = evalEnv)
      out_nfC = nfC(evalEnv$arg1)
  }
  out <- savedOutputs$out
  attributes(out) <- attributes(out_nfR) <- attributes(out_nfC) <- NULL
  try(test_that(paste0("Test of coreRfeature (direct R vs. R nimbleFunction): ", input$name),
                expect_identical(out, out_nfR)))
  try(test_that(paste0("Test of math (direct R vs. C++ nimbleFunction): ", input$name),
                expect_identical(out, out_nfC)))
  # unload DLL as R doesn't like to have too many loaded
  if(.Platform$OS.type != 'windows') nimble:::clearCompiled(nfR) ##dyn.unload(project$cppProjects[[1]]$getSOName())
  invisible(NULL)

}

diagTests <- list(
    list(name = "diag(scalar)", expr = quote(out <- diag(arg1)), args = list(arg1 = quote(double(0))),
         setArgVals = quote({arg1 <- 3}), outputType = double(1)))

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
         setArgVals = quote({arg1 <- matrix(as.numeric(1:4), nrow = 2); arg2 <- as.numeric(10:11)}), outputType = quote(double(1)))
## add some that use a c() in an expression
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

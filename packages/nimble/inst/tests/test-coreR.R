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
  eval(input$setArgVals)
  eval(input$expr)
   if(nArgs == 3) {
      out_nfR = nfR(arg1, arg2, arg3)
      out_nfC = nfC(arg1, arg2, arg3)
  }  
  if(nArgs == 2) {
      out_nfR = nfR(arg1, arg2)
      out_nfC = nfC(arg1, arg2)
  }
  if(nArgs == 1) {
    out_nfR = nfR(arg1)
    out_nfC = nfC(arg1)
  }
  attributes(out) <- attributes(out_nfR) <- attributes(out_nfC) <- NULL
  try(test_that(paste0("Test of coreRfeature (direct R vs. R nimbleFunction): ", input$name),
                expect_identical(out, out_nfR)))
  try(test_that(paste0("Test of math (direct R vs. C++ nimbleFunction): ", input$name),
                expect_identical(out, out_nfC)))
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
         setArgVals = quote({arg1 <- matrix(as.numeric(1:4), nrow = 2); arg2 <- as.numeric(10:11)}), outputType = quote(double(1)))
## add some that use a c() in an expression
)


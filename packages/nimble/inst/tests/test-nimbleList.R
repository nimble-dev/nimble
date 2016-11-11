## Tests for building, compiling, and using nimbleList objects
## There are three distinct ways that a nimbleList can be created for use in a nimbleFunction:
## 1) nimbleListDef is created outside of the nimbleFunction, nimbleListDef$new() used to create new nimbleList in setup code
## 2) nimbleListDef is created outside of the nimbleFunction, nimbleListDef$new() used to create new nimbleList in run code
## 3) nimbleListDef is created in setup code of the nimbleFunction, nimbleListDef$new() used to create new nimbleList in setup code
## 4) nimbleListDef is created in setup code of the nimbleFunction, nimbleListDef$new() used to create new nimbleList in run code

context('nimbleList() tests')


library(nimble)
library(testthat)
nimbleOptions(showCompilerOutput = TRUE)
 nimbleOptions(debugCppLineByLine = FALSE)
# nimbleOptions(debugSizeProcessing = TRUE)
# debug(nimble:::sizeAssignAfterRecursing)
options( warn = 1 )

########
## Test of creating new nimbleList in run code and specifying initial values for that list
## Here, the nlDef is created in the global environment, outside of setup code
########

nlTestFunc1 <- nimbleFunction(
  setup = function(){
    doubleMatrix <- diag(1)
  },
  run = function(){
    doubleScalar <- 1
    doubleVector <- numeric(2, 1)
    newList1 <- testListDef1$new(nlVector = doubleVector, nlMatrix = doubleMatrix)
    newList1$nlScalar <- doubleScalar
    returnType(testListDef1())
    return(newList1)
  }
)

testTypes <- list(vars = c('nlScalar', 'nlVector', 'nlMatrix'), types = c('double(0)', "double(1)", "double(2)"))
testListDef1 <- nimbleList(testTypes)
testInst <- nlTestFunc1()
RnimbleList <- testInst$run()
CtestInst <- compileNimble(testInst, control = list(debug =  F))
CnimbleList <- CtestInst$run()

## test for correct values of R nimbleList
expect_identical(RnimbleList$nlScalar, 1)
expect_identical(RnimbleList$nlVector, c(1, 1))
expect_identical(RnimbleList$nlMatrix, diag(1))
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nlScalar, CnimbleList$nlScalar)
expect_identical(RnimbleList$nlVector, CnimbleList$nlVector)
expect_identical(RnimbleList$nlMatrix, CnimbleList$nlMatrix)
test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })

########
## Test of creating new nimbleList in run code and specifying initial value for that list
## Here, the nlDef is created in setup code.  Inital value is an expression.
########

nlTestFunc2 <- nimbleFunction(
  setup = function(){
    testTypes <- list(vars = c('nlScalar'), types = c('double(0)'))
    testListDef2 <- nimbleList(testTypes)
    doubleMatrix <- diag(1)
  },
  run = function(){
    newList2 <- testListDef2$new(nlScalar = doubleMatrix[1,1]*2)
    returnType(testListDef2())
    return(newList2)
  }
)

testInst <- nlTestFunc2()
RnimbleList <- testInst$run()
ctestInst <- compileNimble(testInst, control = list(debug =  F))
CnimbleList <- ctestInst$run()

## test for correct values of R nimbleList
expect_identical(RnimbleList$nlScalar, 2)
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nlScalar, CnimbleList$nlScalar)
test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })


########
## Test of creating new nimbleList in setup code and specifying initial values for that list.
## Here, the nlDef is created in setup code
########

nlTestFunc3 <- nimbleFunction(
  setup = function(){
    testTypes <- list(vars = c('nlCharacter'), types = c('character(0)'))
    testListDef3 <- nimbleList(testTypes)
    setupList3 <- testListDef3$new(nlCharacter = "hello world")
  },
  run = function(){
    returnType(testListDef3())
    return(setupList3)
  }
)

testInst <- nlTestFunc3()
RnimbleList <- testInst$run()
ctestInst <- compileNimble(testInst)
CnimbleList <- ctestInst$run()

## test for correct values of R nimbleList
expect_identical(RnimbleList$nlCharacter, "hello world")
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nlCharacter, CnimbleList$nlCharacter)
test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })


########
## Test of using a nimbleList as a run function argument 
########


nlTestFunc4 <- nimbleFunction(
  setup = function(){
  },
  run = function(argList4 = testListDef4()){
    argList4$nlMatrix[1,2] <- 10
    returnType(testListDef4())
    return(argList4)
  }
)

testTypes <- list(vars = c('nlMatrix'), types = c('double(2)'))
testListDef4 <- nimble:::nimbleList(testTypes)
testList4 <- testListDef4$new(nlMatrix = diag(2))

testInst <- nlTestFunc4()
RnimbleList <- testInst$run(testList4)
ctestInst <- compileNimble(testInst, control = list(debug =  F))
CnimbleList <- ctestInst$run(testList4)

## test for correct values of R nimbleList
expect_identical(RnimbleList$nlMatrix, matrix(c(1,0,10,1),nrow =  2))
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nlMatrix, CnimbleList$nlMatrix)
## test for identical values of testList4 and CnimbleList
expect_identical(testList4$nlMatrix, CnimbleList$nlMatrix)
test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
})



########
## Test of using a nimbleList as an argument to a function within a nimbleFunction  
########


innerNlTestFunc5 <- nimbleFunction(
  setup = function(){},
  run = function(argList5 = testListDef5()){
    argList5$a <- argList5$a*2
    returnType(testListDef5())
    return(argList5)
  })

nlTestFunc5 <- nimbleFunction(
  setup = function(){
    innerFunc5 <- innerNlTestFunc5()
  },
  run = function(argList5 = testListDef5()){
    argList5$a <- 1
    argList5 <- innerFunc5$run(argList5)
    returnType(testListDef5())
    return(argList5)
  }
)

testTypes <- list(vars = c('a'), types = c('double(0)'))
testListDef5 <- nimble:::nimbleList(testTypes)
testList5 <- testListDef5$new(a = 0)

testInst <- nlTestFunc5()
RnimbleList <- testInst$run(testList5)
ctestInst <- compileNimble(testInst, control = list(debug =  F))
CnimbleList <- ctestInst$run(testList5)

## test for correct values of R nimbleList
expect_identical(RnimbleList$a, 2)
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$a, CnimbleList$a)
## test for identical values of testList5 and CnimbleList
expect_identical(testList5$a, CnimbleList$a)
test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })


########
## Test of using a nimbleList returned from a nimbleFunction as a run argument to another nimble function
########


nlTestFunc6 <- nimbleFunction(
  setup = function(){
  },
  run = function(argList6 = testListDef6()){
    argList6$nlMatrix[1,2] <- argList6$nlMatrix[1,2] + 10
    returnType(testListDef6())
    return(argList6)
  }
)

testTypes <- list(vars = c('nlMatrix'), types = c('double(2)'))
testListDef6 <- nimble:::nimbleList(testTypes)
testList6 <- testListDef4$new(nlMatrix = diag(2))

testInst <- nlTestFunc6()
RnimbleList <- testInst$run(testList6)
RnimbleList <- testInst$run(RnimbleList)
ctestInst <- compileNimble(testInst, control = list(debug =  F))
CnimbleList <- ctestInst$run(testList6)
CnimbleList <- ctestInst$run(CnimbleList)

## test for correct values of R nimbleList
expect_identical(RnimbleList$nlMatrix, matrix(c(1,0,40,1),nrow =  2))
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nlMatrix, CnimbleList$nlMatrix)
##test for identical values of testList6 and CnimbleList
expect_identical(testList6$nlMatrix, CnimbleList$nlMatrix)

test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })


########
## Test of using two nimbleLists as  run function arguments 
########


nlTestFunc7 <- nimbleFunction(
  setup = function(){
    diamat <- diag(2)
  },
  run = function(argList7a = testListDef7(), argList7b = testListDef7()){
    argList7a$nlMatrix <- argList7a$nlMatrix %*% argList7b$nlMatrix
    returnType(testListDef7())
    return(argList7a)
  }
)

testTypes <- list(vars = c('nlMatrix'), types = c('double(2)'))
testListDef7 <- nimble:::nimbleList(testTypes)
testList7a <- testListDef7$new(nlMatrix = matrix(1, nrow = 2, ncol = 2))
testList7b <- testListDef7$new(nlMatrix = matrix(2, nrow = 2, ncol = 2))

testInst <- nlTestFunc7()
RnimbleList <- testInst$run(testList7a, testList7b)
ctestInst <- compileNimble(testInst, control = list(debug =  F))
CnimbleList <- ctestInst$run(testList7a, testList7b)

## test for correct values of R nimbleList
expect_identical(RnimbleList$nlMatrix, matrix(c(16,16,16,16),nrow =  2))
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nlMatrix, CnimbleList$nlMatrix)
## test for identical values of testList7a and CnimbleList
expect_identical(testList7a$nlMatrix, CnimbleList$nlMatrix)
test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })


########
## Test of assigning a new nimbleList object to a nimbleList object that was passed as a run argument
########


nlTestFunc8 <- nimbleFunction(
  setup = function(){
    diamat <- diag(2)
  },
  run = function(argList8 = testListDef8()){
    argList8 <- testListDef8$new(nlMatrix = diamat)
    returnType(testListDef8())
    return(argList8)
  }
)

testTypes <- list(vars = c('nlMatrix'), types = c('double(2)'))
testListDef8 <- nimble:::nimbleList(testTypes)
testList8 <- testListDef8$new(nlMatrix = matrix(1, nrow = 2, ncol = 2))

trace(nimble:::cppNewNimbleList, browser)

testInst <- nlTestFunc8()
RnimbleList <- testInst$run(testList8)
ctestInst <- compileNimble(testInst, control = list(debug =  F))
CnimbleList <- ctestInst$run(testList8)

## test for correct values of R nimbleList
expect_identical(RnimbleList$nlMatrix, diag(2))
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nlMatrix, CnimbleList$nlMatrix)
## test that testList8 doesn't get changed by function
expect_identical(testList8$nlMatrix, matrix(1, nrow = 2, ncol = 2))
test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })



########
## Test of using two nimbleLists as  run function arguments, assigning one list to the other within the function
########

nlTestFunc9 <- nimbleFunction(
  setup = function(){
  },
  run = function(argList9a = testListDef9(), argList9b = testListDef9()){
    argList9a <- argList9b
    returnType(testListDef9())
    return(argList9a)
  }
)

testTypes <- list(vars = c('nlMatrix'), types = c('double(2)'))
testListDef9 <- nimble:::nimbleList(testTypes)
testList9a <- testListDef9$new(nlMatrix = matrix(1, nrow = 2, ncol = 2))
testList9b <- testListDef9$new(nlMatrix = matrix(2, nrow = 2, ncol = 2))

testInst <- nlTestFunc9()
RnimbleList <- testInst$run(testList9a, testList9b)
ctestInst <- compileNimble(testInst, control = list(debug =  F))
CnimbleList <- ctestInst$run(testList9a, testList9b)

## test for correct values of R nimbleList
expect_identical(RnimbleList$nlMatrix, matrix(c(16,16,16,16),nrow =  2))
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nlMatrix, CnimbleList$nlMatrix)
## test for identical values of testList7a and CnimbleList
expect_identical(testList7a$nlMatrix, CnimbleList$nlMatrix)
test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })


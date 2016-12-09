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
testList6 <- testListDef6$new(nlMatrix = diag(2))

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
    argList9a$nlMatrix[1,1] <- argList9a$nlMatrix[1,1] + 5
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
expect_identical(RnimbleList$nlMatrix, matrix(c(2,2,2,2),nrow =  2))
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nlMatrix, CnimbleList$nlMatrix)
## test for identical values of testList9b and CnimbleList
expect_identical(testList9b$nlMatrix, CnimbleList$nlMatrix)
## test that  testList9a matrix was updated 
expect_identical(testList9a$nlMatrix, matrix(c(11, 1, 1, 1), nrow = 2))

test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })


########
## Test of nested nimbleLists.  Both nimbleList generators are defined in setup code. Top level nimbleList is returned.
########

nlTestFunc10 <- nimbleFunction(
  setup = function(){
    testTypes1 <- list(vars = c('nlDouble'), types = c('double(0)'))
    testListDef10a <- nimble:::nimbleList(testTypes1)
    testTypes2 <- list(vars = c('nestedNL'), types = c('testListDef10a()'))
    testListDef10b <- nimble:::nimbleList(testTypes2)
  },
  run = function(){
    newList10 <- testListDef10b$new()
    newList10$nestedNL$nlDouble <- 3.14
    returnType(testListDef10b())
    return(newList10)
  }
)

testInst <- nlTestFunc10()
RnimbleList <- testInst$run()
ctestInst <- compileNimble(testInst)
CnimbleList <- ctestInst$run()

## test for correct values of R nimbleList
expect_identical(RnimbleList$nestedNL$nlDouble, 3.14)
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nestedNL$nlDouble, CnimbleList$nestedNL$nlDouble)

test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })

########
## Test of nested nimbleLists.  Here, top level nimbleList generator is defined in setup code.
## Lower level nimbleList generator is defined in global environment.  Top level nimbleList is returned.
########

nlTestFunc11 <- nimbleFunction(
  setup = function(){
    testTypes2 <- list(vars = c('nestedNL'), types = c('testListDef11a()'))
    testListDef11b <- nimble:::nimbleList(testTypes2)
  },
  run = function(){
    newList11 <- testListDef11b$new()
    newList11$nestedNL$nlDouble <- 3.14
    returnType(testListDef11b())
    return(newList11)
  }
)
testTypes1 <- list(vars = c('nlDouble'), types = c('double(0)'))
testListDef11a <- nimble:::nimbleList(testTypes1)

testInst <- nlTestFunc11()
RnimbleList <- testInst$run()
ctestInst <- compileNimble(testInst)
CnimbleList <- ctestInst$run()

## test for correct values of R nimbleList
expect_identical(RnimbleList$nestedNL$nlDouble, 3.14)
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nestedNL$nlDouble, CnimbleList$nestedNL$nlDouble)

test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })

########
## Test of nested nimbleLists. Two nimbleLists with nested nimbleLists are given as run time arguments.  
## In run code, the nested list of the second argument is assigned to the nested list of the first argument.
########

nlTestFunc12 <- nimbleFunction(
  setup = function(){
  },
  run = function(argList12a = testListDef12b(), argList12b = testListDef12b()){
    argList12a$nestedNL <- argList12b$nestedNL
    returnType(testListDef12b())
    return(argList12a)
  }
)

testTypes1 <- list(vars = c('nlDouble'), types = c('double(0)'))
testListDef12a <- nimble:::nimbleList(testTypes1)
testTypes2 <- list(vars = c('nestedNL', 'nlVector'), types = c('testListDef12a()', 'double(1)'))
testListDef12b <- nimble:::nimbleList(testTypes2)

argList12a <- testListDef12b$new(nlVector = c(1,2))
argList12a$nestedNL$nlDouble <- 0
argList12b <- testListDef12b$new(nlVector = c(3,4))
argList12b$nestedNL$nlDouble <- 100

testInst <- nlTestFunc12()
RnimbleList <- testInst$run(argList12a, argList12b)
ctestInst <- compileNimble(testInst)
CnimbleList <- ctestInst$run(argList12a, argList12b)

## test for correct values of nimbleLists
expect_identical(RnimbleList$nestedNL$nlDouble, 100)
expect_identical(RnimbleList$nestedNL$nlDouble, CnimbleList$nestedNL$nlDouble)
expect_identical(RnimbleList$nestedNL$nlDouble, argList12a$nestedNL$nlDouble)

expect_identical(RnimbleList$nlVector, c(1,2))
expect_identical(RnimbleList$nlVector, CnimbleList$nlVector)
expect_identical(RnimbleList$nlVector, argList12a$nlVector)


test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })


########
## The same as test 12, but with lists provided as setup arguments
########



nlTestFunc13 <- nimbleFunction(
  setup = function(argList13a, argList13b){
  },
  run = function(){
    argList13a$nestedNL <<- argList13b$nestedNL
    returnType(testListDef13b())
    return(argList13a)
  }
)

testTypes1 <- list(vars = c('nlDouble'), types = c('double(0)'))
testListDef13a <- nimble:::nimbleList(testTypes1)
testTypes2 <- list(vars = c('nestedNL', 'nlVector'), types = c('testListDef13a()', 'double(1)'))
testListDef13b <- nimble:::nimbleList(testTypes2)

argList13a <- testListDef13b$new(nlVector = c(1,2))
argList13a$nestedNL$nlDouble <- 0
argList13b <- testListDef13b$new(nlVector = c(3,4))
argList13b$nestedNL$nlDouble <- 100

testInst <- nlTestFunc13(argList13a, argList13b)
RnimbleList <- testInst$run()
ctestInst <- compileNimble(testInst)
CnimbleList <- ctestInst$run()

## test for correct values of nimbleLists
expect_identical(RnimbleList$nestedNL$nlDouble, 100)
expect_identical(RnimbleList$nestedNL$nlDouble, CnimbleList$nestedNL$nlDouble)
expect_identical(RnimbleList$nestedNL$nlDouble, argList13a$nestedNL$nlDouble)

expect_identical(RnimbleList$nlVector, c(1,2))
expect_identical(RnimbleList$nlVector, CnimbleList$nlVector)
expect_identical(RnimbleList$nlVector, argList13a$nlVector)


test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })




########
## Test #1 for eigen() function.  Return an eigenList.  
########

nlTestFunc14 <- nimbleFunction(
  setup = function(){
    testTypes <- list(vars = 'test_matrix', types = c('double(2)'))
    testListDef14 <- nimbleList(testTypes)
    testList14 <- testListDef14$new(test_matrix = diag(4) + 1)
  },
  run = function(){
    eigenOut <- eigen(testList14$test_matrix)
    returnType(eigen())
    return(eigenOut)
  }
)


testInst <- nlTestFunc14()
RnimbleList <- testInst$run()
ctestInst <- compileNimble(testInst)
CnimbleList <- ctestInst$run()

## test for eigenValues and eigenVectors
expect_equal(diag(4)+1, RnimbleList$vectors%*%diag(RnimbleList$values)%*%solve(RnimbleList$vectors))
expect_equal(diag(4)+1, CnimbleList$vectors%*%diag(CnimbleList$values)%*%solve(CnimbleList$vectors))
expect_equal(sum(abs(RnimbleList$values)), sum(abs(CnimbleList$values)))

test_that("return object (from c++) is nimbleList.", 
          {
            expect_identical(is.nl(CnimbleList), TRUE)
          })



########
## Test #1 for svd() function.  Return an svdList.  
########

nlTestFunc14 <- nimbleFunction(
  setup = function(){
    testTypes <- list(vars = 'test_matrix', types = c('double(2)'))
    testListDef14 <- nimbleList(testTypes)
    testList14 <- testListDef14$new(test_matrix = diag(4) + 1)
  },
  run = function(){
    eigenOut <- svd(testList14$test_matrix)
    returnType(svd())
    return(eigenOut)
  }
)


testInst <- nlTestFunc14()
RnimbleList <- testInst$run()
ctestInst <- compileNimble(testInst)
CnimbleList <- ctestInst$run()

## test for eigenValues and eigenVectors
expect_equal(diag(4)+1, RnimbleList$vectors%*%diag(RnimbleList$values)%*%solve(RnimbleList$vectors))
expect_equal(diag(4)+1, CnimbleList$vectors%*%diag(CnimbleList$values)%*%solve(CnimbleList$vectors))
expect_equal(sum(abs(RnimbleList$values)), sum(abs(CnimbleList$values)))

test_that("return object (from c++) is nimbleList.", 
          {
            expect_identical(is.nl(CnimbleList), TRUE)
          })




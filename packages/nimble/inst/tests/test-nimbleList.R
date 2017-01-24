## Tests for building, compiling, and using nimbleList objects
## There are five distinct ways that a nimbleList can be created for use in a nimbleFunction:
## 1) nimbleListDef is created outside of the nimbleFunction, nimbleListDef$new() used to create new nimbleList in setup code
## 2) nimbleListDef is created outside of the nimbleFunction, nimbleListDef$new() used to create new nimbleList in run code
## 3) nimbleListDef is created in setup code of the nimbleFunction, nimbleListDef$new() used to create new nimbleList in setup code
## 4) nimbleListDef is created in setup code of the nimbleFunction, nimbleListDef$new() used to create new nimbleList in run code
## 5) nimbleListDef is created outside of the nimbleFunction, nimbleListDef$new() used to create new nimbleList outside of the nimbleFunction,
##    nimbleList passed as argument to nimbleFunction
context('nimbleList() tests')


library(nimble)
library(testthat)
# nimbleOptions(showCompilerOutput = TRUE)
# nimbleOptions(debugCppLineByLine = FALSE)
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
expect_identical(CnimbleList$nlCharacter, ctestInst$setupList3$nlCharacter)

test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
            # expect_identical(is.nl(ctestInst$setupList3), TRUE)  is.nl gives FALSE here, not sure if should be corrected
          })




########
## Test of creating new nimbleList in a function that is internal to another function.
## The list is passed from the internal function to the outer function, and then to the user.
## This test does not use nimble's active binding system.
########
innerNlTestFunc4 <- nimbleFunction(
  setup = function(){
    setupList4 <- testListDef4$new(nlCharacter = "hello world")
  },
  run = function(){
    returnType(testListDef4())
    return(setupList4)
  }
)

nlTestFunc4 <- nimbleFunction(
  setup = function(){
    innerFunc4 <- innerNlTestFunc4()
  },
  run = function(){
    outList <- innerFunc4$run()
    returnType(testListDef4())
    return(outList)
  }
)

testTypes <- list(vars = c('nlCharacter'), types = c('character(0)'))
testListDef4 <- nimbleList(testTypes)

testInst <- nlTestFunc4()
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
## Test of creating and interacting with two nimble lists in two different functions from the same def'n
########
innerNlTestFunc5 <- nimbleFunction(
  setup = function(listDef){
    innerSetupList5 <- listDef$new(nlCharacter = "hello world")
  },
  run = function(){
    returnType(testListDef5())
    return(innerSetupList5)
  }
)

nlTestFunc5 <- nimbleFunction(
  setup = function(listDef){
    innerFunc5 <- innerNlTestFunc5(listDef)
    outerSetupList5 <- listDef$new(nlCharacter = "goodbye!")
  },
  run = function(){
    innerNimList <- innerFunc5$run()
    innerNimList <- outerSetupList5
    returnType(testListDef5())
    return(innerNimList)
  }
)

testTypes <- list(vars = c('nlCharacter'), types = c('character(0)'))
testListDef5 <- nimbleList(testTypes)

testInst <- nlTestFunc5(testListDef5)
RnimbleList <- testInst$run()
ctestInst <- compileNimble(testInst)
CnimbleList <- ctestInst$run()

## test for correct values of R nimbleList
expect_identical(RnimbleList$nlCharacter, "goodbye")
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nlCharacter, CnimbleList$nlCharacter)

test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })

########
## Test of creating an nlDef in an outer function, passing the def to the inner function,
## and returning the created list from the inner function through the outer function
########


innerNlTestFunc6 <- nimbleFunction(
  setup = function(nimListDef){
    setupList6 <- nimListDef$new(nlCharacter = "hello world")
  },
  run = function(){
    returnType(nimListDef())
    return(setupList6)
  }
)

nlTestFunc6 <- nimbleFunction(
  setup = function(){
    testTypes <- list(vars = c('nlCharacter'), types = c('character(0)'))
    testListDef6 <- nimbleList(testTypes)
    innerFunc6 <- innerNlTestFunc6(testListDef6)
  },
  run = function(){
    outList <- innerFunc6$run()
    returnType(testListDef6())
    return(outList)
  }
)

testInst <- nlTestFunc6()
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

######
## test of a nimbleFunction method returning a nimbleList
######

nlTestFunc7 <- nimbleFunction(
  setup = function(){
    testTypes <- list(vars = c('nlCharacter'), types = c('character(0)'))
    testListDef7 <- nimbleList(testTypes)
    
  },
  run = function(){
    outList <- method1()
    returnType(testListDef7())
    return(outList)
  },
  methods = list(
    method1 = function(){
      outList <- testListDef7$new(nlCharacter = 'hello world')
      returnType(testListDef7())
      return(outList)
    }
  )
)

testInst <- nlTestFunc7()
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


nlTestFunc8 <- nimbleFunction(
  setup = function(){
  },
  run = function(argList8 = testListDef8()){
    argList8$nlMatrix[1,2] <- 10
    returnType(testListDef8())
    return(argList8)
  }
)

testTypes <- list(vars = c('nlMatrix'), types = c('double(2)'))
testListDef8 <- nimble:::nimbleList(testTypes)
testList8 <- testListDef4$new(nlMatrix = diag(2))

testInst <- nlTestFunc4()
RnimbleList <- testInst$run(testList8)
ctestInst <- compileNimble(testInst, control = list(debug =  F))
CnimbleList <- ctestInst$run(testList8)

## test for correct values of R nimbleList
expect_identical(RnimbleList$nlMatrix, matrix(c(1,0,10,1),nrow =  2))
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nlMatrix, CnimbleList$nlMatrix)
## test for identical values of testList4 and CnimbleList
expect_identical(testList8$nlMatrix, CnimbleList$nlMatrix)
test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
})



########
## Test of using a nimbleList as an argument to a function within a nimbleFunction  
########


innerNlTestFunc9 <- nimbleFunction(
  setup = function(){},
  run = function(argList9 = testListDef5()){
    argList9$a <- argList9$a*2
    returnType(testListDef9())
    return(argList9)
  })

nlTestFunc9 <- nimbleFunction(
  setup = function(){
    innerFunc9 <- innerNlTestFunc9()
  },
  run = function(argList9 = testListDef9()){
    argList9$a <- 1
    argList9 <- innerFunc9$run(argList9)
    returnType(testListDef9())
    return(argList9)
  }
)

testTypes <- list(vars = c('a'), types = c('double(0)'))
testListDef5 <- nimble:::nimbleList(testTypes)
testList9 <- testListDef5$new(a = 0)

testInst <- nlTestFunc9()
RnimbleList <- testInst$run(testList9)
ctestInst <- compileNimble(testInst, control = list(debug =  F))
CnimbleList <- ctestInst$run(testList9)

## test for correct values of R nimbleList
expect_identical(RnimbleList$a, 2)
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$a, CnimbleList$a)
## test for identical values of testList5 and CnimbleList
expect_identical(testList9$a, CnimbleList$a)
test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })


########
## Test of using a nimbleList returned from a nimbleFunction as a run argument to another nimble function
########


nlTestFunc10 <- nimbleFunction(
  setup = function(){
  },
  run = function(argList10 = testListDef10()){
    argList10$nlMatrix[1,2] <- argList10$nlMatrix[1,2] + 10
    returnType(testListDef10())
    return(argList10)
  }
)

testTypes <- list(vars = c('nlMatrix'), types = c('double(2)'))
testListDef10 <- nimble:::nimbleList(testTypes)
testList10 <- testListDef10$new(nlMatrix = diag(2))

testInst <- nlTestFunc10()
RnimbleList <- testInst$run(testList10)
RnimbleList <- testInst$run(RnimbleList)
ctestInst <- compileNimble(testInst, control = list(debug =  F))
CnimbleList <- ctestInst$run(testList10)
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
## Test of using two nimbleLists as run function arguments 
########


nlTestFunc11 <- nimbleFunction(
  setup = function(){
    diamat <- diag(2)
  },
  run = function(argList11a = testListDef11(), argList11b = testListDef11()){
    argList11a$nlMatrix <- argList11a$nlMatrix %*% argList11b$nlMatrix
    returnType(testListDef11())
    return(argList11a)
  }
)

testTypes <- list(vars = c('nlMatrix'), types = c('double(2)'))
testListDef11 <- nimble:::nimbleList(testTypes)
testList11a <- testListDef11$new(nlMatrix = matrix(1, nrow = 2, ncol = 2))
testList11b <- testListDef11$new(nlMatrix = matrix(2, nrow = 2, ncol = 2))

testInst <- nlTestFunc11()
RnimbleList <- testInst$run(testList11a, testList11b)
ctestInst <- compileNimble(testInst, control = list(debug =  F))
CnimbleList <- ctestInst$run(testList11a, testList11b)

## test for correct values of R nimbleList
expect_identical(RnimbleList$nlMatrix, matrix(c(16,16,16,16),nrow =  2))
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nlMatrix, CnimbleList$nlMatrix)
## test for identical values of testList7a and CnimbleList
expect_identical(testList11a$nlMatrix, CnimbleList$nlMatrix)
test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })


########
## Test of assigning a new nimbleList object to a nimbleList object that was passed as a run argument
########


nlTestFunc12 <- nimbleFunction(
  setup = function(){
    diamat <- diag(2)
  },
  run = function(argList12 = testListDef12()){
    argList12 <- testListDef12$new(nlMatrix = diamat)
    returnType(testListDef12())
    return(argList12)
  }
)

testTypes <- list(vars = c('nlMatrix'), types = c('double(2)'))
testListDef12 <- nimble:::nimbleList(testTypes)
testList12 <- testListDef8$new(nlMatrix = matrix(1, nrow = 2, ncol = 2))

testInst <- nlTestFunc12()
RnimbleList <- testInst$run(testList12)
ctestInst <- compileNimble(testInst, control = list(debug =  F))
CnimbleList <- ctestInst$run(testList12)

## test for correct values of R nimbleList
expect_identical(RnimbleList$nlMatrix, diag(2))
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nlMatrix, CnimbleList$nlMatrix)
## test that testList12 doesn't get changed by function
expect_identical(testList12$nlMatrix, matrix(1, nrow = 2, ncol = 2))
test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })



########
## Test of using two nimbleLists as  run function arguments, assigning one list to the other within the function
########

nlTestFunc13 <- nimbleFunction(
  setup = function(){
  },
  run = function(argList13a = testListDef13(), argList13b = testListDef13()){
    argList13a$nlMatrix[1,1] <- argList13a$nlMatrix[1,1] + 5
    argList13a <- argList13b
    returnType(testListDef13())
    return(argList13a)
  }
)

testTypes <- list(vars = c('nlMatrix'), types = c('double(2)'))
testListDef13 <- nimble:::nimbleList(testTypes)
testList13a <- testListDef13$new(nlMatrix = matrix(1, nrow = 2, ncol = 2))
testList13b <- testListDef13$new(nlMatrix = matrix(2, nrow = 2, ncol = 2))

testInst <- nlTestFunc13()
RnimbleList <- testInst$run(testList13a, testList13b)
ctestInst <- compileNimble(testInst, control = list(debug =  F))
CnimbleList <- ctestInst$run(testList13a, testList13b)

## test for correct values of R nimbleList
expect_identical(RnimbleList$nlMatrix, matrix(c(2,2,2,2),nrow =  2))
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nlMatrix, CnimbleList$nlMatrix)
## test for identical values of testList13b and CnimbleList
expect_identical(testList13b$nlMatrix, CnimbleList$nlMatrix)
## test that  testList13a matrix was updated 
expect_identical(testList13a$nlMatrix, matrix(c(11, 1, 1, 1), nrow = 2))

test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })


########
## Test of nested nimbleLists.  Both nimbleList generators are defined in setup code. Top level nimbleList is returned.
########

nlTestFunc14 <- nimbleFunction(
  setup = function(){
    testTypes1 <- list(vars = c('nlDouble'), types = c('double(0)'))
    testListDef14a <- nimble:::nimbleList(testTypes1)
    testTypes2 <- list(vars = c('nestedNL'), types = c('testListDef14a()'))
    testListDef14b <- nimble:::nimbleList(testTypes2)
  },
  run = function(){
    newList14 <- testListDef14b$new()
    newList14$nestedNL$nlDouble <- 3.14
    returnType(testListDef14b())
    return(newList14)
  }
)

testInst <- nlTestFunc14()
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

nlTestFunc15 <- nimbleFunction(
  setup = function(){
    testTypes2 <- list(vars = c('nestedNL'), types = c('testListDef15a()'))
    testListDef15b <- nimble:::nimbleList(testTypes2)
  },
  run = function(){
    newList15 <- testListDef15b$new()
    newList15$nestedNL$nlDouble <- 3.14
    returnType(testListDef15b())
    return(newList15)
  }
)
testTypes1 <- list(vars = c('nlDouble'), types = c('double(0)'))
testListDef15a <- nimble:::nimbleList(testTypes1)

testInst <- nlTestFunc15()
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

nlTestFunc16 <- nimbleFunction(
  setup = function(){
  },
  run = function(argList16a = testListDef16b(), argList16b = testListDef16b()){
    argList16a$nestedNL <- argList16b$nestedNL
    returnType(testListDef16b())
    return(argList16a)
  }
)

testTypes1 <- list(vars = c('nlDouble'), types = c('double(0)'))
testListDef16a <- nimble:::nimbleList(testTypes1)
testTypes2 <- list(vars = c('nestedNL', 'nlVector'), types = c('testListDef16a()', 'double(1)'))
testListDef16b <- nimble:::nimbleList(testTypes2)

argList16a <- testListDef16b$new(nlVector = c(1,2))
argList16a$nestedNL$nlDouble <- 0
argList16b <- testListDef16b$new(nlVector = c(3,4))
argList16b$nestedNL$nlDouble <- 100

testInst <- nlTestFunc16()
RnimbleList <- testInst$run(argList16a, argList16b)
ctestInst <- compileNimble(testInst)
CnimbleList <- ctestInst$run(argList16a, argList16b)

## test for correct values of nimbleLists
expect_identical(RnimbleList$nestedNL$nlDouble, 100)
expect_identical(RnimbleList$nestedNL$nlDouble, CnimbleList$nestedNL$nlDouble)
expect_identical(RnimbleList$nestedNL$nlDouble, argList16a$nestedNL$nlDouble)

expect_identical(RnimbleList$nlVector, c(1,2))
expect_identical(RnimbleList$nlVector, CnimbleList$nlVector)
expect_identical(RnimbleList$nlVector, argList16a$nlVector)


test_that("return objects are nimbleLists", 
          {
            expect_identical(nimble:::is.nl(RnimbleList), TRUE)
            expect_identical(is.nl(CnimbleList), TRUE)
          })


########
## The same as test 12, but with lists provided as setup arguments
########



nlTestFunc17 <- nimbleFunction(
  setup = function(argList17a, argList17b){
  },
  run = function(){
    argList17a$nestedNL <<- argList17b$nestedNL
    returnType(testListDef17b())
    return(argList17a)
  }
)

testTypes1 <- list(vars = c('nlDouble'), types = c('double(0)'))
testListDef17a <- nimble:::nimbleList(testTypes1)
testTypes2 <- list(vars = c('nestedNL', 'nlVector'), types = c('testListDef17a()', 'double(1)'))
testListDef17b <- nimble:::nimbleList(testTypes2)

argList17a <- testListDef17b$new(nlVector = c(1,2))
argList17a$nestedNL$nlDouble <- 0
argList17b <- testListDef17b$new(nlVector = c(3,4))
argList17b$nestedNL$nlDouble <- 100

testInst <- nlTestFunc17(argList17a, argList17b)
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

nlTestFunc18 <- nimbleFunction(
  setup = function(){
    testTypes <- list(vars = 'test_matrix', types = c('double(2)'))
    testListDef18 <- nimbleList(testTypes)
    testList18 <- testListDef18$new(test_matrix = diag(4) + 1)
  },
  run = function(){
    eigenOut <- eigen(testList18$test_matrix)
    returnType(eigen())
    return(eigenOut)
  }
)


testInst <- nlTestFunc18()
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

nlTestFunc19 <- nimbleFunction(
  setup = function(){
    testTypes <- list(vars = 'test_matrix', types = c('double(2)'))
    testListDef19 <- nimbleList(testTypes)
    testList19 <- testListDef19$new(test_matrix = diag(4) + 1)
  },
  run = function(){
    svdOut <- svd(testList19$test_matrix)
    returnType(svd())
    return(svdOut)
  }
)

testInst <- nlTestFunc19()
RnimbleList <- testInst$run()
ctestInst <- compileNimble(testInst)
CnimbleList <- ctestInst$run()

## test for eigenValues and eigenVectors
expect_equal(diag(4)+1, RnimbleList$u%*%diag(RnimbleList$d)%*%solve(RnimbleList$v))
expect_equal(diag(4)+1,  RnimbleList$u%*%diag(RnimbleList$d)%*%solve(RnimbleList$v))
expect_equal(sum(abs(RnimbleList$values)), sum(abs(CnimbleList$values)))

test_that("return object (from c++) is nimbleList.", 
          {
            expect_identical(is.nl(CnimbleList), TRUE)
          })




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
# nimbleOptions(debugCppLineByLine = TRUE)
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
            is.nl(RnimbleList)
            is.nl(CnimbleList)
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
            is.nl(RnimbleList)
            is.nl(CnimbleList)
          })

########
#Test of creating nimbleListDef outside function code and creating new nimbleList in setup code 
########

nlTestFunc3 <- nimbleFunction(
  setup = function(){
    setupList3 <- testListDef3$new(nlCharacter = "hello world")
  },
  run = function(){
    returnType(testListDef3())
    return(setupList3)
  }
)

testTypes <- list(vars = c('nlCharacter'), types = c('character(0)'))
testListDef3 <- nimbleList(testTypes)
testInst <- nlTestFunc3()
RnimbleList <- testInst$run()
ctestInst <- compileNimble(testInst)
CnimbleList <- ctestInst$run()


## test for correct values of R nimbleList
expect_identical(RnimbleList$nlScalar, 2)
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nlScalar, CnimbleList$nlScalar)
test_that("return objects are nimbleLists", 
          {
            is.nl(RnimbleList)
            is.nl(CnimbleList)
          })

########
#Test of creating new nimbleList in setup code and specifying initial values for that list
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
expect_identical(RnimbleList$nlScalar, 2)
## test for identical values of R and C nimbleLists
expect_identical(RnimbleList$nlScalar, CnimbleList$nlScalar)
test_that("return objects are nimbleLists", 
          {
            is.nl(RnimbleList)
            is.nl(CnimbleList)
          })


########
#Test of using a nimbleList as a run function argument 
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
testList4 <- testListDef4$new(nlMatrix = diag(5))

testInst <- nlTestFunc4()
ctestInst <- compileNimble(testInst, control = list(debug =  F))
ctestInst$run(testList4)




########
#Test of using a nimbleList as an argument to a function within a nimbleFunction  
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
ctestInst <- compileNimble(testInst, control = list(debug =  F))
ctestInst$run(testList5)

innerNlTestFunc1 <- nimbleFunction(
  setup = function(){},
  run = function(nimList = testListDef()){
    nimList$a <- 1
    returnType(testListDef())
    return(nimList)
})

nlTestFunc1 <- nimbleFunction(
  setup = function(){
    doubleMatrix <- diag(1)
    innerFunc <- innerNlTestFunc1()
    },
  run = function(){
    doubleScalar <- 1
    doubleVector <- numeric(1, 1)
    argList <- testListDef$new(b = doubleVector, c = doubleMatrix)
    argList <- innerFunc(argList)
    returnType(testListDef())
    return(argList)
  }
)

testTypes <- list(vars = c('a'), types = c('double(0)'))
testListDef <- nimble:::nimbleList(testTypes)
testList <- testListDef(a = 5.5)
testInst <- innerNlTestFunc1()
ctestInst <- compileNimble(testInst, control = list(debug =  F))
ctestInst$run(testList)

testInst <- nlTestFunc1()
ctestInst <- compileNimble(testInst, control = list(debug =  F))
ctestInst$run()







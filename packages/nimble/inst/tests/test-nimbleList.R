## Tests for building, compiling, and using nimbleList objects
## These use some of the same internals (accessors), so they are in the same testing file.
## Checks are made internally for uncompiled and compiled cases.  Then uncompiled and compiled outcomes are compared to check that they behaved identically.

context('nimbleList() tests')


library(nimble)
nimbleOptions(showCompilerOutput = TRUE)
# nimbleOptions(debugCppLineByLine = TRUE)
# nimbleOptions(debugSizeProcessing = TRUE)
# debug(nimble:::sizeAssignAfterRecursing)
options( warn = 1 )

########
#Test of creating new nimbleList in run code and specifying initial values for that list
########

innerNlTestFunc1 <- nimbleFunction(
  setup = function(){},
  run = function(nimList = testListDef()){
    nimList$a[1] <- 1
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

testTypes <- list(vars = c('a','b','c'), types = c('double(1)', 'double(1)', 'double(2)'))
testListDef <- nimble:::nimbleList(testTypes)


testInst <- nlTestFunc1()
ctestInst <- compileNimble(testInst, control = list(debug =  F))
ctestInst$run()







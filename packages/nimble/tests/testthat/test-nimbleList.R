## Tests for building, compiling, and using nimbleList objects
## There are five distinct ways that a nimbleList can be created for use in a nimbleFunction:
## 1) nimbleListDef is created outside of the nimbleFunction, nimbleListDef$new() used to create new nimbleList in setup code
## 2) nimbleListDef is created outside of the nimbleFunction, nimbleListDef$new() used to create new nimbleList in run code
## 3) nimbleListDef is created in setup code of the nimbleFunction, nimbleListDef$new() used to create new nimbleList in setup code
## 4) nimbleListDef is created in setup code of the nimbleFunction, nimbleListDef$new() used to create new nimbleList in run code
## 5) nimbleListDef is created outside of the nimbleFunction, nimbleListDef$new() used to create new nimbleList outside of the nimbleFunction,
##    nimbleList passed as argument to nimbleFunction

source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)


context('Testing nimbleLists')

########
## Test of creating new nimbleList in run code and specifying initial values for that list
## Here, the nlDef is created in the global environment, outside of setup code
########
test_that("nimbleList test 1: return objects are nimbleLists", {
    testListDef1 <- nimbleList(nlScalar = double(0), nlVector = double(1), nlMatrix = double(2))
    temporarilyAssignInGlobalEnv(testListDef1)
    
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
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})

########
## Test of creating new nimbleList in run code and specifying initial value for that list
## Here, the nlDef is created in setup code.  Inital value is an expression.
########
test_that("nimbleList test 2: return objects are nimbleLists", {
    nlTestFunc2 <- nimbleFunction(
        setup = function(){
            testListDef2 <- nimbleList(nlScalar = double(0))
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
    CtestInst <- compileNimble(testInst, control = list(debug =  F))
    CnimbleList <- CtestInst$run()
    
    ## test for correct values of R nimbleList
    expect_identical(RnimbleList$nlScalar, 2)
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$nlScalar, CnimbleList$nlScalar)
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})


########
## Test of creating new nimbleList in setup code and specifying initial values for that list.
## Here, the nlDef is created in setup code
########
test_that("nimbleList test 3: return objects are nimbleLists", {
    nlTestFunc3 <- nimbleFunction(
        setup = function(){
            testTypes <- list(nimbleType(name = 'nlCharacter', type = 'character', dim = 0))
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
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    ## test for correct values of R nimbleList
    expect_identical(RnimbleList$nlCharacter, "hello world")
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$nlCharacter, CnimbleList$nlCharacter)
    expect_identical(CnimbleList$nlCharacter, CtestInst$setupList3$nlCharacter)
    
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
            # expect_identical(is.nl(CtestInst$setupList3), TRUE)  is.nl gives FALSE here, not sure if should be corrected
})


########
## Test of using a character vector in a nimbleList.  nimbleList created in setup code
########
test_that("nimbleList test 3a: return objects are nimbleLists", {
    nlTestFunc3a <- nimbleFunction(
        setup = function(){
            testTypes <- list(nimbleType(name = 'nlCharacters', type = 'character', dim = 1))
            testListDef3a <- nimbleList(testTypes)
            setupList3a <- testListDef3a$new(nlCharacters = c("hello", "world"))
        },
        run = function(){
            returnType(testListDef3a())
            return(setupList3a)
        }
    )
    
    testInst <- nlTestFunc3a()
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    ## test for correct values of R nimbleList
    expect_identical(RnimbleList$nlCharacters[1], "hello")
    expect_identical(RnimbleList$nlCharacters[2], "world")
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$nlCharacters, CnimbleList$nlCharacters)
    expect_identical(CnimbleList$nlCharacters, CtestInst$setupList3a$nlCharacters)
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})

########
## Test of using a character vector in a nimbleList.  nimbleList created in run code
########
test_that("nimbleList test 3b: return objects are nimbleLists", {
    nlTestFunc3b <- nimbleFunction(
        setup = function(){
            testTypes <- list(nimbleType(name = 'nlCharacters', type = 'character', dim = 1))
            testListDef3b <- nimbleList(testTypes)
            charVec <- c('hello', 'world')
        },
        run = function(){
            runList3b <- testListDef3b$new(nlCharacters = charVec)
            returnType(testListDef3b())
            return(runList3b)
        }
    )
    
    testInst <- nlTestFunc3b()
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    ## test for correct values of R nimbleList
    expect_identical(RnimbleList$nlCharacters[1], "hello")
    expect_identical(RnimbleList$nlCharacters[2], "world")
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$nlCharacters, CnimbleList$nlCharacters)
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
    ## expect_identical(is.nl(CtestInst$setupList3), TRUE)  is.nl gives FALSE here, not sure if should be corrected
})



########
## Test of creating new nimbleList in a function that is internal to another function.
## The list is passed from the internal function to the outer function, and then to the user.
## This test does not use nimble's active binding system.
########
test_that("nimbleList test 4: return objects are nimbleLists", {
    testListDef4 <- nimbleList(list(nimbleType(name = 'nlCharacter', type = 'character', dim = 0)))
    temporarilyAssignInGlobalEnv(testListDef4)
    
    innerNlTestFunc4 <- nimbleFunction(
        setup = function(){
            setupList4 <- testListDef4$new(nlCharacter = "hello world")
        },
        run = function(){
            returnType(testListDef4())
            return(setupList4)
        }
    )
    temporarilyAssignInGlobalEnv(innerNlTestFunc4)
    
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
    
    testInst <- nlTestFunc4()
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    ## test for correct values of R nimbleList
    expect_identical(RnimbleList$nlCharacter, "hello world")
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$nlCharacter, CnimbleList$nlCharacter)
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})

########
## Test of creating and interacting with two nimble lists in two different functions from the same def'n
########
test_that("nimbleList test 5: return objects are nimbleLists", {
    testListDef5 <- nimbleList(nlCharacter = character(0))
    temporarilyAssignInGlobalEnv(testListDef5)
    
    innerNlTestFunc5 <- nimbleFunction(
        setup = function(listDef){
            innerSetupList5 <- listDef$new(nlCharacter = "hello world")
        },
        run = function(){
            returnType(testListDef5())
            return(innerSetupList5)
        }
    )
    temporarilyAssignInGlobalEnv(innerNlTestFunc5)
    
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
    
    testInst <- nlTestFunc5(testListDef5)
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    ## test for correct values of R nimbleList
    expect_identical(RnimbleList$nlCharacter, "goodbye!")
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$nlCharacter, CnimbleList$nlCharacter)
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})

########
## Test of creating an nlDef in an outer function, passing the def to the inner function,
## and returning the created list from the inner function through the outer function
########
test_that("nimbleList test 6: return objects are nimbleLists", {
    innerNlTestFunc6 <- nimbleFunction(
        setup = function(nimListDef){
            setupList6 <- nimListDef$new(nlCharacter = "hello world")
        },
        run = function(){
            returnType(nimListDef())
            return(setupList6)
        }
    )
    
    temporarilyAssignInGlobalEnv(innerNlTestFunc6)
    nlTestFunc6 <- nimbleFunction(
        setup = function(){
            testTypes <- list(nimbleType('nlCharacter', 'character', 0))
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
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    ## test for correct values of R nimbleList
    expect_identical(RnimbleList$nlCharacter, "hello world")
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$nlCharacter, CnimbleList$nlCharacter)
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})
    
######
## test of a nimbleFunction method returning a nimbleList
######
test_that("nimbleList test 7: return objects are nimbleLists", {
    nlTestFunc7 <- nimbleFunction(
        setup = function(){
            testListDef7 <- nimbleList(list(nimbleType('nlVector', 'double', 1)))
        },
        run = function(){
            outList <- method1()
            outList$nlVector[1] <- method2(outList$nlVector*2)[1]
            returnType(testListDef7())
            return(outList)
        },
        methods = list(
            method1 = function(){
                outList <- testListDef7$new(nlVector = numeric(5, 1))
                returnType(testListDef7())
                return(outList)
            },
            method2 = function(x = double(1)){
                x2 <- x*2
                return(x2)
                returnType(double(1))
            }
        )
    )
    
    testInst <- nlTestFunc7()
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    ## test for correct values of R nimbleList
    expect_identical(CnimbleList$nlVector, c(4,1,1,1,1))
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$nlVector, CnimbleList$nlVector)
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})

test_that("nimbleList test 7a: return objects are nimbleLists", {
    innerNlTestFunc7a <- nimbleFunction(
        setup = function(nlDef){
            nlDef2 <- nimbleList(list(nimbleType('nlDouble', 'double', 0)))
        },
        methods = list(
            method1 = function(){
                doubleReturn <- nlDef2$new(nlDouble = 2 - 1)$nlDouble + 1
                returnType(double(0))
                return(doubleReturn)
            }
        )
    )
    temporarilyAssignInGlobalEnv(innerNlTestFunc7a)
    
    nlTestFunc7a <- nimbleFunction(
        setup = function(){
            testListDef7 <- nimbleList(list(nimbleType('nlVector', 'double', 1)))
            innerFunc7a <- innerNlTestFunc7a(testListDef7) 
        },
        methods = list(
            method1 = function(){
                outList <- testListDef7$new(nlVector = numeric(5, 1))
                outList$nlVector[1] <- innerFunc7a$method1()*2
                returnType(testListDef7())
                return(outList)
            },
            method2 = function(x = double(1)){
                x2 <- x*2
                return(x2)
                returnType(double(1))
            }
        )
    )
    
    testInst <- nlTestFunc7a()
    RnimbleList <- testInst$method1()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$method1()
    
    ## test for correct values of R nimbleList
    expect_identical(CnimbleList$nlVector, c(4,1,1,1,1))
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$nlVector, CnimbleList$nlVector)
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})

test_that("nimbleList test 7b: return objects are nimbleLists", {
    nlTestFunc7b <- nimbleFunction(
        setup = function(){
            testListDef7 <- nimbleList(list(nimbleType('nlDouble', 'double', 0)))
        },
        methods = list(
            method1 = function(){
                outList <- testListDef7$new(nlDouble = sum(numeric(5, 1))+ testListDef7$new(nlDouble = 1)$nlDouble)
                outList$nlDouble <- rep(10,2)[1] + outList$nlDouble
                returnType(testListDef7())
                return(outList)
            }
        )
    )
    
    testInst <- nlTestFunc7b()
    RnimbleList <- testInst$method1()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$method1()
    
    ## test for correct values of R nimbleList
    expect_identical(CnimbleList$nlDouble, 16)
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$nlDouble, CnimbleList$nlDouble)
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})



########
## Test of using a nimbleList as a run function argument 
########

test_that("nimbleList test 8: return objects are nimbleLists", {
    testListDef8 <- nimbleList(nlMatrix = double(2))
    temporarilyAssignInGlobalEnv(testListDef8)
    
    nlTestFunc8 <- nimbleFunction(
        setup = function(){
        },
        run = function(argList8 = testListDef8()){
            argList8$nlMatrix[1,2] <- 10
            returnType(testListDef8())
            return(argList8)
        }
    )
    
    testList8 <- testListDef8$new(nlMatrix = diag(2))
    
    testInst <- nlTestFunc8()
    RnimbleList <- testInst$run(testList8)
    CtestInst <- compileNimble(testInst, control = list(debug =  F))
    CnimbleList <- CtestInst$run(testList8)
    
    ## test for correct values of R nimbleList
    expect_identical(RnimbleList$nlMatrix, matrix(c(1,0,10,1),nrow =  2))
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$nlMatrix, CnimbleList$nlMatrix)
    ## test for identical values of testList4 and CnimbleList
    expect_identical(testList8$nlMatrix, CnimbleList$nlMatrix)
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})



########
## Test of using a nimbleList as an argument to a function within a nimbleFunction  
########

test_that("nimbleList test 9: return objects are nimbleLists", {
    testListDef9 <- nimbleList(a = double(0))
    temporarilyAssignInGlobalEnv(testListDef9)
    
    innerNlTestFunc9 <- nimbleFunction(
        setup = function(){},
        run = function(argList9 = testListDef9()){
            argList9$a <- argList9$a*2
            returnType(testListDef9())
            return(argList9)
        })
    temporarilyAssignInGlobalEnv(innerNlTestFunc9)
    
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
    
    testList9 <- testListDef9$new(a = 0)
    
    testInst <- nlTestFunc9()
    RnimbleList <- testInst$run(testList9)
    CtestInst <- compileNimble(testInst, control = list(debug =  F))
    CnimbleList <- CtestInst$run(testList9)
    
    ## test for correct values of R nimbleList
    expect_identical(RnimbleList$a, 2)
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$a, CnimbleList$a)
    ## test for identical values of testList5 and CnimbleList
    expect_identical(testList9$a, CnimbleList$a)
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})


########
## Test of using a nimbleList returned from a nimbleFunction as a run argument to another nimble function
########

test_that("nimbleList test 10: return objects are nimbleLists", {
    testListDef10 <- nimbleList(nlMatrix = double(2))
    temporarilyAssignInGlobalEnv(testListDef10)
    
    nlTestFunc10 <- nimbleFunction(
        setup = function(){
        },
        run = function(argList10 = testListDef10()){
            argList10$nlMatrix[1,2] <- argList10$nlMatrix[1,2] + 10
            returnType(testListDef10())
            return(argList10)
        }
    )
    
    testList10 <- testListDef10$new(nlMatrix = diag(2))
    
    testInst <- nlTestFunc10()
    RnimbleList <- testInst$run(testList10)
    RnimbleList <- testInst$run(RnimbleList)
    CtestInst <- compileNimble(testInst, control = list(debug =  F))
    CnimbleList <- CtestInst$run(testList10)
    CnimbleList <- CtestInst$run(CnimbleList)
    
    ## test for correct values of R nimbleList
    expect_identical(RnimbleList$nlMatrix, matrix(c(1,0,40,1),nrow =  2))
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$nlMatrix, CnimbleList$nlMatrix)
    ##test for identical values of testList6 and CnimbleList
    expect_identical(testList10$nlMatrix, CnimbleList$nlMatrix)
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})


########
## Test of using two nimbleLists as run function arguments 
########

test_that("nimbleList test 11: return objects are nimbleLists", {
    testTypes <- list(nimbleType('nlMatrix', 'double', 2))
    testListDef11 <- nimble:::nimbleList(testTypes)
    temporarilyAssignInGlobalEnv(testListDef11)
    
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
    
    testList11a <- testListDef11$new(nlMatrix = matrix(1, nrow = 2, ncol = 2))
    testList11b <- testListDef11$new(nlMatrix = matrix(2, nrow = 2, ncol = 2))
    
    testInst <- nlTestFunc11()
    RnimbleList <- testInst$run(testList11a, testList11b)
    CtestInst <- compileNimble(testInst, control = list(debug =  F))
    CnimbleList <- CtestInst$run(testList11a, testList11b)
    
    ## test for correct values of R nimbleList
    expect_identical(RnimbleList$nlMatrix, matrix(c(16,16,16,16),nrow =  2))
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$nlMatrix, CnimbleList$nlMatrix)
    ## test for identical values of testList7a and CnimbleList
    expect_identical(testList11a$nlMatrix, CnimbleList$nlMatrix)
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})


########
## Test of assigning a new nimbleList object to a nimbleList object that was passed as a run argument
########

test_that("nimbleList test 12: return objects are nimbleLists", {
    testListDef12 <- nimble:::nimbleList(nlMatrix = double(2))
    testList12 <- testListDef12$new(nlMatrix = matrix(1, nrow = 2, ncol = 2))
    temporarilyAssignInGlobalEnv(testListDef12)
    
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
    
    testListDef12 <- nimble:::nimbleList(nlMatrix = double(2))
    testList12 <- testListDef12$new(nlMatrix = matrix(1, nrow = 2, ncol = 2))
    
    testInst <- nlTestFunc12()
    RnimbleList <- testInst$run(testList12)
    CtestInst <- compileNimble(testInst, control = list(debug =  F))
    CnimbleList <- CtestInst$run(testList12)
    
    ## test for correct values of R nimbleList
    expect_identical(RnimbleList$nlMatrix, diag(2))
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$nlMatrix, CnimbleList$nlMatrix)
    ## test that testList12 doesn't get changed by function
    expect_identical(testList12$nlMatrix, matrix(1, nrow = 2, ncol = 2))
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})



########
## Test of using two nimbleLists as  run function arguments, assigning one list to the other within the function
########

test_that("nimbleList test 13: return objects are nimbleLists",  {
    testListDef13 <- nimble:::nimbleList(nlMatrix = double(2))
    testList13a <- testListDef13$new(nlMatrix = matrix(1, nrow = 2, ncol = 2))
    testList13b <- testListDef13$new(nlMatrix = matrix(2, nrow = 2, ncol = 2))
    
    temporarilyAssignInGlobalEnv(testListDef13)
    
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
    
    
    testInst <- nlTestFunc13()
    RnimbleList <- testInst$run(testList13a, testList13b)
    CtestInst <- compileNimble(testInst, control = list(debug =  F))
    CnimbleList <- CtestInst$run(testList13a, testList13b)
    
    ## test for correct values of R nimbleList
    expect_identical(RnimbleList$nlMatrix, matrix(c(2,2,2,2),nrow =  2))
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$nlMatrix, CnimbleList$nlMatrix)
    ## test for identical values of testList13b and CnimbleList
    expect_identical(testList13b$nlMatrix, CnimbleList$nlMatrix)
    ## test that  testList13a matrix was updated 
    expect_identical(testList13a$nlMatrix, matrix(c(11, 1, 1, 1), nrow = 2))
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})


########
## Test of nested nimbleLists.  Both nimbleList generators are defined in setup code. Top level nimbleList is returned.
########

test_that("nimbleList test 14: return objects are nimbleLists", {
    nlTestFunc14 <- nimbleFunction(
        setup = function(){
            testListDef14a <- nimbleList(nlDouble = double(0))
            testTypes2 <- list(nimbleType('nestedNL', 'testListDef14a'))
            testListDef14b <- nimbleList(testTypes2)
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
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    ## test for correct values of R nimbleList
    expect_identical(RnimbleList$nestedNL$nlDouble, 3.14)
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$nestedNL$nlDouble, CnimbleList$nestedNL$nlDouble)
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})

########
## Test of nested nimbleLists.  Here, top level nimbleList generator is defined in setup code.
## Lower level nimbleList generator is defined in global environment.  Top level nimbleList is returned.
########

test_that("nimbleList test 15: return objects are nimbleLists", {
    testListDef15a <- nimbleList(nlDouble = double(0))
    temporarilyAssignInGlobalEnv(testListDef15a)
    
    nlTestFunc15 <- nimbleFunction(
        setup = function(){
            testListDef15b <- nimbleList(nestedNL = testListDef15a())
        },
        run = function(){
            newList15 <- testListDef15b$new()
            newList15$nestedNL$nlDouble <- 3.14
            returnType(testListDef15b())
            return(newList15)
        }
    )
    
    testInst <- nlTestFunc15()
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    ## test for correct values of R nimbleList
    expect_identical(RnimbleList$nestedNL$nlDouble, 3.14)
    ## test for identical values of R and C nimbleLists
    expect_identical(RnimbleList$nestedNL$nlDouble, CnimbleList$nestedNL$nlDouble)
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})

########
## Test of nested nimbleLists. Two nimbleLists with nested nimbleLists are given as run time arguments.  
## In run code, the nested list of the second argument is assigned to the nested list of the first argument.
########

test_that("nimbleList test 16: return objects are nimbleLists", {
    testTypes1 <- list(nimbleType('nlDouble', 'double', 0))
    testListDef16a <- nimbleList(testTypes1)
    testListDef16b <- nimbleList(nestedNL = testListDef16a(), nlVector = double(1))
    temporarilyAssignInGlobalEnv(testListDef16a)
    temporarilyAssignInGlobalEnv(testListDef16b)
    
    nlTestFunc16 <- nimbleFunction(
        setup = function(){
        },
        run = function(argList16a = testListDef16b(), argList16b = testListDef16b()){
            argList16a$nestedNL <- argList16b$nestedNL
            returnType(testListDef16b())
            return(argList16a)
        }
    )
    
    
    argList16a <- testListDef16b$new(nlVector = c(1,2))
    argList16a$nestedNL$nlDouble <- 0
    argList16b <- testListDef16b$new(nlVector = c(3,4))
    argList16b$nestedNL$nlDouble <- 200
    
    testInst <- nlTestFunc16()
    RnimbleList <- testInst$run(argList16a, argList16b)
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run(argList16a, argList16b)
    
    ## test for correct values of nimbleLists
    expect_identical(RnimbleList$nestedNL$nlDouble, 200)
    expect_identical(RnimbleList$nestedNL$nlDouble, CnimbleList$nestedNL$nlDouble)
    expect_identical(RnimbleList$nestedNL$nlDouble, argList16a$nestedNL$nlDouble)
    
    expect_identical(RnimbleList$nlVector, c(1,2))
    expect_identical(RnimbleList$nlVector, CnimbleList$nlVector)
    expect_identical(RnimbleList$nlVector, argList16a$nlVector)
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)

    ## now rerun function with different second argument and see if nl's are updated appropriately
    
    argList16b <- testListDef16b$new(nlVector = c(3,4))
    argList16b$nestedNL$nlDouble <- 500
    
    CnimbleList <- CtestInst$run(argList16a, argList16b)
    
    ## test for correct values of nimbleLists
    expect_identical(RnimbleList$nestedNL$nlDouble, 500)
    expect_identical(RnimbleList$nestedNL$nlDouble, CnimbleList$nestedNL$nlDouble)
    expect_identical(RnimbleList$nestedNL$nlDouble, argList16a$nestedNL$nlDouble)
    
    expect_identical(RnimbleList$nlVector, c(1,2))
    expect_identical(RnimbleList$nlVector, CnimbleList$nlVector)
    expect_identical(RnimbleList$nlVector, argList16a$nlVector)
    
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})



########
## The same as test 16, but with lists provided as setup arguments
########

test_that("nimbleList test 17: return objects are nimbleLists", {
    testListDef17a <- nimbleList(nlDouble = double(0))
    testListDef17b <- nimbleList(nestedNL = testListDef17a(), nlVector = double(1))
    temporarilyAssignInGlobalEnv(testListDef17a)
    temporarilyAssignInGlobalEnv(testListDef17b)
    
    nlTestFunc17 <- nimbleFunction(
        setup = function(argList17a, argList17b){
        },
        run = function(){
            argList17a$nestedNL <<- argList17b$nestedNL
            returnType(testListDef17b())
            return(argList17a)
        }
    )
    
    argList17a <- testListDef17b$new(nlVector = c(1,2))
    argList17a$nestedNL$nlDouble <- 0
    argList17b <- testListDef17b$new(nlVector = c(3,4))
    argList17b$nestedNL$nlDouble <- 100
    
    testInst <- nlTestFunc17(argList17a, argList17b)
    RnimbleList <- testInst$run()
    ## note - we expect warning here, since when testInst was run, the nestedNL's of the two argLists became the same, 
    ## thus will be added twice in compilation.
    expect_warning(CtestInst <- compileNimble(testInst),
                   "Adding a specialized nimbleList to a project to which it already belongs")
    CnimbleList <- CtestInst$run()
    
    ## test for correct values of nimbleLists
    expect_identical(RnimbleList$nestedNL$nlDouble, 100)
    expect_identical(RnimbleList$nestedNL$nlDouble, CnimbleList$nestedNL$nlDouble)
    expect_identical(RnimbleList$nestedNL$nlDouble, argList17a$nestedNL$nlDouble)
    
    expect_identical(RnimbleList$nlVector, c(1,2))
    expect_identical(RnimbleList$nlVector, CnimbleList$nlVector)
    expect_identical(RnimbleList$nlVector, argList17a$nlVector)
    
    expect_identical(nimble:::is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})


########
## Test of nested nimbleLists.  Here, the inner list is returned - goes through a slightly different size
## processing flow than the above examples, so important to test
########

test_that("nimbleList test 18: return objects are nimbleLists", {
    testListDef18a <- nimbleList(nlDouble = double(0))
    temporarilyAssignInGlobalEnv(testListDef18a)
    
    nlTestFunc18 <- nimbleFunction(
        setup = function(argList){
            testListDef18b <- nimbleList(nestedNL = argList())
        },
        run = function(){
            newList18 <- testListDef18b$new()
            newList18$nestedNL$nlDouble <- 3.14
            outList <- newList18$nestedNL
            returnType(argList())
            return(outList)
        }
    )
    testInst <- nlTestFunc18(testListDef18a)
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    expect_identical(RnimbleList$nlDouble, 3.14)
    expect_identical(RnimbleList$nlDouble, CnimbleList$nlDouble)
    
    expect_identical(is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})



########
## Test of nested nimbleLists. List is created in setup code and nested list is returned from run code.
########

test_that("nimbleList test 19: return objects are nimbleLists", {
    nlTestFunc19 <- nimbleFunction(
        setup = function(){
            testListDef19a <- nimbleList(nlDouble = double(0))
            testListDef19b <- nimbleList(nestedNL = testListDef19a())
            testList19 <- testListDef19b$new()
        },
        run = function(){
            testList19$nestedNL$nlDouble <<- 3.14
            returnType(testListDef19a())
            return(testList19$nestedNL)
        }
    )
    
    testInst <- nlTestFunc19()
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    expect_identical(RnimbleList$nlDouble, 3.14)
    expect_identical(RnimbleList$nlDouble, CnimbleList$nlDouble)
    expect_identical(CtestInst$testList19$nestedNL$nlDouble, 3.14)

    expect_identical(is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})



########
## Test of nested nimbleLists. Both lists are created in run code and nested list is returned from run code.
########

test_that("nimbleList test 20: return objects are nimbleLists", {
    nlTestFunc20 <- nimbleFunction(
        setup = function(){
            testListDef20a <- nimbleList(nlDouble = double(0))
            testListDef20b <- nimbleList(nestedNL = testListDef20a())
        },
        run = function(){
            testList20a <- testListDef20a$new(nlDouble = 3.14)
            testList20 <- testListDef20b$new(nestedNL = testList20a)
            returnType(testListDef20a())
            return(testList20$nestedNL)
        }
    )
    
    testInst <- nlTestFunc20()
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    expect_identical(RnimbleList$nlDouble, 3.14)
    expect_identical(RnimbleList$nlDouble, CnimbleList$nlDouble)
    
    expect_identical(is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})


########
## Test of double-nested nimbleLists. List is created in setup code and nested list is returned from run code.
########

test_that("nimbleList test 21: return objects are nimbleLists", {
    nlTestFunc21 <- nimbleFunction(
        setup = function(){
            testListDef21a <- nimbleList(nlDouble = double(0))
            testListDef21b <- nimbleList(nestedNL = testListDef21a())
            testListDef21c <- nimbleList(nestedNL = testListDef21b())
            testList21 <- testListDef21c$new()
        },
        run = function(){
            testList21$nestedNL$nestedNL$nlDouble <<- 3.14
            returnType(testListDef21a())
            return(testList21$nestedNL$nestedNL)
        }
    )
    
    testInst <- nlTestFunc21()
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    expect_identical(RnimbleList$nlDouble, 3.14)
    expect_identical(RnimbleList$nlDouble, CnimbleList$nlDouble)
    
    expect_identical(is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})


test_that("nimbleList test 21abc: return objects are nimbleLists", {
    nlTestFunc21a <- nimbleFunction(
        setup = function(){
            testListDef21a <- nimbleList(nlDouble = double(0))
            testListDef21b <- nimbleList(nestedNL = testListDef21a())
            testListDef21c <- nimbleList(nestedNL = testListDef21b())
        },
        run = function(){
            testList21a <- testListDef21c$new(nestedNL = testListDef21b$new(nestedNL = 
                                              testListDef21a$new(nlDouble =
                                                                     testListDef21a$new(nlDouble = 
                                                                                            sum(rep(1,5)))$nlDouble)))
            returnType(testListDef21a())
            return(testList21a$nestedNL$nestedNL)
        }
    )

    
    testInst <- nlTestFunc21a()
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    expect_identical(RnimbleList$nlDouble, 5)
    expect_identical(RnimbleList$nlDouble, CnimbleList$nlDouble)
    expect_identical(is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})



########
## Test of double-nested nimbleLists. Lists are created in run code using nested $new() calls.  One nlDef is created outside of setup.
########

test_that("nimbleList test 22: return objects are nimbleLists", {
    testListDef22a <- nimbleList(nlDouble = double(0))
    temporarilyAssignInGlobalEnv(testListDef22a)
    
    nlTestFunc22 <- nimbleFunction(
        setup = function(){
            testListDef22b <- nimbleList(nestedNL = testListDef22a())
            testListDef22c <- nimbleList(nestedNL = testListDef22b())
        },
        run = function(){
            testList22 <- testListDef22c$new(nestedNL = testListDef22b$new(nestedNL = testListDef22a$new(nlDouble = 3.14)))
            returnType(testListDef22a())
            return(testList22$nestedNL$nestedNL)
        }
    )
    
    testInst <- nlTestFunc22()
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    expect_identical(RnimbleList$nlDouble, 3.14)
    expect_identical(RnimbleList$nlDouble, CnimbleList$nlDouble)
    
    expect_identical(is.nl(RnimbleList), TRUE)
    expect_identical(is.nl(CnimbleList), TRUE)
})




########
## Test of using nlDef$new()$nlVar in an expression 
########

test_that("nimbleList Test 23", {
    testListDef23 <- nimbleList(nlDouble = double(0))
    temporarilyAssignInGlobalEnv(testListDef23)
    
    nlTestFunc23 <- nimbleFunction(
        setup = function(){},
        run = function(){
            piDigits <- testListDef23$new(nlDouble = 2.14)$nlDouble + 1
            returnType(double())
            return(piDigits)
        }
    )
    
    testInst <- nlTestFunc23()
    RpiDigits <- testInst$run()
    CtestInst <- compileNimble(testInst, control = list(debug = F))
    CpiDigits <- CtestInst$run()
    
    expect_identical(RpiDigits, 3.14)
    expect_identical(RpiDigits, CpiDigits)
})

########
## Test #1 for eigen() function.  Return an eigenList.
########

test_that("nimbleList test 24: return object (from C++) is nimbleList", {
    nlTestFunc24 <- nimbleFunction(
        setup = function(){
            testListDef24 <- nimbleList(test_matrix = double(2))
            testList24 <- testListDef24$new(test_matrix = diag(4) + 1)
        },
        run = function(){
            eigenOut <- eigen(testList24$test_matrix)
            returnType(eigenNimbleList())
            return(eigenOut)
        }
    )
    
    
    testInst <- nlTestFunc24()
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    ## test for eigenValues and eigenVectors
    expect_equal(diag(4)+1, RnimbleList$vectors%*%diag(RnimbleList$values)%*%solve(RnimbleList$vectors))
    expect_equal(CnimbleList$vectors, RnimbleList$vectors)
    expect_equal(CnimbleList$values, RnimbleList$values)
    expect_identical(is.nl(CnimbleList), TRUE)
})



########
## Test #1 for svd() function.  Return an svdList.
########


test_that("nimbleList test 25: return object (from C++) is nimbleList", {
    nlTestFunc25 <- nimbleFunction(
        setup = function(){
            testListDef25 <- nimbleList(test_matrix = double(2))
            testList25 <- testListDef25$new(test_matrix = diag(4) + 1)
        },
        run = function(){
            svdOut <- svd(testList25$test_matrix, 'thin')
            returnType(svdNimbleList())
            return(svdOut)
        }
    )
    
    testInst <- nlTestFunc25()
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    ## test for singluar values
    expect_equal(diag(4)+1, RnimbleList$u%*%diag(RnimbleList$d)%*%solve(RnimbleList$v))
    expect_equal(RnimbleList$u, CnimbleList$u)
    expect_equal(RnimbleList$d, CnimbleList$d)
    expect_equal(RnimbleList$v, CnimbleList$v)

    expect_identical(is.nl(CnimbleList), TRUE)
})

########
## Test #2 for eigen() function.  Use eigen() to specify a nl element.
########
test_that("nimbleList test 26: return object (from C++) is nimbleList", {
    nlTestFunc26 <- nimbleFunction(
        setup = function(){
            testListDef26 <- nimbleList(testEigen = eigenNimbleList())
            testList26 <- testListDef26$new()
            testMat <- diag(2)
        },
        run = function(){
            eigenOut <- eigen(testMat)
            testList26$testEigen <<- eigenOut
            returnType(testListDef26())
            return(testList26)
        }
    )
    
    testInst <- nlTestFunc26()
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    ## test for eigenValues and eigenVectors
    expect_equal(diag(2), RnimbleList$testEigen$vectors%*%diag(RnimbleList$testEigen$values)%*%solve(RnimbleList$testEigen$vectors))
    expect_equal(CnimbleList$testEigen$vectors, RnimbleList$testEigen$vectors)
    expect_equal(CnimbleList$testEigen$values, RnimbleList$testEigen$values)
    expect_identical(is.nl(CnimbleList), TRUE)
})


########
## Test #2 for svd() function.  Use nimSvd() to specify a nl element.
########

test_that("nimbleList test 27: return object (from C++) is nimbleList", {
    nlTestFunc27 <- nimbleFunction(
        setup = function(){
            testListDef27 <- nimbleList(testSvd = svdNimbleList())
            testList27 <- testListDef27$new()
            testMat <- diag(2)
        },
        run = function(){
            svdOut <- svd(testMat)
            testList27$testSvd <<- svdOut
            returnType(svdNimbleList())
            return(testList27$testSvd)
        }
    )
    
    testInst <- nlTestFunc27()
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    expect_equal(diag(2), RnimbleList$u%*%diag(RnimbleList$d)%*%solve(RnimbleList$v))
    expect_equal(RnimbleList$u, CnimbleList$u)
    expect_equal(RnimbleList$d, CnimbleList$d)
    expect_equal(RnimbleList$v, CnimbleList$v)
    expect_equal(testInst$testList27$testSvd$d, CtestInst$testList27$testSvd$d)
    expect_equal(testInst$testList27$testSvd$u, CtestInst$testList27$testSvd$u)
    expect_equal(testInst$testList27$testSvd$v, CtestInst$testList27$testSvd$v)
    
    expect_identical(is.nl(CnimbleList), TRUE)
})

########
## Test #3 for eigen() function.  Use eigen() with a non-symmetric matrix.
########


test_that("nimbleList test 28: return object (from C++) is nimbleList", {
    nlTestFunc28 <- nimbleFunction(
        setup = function(){
            testMat <- matrix(c(1.0, 2.0, 3.0, 4.0), nrow = 2)
        },
        run = function(){
            eigenOut <- eigen(testMat)
            returnType(eigenNimbleList())
            return(eigenOut)
        }
    )
    
    testInst <- nlTestFunc28()
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    ## test for eigenValues and eigenVectors
    expect_equal(matrix(c(1.0, 2.0, 3.0, 4.0), nrow = 2), RnimbleList$vectors%*%diag(RnimbleList$values)%*%solve(RnimbleList$vectors))
    expect_equal(CnimbleList$vectors, RnimbleList$vectors)
    expect_equal(CnimbleList$values, RnimbleList$values)
    expect_identical(is.nl(CnimbleList), TRUE)
})

########
## Test #3 for svd() function.  Use svd() with a non-symmetric matrix.
########

test_that("nimbleList Test 29", {
    nlTestFunc29 <- nimbleFunction(
        setup = function(){
            testMat <- matrix(c(1.0, 2.0, 3.0, 4.0), nrow = 2)
        },
        run = function(){
            svdOut <- svd(testMat)
            returnType(svdNimbleList())
            return(svdOut)
        }
    )
    
    testInst <- nlTestFunc29()
    RnimbleList <- testInst$run()
    CtestInst <- compileNimble(testInst)
    CnimbleList <- CtestInst$run()
    
    expect_equal(matrix(c(1.0, 2.0, 3.0, 4.0), nrow = 2), RnimbleList$u%*%diag(RnimbleList$d)%*%solve(RnimbleList$v))
    expect_equal(RnimbleList$u, CnimbleList$u)
    expect_equal(RnimbleList$d, CnimbleList$d)
    expect_equal(RnimbleList$v, CnimbleList$v)
})

########
## Testing for use of eigen() in BUGS models
########

test_that("nimEigen Model Test 1", {
    modelCode1 <- nimbleCode({
        meanVec[1:5] <- eigen(constMat[,])$values[1:5]
        covMat[1:5,1:5] <- eigen(constMat[1:5,])$vectors[1:5,1:5]%*%t(eigen(constMat[1:5,1:5])$vectors[,])
        y[] ~ dmnorm(meanVec[], cov = covMat[1:5,1:5])
    })
    set.seed(1)    
    data <- list(y = rnorm(5, 0, 1))
    constants <- list(constMat = diag(5) + 1)
    dims <- list(y = c(5), meanVec = c(5), constMat = c(5,5))
    
    Rmodel1 <- nimbleModel(modelCode1, data = data, dimensions = dims, constants = constants)
    Cmodel1 <- compileNimble(Rmodel1)
    expect_silent(Rmodel1$simulate())
    expect_silent(Cmodel1$simulate())
    expect_identical(Rmodel1$calculate(), Cmodel1$calculate())
})

  
########
## Testing for use of svd() in BUGS models
########

test_that("nimSvd Model Test 1", {
    modelCode2 <- nimbleCode({
        meanVec[1:5] <- svd(constMat[,]%*%constMat[,])$d[1:5]
        covMat[1:5,1:5] <- svd(constMat[1:5,])$u[1:5,1:5]
        y[] ~ dmnorm(meanVec[], cov = covMat[1:5,1:5])
    })
    set.seed(1)    
    data <- list(y = rnorm(5, 0, 1))
    constants <- list(constMat = diag(5))
    dims <- list(y = c(5), meanVec = c(5))
    
    Rmodel2 <- nimbleModel(modelCode2, data = data, dimensions = dims, constants = constants)
    Cmodel2 <- compileNimble(Rmodel2)

    expect_silent(Rmodel2$simulate())
    expect_silent(Cmodel2$simulate())
    expect_identical(Rmodel2$calculate(), Cmodel2$calculate())
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)

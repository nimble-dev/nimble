source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

# Might want to make a test_dslCheck function.

RwarnLevel <- options('warn')$warn
options(warn = 1)

context('Testing of nimbleFunction (DSL) code checking')

test_that("Test of DSL check of valid RCfunction", expect_silent(
test1 <- nimbleFunction(
    run = function(x = double(0)) {
        returnType(double())
        y = cos(x)
        tmp <- c(seq(0, 1, by = 0.1), rep(4, 5))
        tmp2 <- logical(5)
        tmp3 <- numeric(3)
        tmp4 <- identityMatrix(3)
        return(y)
    }
    )
))

test_that("Test of DSL check of valid nimbleFunction", expect_silent(
test2 <- nimbleFunction(
    setup = function(model, node) {
        calcNodes <- model$getDependencies(node)
    },
    run = function(x = double(0)) {
        returnType(double())
        y = model$calculate(calcNodes)
        return(y)
    }
    )
))

test_that("Test of DSL check of invalid RCfunction with R code present", expect_message(
test3 <- nimbleFunction(
    run = function(x = double(1), y = double(1)) {
        returnType(double(1))
        out <- lm(y ~ x)
        return(out$coefficients)
    }
    ), "Detected possible use of R functions"
))


test_that("Test of DSL check of invalid RCfunction with check turned off", expect_silent(
test3a <- nimbleFunction(
    run = function(x = double(1), y = double(1)) {
        returnType(double(1))
        out <- lm(y ~ x)
        return(out$coefficients)
    }, check = FALSE
    )
))

# we'd like to catch this but 'aa' could be a nf and don't want to catch that; at the moment the only thing we readily have access to at time of checking is the names of objects in setup, not their types
test_that("Test of DSL check of invalid nimbleFunction with R code present", expect_failure(expect_warning(
test4 <- nimbleFunction(
    setup = function() {
        aa <- function() { return(0) }
    },
    run = function(x = double(0)) {
        returnType(double(0))
        tmp <- aa()
        return(tmp)
    }
)
)))

test_that("Test of DSL check of valid nimbleFunction with other nimbleFunctions present", expect_silent(
test5 <- nimbleFunction(
    setup = function(model, target) {
        calcNodes <- model$getDependencies(target)
        my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
    },
    run = function(par = double(1)) {
        returnType(double(0))
        values(model, target) <<- par
        ans <- model$calculate(calcNodes)
        tmp = my_decideAndJump(0,0,0,0)
        return(tmp)
    }
)
))

test_that("Test of DSL check of valid nimbleFunction with other methods present", expect_silent(
test6 <- nimbleFunction(
    setup = function(model, target) {
        calcNodes <- model$getDependencies(target)
    },
    run = function(par = double(1)) {
        returnType(double(0))
        values(model, target) <<- par
        ans <- model$calculate(calcNodes)
        tmp = auxil()
        return(tmp)
    }, methods = list(
           auxil = function() {
               returnType(double(0))
               return(3)
           })
)
))

test_that("Test of DSL check of valid nimbleFunction using undefined RC function", expect_message(
test7 <- nimbleFunction(
    setup = function(model, target) {
        calcNodes <- model$getDependencies(target)
    },
    run = function(par = double(1)) {
        returnType(double(0))
        values(model, target) <<- par
        ans <- model$calculate(calcNodes)
        tmp = myRCfun()
        return(tmp)
    }), "For this nimbleFunction to compile"
))

test_that("Test of DSL check of valid nimbleFunction using RC function from setup", expect_silent(
test8 <- nimbleFunction(
    setup = function(model, target) {
        calcNodes <- model$getDependencies(target)
        myRCfun <- nimbleFunction(run = function() {})
    },
    run = function(par = double(1)) {
        returnType(double(0))
        values(model, target) <<- par
        ans <- model$calculate(calcNodes)
        tmp = myRCfun()
        return(tmp)
    })
))

myRCfun <- nimbleFunction(run = function() {})
assign('myRCfun', myRCfun, envir = .GlobalEnv)

test_that("Test of DSL check of valid nimbleFunction using RC function from R", expect_silent( 
test9 <- nimbleFunction(
    setup = function(model, target) {
        calcNodes <- model$getDependencies(target)
    },
    run = function(par = double(1)) {
        returnType(double(0))
        values(model, target) <<- par
        ans <- model$calculate(calcNodes)
        tmp = myRCfun()
        return(tmp)
    })
))

test_that("Test of DSL check of valid nimbleFunction with undefined nimbleList", expect_message(
test10 <- nimbleFunction(
  setup = function(){    
  },
  run = function(){
    outList <- method1()
    returnType(testListDef())
    return(outList)
  },
  methods = list(
    method1 = function(){
      outList <- testListDef$new(nlCharacter = 'hello world')
       returnType(testListDef())
      return(outList)
    }
  )
), "For this nimbleFunction to compile"
))

test_that("Test of DSL check of valid nimbleFunction with nimbleList in setup", expect_silent(
test11 <- nimbleFunction(
  setup = function(){
      testListDef <- nimbleList(nlScalar = double(0), nlVector = double(1), nlMatrix = double(2))
  },
  run = function(){
    outList <- method1()
    returnType(testListDef())
     return(outList)
  },
  methods = list(
    method1 = function(){
      outList <- testListDef$new(nlCharacter = 'hello world')
       returnType(testListDef())
      return(outList)
    }
  )
)
))

testListDef <- nimbleList(nlScalar = double(0), nlVector = double(1), nlMatrix = double(2))
assign('testListDef', testListDef, envir = .GlobalEnv)
test_that("Test of DSL check of valid nimbleFunction with nimbleList in R", expect_silent(
test12 <- nimbleFunction(
  setup = function(){
  },
  run = function(){
    outList <- method1()
    returnType(testListDef())
     return(outList)
  },
  methods = list(
    method1 = function(){
      outList <- testListDef$new(nlCharacter = 'hello world')
       returnType(testListDef())
      return(outList)
    }
  )
  )))

## Actually not any different than test7.
test_that("Test of DSL check of valid nimbleFunction with call to undefined nimbleFunction", expect_message(
test13 <- nimbleFunction(
    setup = function(model, target) {
        calcNodes <- model$getDependencies(target)
    },
    run = function(par = double(1)) {
        returnType(double(0))
        values(model, target) <<- par
        ans <- model$calculate(calcNodes)
        tmp = mynf()
        return(tmp)
    }), "For this nimbleFunction to compile"
))


test_that("Test of DSL check of valid nimbleFunction using nimbleFunction from setup", expect_silent({
    nfGen <- nimbleFunction(
    setup = function() {},
    run = function() {}
    );
test14 <- nimbleFunction(
    setup = function(model, target) {
        calcNodes <- model$getDependencies(target)
        mynf <- nfGen()
    },
    run = function(par = double(1)) {
        returnType(double(0))
        values(model, target) <<- par
        ans <- model$calculate(calcNodes)
        tmp = mynf()
        return(tmp)
    })
}
))

test_that("Test of DSL check of valid nimbleFunction with calls to undefined nimbleFunctions and list elements", {
    nf_base= nimbleFunctionVirtual(
        name = 'nfbase',
        run = function() {returnType(double(0))},
        methods = list(
            foo = function() { }
        )
    )

    nfextend = nimbleFunction(
        contains = nf_base,
        setup = function() {},
        run = function() {
            returnType(double(0))
            return(3)
        },
        methods = list(
            foo2 = function() {},
            foo3 = function() {}
        )
    )

    expect_silent(
        nf_more <- nimbleFunction(
            setup = function() {
                fxns = nimbleFunctionList(nf_base)
                fxns[[1]] <- nf_extend()
            },
            run = function() {
                returnType(double(0))
                fxns[[1]]$foo()
                fxns[[1]]$foo2()
                i <- 1
                fxns[[i]]$foo3()
                out = fxns[[1]]$run()
                fxns[[i]]$foo4()
                return(out)
            }))
})

test_that("Test of detection of nf methods in findMethods", {
    code <- quote({
        m$cc1
        m$cc2()
        m$cc3(a)
        m[[1]]$cc4
        m[[i]]$cc5
        m[[1]]$cc6()
        m[[i]]$cc7()
        m[[x+3]]$cc8
        m[[x+3]]$cc9()
        m[[1]]$cc10(mm$cc11)
        m[[1]]$cc12(mm$cc13())
        m[[mm$cc14]]$cc15()
        m[[mm$cc16()]]$cc17()
    })
    expr <- nimble:::RparseTree2ExprClasses(code)
    expect_identical(nimble:::findMethodsInExprClass(expr), c('cc2', 'cc3', 'cc6', 'cc7', 'cc9', 'cc10', 'cc12', 'cc13', 'cc15', 'cc17', 'cc16'))
})

test_that("Handling of negative indexing in nimbleFunction code", {
    test = nimbleFunction(
        run = function(x=double(2)) {
            nimPrint(sum(x[1:3, 6:1]))
        }
    )
    expect_error(cTest <- compileNimble(test), "negative indexing")

    test = nimbleFunction(
        run = function(x=double(2)) {
            for(i in 3:1) {}
        }
    )
    expect_error(cTest <- compileNimble(test), "negative indexing")

    ## We don't check run-time negative indexing.
    test = nimbleFunction(
        run = function(x=double(2), n1 = double(0), n2 = double(0)) {
            nimPrint(sum(x[1:3, n1:n2]))
        }
    )
    cTest <- compileNimble(test)
    expect_failure(expect_message(cTest(matrix(rnorm(16), 4, 4), 3, 2),
                   "Run-time negative indexing error"))
    
})

test_that("Test of error trapping for reduction functions", {
    nf <- nimbleFunction(
        run = function(x = double(0)) {
            returnType(double())
            y <- mean(x)
            return(0)
        }
    ) 
    expect_error(cnf <- compileNimble(nf), "does not support reduction")
    
    nf <- nimbleFunction(
        run = function(x = double(1)) {
            returnType(double())
            y <- mean(x)
            return(0)
        }
    ) 
    cnf <- compileNimble(nf)
})


test_that("nimbleFunctionVirtual works with abstract and non-abstract methods", {
  # methods in nimbleFunctionVirtual are by default *abstract*, which means they
  # are declared in the base class but are not defined. (In C++ they have "=0" at the
  # end of the declaration.)  Hence they can't be called.
  # This means that any derived class *must* provide a definition to become a
  # fully defined (non-abstract) class and thus be usable.
  #
  # In order to provide nimbleFunctionVirtual methods that are not abstract,
  # make an entry in the methodControl list setting abstract=FALSE,
  # as shown in baseNF below. 
  # non-abstract methods are limited to default do-nothing ({}) behavior.
  # This means that to do anything interesting, a method must still be defined in
  # a derived class.  But the do-nothing behavior allows the method to be
  # simply ignored in a derived class that doesn't want to define it.
  # This is a bare-bones feature that could be expanded in the future (e.g.
  # by a meaningful base class method definition instead of simply {}).
  baseNF <- nimbleFunctionVirtual(
    run = function(x = double()) {returnType(double())},
    methods = list(
      foo = function(x = double(1)) {returnType(double(1))},
      hw = function(x = double(1)) {a <- 1; returnType(void())}    ),
    methodControl = list(hw = list(required = FALSE))
  )

  dNFa <- nimbleFunction(
    contains = baseNF,
    setup = TRUE,
    run = function(x = double()) {
      return(x + 1)
      returnType(double())
    },
    methods = list(
      foo = function(x = double(1)) {
        return(x + 1)
        returnType(double(1))
      },
      hw = function(x = double(1)) {
        cat("hello world from a dNFa object.\n")
      }
    )
  )

  dNFb <- nimbleFunction(
    contains = baseNF,
    setup = TRUE,
    run = function(x = double()) {
      return(x + 2)
      returnType(double())
    },
    methods = list(
      foo = function(x = double(1)) {
        return(x + 2)
        returnType(double(1))
      } # NO how method provided.
    )
  )

  ## message capturing in testthat is a pain
  ## Check for warnings: foo was required but not provided
  res <- capture_messages(    dNFc <- nimbleFunction(
                                contains = baseNF,
                                setup = TRUE,
                                run = function(x = double()) {
                                  return(x + 2)
                                  returnType(double())
                                }
                              ) )
  expect_true(grepl("[Warning]", res))

  ## Check for no warning: run is required but defaults to function(){}, so no warning is issued
  res <- capture_messages(
    dNFc <- nimbleFunction(
                           contains = baseNF,
                           setup = TRUE,
                           methods = list(
                             foo = function(x = double(1)) {
                               return(x + 2)
                               returnType(double(1))
                             } # NO hw method provided, but not required
                           )
    )
  )
  expect_identical(res, character())

  useBase <- nimbleFunction(
    setup = function() {
      NFlist <- nimbleFunctionList(baseNF)
      NFlist[[1]] <- dNFa()
      NFlist[[2]] <- dNFb()
    },
    run = function(x = double(1)) {
      for(i in 1:2) {
        x[1] <- NFlist[[i]]$run(x[1])
        x <- NFlist[[i]]$foo(x)
        NFlist[[i]]$hw(x)
      }
      return(x)
      returnType(double(1))
    }
  )

  useBase1 <- useBase()
  CuseBase1 <- compileNimble(useBase1)
  res <- capture.output(ans <- CuseBase1$run( c(10, 100) ))
  expect_identical(ans, c(16, 103))
  expect_true(grepl("hello world", res)) 
})

test_that("Missing variables in nf code cause error (not warning)", {

    test1 <- nimbleFunction(
    run = function() {
        a  <- 7
        c <- a+b  # warning about 'b' not yet created
        return(c)
        returnType(double(0))
    })
    expect_error(compileNimble(test1), "is not available")

    test2 <- nimbleFunction(
    run = function() {
        a  <- 7
        return(a + b)
        returnType(double(0))
    })
    expect_error(compileNimble(test1), "is not available")
})


options(warn = RwarnLevel)


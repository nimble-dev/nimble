source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

test_that("nimbleFunction use of model variables and nodes works",
{
    code <- nimbleCode({
        a ~ dnorm(0, 1)
        for(i in 1:5) {
            b[i] ~ dnorm(0, 1)
            for(j in 1:3)
                c[i,j] ~ dnorm(0, 1)
        }
    })
    constants <- list()
    data <- list()
    inits <- list(a = 1, b=1:5, c = array(1:30, c(5,6)))

    Rmodel <- nimbleModel(code, constants, data, inits)

    nfDef <- nimbleFunction(
        setup = function(model) {
            aNode <- 'a'                ## CORRECT
            bNode <- 'b'
            cNode <- 'c'
            bNode2 <- 'b[2]'
            cNode34 <- 'c[3, 4]'
        },
        run = function(aCheck = double(),
                       bCheck = double(1),
                       cCheck = double(2)) {
            ok <- TRUE

            bNode2v <- model[[bNode2]]
            cNode34v <- model[[cNode34]]
            bNode2v2 <- model[['b[2]']]
            cNode34v2 <- model[['c[3, 4]']]

            ok <- ok & bNode2v == bCheck[2]
            ok <- ok & bNode2v == bNode2v2
            ok <- ok & cNode34v == cCheck[3, 4]
            ok <- ok & cNode34v2 == cNode34v
            
            a1 <- model[['a']]
            ok <- ok & a1 == aCheck
            b1 <- model[['b']]
            ok <- ok & sum(dim(b1) != dim(bCheck)) == 0
            ok <- ok & sum(b1 != bCheck) == 0
            c1 <- model[['c']]
            ok <- ok & sum(dim(c1) != dim(cCheck)) == 0
            ok <- ok & sum(c1 != cCheck) == 0

            ## This isn't a very creative check of assignment
            ## but at least it should compile.  If the value
            ## assignment fails, it should be revealed in the next
            ## checks.
            model[['a']] <<- a1 + 1
            model[['b']] <<- b1 + 1
            model[['c']] <<- c1 + 1
            aCheck <- aCheck + 1
            bCheck <- bCheck + 1
            cCheck <- cCheck + 1
            
            a2 <- model$a
            ok <- ok & sum(a2 != aCheck) == 0
            b2 <- model$b
            ok <- ok & sum(dim(b2) != dim(bCheck)) == 0
            ok <- ok & sum(b2 != bCheck) == 0
            c2 <- model$c
            ok <- ok & sum(dim(c2) != dim(cCheck)) == 0
            ok <- ok & sum(c2 != cCheck) == 0

            model$a <<- a2 + 1
            model$b <<- b2 + 1
            model$c <<- c2 + 1
            aCheck <- aCheck + 1
            bCheck <- bCheck + 1
            cCheck <- cCheck + 1
            bSubCheck <- bCheck[2:4]
            cSubCheck <- cCheck[2:4, 3:5]

            ## model$a is a vector, so model$a[1] is allowed
            a3 <- model$a[1]
            ok <- ok & a3 == aCheck
            b3 <- model$b[2:4]
            ok <- ok & sum(dim(b3) != dim(bSubCheck)) == 0
            ok <- ok & sum(b3 != bSubCheck) == 0
            c3 <- model$c[2:4 , 3:5]
            ok <- ok & sum(dim(c3) != dim(cSubCheck)) == 0
            ok <- ok & sum(c3 != cSubCheck) == 0

            model$a[1] <<- a3 + 1
            model$b[2:4] <<- bSubCheck + 1
            model$c[2:4, 3:5] <<- cSubCheck + 1
            aCheck <- aCheck + 1
            bCheck[2:4] <- bSubCheck + 1
            cCheck[2:4, 3:5] <- cSubCheck + 1
            bSubCheck <- bCheck[2:4]
            cSubCheck <- cCheck[2:4, 3:5]
            
            ## model[['a']] is a scalar, so model[['a']][1] is not allowed
            ##            a4 <- model[['a']][1]
            ##            ok <- ok & a4 == aCheck
            b4 <- model[['b']][2:4]
            ok <- ok & sum(dim(b4) != dim(bSubCheck)) == 0
            ok <- ok & sum(b4 != bSubCheck) == 0
            c4 <- model[['c']][2:4 , 3:5]
            ok <- ok & sum(dim(c4) != dim(cSubCheck)) == 0
            ok <- ok & sum(c4 != cSubCheck) == 0

            model[['b']][2:4] <<- bSubCheck + 1
            model[['c']][2:4, 3:5] <<- cSubCheck + 1
            bCheck[2:4] <- bSubCheck + 1
            cCheck[2:4, 3:5] <- cSubCheck + 1
            bSubCheck <- bCheck[2:4]
            cSubCheck <- cCheck[2:4, 3:5]
            
            a5 <- model[[aNode]]
            ok <- ok & a5 == aCheck
            b5 <- model[[bNode]]
            ok <- ok & sum(dim(b5) != dim(bCheck)) == 0
            ok <- ok & sum(b5 != bCheck) == 0
            c5 <- model[[cNode]]
            ok <- ok & sum(dim(c5) != dim(cCheck)) == 0
            ok <- ok & sum(c5 != cCheck) == 0

            model[[aNode]] <<- a5 + 1
            model[[bNode]] <<- b5 + 1
            model[[cNode]] <<- c5 + 1
            aCheck <- aCheck + 1
            bCheck <- bCheck + 1
            cCheck <- cCheck + 1
            bSubCheck <- bCheck[2:4]
            cSubCheck <- cCheck[2:4, 3:5]

            ## model[[aNode]] is a scalar, so model[[aNode]][1] is not allowed
            ## a6 <- model[[aNode]][1]
            ## ok <- ok & a6 == aCheck
            b6 <- model[[bNode]][2:4]
            ok <- ok & sum(dim(b6) != dim(bSubCheck)) == 0
            ok <- ok & sum(b6 != bSubCheck) == 0
            c6 <- model[[cNode]][2:4 , 3:5]
            ok <- ok & sum(dim(c6) != dim(cSubCheck)) == 0
            ok <- ok & sum(c6 != cSubCheck) == 0

            model[[bNode]][2:4] <<- bSubCheck + 1
            model[[cNode]][2:4, 3:5] <<- cSubCheck + 1
            
            return(ok)
            returnType(logical())
        }
    )
    
    Rnf <- nfDef(Rmodel)
    uncompiledAns <- Rnf$run(Rmodel$a, Rmodel$b, Rmodel$c)
    expect_true(uncompiledAns)
    compiled <- compileNimble(Rmodel, Rnf)
    compiledAns <- compiled$Rnf$run(Rmodel$a, Rmodel$b, Rmodel$c)
    expect_true(compiledAns)    
})

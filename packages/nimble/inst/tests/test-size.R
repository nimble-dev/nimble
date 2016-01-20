library(nimble, lib.loc='/tmp/nim05_size')
library(testthat)

#source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
#context("Testing of size/dimension checks in NIMBLE code")


fit_model <- function(input, check = TRUE, useConst = FALSE) {
    if(useConst) 
            m <- nimbleModel(code = input$expr, data = input$data, constants = input$inits, check = check)
    else m <- nimbleModel(code = input$expr, data = input$data, inits = input$inits, check = check)
}


test_size <- function(input, verbose = TRUE) {
    errorMsg <- paste0(ifelse(input$knownProblem, "KNOWN ISSUE: ", ""), "Result does not match ", input$expectPass)
    if(verbose) cat("### Testing", input$name, " with RHS variable ###\n")
    result <- try(
        m <- nimbleModel(code = input$expr, data = input$data, inits = input$inits)
    )    
    try(test_that(paste0("Test of size/dimension check: ", input$name),
                  expect_that(!is(result, "try-error"), equals(input$expectPass),
                              errorMsg)))
    if(!is(result, "try-error")) {
        result <- try(
            { calculate(m); out <- calculate(m)} )          
        try(test_that(paste0("Test of size/dimension check: ", input$name),
                      expect_that(!is(result, "try-error"), equals(input$expectPass),
                                  errorMsg)))
    }
    if(verbose) cat("### Testing", input$name, "with RHS constant ###\n")
    if(!is.null(input$expectPassWithConst)) input$expectPass <- input$expectPassWithConst
    result <- try(
        m <- nimbleModel(code = input$expr, data = input$data, constants = input$inits)
    )
    try(test_that(paste0("Test of size/dimension check: ", input$name),
                  expect_that(!is(result, "try-error"), equals(input$expectPass),
                              errorMsg)))
    if(!is(result, "try-error")) {
        result <- try(
            { calculate(m); out <- calculate(m)} )          
        try(test_that(paste0("Test of size/dimension check: ", input$name),
                      expect_that(!is(result, "try-error"), equals(input$expectPass),
                                  errorMsg)))
    }
    invisible(NULL)
}
    
### basic tests of scalar distribution


testsScalar <- list(
    list(name = 'scalar stochastic, basic', expectPass = TRUE,
         expr = quote({y ~ dnorm(mu, sd = sig)}), 
         inits = list(mu = 0, sig = 1) ),
    list(name = 'scalar stochastic, parameter expression', expectPass = TRUE,
         expr = quote({y ~ dnorm(mu1 + mu2, sd = sig)}), 
         inits = list(mu1 = 0, mu2 = 0, sig = 1) ),
    list(name = 'scalar stochastic, parameter expression non-scalar, no indices', expectPass = TRUE,
         expr = quote({y ~ dnorm(mu1%*%mu2, sd = sig)}), 
         inits = list(mu1 = c(1,1), mu2 = c(1,1), sig = 1) ),
    list(name = 'scalar stochastic, parameter expression non-scalar', expectPass = TRUE,
         expr = quote({y ~ dnorm(mu1[1:2]%*%mu2[1:2], sd = sig)}), 
         inits = list(mu1 = c(1,1), mu2 = c(1,1), sig = 1) ),
    list(name = 'scalar stochastic, non-scalar parameter', expectPass = FALSE,
         expr = quote({y ~ dnorm(mu1[1:2, 1:2]%*%mu2[1:2], sd = sig)}), 
         inits = list(mu1 = matrix(1:4, 2), mu2 = c(1,1), sig = 1) ),
    list(name = 'scalar stochastic, non-scalar value', expectPass = FALSE,
         expr = quote({y[1:2] ~ dnorm(mu, sd = sig)}), 
         inits = list(mu = 1, sig = 1) ),
    list(name = 'scalar stochastic, non-scalar value and parameter', expectPass = FALSE,
         expr = quote({y[1:2] ~ dnorm(mu1[1:2, 1:2]%*%mu2[1:2], sd = sig)}), 
         inits = list(mu1 = matrix(1:4, 2), mu2 = c(1,1), sig = 1) ),

    list(name = 'scalar stochastic, scalar within multivar variable', expectPass = FALSE,
         expr = quote({for(i in 1:1)
                           for(j in 1:1)
                               y[i,j] ~ dnorm(mu[i,j], sd = sig)}), 
         inits = list(mu = 0, sig = 1) ),
    # ERRORS in R but not C nodeFunction operation, passes 

    list(name = 'scalar stochastic, scalar within multivar variable 2', expectPass = TRUE,
         expr = quote({for(i in 1:1)
                           for(j in 1:1)
                               y[i,j] ~ dnorm(mu[i,j], sd = sig)}), 
         inits = list(mu = matrix(0, 1,1), sig = 1) ),
    # ERRORS when mu is constant (inconsistent dimensionality in addMissingIndexRecurse), gives warning when mu is RHS variable

    list(name = 'scalar stochastic, scalar within multivar variable 3', expectPass = TRUE,
         expr = quote({for(i in 1:1)
                           for(j in 1:1)
                               y[i,j] ~ dnorm(mu[i,j], sd = sig)}), 
         inits = list(mu = matrix(0, 2, 2), sig = 1) )
    # gives warning when mu is constant
)



# sapply(testsScalar, test_size)

testsMultivarParam <- list(
    list(name = 'mv param stochastic, basic', expectPass = TRUE,
         expr = quote({y ~ dcat(p[1:3])}), 
         inits = list(p = rep(1,3)/3) ),
    list(name = 'mv param stochastic, no indices', expectPass = FALSE,
         expr = quote({y ~ dcat(p)}), 
         inits = list(p = rep(1,3)/3) ),
    # with p as constant ERRORS in compileNimble(); not caught in size check, but replaceConstantsRecurse warning regarding dimensionality is given
    
    list(name = 'mv param stochastic, scalar param', expectPass = FALSE,
         expr = quote({y ~ dcat(p[1])}), 
         inits = list(p = 1) ),
    # with p as constant ERRORS while defining model with message regarding inconsistent dimensionality, caught in size check
    
    list(name = 'mv param stochastic, scalar param, no indices', expectPass = FALSE,
         expr = quote({y ~ dcat(p)}), 
         inits = list(p = 1) ),
    list(name = 'mv param stochastic, parameter expression', expectPass = TRUE,
         expr = quote({y ~ dcat(p1[1:3] + p2[1:3])}), 
         inits = list(p1 = rep(1, 3)/6, p2 = rep(1,3)/6 )),
    list(name = 'mv param stochastic, parameter expression, matrices', expectPass = FALSE,
         expr = quote({y ~ dcat(p1[1:3,1:3] %*% p2[1:3])}), 
         inits = list(p1 = matrix(rep(1, 9)/9, 3), p2 = rep(1,3))),
    # parameter is one-column matrix so fails check (and fails in compileNimble)

    ## list(name = 'mv param stochastic, parameter expression, matrices subindexing', expectPass = FALSE,
    ##     expr = quote({y ~ dcat((p1[1:3,1:3] %*% p2[1:3])[1:3,1])}), 
    ##     inits = list(p1 = matrix(rep(1, 9)/9, 3), p2 = rep(1,3))),
    # ERRORS OUT in genVarInfo3() - Error in BUGSdecl$symbolicParentNodes[[iV]] : subscript out of bounds

    list(name = 'mv param stochastic, non-scalar value', expectPass = FALSE,
         expr = quote({y[1:2] ~ dcat(p[1:3])}), 
         inits = list(p = rep(1,3)/3)),

    list(name = 'mv param stochastic, scalar value within multivar variable ', expectPass = TRUE,
         expr = quote({for(i in 1:1)
                           for(j in 1:1)
                               y[i,j] ~ dcat(p[1:3])}),
         inits = list(p = rep(1,3)/3)), 
    list(name = 'mv param stochastic, vector variable within matrix variable ', expectPass = TRUE,
         expr = quote({for(i in 1:1)
                               y[i] ~ dcat(p[1:3, i])}),
         inits = list(p = matrix(rep(1,3)/3, ncol = 1)))
)

# sapply(testsMultivarParam, test_size)

testsDeterm <- list(
    list(name = 'deterministic, basic', expectPass = TRUE,
         expr = quote({y <- a + b}),
         inits = list(a = 3, b = 3 )),
    list(name = 'deterministic, basic with node in variable, incorrect init size', expectPass = FALSE,
         expr = quote({for(i in 1:1)
                           for(j in 1:1)
                               y[i,j] <- a[i,j] + b}),
         inits = list(a = 3, b = 3 )),
    
    list(name = 'deterministic, basic with node in variable', expectPass = TRUE,
         expr = quote({for(i in 1:1)
                           for(j in 1:1)
                               y[i,j] <- a[i,j] + b}),
         inits = list(a = matrix(3, 1, 1), b = 3 )),
    # ERRORS when a is constant

    list(name = 'deterministic, basic with node in variable', expectPass = TRUE,
         expr = quote({for(i in 1:2)
                           for(j in 1:2)
                               y[i,j] <- a[i,j] + b}),
         inits = list(a = matrix(3, 2, 2), b = 3 )),

    list(name = 'deterministic, basic with node in variable', expectPass = TRUE,
         expr = quote({for(i in 1:1)
                           for(j in 1:1)
                               y[i,j] <- a[i,j] + b}),
         inits = list(a = matrix(3, 2, 2), b = 3 )),
    # gives warning when a is constant "In cbind(edgesFrom, edgesTo)"
    
    list(name = 'deterministic, non-scalar expression, no indices', expectPass = FALSE,
         expr = quote({y <- a %*% b}),
         inits = list(a = rep(1,2), b = rep(1,2) )),
    # doesn't catch that this will ERROR at compileNimble()
    
    list(name = 'deterministic, non-scalar expression', expectPass = TRUE,
         expr = quote({y <- a[1:2] %*% b[1:2]}),
         inits = list(a = rep(1,2), b = rep(1,2) )),

    list(name = 'deterministic, non-scalar expression, dimension mismatch', expectPass = FALSE,
         expr = quote({y <- a[1:2,1:2] %*% b[1:2]}),
         inits = list(a = diag(rep(1,2)), b = rep(1,2) )),

    # LHS vec
    list(name = 'deterministic, vector value, dimension mismatch', expectPass = FALSE,
         expr = quote({y[1:2] <- a + b}),
         inits = list(a = 3, b = 3)),
    
    list(name = 'deterministic, vector value, missing indices', expectPass = FALSE,
         expr = quote({y[1:2] <- a + b}),
         inits = list(a = c(0,0), b = c(0,0))),
    # test_size doesn't catch error

    list(name = 'deterministic, basic vector', expectPass = TRUE,
         expr = quote({y[1:2] <- a[1:2,1:2] %*% b[1:2]}),
         inits = list(a = diag(rep(1,2)), b = rep(1,2) )),
    
    list(name = 'deterministic, basic vector, missing indices', expectPass = FALSE,
         expr = quote({y[1:2] <- a %*% b[1:2]}),
         inits = list(a = diag(rep(1,2)), b = rep(1,2) )),
    # test_size doesn't catch error for RHS variable case
    
    list(name = 'deterministic, basic vector, dimension mismatch', expectPass = FALSE,
         expr = quote({y[1:2] <- a %*% b}),
         inits = list(a = 3, b = 2 )),
    list(name = 'deterministic, basic vector, RHS dimension mismatch', expectPass = FALSE,
         expr = quote({y[1:2] <- a[1:2] + b[1:2, 1:2]}),
         inits = list(a = rep(1,2), b = diag(rep(1,2)) )),
    list(name = 'deterministic, basic vector, size mismatch', expectPass = FALSE,
         expr = quote({y[1:3] <- a[1:2,1:2] %*% b[1:2]}),
         inits = list(a = diag(rep(1,2)), b = rep(1,2) )),
    list(name = 'deterministic, basic vector, size mismatch 2', expectPass = FALSE,
         expr = quote({y[1:2] <- a[1:3,1:3] %*% b[1:3]}),
         inits = list(a = diag(rep(1,3)), b = rep(1,3) )),
    list(name = 'deterministic, basic vector dimension mismatch', expectPass = FALSE,
         expr = quote({y[1:3] <- a[1:3,1:3] %*% b[1:3,1:3]}),
         inits = list(a = diag(rep(1,3)), b = diag(rep(1,3)) )),
    list(name = 'deterministic, basic matrix', expectPass = TRUE,
         expr = quote({y[1:3, 1:3] <- a[1:3,1:3] %*% b[1:3,1:3]}),
         inits = list(a = diag(rep(1,3)), b = diag(rep(1,3)) )),

    # nested nodes
    list(name = 'deterministic, nodes within multivar variables 1', expectPass = TRUE,
         expr = quote({
             for(i in 1:1)
                 y[1:2, i] <- a[1:2,1:2] %*% b[1:2, i]}),
         inits = list(a = diag(rep(1,2)), b = diag(rep(1,2)) )),
    list(name = 'deterministic, nodes within multivar variables, size mismatch', expectPass = FALSE,
         expr = quote({
             for(i in 1:1)
                 y[1:3] <- a[1:2,1:2] %*% b[1:2, i]}),
         inits = list(a = diag(rep(1,2)), b = rep(1,2) )),
    list(name = 'deterministic, nodes within multivar variables, dim mismatch', expectPass = FALSE,
         expr = quote({
             for(i in 1:1)
                 y[1:2, i] <- a[1:2,1:2] %*% b[1:2, 1:2]}),
         inits = list(a = diag(rep(1,2)), b = diag(rep(1,2)) )),
    list(name = 'deterministic, nodes within multivar variables, input wrong dimension', expectPass = FALSE,
         expr = quote({
             for(i in 1:1)
                 y[1:2, i] <- a[1:2,1:2] %*% b[1:2, 1:2]}),
         inits = list(a = diag(rep(1,2)), b = rep(1,2) ))
)

# sapply(testsDeterm, test_size)

#HERE

testsMultivar <- list(
)

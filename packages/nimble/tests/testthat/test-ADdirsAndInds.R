# Tests of controlling nimDerivs behavior with
# combinations of wrt, outInds, inDir, and/or outDir.
source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)

test_that('Derivatives with double-taping work',
  {
    ADfun1 <- nimbleFunction(
      setup = function(){},
      run = function(){},
      methods = list(
        jacobianRun = function(x = double(1)) {
          out <- derivs(zo(x), wrt = 1:2, order = 1)
          ans <- nimNumeric(length = 6, value = out$jacobian)
          returnType(double(1))
          return(ans)
        },
        zo = function(x = double(1)) {
          ans <- c(x[1]^3, x[2]^3, (x[1]^2*x[2]^2))
          return(ans)
          returnType(double(1))
        },
        metaDerivsRun = function(x = double(1),
                                 order = double(1)) {
          out <- derivs(jacobianRun(x), wrt = 1:2,
                        order = order)
          returnType(ADNimbleList())
          return(out)
        }
      ), buildDerivs = c('jacobianRun', 'zo')
    )

    ADfunInst <- ADfun1()
    xRec <- c(2.2, 3.3)
    x <- c(1.6, 2.8)
    Rderivs <- ADfunInst$jacobianRun(x)
    Rderivs <- ADfunInst$metaDerivsRun(x, order = 1)
    temporarilyAssignInGlobalEnv(ADfunInst)
    cADfunInst <- compileNimble(ADfunInst)
    cADfunInst$metaDerivsRun(xRec, 1)
    cderivs <- cADfunInst$metaDerivsRun(x, 1)
    expect_equal(cderivs$value, Rderivs$value)
    expect_equal(cderivs$jacobian, Rderivs$jacobian, tolerance = 1e-5)
    expect_equal(cderivs$hessian, Rderivs$hessian)
  }
)


test_that('Derivatives with double-taping, indices, and directions work',
  {
    ADfun1 <- nimbleFunction(
      setup = function(){},
      run = function(x = double(1)) {
        ans <- c(x[1]^3, x[2]^3, (x[1]^2*x[2]^2))
        return(ans)
        returnType(double(1))
      },
      methods = list(
        derivsRun = function(x = double(1),
                             wrt = double(1),
                             outInds = double(1),
                             inDir = double(1),
                             outDir = double(1),
                             order = double(1)) {
          ans <- derivs(run(x), wrt = wrt, outInds = outInds,
                        inDir = inDir, outDir = outDir, order = order)
          return(ans)
          returnType(ADNimbleList())
        },
        jacobianRun = function(x = double(1)) {
          out <- derivs(run(x), wrt = 1:2, order = 1)
          ans <- nimNumeric(length = 6, value = out$jacobian)
          returnType(double(1))
          return(ans)
        },

        metaDerivsRun = function(x = double(1),
                                 order = double(1)) {
          out <- derivs(jacobianRun(x), wrt = 1:2,
                        order = order)
          returnType(ADNimbleList())
          return(out)
        }
      ), buildDerivs = c('jacobianRun', 'run')
    )

    ADfunInst <- ADfun1()
    xRec <- c(2.2, 3.3)
    x <- c(1.6, 2.8)
    Rderivs <- ADfunInst$jacobianRun(x)
    Rderivs <- ADfunInst$metaDerivsRun(x, order = 1)
    temporarilyAssignInGlobalEnv(ADfunInst)
    cADfunInst <- compileNimble(ADfunInst)

    ## record
    cADfunInst$derivsRun(x = xRec,
                         wrt = 1:2,
                         outInds = numeric(),
                         inDir = numeric(),
                         outDir = numeric(),
                         order = 1:2)

    ## Future args:
    wrt <- list(1:2, 1, 2, c(2, 1), c(2, 1, 2))
    outInds <- list(1:3, 1, 2, 3, c(3, 1), c(3, 1, 3))
    inDirs <- list(c(-0.3, -0.7))
    outDirs <- list(c(-2, -1, -3))
    orders <- list(0, 1, 2, c(0, 1), c(0, 2), c(1, 2), 0:2)

    correct <- cADfunInst$derivsRun(x = x,
                         wrt = numeric(),
                         outInds = numeric(),
                         inDir = numeric(),
                         outDir = numeric(),
                         order = 0:2)

    # Here is a comprehensive testing function

    check_over_orders <- function(x, orders, wrt=numeric(),
                                  outInds=numeric(), inDir=numeric(),
                                  outDir=numeric()) {
      n <- length(x)
      for(ord in orders) {
        this_x <- x + rnorm(length(x), 0, 0.1)
        check <- cADfunInst$derivsRun(x = this_x,
                                      wrt = wrt,
                                      outInds = outInds,
                                      inDir = inDir,
                                      outDir = outDir,
                                      order = ord)
        correct <- cADfunInst$derivsRun(x = this_x,
                                        wrt = numeric(),
                                        outInds = numeric(),
                                        inDir = numeric(),
                                        outDir = numeric(),
                                        order = 0:2)
        m <- nrow(correct$jacobian)
        if(0 %in% ord) {
          correct_value <- correct$value
          if(length(outInds))
            correct_value <- correct_value[outInds,drop=FALSE]
          expect_equal(check$value, correct_value)
        }
        if(1 %in% ord) {
          if(length(inDir) > 0) {
            correct_jac <- correct$jacobian %*% matrix(inDir, ncol=1)
            if(length(outInds) > 0)
              correct_jac <- correct_jac[outInds,, drop=FALSE]
          } else {
            if(length(outDir) > 0 && !(2 %in% ord))
              correct_jac <- matrix(outDir, nrow=1) %*% correct$jacobian
            else {
              correct_jac <- correct$jacobian
              if(length(outInds) > 0)
                correct_jac <- correct_jac[outInds,, drop=FALSE]
            }
            if(length(wrt) > 0)
              correct_jac <- correct_jac[,wrt,drop=FALSE]
          }
          expect_equal(check$jacobian,
                       correct_jac)
        }
        if(2 %in% ord) {
          if(length(inDir) && length(outDir)) {
            correct_hess <- inDir %*%
              apply(correct$hessian, 1, \(H) H %*% outDir) |>
              array(dim=c(1, n, 1))
            if(length(wrt))
              correct_hess <- correct_hess[,wrt,, drop=FALSE]
          } else if(length(inDir) && !length(outDir)) {
            correct_hess <- correct$hessian |>
              apply(1, \(x) inDir %*% x) |> t() |>  array(dim=c(1, n, m))
            if(length(wrt))
              correct_hess <- correct_hess[,wrt,, drop=FALSE]
            if(length(outInds))
              correct_hess <- correct_hess[,, outInds, drop=FALSE]
          } else if(!length(inDir) && length(outDir)) {
            correct_hess <- apply(correct$hessian, 1, \(H) H %*% outDir) |>
              array(dim=c(n,n,1))
            if(length(wrt))
              correct_hess <- correct_hess[wrt, wrt, , drop=FALSE]
          } else {
            correct_hess <- correct$hessian
            if(length(wrt))
              correct_hess <- correct_hess[wrt, wrt, , drop=FALSE]
            if(length(outInds))
              correct_hess <- correct_hess[,, outInds, drop=FALSE]
          }
          expect_equal(check$hessian,
                       correct_hess)
        }
      }
    }

    ## Next we use the above with all combinations of inputs
    ## This is rather verbose. If we included a case "numeric()"
    ## in each of the input lists above, then we could simply use
    ## one big quad-nested loop.
    ## However, in building up, it was useful to work in
    ## systematic combinations as done below.
    ## We could consolidate this in the future.
    ## But they run fast, as they use a single compiled object.

    ## wrt alone
    for(this_wrt in wrt) {
      check_over_orders(x=x,orders,wrt=this_wrt)
    }
    ## outInds alone
    for(this_outInds in outInds) {
      check_over_orders(x = x, orders, outInds = this_outInds)
    }
    ## wrt and outInds together -- the two kinds of indexing
    for(this_wrt in wrt) {
      for(this_outInds in outInds) {
        check_over_orders(x=x, orders=orders, wrt=this_wrt, outInds=this_outInds)
      }
    }
    ## inDir alone
    for(this_inDir in inDirs) {
      check_over_orders(x=x, orders=orders, inDir = this_inDir)
    }
    ## inDir and wrt -- the two kinds of "x" controls
    for(this_wrt in wrt) {
      for(this_inDir in inDirs) {
        check_over_orders(x=x, orders=orders, wrt = this_wrt, inDir = this_inDir)
      }
    }
    ## outDir alone
    for(this_outDir in outDirs) {
      check_over_orders(x=x, orders=orders, outDir = this_outDir)
    }
    ## outDir and outInds -- the two kinds of "y" controls
    for(this_outInds in outInds) {
      for(this_outDir in outDirs) {
        check_over_orders(x=x, orders=orders, outInds = this_outInds, outDir = this_outDir)
      }
    }
    ## inDir and outDir -- the two kinds of directions
    for(this_inDir in inDirs) {
      for(this_outDir in outDirs) {
        check_over_orders(x=x, orders=orders, outDir = this_outDir, inDir=this_inDir)
      }
    }

    ## wrt, inDir, and outDir
    for(this_wrt in wrt) {
      for(this_inDir in inDirs) {
        for(this_outDir in outDirs) {
          check_over_orders(x=x, orders=orders, wrt=this_wrt,
                            outDir = this_outDir, inDir=this_inDir)
        }
      }
    }

    ## wrt and outDir
    for(this_wrt in wrt) {
      for(this_outDir in outDirs) {
        check_over_orders(x=x, orders=orders, wrt=this_wrt, outDir = this_outDir)
      }
    }

    ## outInds, inDir, and outDir
    for(this_outInds in outInds) {
      for(this_inDir in inDirs) {
        for(this_outDir in outDirs) {
          check_over_orders(x=x, orders=orders,
                            outInds = this_outInds, outDir = this_outDir, inDir=this_inDir)
        }
      }
    }

    ## outInds and  inDir
    for(this_outInds in outInds) {
      for(this_inDir in inDirs) {
        check_over_orders(x=x, orders=orders,
                          outInds = this_outInds, inDir=this_inDir)
      }
    }

    ## all 4 controls
    for(this_wrt in wrt) {
      for(this_outInds in outInds) {
        for(this_inDir in inDirs) {
          for(this_outDir in outDirs) {
            check_over_orders(x=x, orders=orders, wrt=this_wrt,
                              outInds = this_outInds, outDir = this_outDir, inDir=this_inDir)
          }
        }
      }
    }


    ## ## What follows is some initial testing for manually inspecting cases.
    ## ## It may be useful if anything breaks to be able to run through these.
    ## foo <- function(t, x_, inDir) cADfunInst$run(x_ + t*inDir)
    ## ######################
    ## ## play with order 0 #
    ## ######################
    ## ## both wrt         ##
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 0)
    ## ## wrt flipped, all out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = c(2, 1),
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 0)

    ## ## single wrt, all out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1,
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 0)

    ## ## single wrt, all out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 2,
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 0)
    ## ## all wrt, single out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = 1,
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 0)

    ## ## all wrt, single out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = 2,
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 0)

    ## ## all wrt, single out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = 3,
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 0)


    ## ## all wrt, permuted out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = c(3, 1),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 0)

    ## ## single wrt, permuted out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 2,
    ##                      outInds = c(3, 1),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 0)

    ## # permuted wrt and permuted out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = c(2, 1),
    ##                      outInds = c(3, 1),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 0)

    ## ## inDir, all out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = numeric(),
    ##                      outInds = numeric(),
    ##                      inDir = c(-.3, -.7),
    ##                      outDir = numeric(),
    ##                      order = 0)

    ## ## outDir, all in
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = c(-2, -1, -3),
    ##                      order = 0)


    ## #####################
    ## ## play with order 1
    ## #####################
    ## ## both wrt
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1)
    ## ## wrt flipped, all out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = c(2, 1),
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1)

    ## ## single wrt, all out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1,
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1)

    ## ## single wrt, all out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 2,
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1)
    ## ## all wrt, single out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = 1,
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1)

    ## ## all wrt, single out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = 2,
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1)

    ## ## all wrt, single out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = 3,
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1)


    ## ## all wrt, permuted out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = c(3, 1),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1)

    ## ## single wrt, permuted out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 2,
    ##                      outInds = c(3, 1),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1)

    ## # permuted wrt and permuted out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = c(2, 1),
    ##                      outInds = c(3, 1),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1)

    ## ## inDir, all out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = numeric(),
    ##                      outInds = numeric(),
    ##                      inDir = c(-.3, -.7),
    ##                      outDir = numeric(),
    ##                      order = 1)

    ## ## outDir, all in
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = c(-2, -1, -3),
    ##                      order = 1)


    ## #####################
    ## ## play with order 2
    ## #####################
    ## ## both wrt
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 2)
    ## ## wrt flipped, all out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = c(2, 1),
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 2)

    ## ## single wrt, all out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1,
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 2)

    ## ## single wrt, all out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 2,
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 2)

    ## ## all wrt, single out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = 1,
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 2)

    ## ## all wrt, single out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = 2,
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 2)

    ## ## all wrt, single out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = 3,
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 2)


    ## ## all wrt, permuted out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = c(3, 1),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 2)

    ## ## single wrt, permuted out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 2,
    ##                      outInds = c(3, 1),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 2)

    ## # permuted wrt and permuted out
    ## cADfunInst$Derivsrun(x = x,
    ##                      wrt = c(2, 1),
    ##                      outInds = c(3, 1),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 2)

    ## ## inDir, all out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = numeric(),
    ##                      outInds = numeric(),
    ##                      inDir = c(-.3, -.7),
    ##                      outDir = numeric(),
    ##                      order = 2)

    ## ## outDir, all in
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = c(-2, -1, -3),
    ##                      order = 2)

    ## #    apply(correct$hessian, 1, \(H) H %*% c(-2, -1, -3))

    ## ## inDir and outDir
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = numeric(),
    ##                      outInds = numeric(),
    ##                      inDir = c(-.3, -.7),
    ##                      outDir = c(-2, -1, -3),
    ##                      order = 2)

    ## #c(-.3, -.7) %*% apply(correct$hessian, 1, \(H) H %*% c(-2, -1, -3))


    ## #####################
    ## ## play with orders 1 & 2
    ## #####################
    ## ## both wrt
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1:2)
    ## ## wrt flipped, all out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = c(2, 1),
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1:2)

    ## ## single wrt, all out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1,
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1:2)

    ## ## single wrt, all out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 2,
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1:2)

    ## ## all wrt, single out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = 1,
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1:2)

    ## ## all wrt, single out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = 2,
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1:2)

    ## ## all wrt, single out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = 3,
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1:2)


    ## ## all wrt, permuted out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = c(3, 1),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1:2)

    ## ## single wrt, permuted out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 2,
    ##                      outInds = c(3, 1),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1:2)

    ## # permuted wrt and permuted out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = c(2, 1),
    ##                      outInds = c(3, 1),
    ##                      inDir = numeric(),
    ##                      outDir = numeric(),
    ##                      order = 1:2)

    ## ## inDir, all out
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = numeric(),
    ##                      outInds = numeric(),
    ##                      inDir = c(-.3, -.7),
    ##                      outDir = numeric(),
    ##                      order = 1:2)

    ## ## outDir, all in
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = 1:2,
    ##                      outInds = numeric(),
    ##                      inDir = numeric(),
    ##                      outDir = c(-2, -1, -3),
    ##                      order = 1:2)

    ## #    apply(correct$hessian, 1, \(H) H %*% c(-2, -1, -3))

    ## ## inDir and outDir
    ## cADfunInst$derivsRun(x = x,
    ##                      wrt = numeric(),
    ##                      outInds = numeric(),
    ##                      inDir = c(-.3, -.7),
    ##                      outDir = c(-2, -1, -3),
    ##                      order = 1:2)

    ## #c(-.3, -.7) %*% apply(correct$hessian, 1, \(H) H %*% c(-2, -1, -3))
  }
)

## test_that('Derivatives with double-taping, indices, and directions work',
##   {
##     ADfun1 <- nimbleFunction(
##       setup = function(){},
##       run = function(x = double(1)) {
##         ans <- c(x[1]^3, x[2]^3, (x[1]^2*x[2]^2))
##         return(ans)
##         returnType(double(1))
##       },
##       methods = list(
##         derivsRun = function(x = double(1),
##                              wrt = double(1),
##                              outInds = double(1),
##                              inDir = double(1),
##                              outDir = double(1),
##                              order = double(1)) {
##           ans <- derivs(run(x), wrt = wrt, outInds = outInds,
##                         inDir = inDir, outDir = outDir, order = order)
##           return(ans)
##           returnType(ADNimbleList())
##         },
##         jacobianRun = function(x = double(1)) {
##           out <- derivs(run(x), wrt = 1:2, order = 1)
##           ans <- nimNumeric(length = 6, value = out$jacobian)
##           returnType(double(1))
##           return(ans)
##         },

##         metaDerivsRun = function(x = double(1),
##                                  order = double(1)) {
##           out <- derivs(jacobianRun(x), wrt = 1:2,
##                         order = order)
##           returnType(ADNimbleList())
##           return(out)
##         }
##       ), buildDerivs = c('jacobianRun', 'run')
##     )

##     ADfunInst <- ADfun1()
##     xRec <- c(2.2, 3.3)
##     x <- c(1.6, 2.8)
##     Rderivs <- ADfunInst$jacobianRun(x)
##     Rderivs <- ADfunInst$metaDerivsRun(x, order = 1)
##     temporarilyAssignInGlobalEnv(ADfunInst)
##     cADfunInst <- compileNimble(ADfunInst)

##     ## record
##     cADfunInst$derivsRun(x = xRec,
##                          wrt = 1:2,
##                          outInds = numeric(),
##                          inDir = numeric(),
##                          outDir = numeric(),
##                          order = 1:2)

##     ## Future args:
##     wrt <- list(1:2, 1, 2, c(2, 1), c(2, 1, 2))
##     outInds <- list(1:3, 1, 2, 3, c(3, 1), c(3, 1, 3))
##     inDirs <- list(c(-0.3, -0.7))
##     outDirs <- list(c(-2, -1, -3))
##     orders <- list(0, 1, 2, c(0, 1), c(0, 2), c(1, 2), 0:2)

##     cADfunInst$derivsRun(x = xRec,
##                          wrt = 1:2,
##                          outInds = numeric(),
##                          inDir = rep(inDirs[[1]], 4),
##                          outDir = numeric(),
##                          order = 1:2)

##   })

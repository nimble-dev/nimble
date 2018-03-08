
######################################################################
######################################################################


## Testing bnp distributions

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
context('Testing NIMBLE distributions')

## dCRP nimble function

x <- c(1,1,3,1,1,3)
conc <- 1

truth <- (conc/(conc+1-1))*(1/(conc+2-1))*(conc/(conc+3-1))*
  (1/(conc+4-1))*(1/(conc+5-1))*(1/(conc+6-1))

try(test_that("dCRP nimble function calculates density correctly: ",
              expect_equal(dCRP(x, conc, log=FALSE),
                           truth,
                           info = paste0("incorrect dCRP nimble function calculation"))))

## test use through log scale
ltruth <- log((conc/(conc+1-1))*(1/(conc+2-1))*(conc/(conc+3-1))*
  (1/(conc+4-1))*(1/(conc+5-1))*(1/(conc+6-1)))

try(test_that("dCRP nimble function calculates density correctly in log scale: ",
              expect_equal(dCRP(x, conc, log=TRUE),
                           ltruth,
                           info = paste0("incorrect dCRP nimble function calculation in log scale"))))


## test use through compile nimble function
cdCRP <- compileNimble(dCRP)

try(test_that("Test that dCRP works correctly in compiled nimble function: ",
              expect_equal(cdCRP(x, conc), (truth), 
                           info = paste0("incorrect dCRP value in compiled nimble function"))))


## test use through compile nimble function in log scale

try(test_that("Test that dCRP works correctly in compiled nimble function in log scale: ",
              expect_equal(cdCRP(x, conc, log=TRUE), (ltruth), 
                           info = paste0("incorrect dCRP value in compiled nimble function  in log scale"))))

## test use in model
CRP_code <- nimbleCode({
  x[1:6] ~ dCRP(conc)
})

Consts <- list(conc = 1)
Inits <- list(x = c(1,1,3,1,1,3))
CRP_model <- nimbleModel(CRP_code, data=Inits, constants=Consts)

CRP_model$x <- x

try(test_that("Test that dCRP calculation is correct in model likelihood calculation: ",
              expect_equal(exp(CRP_model$calculate()), (truth),
                           info = paste0("incorrect likelihood value for dCRP"))))

c_CRP_model <- compileNimble(CRP_model)
c_CRP_model$x
try(test_that("Test that dCRP (compiled) calculation is correct in model likelihood calculation: ",
              expect_equal(exp(c_CRP_model$calculate()), (truth),
                           info = paste0("incorrect likelihood value for dCRP (compiled)"))))


## random sampling
set.seed(0)
n <- 6
r_samps <- t(replicate(10000, rCRP(n = n, conc)))
# K is the number of unique components in x of length 6
true_EK <- sum(conc/(conc+1:n-1))


try(test_that("Test that random samples (rCRP) have correct mean of K: ",
              expect_equal(mean(apply(r_samps, 1, function(x)length(unique(x)))), true_EK, 
                           tol = 0.01,
                           info = "Difference in expected mean of K exceeds tolerance")))

## sampling via `simulate`: does not work!
#set.seed(0)
#simul_samp <- function(model) {
#  model$simulate()
#  return(model$x)
#}

#simul_samps <- t(replicate(10000, simul_samp(c_CRP_model)))

#c_CRP_model$simulate()

#try(test_that("Test that random samples (simulate) have correct mean of K: ",
#              expect_equal(mean(apply(r_samps, 1, function(x)length(unique(x)))), true_EK, 
#                           tol = 0.01,
#                           info = "Difference in expected mean of K exceeds tolerance")))


# check parameters: conc>0; in sampling check returned values




## Test stick_breaking nimble function:



source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of compareMCMCs")
iCase <- 1

code1 <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dnorm(a, 1)
    for(i in 1:3) y[i] ~ dnorm(b, 1)
})

input1 <- list(code = code1, data = list(y = 1:3), inits = list(a = 0.5))

test <- try(results1 <- compareMCMCs(input1, MCMCs = c('nimble')))
expect_true(!inherits(test, 'try-error'), paste('case', iCase))
iCase <- iCase + 1

## This will generate a warning -- sensibly -- which is ok.
test <- try(make_MCMC_comparison_pages(results1, 'model1'))
expect_true(!inherits(test, 'try-error'), paste('case', iCase))
iCase <- iCase + 1

## two MCMC methods
test <- try(results2 <- compareMCMCs(input1, MCMCs = c('nimble','noConj')))
expect_true(!inherits(test, 'try-error'), paste('case', iCase))
iCase <- iCase + 1

test <- try(make_MCMC_comparison_pages(results2, 'model1b'))
expect_true(!inherits(test, 'try-error'), paste('case', iCase))
iCase <- iCase + 1

## rename a result
test <- try(results2[[1]] <- rename_MCMC_comparison_method('nimble', 'another nimble', results2[[1]]))
expect_true(!inherits(test, 'try-error'), paste('case', iCase))
iCase <- iCase + 1

test <- try(results3 <- combine_MCMC_comparison_results(results1[[1]], results2[[1]], name = 'combined results'))
expect_true(!inherits(test, 'try-error'), paste('case', iCase))
iCase <- iCase + 1

test <- try(make_MCMC_comparison_pages(results3, 'model1c'))
expect_true(!inherits(test, 'try-error'), paste('case', iCase))
iCase <- iCase + 1

#### Now repeat in case of multiple parameters
code1 <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dnorm(a, 1)
    sigma ~ dunif(0.5, 1)
    for(i in 1:3) y[i] ~ dnorm(b, sd = sigma)
})

input1 <- list(code = code1, data = list(y = 1:3), inits = list(a = 0.5))
test <- try(results1 <- compareMCMCs(input1, MCMCs = c('nimble')))
## This will generate a warning -- sensibly -- which is ok.
expect_true(!inherits(test, 'try-error'), paste('case', iCase))
iCase <- iCase + 1

test <- try(make_MCMC_comparison_pages(results1, 'model1'))
expect_true(!inherits(test, 'try-error'), paste('case', iCase))
iCase <- iCase + 1

## two MCMC methods
test <- try(results2 <- compareMCMCs(input1, MCMCs = c('nimble','noConj')))
expect_true(!inherits(test, 'try-error'), paste('case', iCase))
iCase <- iCase + 1

test <- try(make_MCMC_comparison_pages(results2, 'model1b'))
expect_true(!inherits(test, 'try-error'), paste('case', iCase))
iCase <- iCase + 1

## rename a result
test <- try(results2[[1]] <- rename_MCMC_comparison_method('nimble', 'another nimble', results2[[1]]))
expect_true(!inherits(test, 'try-error'), paste('case', iCase))
iCase <- iCase + 1

test <- try(results3 <- combine_MCMC_comparison_results(results1[[1]], results2[[1]], name = 'combined results'))
expect_true(!inherits(test, 'try-error'), paste('case', iCase))
iCase <- iCase + 1

test <- try(make_MCMC_comparison_pages(results3, 'model1c'))
expect_true(!inherits(test, 'try-error'), paste('case', iCase))
iCase <- iCase + 1

test <- try({
    mfiles <- list.files("model1", full.names = TRUE)
    file.remove(mfiles)
    file.remove('model1')
    mfiles <- list.files("model1b", full.names = TRUE)
    file.remove(mfiles)
    file.remove('model1b')
    mfiles <- list.files("model1c", full.names = TRUE)
    file.remove(mfiles)
    file.remove('model1c')
})

expect_true(!inherits(test, 'try-error'), paste('clearing files: case', iCase))

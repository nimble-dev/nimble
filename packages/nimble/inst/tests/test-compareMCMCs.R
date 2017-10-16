source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)

context("Testing of compareMCMCs")

code1 <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dnorm(a, 1)
    for(i in 1:3) y[i] ~ dnorm(b, 1)
})

input1 <- list(code = code1, data = list(y = 1:3), inits = list(a = 0.5))

test_that("Basic MCMC comparison", {
    expect_output(results1 <- compareMCMCs(input1, MCMCs = c('nimble')))
    expect_message(make_MCMC_comparison_pages(results1, 'model1'), "geom_path")
    expect_output(results2 <- compareMCMCs(input1, MCMCs = c('nimble','noConj')))
    expect_silent(make_MCMC_comparison_pages(results2, 'model1b'))
    expect_silent(results2[[1]] <- rename_MCMC_comparison_method('nimble', 'another nimble', results2[[1]]))
    expect_silent(results3 <- combine_MCMC_comparison_results(results1[[1]], results2[[1]], name = 'combined results'))
    expect_silent(make_MCMC_comparison_pages(results3, 'model1c'))
})



#### Now repeat in case of multiple parameters
code1 <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dnorm(a, 1)
    sigma ~ dunif(0.5, 1)
    for(i in 1:3) y[i] ~ dnorm(b, sd = sigma)
})

input1 <- list(code = code1, data = list(y = 1:3), inits = list(a = 0.5))

test_that("Basic MCMC comparison, multiple parameters", {
    expect_output(results1 <- compareMCMCs(input1, MCMCs = c('nimble')))
    expect_message(make_MCMC_comparison_pages(results1, 'model1'), "geom_path")
    expect_output(results2 <- compareMCMCs(input1, MCMCs = c('nimble','noConj')))
    expect_silent(make_MCMC_comparison_pages(results2, 'model1b'))
    expect_silent(results2[[1]] <- rename_MCMC_comparison_method('nimble', 'another nimble', results2[[1]]))
    expect_silent(results3 <- combine_MCMC_comparison_results(results1[[1]], results2[[1]], name = 'combined results'))
    expect_silent(make_MCMC_comparison_pages(results3, 'model1c'))
    mfiles <- list.files("model1", full.names = TRUE)
    expect_identical(file.remove(mfiles), rep(TRUE, 5))
    expect_identical(file.remove('model1'), TRUE)
    mfiles <- list.files("model1b", full.names = TRUE)
    expect_identical(file.remove(mfiles), rep(TRUE, 5))
    expect_identical(file.remove('model1b'), TRUE)
    mfiles <- list.files("model1c", full.names = TRUE)
    expect_identical(file.remove(mfiles), rep(TRUE, 5))
    expect_identical(file.remove('model1c'), TRUE)
})



options(warn = RwarnLevel)

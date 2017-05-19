#!/bin/sh -ex

# This runner profiles the tests so that we can prioritize for faster failures.
RUNNER="/usr/bin/time -v Rscript test-files.R"

$RUNNER packages/nimble/inst/tests/test-types.R  # takes 0.94 seconds
$RUNNER packages/nimble/inst/tests/test-checkDSL.R  # takes 1.78 seconds
$RUNNER packages/nimble/inst/tests/test-misc.R  # takes 1.89 seconds
$RUNNER packages/nimble/inst/tests/test-setData.R  # takes 2.2 seconds
$RUNNER packages/nimble/inst/tests/test-getDependencies.R  # takes 3.87 seconds
$RUNNER packages/nimble/inst/tests/test-errors.R  # takes 4.03 seconds
$RUNNER packages/nimble/inst/tests/test-modelValuesInterface.R  # takes 5.33 seconds
$RUNNER packages/nimble/inst/tests/test-optim.R  # takes 9.52 seconds
$RUNNER packages/nimble/inst/tests/test-nimbleFunctionInterfaces.R  # takes 11.08 seconds
$RUNNER packages/nimble/inst/tests/test-dsl_dists.R  # takes 13.25 seconds
$RUNNER packages/nimble/inst/tests/test-getDistributionInfo.R  # takes 14.15 seconds
$RUNNER packages/nimble/inst/tests/test-size.R  # takes 35.44 seconds
$RUNNER packages/nimble/inst/tests/test-declare.R  # takes 44.82 seconds
$RUNNER packages/nimble/inst/tests/test-nimbleList_RCfun.R  # takes 53.15 seconds
$RUNNER packages/nimble/inst/tests/test-mcem.R  # takes 65.57 seconds
$RUNNER packages/nimble/inst/tests/test-getBound.R  # takes 88.2 seconds
$RUNNER packages/nimble/inst/tests/test-compareMCMCs.R  # takes 101.4 seconds
$RUNNER packages/nimble/inst/tests/test-distributions.R  # takes 119 seconds
$RUNNER packages/nimble/inst/tests/test-nimbleList.R  # takes 128.91 seconds
$RUNNER packages/nimble/inst/tests/test-user.R  # takes 134.32 seconds
$RUNNER packages/nimble/inst/tests/test-allocation.R  # takes 200.52 seconds
$RUNNER packages/nimble/inst/tests/test-math.R  # takes 217.93 seconds
$RUNNER packages/nimble/inst/tests/test-coreR.R  # takes 273.36 seconds
$RUNNER packages/nimble/inst/tests/test-getParam.R  # takes 278.14 seconds
$RUNNER packages/nimble/inst/tests/test-copy.R  # takes 521.16 seconds
$RUNNER packages/nimble/inst/tests/test-filtering.R  # takes 565.14 seconds
$RUNNER packages/nimble/inst/tests/test-trunc.R  # takes 579.12 seconds
$RUNNER packages/nimble/inst/tests/test-models.R  # takes 594.9 seconds
$RUNNER packages/nimble/inst/tests/test-numericTypes.R  # takes 1220.14 seconds
$RUNNER packages/nimble/inst/tests/test-mcmc.R  # takes 1684.73 seconds

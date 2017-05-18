#!/bin/sh -ex

# This runner profiles the tests so that we can prioritize for faster failures.
RUNNER="/usr/bin/time -v Rscript"

$RUNNER packages/nimble/inst/tests/test-allocation.R
$RUNNER packages/nimble/inst/tests/test-checkDSL.R
$RUNNER packages/nimble/inst/tests/test-compareMCMCs.R
$RUNNER packages/nimble/inst/tests/test-copy.R
$RUNNER packages/nimble/inst/tests/test-coreR.R
$RUNNER packages/nimble/inst/tests/test-declare.R
$RUNNER packages/nimble/inst/tests/test-distributions.R
$RUNNER packages/nimble/inst/tests/test-dsl_dists.R
$RUNNER packages/nimble/inst/tests/test-errors.R
$RUNNER packages/nimble/inst/tests/test-filtering.R
$RUNNER packages/nimble/inst/tests/test-getBound.R
$RUNNER packages/nimble/inst/tests/test-getDependencies.R
$RUNNER packages/nimble/inst/tests/test-getDistributionInfo.R
$RUNNER packages/nimble/inst/tests/test-getParam.R
$RUNNER packages/nimble/inst/tests/test-math.R
$RUNNER packages/nimble/inst/tests/test-mcem.R
$RUNNER packages/nimble/inst/tests/test-mcmc.R
$RUNNER packages/nimble/inst/tests/test-misc.R
$RUNNER packages/nimble/inst/tests/test-modelValuesInterface.R
$RUNNER packages/nimble/inst/tests/test-models.R
$RUNNER packages/nimble/inst/tests/test-nimbleFunctionInterfaces.R
$RUNNER packages/nimble/inst/tests/test-nimbleList.R
$RUNNER packages/nimble/inst/tests/test-nimbleList_RCfun.R
$RUNNER packages/nimble/inst/tests/test-numericTypes.R
$RUNNER packages/nimble/inst/tests/test-optim.R
$RUNNER packages/nimble/inst/tests/test-setData.R
$RUNNER packages/nimble/inst/tests/test-size.R
$RUNNER packages/nimble/inst/tests/test-trunc.R
$RUNNER packages/nimble/inst/tests/test-types.R
$RUNNER packages/nimble/inst/tests/test-user.R

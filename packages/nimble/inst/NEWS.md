#                            CHANGES IN VERSION 1.1.0 (January 2024) 

## USER LEVEL CHANGES

- Enhance use of AD in models:

    -- Allow use of stochastic indexing in models with AD (PR #1389).

    -- Allow use of AD with CAR models, which enables use of HMC on nodes specified
    to have `dcar_normal` or `dcar_proper` distributions (PR #1390).
    
    -- Allow distributions and functions (whether user-defined or built-in) that 
    lack AD support to be used and compiled in AD-enabled models, including
    truncated distributions, `dcat`, `dinterval`, and `dconstraint` (PR #1389).
    
- Add `nimIntegrate`, providing functionality analogous to R's `integrate`,
  to the NIMBLE language, providing one-dimensional numerical integration via
  adaptive quadrature (PR #1357).

- Add `prior_samples` MCMC sampler, which uses an existing set of
  numerical samples to define the prior distribution of model node(s).

- `configureMCMC` will no longer assign samplers to data nodes, even if
  the `nodes` argument includes data nodes (PR #1407).
  
- Add new argument `allowData` to the `addSampler` method of MCMC
  configuration objects, with default value `FALSE`.  When `TRUE`,
  samplers can be assigned to operate on data nodes (PR #1407).
  
- Handle predictive node dependencies of `dCRP` cluster indicators
  by erroring out, but letting users know that setting 
  `MCMCusePredictiveDependenciesInCalculations` to `TRUE` will make MCMC 
  possible (PR #1402).
    
- Rename `expandTarget` and `returnScalarComponents` arguments of `addSampler`
  to be `targetByNode` and `multivariateNodesAsScalars`, respectively, to 
  improve clarity. Also clarify help info for `addSampler` (PR #1375).
  
- Stop instead of warning when encountering run-time size errors, to avoid
  the possibility of incorrect results (PR #1401).
  
- Add a warning to MCMC sampling of `CAR_normal` nodes when `zero_mean=1`
  and the on-the-fly centering causes an invalid model state (PR #1400).
  
- Cleanly error out when a variable name in a model would conflict with 
  C++ keywords (PR #1382).
  
- Warn if indexing in model code uses variables from the user environment
  rather than from `constants` (PR #1380).
  
- Error out rather than warning when there is a missing variable in a 
  nimbleFunction (PR #1379).
  
- Allow users to provide list of lists for `inits` in `nimbleModel` (PR #1376).
 
- Add new control list option, `maxDimCovHistory` to `RW_block` sampler to 
  specify maximum dimension for saving proposal covariance history.
  
- Change argument of `besselK` in manual table to be `x` not `k`.

- Do not allow elements of a `nimbleList` to be named `name`, `predefined`,
  or `where` (PR #1384).

- Better notify users when automatically generated `r` function for a user-
  defined distribution has been removed (PR #1374).

- Error out if there is a default first (`x`) argument in a user-defined
  distribution (PR #1371).
  
- Improve error trapping related to arguments of `getParam` (PR #1370). 

- Add information in `nimbleSMC` and `nimbleHMC` samplers to roxygen and manual.

- Improve documentation of omitting variables from AD derivative tracking.

- Update installation guidance in manual.

## DEVELOPER LEVEL CHANGES

- Update to Eigen 3.4.0 but comment out various pragmas in
  `DisableStupidWarnings.h` preventing CRAN checks from passing (PR #1406).

- Improve efficiency of `mcmc_determineCalcAndCopyNodes` by avoiding repeated
  calls (PR #1333).

- Inline various C++ functions to improve efficiency (PR #1349).

- Add some functionality to support model macros (PR #1361).

- Make change to `nimble-package` documentation to use `"_PACKAGE"`
  instead of `@docType` per CRAN request (issue #1359).
  
- Clean up C++ warnings related to unused variables and typedefs (PR #1408).

- Use `--pre-clean` when invoking `R CMD SHLIB` to do C++ compilation to avoid
  sporadic test failures in WAIC tests, seemingly caused by strange timestamp-
  related behavior of `make` (PR #1393).
  
- Fix case of `||` used instead of `|` causing compilation error on M2 Mac 
  (PR #1392). 
  
- Fix some C++ warnings flagged by CRAN (PR #1386).

- Systematically use `getNimbleOption` in NIMBLE code base.

- Update versions of GitHub Actions canned actions.

## BUG FIXES

- Fix `is.na` and `is.nan` in the NIMBLE DSL to behave in a vectorized fashion
  and ensure they mimic behavior in R apart from logical or integer inputs to
  `is.na` (PR #1394).

- Remove the `RW_multinomial` MCMC sampler, which was found to generate incorrect
  posterior results.  A corrected version of this sampler may be
  re-introduced into the package, depending on user interest.
  
- Fix a bug in conjugacy checking in a case of subsets of multivariate nodes
  (PR #1331).
  
- Fix a recycling rule bug involving `pow_int` with a matrix input (PR #1396).

- Fix name-mangling problem that prevented use of periods in names of packages
  using nimble (PR #1383).
  
- Fix handling of `NULL` in `getNimbleOption`. 
  
- Avoid using wrapped sampler on cluster node parameters when using `dCRP`
  and assigning a joint sampler (e.g., HMC) to the parameters of the cluster 
  nodes (PR #1404).
  
- Avoid spurious warning about missing nimbleFunction when using `nimOptim`
  (PR #1378).
  
- Cleanly error out when an undefined function is used in the `return` statement
  of a nimbleFunction (in particular, nested nimbleFunctions) (PR #1381).
  
- Fix handling of deregistered distribution by removing auto-generated 'r'
  function (PR #1377).
  
- Correct a corner case of reading from a BUGS file (PR #1369).
  
- Fix incorrect eigenization code related to `besselk` that did not actually 
  seem to affect behavior (PR #1385).
  
#                            CHANGES IN VERSION 1.0.1 (June 2023) 

## USER LEVEL CHANGES

- Fix bug (introduced in v. 1.0.0) causing incorrect setting of 'data'
flag for models with variables containing a mix of data nodes and nodes
appearing only on the right-hand side of expressions, for cases where not
all elements of the variable are defined (such as capture-recapture
models). This also addresses an unanticipated change in behavior in 
initializing right-hand side only nodes from the `data` argument (PR #1328).

- Improved error trapping in MCMC sampler for categorical distributions 
(PR #1325).


#                            CHANGES IN VERSION 1.0.0 (May 2023) 

This release introduces tools for automatic differentiation in NIMBLE.

## USER LEVEL CHANGES

- Add tools for automatic differentiation. Functionality includes:

    -- The ability to differentiate nimbleFunctions with respect
    to input arguments in a flexible way.
    
    -- The ability to differentiate `model$calculate()` calls
    with respect to model node elements.
    
    -- Functionality that enables Hamiltonian Monte Carlo (provided in the
    `nimbleHMC` package).
    
    -- A parameter transformation system to work in unconstrained parameter
    spaces when model parameters have constrained domains.
        
- Add Laplace approximation algorithm, allowing a user to approximate
the marginal likelihood (integration/marginalizing over random effects/
latent process values) and find the MLE.

- Improve aspects of WAIC (PR #1256):

    -- Provide aggregate and per-chain WAIC when `perChainWAIC` is `TRUE`.

    -- Report number of `pWAIC` values greater than 0.4.
    
    -- Warn user if WAIC is enabled but not set to `TRUE` in `runMCMC`.
    
    -- Improve clarity in `help(WAIC)`.
    
- Reduce the default adaptation interval for the `AF_slice` sampler to 200.

- Improved control of `addSampler` method of MCMC
configuration objects (PR #1293):

    -- `expandTarget` argument controls whether nodes specified in
       `target` undergo expansion via `expandNodeNames`, adding a
       separate sampler instance for each resulting node.

    -- When `expandTarget = TRUE`, the `scalarComponents` argument is
       passed as the `returnScalarComponents` argument to `expandNodeNames`.

    -- Removed `nodes` argument as redundant.  Target nodes are
       uniquely specified using the `target` argument.

- Improve handling of resetting of scale when `adaptive=FALSE` (PR #1275).

- Allow `dcat` to have a `prob` vector of length 1 (PR #1251).

- Add error trapping for MCMC `thin` argument not positive and integer-valued
(PR #1250).

- Enhance discussion of MCMC initialization in the User Manual (issue #1247).

- Enhance discussion of indexing in loops and variables in User Manual,
including differences from JAGS.

- Rework language in `help(buildMCMC)` to improve clarity.

- Add information on predefined nimbleLists to the User Manual.

- Improve error message when user-defined distribution has incorrect dimension
for `x`.

- Improve error message when mistakenly using `T()` in deterministic rather than
  stochastic model declaration.

- Improve error-trapping when a loop index is incorrectly used in indexing
a block of a variable (PR #1289).

- Error trap use of reduction operations on scalars (PR #1281).

- Error trap some cases of incorrect return type in user-defined simulation
('r') functions (PR #1280).

- Error trap case where a variable name is also used as a loop index variable
(PR #1278).

- Error trap setting deterministic nodes as data in `setData` (PR #1269).

- Better error trap of wrong size/dimensions of in `setInits` (PR #1260).

- Export `clearCompiled`, allowing users to clear compiled objects and unload
the shared library produced during nimble's compilation process (on non-Windows
operating systems).

## BUG FIXES

- Fix bug producing integer overflow in `getDependencyPathCountOneNode`
and improve efficiency of checking maximum number of paths in conjugacy
checking (PR #1322).

- Fix long-standing memory leak in `dinvwish_chol` (PR #1320).

- Fix bug giving incorrect results for `runCrossValidate`, when using the 
default MSE loss function and any user-defined loss functions, 
but not the 'predictive' loss function. Also, calculate loss when using
the predictive loss using the compiled model rather than the uncompiled model.
(PRs #1298, #1299).

- Fix bug in conjugacy checking that caused incorrect identification of
conjugate relationships with subsets and supersets of multivariate nodes,
as well as inconsistent slices of multivariate nodes (PR #1290).

- Fix `runMCMC` to reset the WAIC state when rerun (PR #1256).

- Fix reset mechanism for WAIC so that one can get WAIC from all the samples
in a chain that is created by multiple calls to the MCMC `$run` method even
when `getWAIC` is called in the middle (PR #1310).

- Fix bug affecting `d=2` case in uncompiled `RW_lkj_corr_cholesky` sampler
(PR #1273).

- Fix aliasing issue arising in assignment operations such as 
`nf$a <- nf$a[1:3]` (PR 1301).

- Fix assignment into a block of a model variable (PR #1300).

- Fix compilation error when using `nimSeq` with integers (PR #1282).

- Fix scoping of argument evaluation for `replaceSamplers` method (PR #1287).

- Fix a corner case where particular indexing in model code prevents model
building (PR #1279).

- Fix handling of right-hand side only nodes so not flagged as deterministic
(issue #1269).

- Error trap and fix a bug in handling particular cases in list-like subsetting
(PR #1259).

- Make optimDefaultControlList handle default values completely and with
  consistent compiled and uncompiled behavior.
  
- make nimOptim respect parscale control parameter.

- Update manual installation chapter and `INSTALL` links to gfortran on MacOS.

## DEVELOPER LEVEL CHANGES

- Take `nimbleCppADbaseClass` out of `libnimble.a` and instead make a
session-specific .o file from it to address AD-related crashes on some OSes. (PR
#1318)

- Remove `nimOptim_model` functionality. This fixes build warning on MacOS
  involving `nimOptim` (PR #1276).

- Update `Eigen` copyright dates.

- Fix invisible return of `NULL` in nimbleFunctions without explicit return
statements (PR #1254). 

- Fix inconsistent use of offset in internal conversion functions (PR #1277).

- Remove references to unsupported `trace` (issue #1262).


#                            CHANGES IN VERSION 0.13.2 (May 2023) 

## DEVELOPER LEVEL CHANGES

- Remove code triggering Windows warnings about writing bytes into region of size 0 by modifying `setLength`.

- Remove use of C++11 standard, per CRAN requirements. As part of this, replace use of `std::ptr_fun` (PR #1292).


#                            CHANGES IN VERSION 0.13.1 (December 2022) 

## BUG FIXES

- Fix bug in MCMC sampler inclusion/exclusion of predictive nodes from
  target node dependencies (PR #1248).

#                            CHANGES IN VERSION 0.13.0 (November 2022)

## USER LEVEL CHANGES

- Thoroughly revamp handling of predictive nodes in MCMC sampling. If MCMC 
results identical to previous versions of NIMBLE are needed in models with 
posterior predictive nodes, set 
`nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)` 
and `nimbleOptions(MCMCorderPosteriorPredictiveSamplersLast = FALSE)`.

    -- MCMC samplers, by default, will now exclude predictive dependencies 
    from internal sampler calculations.  This can be reverted to the old behavior
    of including predictive dependencies in calculations using 
    `nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)`.
    
    -- At the time of `buildMCMC`, all `posterior_predictive` samplers are
    automatically reordered to operate last among all samplers. Doing so, 
    posterior predictive samples are generated conditional on the other values
    in the MCMC sample.  This reordering can be disabled using 
    `nimbleOptions(MCMCorderPosteriorPredictiveSamplersLast = FALSE)` 
    (but doing so without also setting 
    `nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)` could
    result in samples that are invalid in terms of the joint posterior 
    distribution (but with valid samples marginally).
    
    -- Removal of the `posterior_predictive_branch` sampler.  Filling the same
    role, the `posterior_predictive` sampler now updates all nodes downstream 
    of its `target` node.  Assignment of the `posterior_predictive` sampler 
    happens automatically during MCMC configuration, unless 
    `nimbleOptions(MCMCusePosteriorPredictiveSampler = FALSE)`.
    
    -- Automatic determination of "predictive" model nodes, which are all
    stochastic non-data nodes that have no data nodes anywhere in their 
    downstream dependencies. Tracking of predictive nodes is done 
    automatically, but maybe be disabled using 
    `nimbleOptions(determinePredictiveNodesInModel = FALSE)`.

    -- New arguments `includePredictive` (default value `TRUE`) and 
    `predictiveOnly` (default value `FALSE`), for both the `getNodeNames`
    and the `getDependencies` methods of model objects.  These specify whether
    any predictive nodes are included in the results, and whether only 
    predictive nodes are included, respectively.
    
    -- The MCMC configuration object will issue a warning message if there are
    stochastic non-data nodes which will not undergo MCMC sampling.  This 
    warning can be disabled using 
    `nimbleOptions(MCMCwarnUnsampledStochasticNodes = FALSE)`.

- Add option to WAIC system (via `controlWAIC`) to allow additional burnin (in 
addition to standard MCMC burnin) before calculating online WAIC, thereby 
allowing inspection of initial samples without forcing them to be used for WAIC
(PR #1244).

- For MCMC configuration `addSampler` method, change name of the 
`scalarComponents` argument to `expandComponents` (PR #1215).

- Add new `default` argument for the `addSampler` method of MCMC configuration
objects.  When `default = TRUE`, default samplers (conjugate, or otherwise) will
be added to the specified nodes.  The addition of this argument provides an 
entry point to the logic of default sampler determination and assignment, 
without creating a new MCMC configuration object (PR #1215).

- Add new `nodes` argument for the `addSampler` method of MCMC configuration
objects.  Nodes specified in `nodes` automatically undergo expansion according
to `expandNodeNames` prior to sampler assignment, allowing for easier assignment
of samplers to multiple nodes (PR #1215).

- `rcar_normal` issues an informative error message when invoked from the R 
command line (PR #1243).

- Warn users of unused constants during model building (PR #1242).

- Add `replaceSamplers` method to MCMC configuration objects to simplify 
modifying how a node is sampled (PR #1222).

- Convert `NEWS` to Markdown format for proper rendering in browser 
(issue #1231).

- Indicate model code that produces warnings about unknown nimbleFunctions
(issue #370).

## BUG FIXES

- Avoid error occurring when a model variable name starts with "logProb"
(PR #1240).

- Avoid error occurring when a model variable is named "i" (PR #1239).

- Prevent infinite recursion in particular cases in conjugacy checking 
(PR #1228).

- Fix bug in simulating from `dcar_normal` nodes when multiple nodes passed to
simulate (issue #1238).

- Fix error message about duplicate node declarations (PR #1233).

- Fix another issue with long variable names (PR #1217).

- Fix warning related to `dataNodes` in WAIC. 

## DEVELOPER LEVEL CHANGES

- Remove use of bitwise `|` and `&` operators in C++ code, per CRAN request.

- Refactor `nimbleMCMC` to pull out model creation (PR #1223).

- Fix an issue with nested `nimbleList`s on MacOS (PR #1213).

#                            CHANGES IN VERSION 0.12.2 (February 2022)

## BUG FIXES

- Fix bug in Bayesian nonparametrics (BNP) functionality that gives incorrect 
MCMC results with the dCRP distribution when the parameters of the mixture 
components (i.e., the clusters) have hyperparameters (i.e., the base measure 
parameters) that are unknown and sampled during the MCMC (PR #1202).

- Error trap a similar case to the bug just above that would cause incorrect 
reversible jump MCMC sampling when a parameter being moved into/out of the 
model has an unknown hyperparameter (PR #1193). 

- Fix a bug preventing setup of conjugate sampler for dwishart or dinvwishart nodes 
when using dynamic indexing (PR #1200).

- Fix a bug preventing use of truncation bounds specified via `data` or `constants`
(PR #1209).

- Fix a bug preventing MCMC sampling with the LKJ prior for 2x2 matrices 
(PR #1195).

- Fix a bug in the ordering used by the setSamplers method of MCMC configuration 
objects, occuring when the samplers are specified as a vector of node or 
variable names (PR #1178).

- Fix a bug in `runCrossValidate` affecting extraction of multivariate nodes 
(PR #1206).

- Fix a bug producing incorrect subset assignment into logical vectors in 
nimbleFunction code (PR #1205).

- Fix a bug preventing use of `nimbleExternalCall` with a constant expression 
(PR #1204).

- Fix a bug preventing use of recursion in nimbleFunctions without setup code 
(PR #1203).

- Fix a bug in handling `nimSeq` default `by` value (PR #1199).

- Fix a bug in name mangling affecting certain long names of lifted nodes 
(PR #1192).

- Fix access to member data more than two dimensions in a nested nimbleFunction
(PR #1046).

## DEVELOPER LEVEL CHANGES

- Clean up message verbosity handling (PR #1191).


#                            CHANGES IN VERSION 0.12.1 (October 2021)

## BUG FIXES

- Fix a bug introduced in conjugacy processing in version 0.11.0 (PR #1087)
that causes incorrect MCMC sampling only in specific cases.  The impacted cases
have terms of the form "a[i] + x[i] * beta" (or more simply "x[i] * beta" or
"a[i] + beta"), with beta subject to conjugate sampling and either (i) 'x'
provided via NIMBLE's constants argument and x[1] == 1 or (ii) 'a' provided
via NIMBLE's constants argument and a[1] == 0 (PR #1172).

- Reversible jump MCMC system now detects non-constant hyperparameters
of reversible jump target nodes, and will not operate on such nodes.

## DEVELOPER LEVEL CHANGES

- Use USE_FC_LEN_T-related syntax in F77 calls in dists.cpp per CRAN request
(PR #1173).


#                            CHANGES IN VERSION 0.12.0 (September 2021)

## USER LEVEL CHANGES

- Completely revamp WAIC in NIMBLE, creating an online version that does not
require any particular variable monitors. The new WAIC can calculate conditional
or marginal WAIC and can group data nodes into joint likelihood terms if
desired. In addition there is a new `calculateWAIC` function that will calculate
the basic conditional WAIC on output from MCMC output (either using an MCMC
object or a matrix of samples), without having to enable the WAIC when creating
the MCMC (PR #1136).

- Add LKJ distribution, useful for prior distributions for correlation
matrices, including RW samplers executed on an unconstrained transformed
parameter space, and assigned by default during MCMC configuration (issue #993).

- Improve formatting of standard logging messages (PR 1162).

- Check for negative indexes in nimbleFunction code and error out
during compilation (PR #1166).

- Update handling of missing indexing on the output of functions
(e.g., `myfun()[,1]`) in model code to allow missing indexing (PR #1156)
and to error trap if missing indexing prevents conjugacy checking (PR #1155).

- Improve error message when a user tries to use a dynamic vector of indexes
(PR #1151).

- Add error trapping for use of `[[` in model code (PR #1148).

- Clarify documentation of `getParents` arguments (issue #1143).

## BUG FIXES

- Fix an error in the sampler for the proper CAR distribution (issue #1157)
that gives incorrect MCMC results when the mean of the proper CAR is not the
same value for all locations, e.g., when embedding covariate effects directly
in the `mu` parameter of the `dcar_proper` distribution (PR #1158).

- Fix `isData` to return TRUE whenever any elements of a multivariate data node
are flagged as data. As a result, attempting to carry out MCMC on the non-data
elements will now fail. Formerly if only some elements were flagged as data,
`isData` would only check the first element, potentially leading to elements
flagged as data being overwritten (PR #1165).

- Error trap cases where BNP model had differing number of dependent
stochastic nodes (e.g., observations) or dependent deterministic nodes
per group of elements clustered jointly. Previously we were not error trapping
this, and incorrect MCMC results would be obtained (PR #1152).

- Fix casting upon assignment for vectors in nimbleFunctions, avoiding
compilation errors (PR #1168).

- Fix model processing errors involving long lines of model code related to
`deparse` producing multiple lines of output (PR #1167).

- Ensure that variables in uncompiled modelValues are in the same order as
those in compiled modelValues (PR #1164).

- Update `stick_breaking` to work if given only one probability (PR #1161).

- Robustly handle cases where a user specifies no monitors for an MCMC
(PR #1140).

- Avoid conflicts between parallel worker processes in parallelized
cross-validation (PR #1147).

- Fix `dcat` and `rcat` to return `NaN` and `NA`, respectively, when the
probability vector includes a negative value (issue #1146) as well
as having `rmulti` return `NA` values when given a negative probability.

- Fix `RW_llfunction` sampler example in manual (issue #1135).


## DEVELOPER LEVEL CHANGES

- Use `getParents` from the model API in Bayesian nonparametric functionality
in place of the former bespoke `getParents` written just for BNP (PR #1163).

- Remove legacy S macros PROBLEM, MESSAGE, ERROR, WARN per CRAN request
(PR #1150).

- Pass `where` argument into `nf_checkDSLcode` to improve experience of
package development for packges depending on `nimble` (PR #1124).

- Update code for generating predefined C++ code for nimbleLists (issue #1142).

- Update our package preparation workflow to work with recent versions of
`roxygen` (NCT issue 310).

- Update `configure.ac` (and resulting `configure`) via `autoupdate` per
CRAN request.

- Fix testing issues with `expect_equal` where tests acted as if the checking
used absolute tolerances rather than the actual relative tolerances, as well 
issue with `expect_lt` applied to vectors, all of which would have
prevented our testing from detecting errors.


#                            CHANGES IN VERSION 0.11.1 (May 2021)

## USER LEVEL CHANGES

- Add information about categorical sampler and univariate version of ESS
   sampler to `help(samplers)`.

## BUG FIXES

- Fix to the `posterior_predictive_branch` MCMC sampler, to update
   the log-probabilities of the sampled posterior predictive nodes (PR #1127).

#                            CHANGES IN VERSION 0.11.0 (April 2021)

## USER LEVEL CHANGES

- Add new `posterior_predictive_branch` MCMC sampler, which is
   automatically assigned to trailing dependency node networks of entirely
   non-data nodes (jointly posterior predictive branches).  This sampler
   simulates jointly from the predictive distribtion of these posterior
   predictive node branches, and is designed to improve MCMC mixing of the
   branch, and consequently of the entire model (PR #1086).

- Allow use of elliptical slice sampler for univariate nodes, which can be
   useful in multimodal problems (PR #1109).

- Add a `getParents` method to the model API, allowing one to determine parent
   nodes, analogous to use of `getDependencies` to determine child nodes
   (PR #1094).

- Add `getConditionallyIndependentSets` method (not yet documented) to the model
   API, allowing one to determine nodes that are conditionally independent of
   each other given parent nodes (PR #1094).

- Improve efficiency of conjugate samplers by avoiding unneeded calculations
   when a conjugate relationship does not involve shifting or scaling (PR #1087).

- Allow use of `nimNumeric`, `nimMatrix`, `nimArray` in model code (PR #1096).

- Add progress bar to `getSamplesDPmeasure` (NCT issue 110).

- Allow model definition using `if` without `else`, fixing a longstanding
   oversight (PR #1104).

- Improve warning when multiple nodes provided to `getParam` (PR #1118).

- Check during model building for unnamed elements of `data` and `inits`
   (PR #1117).

- Remove error trapping to prevent use of variables in defining node names,
   such as `getDependencies('y[idx]')` as this is hard to check robustly and
   efficiently (PR #1122).

- Improve error messages when reporting `getParam` cannot calculate a parameter
   when checking a model (PR #1112).

- Error trap cases where model nodes are defined in two different declarations,
   adding check for overlapping multivariate nodes (PR #1110).

- Improve error trapping of mis-formed stochastic declarations in models
   (PR #1106).

- Increase maximum length of compiler output when using
   `compileNimble(..., showCompilerOutput = TRUE)` (NCT issue 205).

- Point to parallelization example on r-nimble.org in relevant places of manual.

## BUG FIXES

- Fix a bug (issue #1091) causing incorrect node names when having more than
   100,000 elements in a vector node or in a dimension of a multi-dimensional
   node (PR #1092).

- Fix `getNodeNames` to return no nodes when `latentOnly` is `TRUE` and model
   contains no latent nodes (PR #1116).

- Fix checking for unknown nimbleFunction methods and improve related error
   trapping (PRs #1107, #1105).

## DEVELOPER LEVEL CHANGES

- Update our testing code/infrastructure to use latest testthat API (PR #1090).

- Shift internal code to use `model$calculate(...)` style rather than
   `calculate(model, ...)` style for various node functions (PR #1114).

- Clean up commented out code (PR #1098) and remove unused test files
   (PR #1097).

- Update to a newer (but not latest) version of Eigen to suppress some compiler
   warnings (PR #1093).


#                            CHANGES IN VERSION 0.10.1 (November 2020)
                            
## USER LEVEL CHANGES

- Add `round` argument to `samplesSummary` (PR #1077).

- `samplesSummary` function (and also `runMCMC(..., summary = TRUE)`) was made
   to be robust against non-valid values in posterior samples array (PR #1075).

## BUG FIXES

- Fix `makeParamInfo` when there is only one declID involved to address a bug
   affecting usage of `getParam`. This bug was introduced in version 0.10.0 when
   reducing memory use of `getParam` (PR #1016). This fixes incorrect behavior
   of conjugate samplers (because of incorrect inputs from `getParam`) under
   certain model structures, in particular state-space style models (PR #1080). 

- Prevent usage of marginal version of WAIC (i.e., when not monitoring all
   direct stochastic parents of data nodes); use of marginal version of WAIC in
   previous versions gave incorrect results (PR #1083).

## DEVELOPER LEVEL CHANGES

- Deprecate `samplerAssignmnentRules` system (PR #1078).

- Deprecate `autoBlock` MCMC option (PR #1079).

#                            CHANGES IN VERSION 0.10.0 (October 2020)

## USER LEVEL CHANGES

- Greatly extend BNP functionality with the CRP (Chinese restaurant process)
   distribution by allowing multiple observations to be grouped together (e.g.,
   for longitudinal or time series data) without requiring they be specified
   as a multivariate node (PR #1033).

- Add a variety of conjugate cases to BNP conjugate samplers (PR #1033).

- Greatly improve efficiency of model and MCMC building and configuration for
   BNP-based models with CRP components (PR #1033).

- Move all sequential Monte Carlo (SMC; aka particle filtering) methods to
   new package `nimbleSMC`, including various particle-filter-based MCMC
   samplers.

- Prevent use of variables in indexes of nodes, such as `y[idx]`, which was
   incorrectly being evaluated based on R scoping rules (PR #1064).

- Allow use of `logdet` in model code (Issue #1018).

- New `resetMV` argument available to `mcmc$run` method. In combination
   with `reset = FALSE`, specifying `resetMV = TRUE` will continue the current
   run of the MCMC, but discard any previously-collected samples
   (PR #1051; thanks to 'DJRP').

- New methods `setMonitors` and `setMonitors2` added for MCMC configuration
   objects.  These methods replace the current set of monitors (or monitors2)
   with the specified variables (PR #1061).

- Add `as.list` method for modelValues objects (PR #1060).

- Update `getSamplesDPmeasure` function to improve efficiency and reduce
   output size; output is now a list of matrices (PR #1059).

- Add `dimensions` argument to `nimbleMCMC` (PR #1058).

- Add method `getWidthHistory` to slice sampler to retrieve sampling history
   information (PR #1057; thanks to 'rpatin').

- Various improvements to the manual. 

## BUG FIXES

- Fix a bug in k-fold cross-validation routine (`runCrossValidate`), where
   the merging of MCMC sampler configurations was done incorrectly and causing
   incorrect results (PR #1068).

- Fix bug giving incorrect `dwish` density when using non-default S
   parameterization (PR #1017).

- Fix incorrect NaN eigenvalues in singular normalized adjacency matrices
   under `dcar_normal` (PR #1019).

- Update all MCMC sampler functions to use a new syntax for control list
   element extraction, which prevents a possible bug caused by R's partial
   matching of list names (PR #1065).

- Define auto-generated simulation ('r') functions for user-defined
   distributions in the global environment to avoid scoping issues (PR #1063).

- Update user-defined distribution processing so user-defined distributions
   can be defined inside functions (PR #1063).

- Fix bug preventing use of `dirName` argument to `compileNimble` (PR #1062).

- Fix a bug preventing model building when there are overly long names of
   model variables resulting from long deterministic expressions in model code
   (PR #1069).

- Fix `buildMCEM` so it works with a compiled model as argument (PR #1028).

- Fix `dmvt` so default unnamed parameters work (PR #1027).

- Fix error in model building in corner case where
   `makeVertexNamesFromIndexArray2` made a simplifying assumption to conclude
   a block of nodes was contiguous (PR #1026).

- Fix bug in `nimbleRcall` causing run-time warnings when `returnType` is void
   (PR #1013). 

## DEVELOPER LEVEL CHANGES

- Improve efficiency of `getParam` implementation, which improves speed for
   MCMC compilation (PR #1016). 

- Improve MCMC sampling efficiency by not copying data nodes, only data node
   logProbs, during sampler execution for various samplers (PR #1040).

- Update Travis testing to use R 4.0.

- Remove `GID_map` internal to modelValues (PR #1032).

- Remove deprecated function `getLoadingNamespace` and (deprecated) use of
   `where=getLoadingNamespace`. Also improve handling of environments set up
   by `nimbleFunction` to make it easier to write packages depending on NIMBLE
   (PRs #1029, 1011).

- Force intermediates of index range expressions to be of type 'int' for use
   in AD (PR #1024).



#                            CHANGES IN VERSION 0.9.1 (May 2020)

## USER LEVEL CHANGES

- Switched from use of `system` to `system2` to avoid problems with installation
   under R 4.0 on Windows (PR #1003).

- Modify various adaptive MCMC samplers so the exponent controlling the scale
   decay of the adaptation is adjustable by user (rather than hard-coded at 0.8
   (PR #981).

- Allow `pmin` and `pmax` to be used in models (PR #982).

- Add documentation for `is.na`,`is.nan`,`any`,`all` (PR #988)

- Add system option `MCMCuseConugacy` to control whether conjugate samplers are
   used (PR #998).

- Adds checks for `niter`, `nburnin` in the `mcmc$run` method (PR #980).

- Modify print handling in `addSampler` and `configureMCMC` (PRs #986, 989).

- Improve handling of NA values in `dCRP` to avoid error messages when building
   models (PR #994).

- Avoid monitoring top-level data nodes in models (PR #1006).

## BUG FIXES

- Modify MCMC `autoBlock` routine to only group Wishart, Inverse-Wishart, and
   Dirichlet nodes with themselves, to avoid violating the constraints of those
   nodes (PR #999).

- Fix incorrect error message from `warnRHSonlyDynIdx` when variable appears
   multiple times on right-hand side of a model expression (PR #997).

- Fix `checkDistributionFunctions` to respect default `nDim=0` when extracting
   first argument, to avoid error when dimension not specified in user-defined
   distributions (PR #992).

- Fix `print` option of `addSampler` (PR #986).

- Improve handling of cases where indexing goes beyond extent of variable in
   `expandNodeNames` and related queries of model structure (PR #977).

## DEVELOPER LEVEL CHANGES

- Use `inherits` rather than testing for equality of `class(object)` (PR #988).

#                            CHANGES IN VERSION 0.9.0 (December 2019)

## USER LEVEL CHANGES

- Added iterated filtering 2 (IF2) algorithm for estimating parameters by
   maximum likelihood in models for which SMC algorithms (particle filters)
   can be used. 

- Added dmnorm-dmnorm conjugacy for BNP mixture models using the dCRP
   distribution (PR #936).

- Added detection of normal-inverse-gamma conjugacy in BNP mixture models using
   the dCRP distribution in the case when intermediate nodes are present
   (PR #944).

- Cleaned up handling of linear predictors in regression-style models to avoid
   an incorrect warning when indexing the result of matrix multiplication
   (PR #929).

- Modified handling of NAs in model nodes under truncation or dynamic indexing
   to avoid printing warnings when calculating the uncompiled model (PR #920).

- Improved speed of 'configureMCMC' (PR #972 and PR #974).

- Added the scalarComponents argument to addSampler method of MCMC configuration
   objects, which enables adding independent univariate samplers to all scalar
   components of a specified target node or variable.

- removeSamplers, setSamplers, and printSamplers methods of MCMC configuration
   object now also accept node names (or sampler indices) through the "..."
   argument.

- Fixed a corner case of pathological behavior in the reflection version of the
   RW sampler when there are both lower and upper bounds and the posterior is
   poorly informed, in which case the proposal scale could grow to be very large
   (PR #948).

- Improved output from printSamplers, including indicating if conjugacy is
   actually used in BNP mixture models using the dCRP distribution.

- Added a warning when user specifies a model with non-dynamic indexes specified
   in 'data' or 'inits' rather than the more computationally-efficient placement
   in 'constants' (PR #963).

- Added a warning when 'inits' argument to nimbleModel is not a list.

- Modified file path handling in Windows that should work on more file-system
   configurations (PR #957).

- Added tip for dealing with installation problems and added information on how
   to properly parallelize nimble functionality to the user manual.

## BUG FIXES

- Remove 'pfResample' option (not the default) for particle MCMC (PMCMC)
   samplers as resampling is not justified by the theory of PMCMC samplers.

- Fixed bug concerning incorrect dependency tracking that affected particle
   filtering for models with deterministic nodes between multiple latent states
   of the same time step.

- Fixed bug concerning tracking re-weighting and resampling in auxiliary
   particle filter in the default 'saveAll=FALSE' case.

- Modified bootstrap filter to always resample in last time step (even if
   'threshold' control argument is less than 1) so that 'mvEWsamples' can always
   be considered equally-weighted.

- Modified auxiliary particle filter to omit multiplication by p(x_t+1 | x_t)
   in look-ahead. This now follows Pitt and Shephard (1999).

- Fixed incorrect handling of rate parameter in ddexp distribution for
   nimbleFunctions and direct use from R that was causing rate to be ignored and
   default scale value to be used. Use of rate parameter in ddexp in models was
   not affected by this issue (PR #950).

- Fixed bug causing incorrect sampling in BNP mixture models with beta-binomial
   or beta-negative-binomial conjugacy (PR #922).

- Fixed bug that prevented recycling rule for distributions ddexp, dexp_nimble,
   dt_nonstandard, and dinvgamma (and their 'p' and 'q' versions) when one or
   more parameter vectors were longer than the primary argument (PR #954).

- Fixed bug preventing recycling rule for distributions dt and dt_nonstandard
   (PR #956).

- Fixed bug in MCMC crossLevel sampler affecting models with non-identity link
   functions, so that underlying conjugate sampler functions no longer modify the
   target node value (PR #925).

- Fixed a bug in corner case of MCMC conjugacy checking system, when linear
   scale of target node is identically equal to zero.

- Fixed a bug in corner case of MCMC conjugacy checking for BNP models -
   removing conjugacy detection in normal-inverse-gamma case where the data
   (i.e., dependent node) variance is scaled.

- Fixed a bug preventing use of 'data' as the name of a data node in the model
   (PR #946).

- Fixed a bug preventing use of 'nimbleModel' via 'do.call' (PR #969).

- Fixed missing export of 'any_na' and 'any_nan' for use in nimbleFunctions.

- Fixed conjugacy checking to detect conjugacy in more complicated cases
   involving linear predictors (PR #958).

## DEVELOPER LEVEL CHANGES

- Add rudimentary support for multiple inheritance (PR #943).

- Match int and unsigned int to silence Windows compiler warnings (PR #914).

- Modified some naming involved in reversible jump MCMC for variable selection.
   Added error trapping if the node being considered has a multivariate prior
   (PR #964).


#                            CHANGES IN VERSION 0.8.0 (June 2019)

## USER LEVEL CHANGES

- Added reversible jump MCMC for variable selection via configureRJ().

- Greatly improved speed of MCMC sampling for Bayesian nonparametric models
   with the dCRP distribution by not sampling parameters of empty clusters.
   (PR #855)

- Improved speed of MCMC configuration (i.e., configureMCMC()) is available
   by setting nimbleOptions(oldConjugacyChecking = FALSE) and/or
   nimbleOptions(useNewConfigureMCMC = TRUE).  These new features are considered
   to be in beta testing.  Default values for these nimbleOptions give old
   behavior. (PR #910, PR #896)

- Removed compareMCMCs() and MCMCsuite(), which are being re-written for
   release in a separate package available currently at
   https://github.com/nimble-dev/compareMCMCs.

- Added printing option for MCMC configuration objects:
   conf$printSamplers(byType = TRUE), will display the node names being sampled,
   grouped together by the sampling algorithm acting on them (e.g., RW sampler).
   (PR #901)

- Added options for control of MCMC behavior: MCMCmonitorAllSampledNodes to
   monitor all sampled nodes and multivariateNodesAsScalars for sampling scalars
   within multivariate nodes. (PR #895, #893)

- Added detection of dCRP clustering of regression coefficients by detecting
   conjugacy when a linear combination is present.

- Added support for recycling-rule in compiled code for dlogis, rlogis,
   qlogis, and plogis. (PR #911)

## BUG FIXES

- Fixed bug in conjugacy checking of CAR model structures. (PR #871)

- model setInits method now checks for unnamed list elements. (PR #880)

- Allow use of TRUE/FALSE for non-parameter arguments in functions in model
   code. (PR #904)

- Correctly manage the RNG state when calling from C++ to R via a nimbleRcall.
   (PR # 897)

- Handle issues when numerical underflow causes exact zeros as proposals in
   RW_dirichlet sampler. (PR #885)

- Fixed a bug in conjugacy detection when models use sum(x[1:n]). (PR #890)

- Fixed bugs in copying in certain cases when the LHS and RHS involve the same
   variable. (PR #886)

- Fixed inefficient handling in configureMCMC of some models with gamma
   distributions. (PR #896)

## DEVELOPER LEVEL CHANGES

- Improved scoping of nimbleFunctions without setup code and nimbleRcall
   and nimbleExternalCall to avoid issues in packages that depend on NIMBLE.
   (PR #889)

- Added a check to detect and error out when clustering deterministic nodes
   in dCRP-based models. (PR #906)

- Added a check to detect and error out if multiple CRP indexes used in an
   expression. (PR #892). 

- Added better error messages for when the slice sampler reaches maximum number
   of contractions (PR #903) and when MCEM optimization calculates a non-valid
   log-density (PR #887).

- Added some links and tweaked documentation for various functions.

- Improved speed in some uses of getNodeNames(). (PR #902)

- Added error trap for a case of incorrect syntax accessing modelValues.
   (PR #899)

- Added various other error traps. (PR #884, #883)

- Removed unused support for running via Tensorflow rather than generating
   and compiling C++. (PR # 882)

- Only calculate eigenvalues in CAR handling to avoid slow (and unneeded)
   eigen computations of eigenvectors. (PR #872)

#                            CHANGES IN VERSION 0.7.1 (March 2019)

## USER LEVEL CHANGES

- Add support for 6-dimensional arrays in model code and in nimbleFunctions.

- Allow use of besselK in model code.

- Add normal-normal conjugacy detection in multivariate regression structures
   using inprod, sum, and matrix multiplication.

- Allow use of 'x' in addition to 'value' in 'types' argument to
   registerDistributions.

## BUG FIXES

- Fix bug in MCMC sampling of dCRP nodes in non-conjugate situations
   introduced in Version 0.7.0. (Issue #859)

- Fix bug in findClusterNodes for dCRP nodes that was not detecting certain
   cases of clustering of multiple parameters in a model code statement and
   thereby not properly handling MCMC configuration for such nodes. (PR #861)

- Avoid protection stack overflow in working with large models, a problem
   introduced in 0.7.0 due to modifications of how PROTECT is used C++ code
   in response to new CRAN checks. (Issue #852)

- Fix bug in corner case of model initialization. (Issue #857)

## DEVELOPER LEVEL CHANGES

- Per CRAN request, remove use of hard-coded paths. One case was linking
   to libnimble.so from nimble.so - everything is now compiled into nimble.so.
   Second case was determination of path to libnimble.so at install time;
   this is now done at run-time. (Issue #858)

#                            CHANGES IN VERSION 0.7.0 (January 2019)

(See also changes below in Version 0.6.13 as Version 0.6.13 existed only
briefly on CRAN because we needed to fix a few minor packaging issues
raised by CRAN.)

## USER LEVEL CHANGES

- Greatly reduced time for setting up model initialization in buildMCMC, fixing
   some computational slowdowns introduced in version 0.6.13's fix of a
   shortcoming in model initialization.

## DEVELOPER LEVEL CHANGES

- Fixed some uses of PROTECT in C++ code, flagged by R's rchk.

#                           CHANGES IN VERSION 0.6.13 (January 2019)

## USER LEVEL CHANGES

- Various changes to greatly improve efficiency of sampling for Bayesian nonparametric (BNP) mixture models using the dCRP distribution.

- Changes to samplers RW, RW_block and categorical that could yield more efficient execution in cases where prior or latent-state distributions have bounded domains (i.e., boundaries that define valid values).

- Added double exponential (Laplace) distribution.

- New "RW_wishart" MCMC sampler, for sampling non-conjugate Wishart and inverse-Wishart nodes.

- Added normal-inverse gamma conjugacy for BNP mixture models using the dCRP distribution.

- getSamplesDPmeasure now works with multivariate cluster parameters and multiple cluster parameters.

- getSamplesDPmeasure now takes an optional argument for the error level that determines the truncation level of the random measure and now returns a list containing the samples and truncation level.

- addSampler method of MCMC configuration objects now accepts ... (dot dot dot) arguments, and uses them as control list elements for the sampling algorithm.

- 5-dimensional arrays now allowed in models.

- Added warning when a node is defined multiple times in model code.

- Revamped handling of showCompilerOutput option in compileNimble so that it is more clear to users to look at the output.

- Updated documentation of WAIC to make more clear that calculation depends on what is monitored.

- Added warning message that Liu-West filter often doesn't work well.

- We now provide an HTML version of the manual at r-nimble.org.

- Allow negative weights in dcar_normal and negative values in C argument to dcar_proper.

## BUG FIXES

- Fixed bug that produced incorrect WAIC calculations when using multiple chains for models with at least one non-scalar monitored variable. Also improved robustness of WAIC calculation.

- Fixed bugs in conjugate samplers for CRP distribution: CRP_conjugate_dgamma_dnorm, CRP_conjugate_dbeta_dbin, CRP_conjugate_dbeta_dnegbin, CRP_conjugate_dgamma_dinvgamma, CRP_conjugate_ddirch_dmulti.

- Fixed bug in categorical sampler that causes value of 1 to be sampled every time when all log probabilities underflow upon exponentiation.

- Improved handling and error-trapping of less common model structures (in particular less common styles of indexing) for models using the dCRP distribution.

- Fixed error in computing truncation level of G in getSamplesDPmeasure function

- Fixed bug that was obviating comparison with gold files in testing system.

- More sophisticated model initialization routine for correctly initializing complex state-space models.

- Fixed a problem with dimension handling in matrix2VecNimArr.

- Fixed issue with checking for valid dynamic index values.

## DEVELOPER LEVEL CHANGES

- Various additional tests for Bayesian nonparametric mixture models.

- Improved MCMC efficiency in RW and RW_block samplers by first checking prior and rejecting without further computation if log prior is -Inf or NA/NaN.

- Modify use of exists() to avoid looking outside local context in some places.

- Improved various error messages and error trapping cases.


#                            CHANGES IN VERSION 0.6-12 (July 2018)

## USER LEVEL CHANGES

- New option for printing MCMC samplers of particular type(s): conf$printSamplers(type = "conjugate"), for example.

- Burnin is now handled natively by mcmc$run method (rather than as a post-processing step in runMCMC).

## BUG FIXES

- Corrected calculation of weights in bootstrap particle filter (calculation was omitting previous weights when particles were not resampled). Also clarified help information regarding equally-weighted samples when resampling not done in an iteration.

- Add return value to auto-generated 'r' function for user-defined distribution to avoid tripping error trap introduced in version 0.6-11.

- Avoid checking for ragged arrays when building models; this fixes an overly zealous check introduced by accident in version 0.6-11.

- configureMCMC no longer assigns a RW_block sampler to non-conjugate inverse-Wishart nodes.

## DEVELOPER LEVEL CHANGES

- Shifted to bookdown/Rmarkdown format for user manual.

- Updated configure.ac to use CXX11 not CXX1X and to avoid calling AC_PROG_CXX before CXX and related variables set based on R configuration, per CRAN request.


#                            CHANGES IN VERSION 0.6-11 (June 2018)

## USER LEVEL CHANGES

- Bayesian nonparametric mixture models can now be used in BUGS code, in particular Chinese Restaurant process (using the 'dCRP' distribution) and stick-breaking (using the 'stick_breaking' function) representations of Dirichlet process models. This allows mixture models with an unknown number of components. NIMBLE's default MCMC configuration will assign specialized samplers to relevant nodes of the model. For this release, this is a beta (experimental) feature.

- Four new resampling methods are available for use in the auxiliary and bootstrap filters.  Resampling methods can be specified by the resamplingMethod control list argument to buildBootstrapFilter and buildAuxiliaryFilter.

- User-defined filtering algorithms can now be provided to the RW_PF and RW_PF_block samplers.

- MCMC thinning intervals (thin and thin2) can now be modified at MCMC runtime.  These may be passed as arguments to either mcmc$run or to runMCMC.

- Both nimbleMCMC and runMCMC functions now drop burnin samples as "pre-thinning", if a thinning interval is specified.

- Increased functionality for the setSeed argument in nimbleMCMC and runMCMC functions.

- New functionality in MCMC, to specify the order in which sampler functions are executed.  This can allow for samplers being repeatedly executed, interleaved, or omitted.  This ordering is specified as part of the MCMC configuration.

- Invalid dynamic indexes are reported and now result in NaN values when calling calculate or simulate but no longer cause execution to error out. This allows MCMC sampling to continue when an invalid index is proposed and rejected.

- The model is fully (re-)initialized between MCMC runs in MCMCsuite(), so that every MCMC method will start from identical conditions.  This was usually but not always the case previously.

## BUG FIXES

- Slice, AFSS, and ESS samplers now have a maximum number of contractions to avoid infinite loops under unusual circumstances.

- Fixed error in generating getParam nodeFunctions with user-supplied distributions using integer types.

- In multivariate distribution arguments, general expressions continue to be disallowed, but expressions inside indexing brackets are now allowed. This allows, e.g., hidden Markov models implemented with a transition matrix and dcat to describe the time-dependence of latent states.

- Fixed bug involving getDependencies with downstream = TRUE.

- Fixed bug with using setSize in uncompiled nimbleFunction.

- Fixed bug in AF_slice sampler causing incorrect dependency updating that affected state-space type models with intermediate deterministic nodes.

- Fixed issues with setInits in complicated cases.

## DEVELOPER LEVEL CHANGES

- Updated some testing to function better with version 2.0.0 of testthat.

- Now allow "to" to be a CmodelValues in R version of nimCopy.

- Added better error trapping when an array is used with an index before it has been created in a nimbleFunction.

- Fixed infinite recursion issue.

- Better error message when invalid parameter name supplied to getParam.

- Added a check that there is a return statement when a returnType is provided.

- Improved error trapping for invalid indexing in BUGS models.


#                           CHANGES IN VERSION 0.6-10 (March 2018)

## USER LEVEL CHANGES

- Data can now be provided as a numeric data frame rather than a matrix.

- To run WAIC, a user now must set 'enableWAIC' to 'TRUE', either in NIMBLE's options or as an argument to buildMCMC().

- If 'enableWAIC' is 'TRUE', buildMCMC() will now check to make sure that the nodes monitored by the MCMC algorithm will lead to a valid WAIC calculation.

- Deprecated use of identityMatrix() in favor of diag().

- Some steps of model and algorithm building and compilation are faster.

- Compiled execution with multivariate distributions or function arguments may be faster.

## BUG FIXES

- Fixed bug related to handling of MCMC history that caused MCMC to stop.

- Fixed a bug in WAIC calculation where downstream nodes were not being re-simulated for certain combinations of models and monitored variables.

## DEVELOPER LEVEL CHANGES

- Made non-scalar arguments in model distributions, like x[1:10] (i.e., those without gaps in indexing), avoid copies for faster C++ calculations.

- Moved parse step of parseEvalNumericManyList() to C++ for faster model processing.

- Sped up makeParamInfo() and related functions/methods to speed up compilation, particularly in models with dynamic indexing.

- Sped up checkForSelfParents() for faster model processing.

- Clarified some error and warning messages when user attempts to dynamically index a constant and when components of user-defined distributions are missing.

- Removed extraneous break/next statements flagged by CRAN.

#                           CHANGES IN VERSION 0.6-9 (January 2018)

## USER LEVEL CHANGES

- Dimensions will now be determined from either 'inits' or 'data' if not otherwise available. 

- Can now specify "nBootReps = NA" in the runCrossValidate() function, which will prevent the Monte Carlo error from being calculated.

- runCrossValidate() now returns the averaged loss over all k folds, instead of the summed loss.

- Added besselK function to DSL.
						   
## BUG FIXES

- Fixed bug that prevented use of dynamic index combined with an index specified by a constant because of duplicated dynamicIndex nodeFunctionNames.

- Added error trapping for dynamic indexing of constants.

- Fixed a bug where calling nimbleModel() on a BUGS model with a node that was its own parent caused R to crash.

- Fixed a bug in runCrossValidate(), where for certain models using the "random" foldFunction would cause an error.

- Fixed a name conflict that caused an error when user-defined distributions had parameters named 'lower' or 'upper'.

- Fixed WAIC bug in getting dependencies of logProb variables.

## DEVELOPER LEVEL CHANGES

- Fixed configure error that caused issues with -fpic compiler flag.

- Removed use of "include" directive in Makefile on a path derived from R_HOME variable at request of CRAN.

- Removed pragma regarding -Wno-ignored-attributes from NimArrBase.h to satisfy CRAN policy.


#		           CHANGES IN VERSION 0.6-8 and 0.6-7 (November 2017)

Note: versions 0.6-8 and 0.6-7 are the same except for a non-user-facing change to the withNimbleOptions function example to pass CRAN checks.

## USER LEVEL CHANGES

- Addition of dcar_proper distribution, the proper Gaussian conditional autoregressive (CAR) distribution, and MCMC sampling support.

- Addition of nimbleMCMC function, providing the most direct one-line invocation of NIMBLE's MCMC engine.

- Addition of runCrossValidate function, which conducts k-fold cross-validation of NIMBLE models fit by MCMC.

- Improved efficiency of MCMC conjugate samplers.

- Added dcat-ddirch conjugacy.

- Added error-trapping for MCMC where there is a Wishart node without conjugacy.

- Added warning that RW block sampler may behave poorly with default initial proposal covariance if elements are on very different scales.

- Modified a variety of warning messages. 

- the nimEigen() function now has a symmetric argument, which can be set to TRUE if a matrix is guaranteed to be symmetric (default = FALSE).

- nimbleModel no longer outputs the names of all uninitialized variables in a model.  Instead, users are directed to use the new $initializeInfo() method if uninitialized variables are detected. 

- Sped up WAIC calculation for large models.

## BUG FIXES

- Fixed MCMC autoBlocking procedure, to work with new sampler default values system.

- Fixed MCMC conjugacy system, detection of indexed node names vs. indexed expressions.

- Fixed bug in MCEM where bounds for multivariate top-level parameters were not stored correctly.

- MCEM outputs a warning message if the provided model has discrete top-level parameters.

- Fixed a bug in conjugacy processing when there were multiple dynamically indexed nodes in an expression.

- Fixed a bug in handling of negative default values.

- Fixed a bug in compiling multiple models in the same session.

- Fixed a bug with step() for vectors.

- Fixed handling of non-trivial conditions for a while() loop.

- Results from R's eigen() and NIMBLE's eigen() functions should now match for non-symmetric matrices.

## DEVELOPER LEVEL CHANGES

- Added new method to generated nimbleFunction C++ classes to copy some member data from R to C++ as a batch, making this step faster.

- Replaced reference classes with R6 classes for exprClass.

- Replaced reference class with S3 class for nodeFunctionVector.

- Added experimental system for expansion of model macros.

- Removed -g debugger flags from CppCode/libnimble.a in unix-alike systems, thereby greatly reducing library size. 

- Removed various legacy browser() statements.

- Cleaned up testing system to make more consistent across test files and in handling of known failures.

- Flag NIMBLE as requiring C++11, mostly for future purposes.

- Fix compilation bug when nimbleFunction has a "." in its name.

#                           CHANGES IN VERSION 0.6-6 (July 2017)

## USER LEVEL CHANGES

- One can now use dynamic indexes in BUGS code - indexes of a variable no longer have to be constants but can be other nodes or functions of other nodes. This allows mixture models with unknown membership. For this release, this is a beta feature and needs to be enabled with nimbleOptions(allowDynamicIndexing = TRUE).

- The Intrinsic Gaussian conditional autoregressive (ICAR) distribution can now be used in BUGS code using the dcar_normal distribution, which behaves similarly to BUGS' car.normal distribution. Specialized MCMC samplers will be assigned to nodes that have this distribution. 

- Enabled use of optim (equivalently, nimOptim) in nimbleFunctions.

- One can now calculate WAIC for model selection using the calculateWAIC method for MCMC objects.

- There is a new nimbleExternalCall function that allows one to call separately compiled code from a nimbleFunction or a model (this feature is experimental for now).

- There is a new nimbleRCall function that allows one to call an R function from a compiled nimbleFunction or model (this feature is experimental for now).

- Enabled use of next in nimbleFunctions.

- A resetFunctions flag have been added to buildMCEM.

- Improved error management from compileNimble

- NaN is now handled correctly as a constant by compiler

- Variable names with "." are now handled correctly by the compiler.

- pi is now handled as a constant or variable in a nimbleFunction (but cannot be used as a constant in model code).

- There is a new returnESS method added for the bootstrap filter and auxiliary particle filter.

## DEVELOPER LEVEL CHANGES

- MCMC sampler control list default values are now internal to sampling algorithms, rather than a system-level option.

- Improvements added to testing system, including safeguarding against silent failures, use of gold files, batching of compiler tests, more expected failures, conversion to running all tests using test_package, and parallelization on Travis-CI.

- Names assigned to package-defined nimbleFunctions to improve C++ readability.

- Use of R registrations for C functions was cleaned up.

- Lots of vestigial code was cleaned up.

- Un-referenced scalar arguments to rankSample.

- We now use -fPIC rather than -fpic in most cases.

- Updated the version of Eigen that we use to version 3.3.4.

## BUG FIXES

- Several bugs fixed in the automated factor slice sampler (AF_slice sampler), having to do with correctly resetting the internal member variables when an MCMC is run multiple times

- Fixed harmless bug in managing split vertex numbering.

- Fixed bugs in test-mcmc that resulted in ignoring comparisons against known results for node names with brackets or underscores.

- Allow sum of length-zero vectors to compile correctly.

#                           CHANGES IN VERSION 0.6-5 (June 2017)

## USER LEVEL CHANGES

- The functions returned from a call to buildMCEM() now include an estimateCov() method, which can be used to estimate the asymptotic covariance of model parameters at their MLE values.  

- nimbleLists can now be used in nimbleFunctions that do not have setup code.

- Enabled use of c(), rep(), seq(), diag() and ':' in BUGS code.

- Added dflat and dhalfflat for improper uniform (prior) distributions on the real line and positive real line; includes conjugacies for dnorm dependents (in mean) with dflat distribution and dnorm dependents (in sd) with dhalfflat distribution.

- Allow compilation of cases like model$getParam(nodes[i], 'mean'), even if nodes is empty, in a nimbleFunction instance (in which case the code would not execute, but it must be able to compile).

- Added inverse-wishart distribution to NIMBLE's provided distributions.

- Improved a variety of error-trapping.

- New nimbleOptions()$verboseErrors, default FALSE.  If TRUE, the call stack will be output when an error is trapped.

## DEVELOPER LEVEL CHANGES

- Added error trapping to report useful error message when the type of the returned object does not match the declared returnType of a nimbleFunction.

- Added onAttach message to point users to website and manual.

- Updated script for creating roxygen-based documentation.

- Added more complete handling of split vertices during graph construction and model-querying methods.

## BUG FIXES

- Fixed underflow error in MCMC binary sampler.

- Fixed bug preventing use of nimC() with an argument like model[[node]].

- Fixed model$getDependencies(nodes) to return only what is needed (and no extra) when nodes contains a subset of elements of a non-scalar node.  e.g. model$getDependencies('x[1]') if BUGS code has 'x[1:5] ~ dmnorm(...)'.

- Fixed bug in model definition processing where some (rare) syntax could lead to very inefficient computation.

- Fixed bug when multiple nimbleFunctions with multiple levels of nesting led to incomplete inclusion of multiple .o files in compilation.

- Fixed compilation of cases where a scalar is computed from vectors and then used in scalar operations, e.g. inprod(X, Y) + a.

- Fixed bug where return(rnorm(2)) and some other cases of a return argument that is an expression returning a vector would fail to compile.

- Fixed bug where identical BUGS declarations with expression arguments in different for-loops could yield incorrect models.

- Fixed bug where a LHS node split in the middle by a RHS usage would fail during model building.

- Fixed bug in compiling code that constructs logical scalars and vectors.

- Fixed bug preventing installation on Solaris.

- Fixed issue with PROTECT flagged by rchk.

- Fixed C++ array deletion issue in two multivariate densities.

#                           CHANGES IN VERSION 0.6-4 (April 2017)

## USER LEVEL CHANGES

- The user manual has been heavily revised and reorganized, with material on NIMBLE programming now more carefully organized.

- Versions of R functions c(), seq(), `:`, rep(), diag(), diag()<-, and which() are now allowed in nimbleFunction run code and in BUGS code.

- Functions numeric(), integer(), logical(), matrix(), and array() can now take non-scalar value arguments and will populate the newly created object by using the contents of value sequentially.

- Distribution functions (d, p, r and q) now follow R's recycling rule, allowing non-scalar arguments and return values.  Unlike R, return values will never have dimension > 1.  Hence matrix arguments will result in a vector return value.

- Logical vectors and operators are now supported.

- Indexing of vectors and matrices can now use arbitrary numeric (with integer contents) and logical vectors.

- Added inverse-gamma distribution to NIMBLE's provided distributions.

- dim() now returns an integer vector.  It is not limited to uses like dim(x)[i].

- Updates to MCMC configuration API:  Member methods conf$addMonitors() and conf$addMonitors2() now accept argument ..., allowing multiple monitors to be added using the syntax conf$addMonitors('x', 'y', 'z').  Introduced a new member method, conf$printMonitors(), which nicely prints the current MCMC monitors.  This takes the place of the former method conf$getMonitors().  Now, conf$getMonitors() returns a character vector of the MCMC monitors, and the new method conf$getMonitors2() returns a character vector of the MCMC monitors2 field.

- Added options to the printSamplers() method of MCMC configuration objects that control the level of detail displayed.  These allow control for displaying default values, displaying values of non-scalar elements, and displaying the dependency lists of conjugate samplers.

- Indexing of nodes as an argument to values() now allowed, so syntax such as `x <- model$values(nodes[i])` or `values(model, nodes[i]) <- x` is supported.

- More general indexing of vectors of node names now allowed in calculate(), simulate(), calculateDiff(), getLogProb() and values().  

- The adaptive factor slice sampler can now be used in NIMBLE's MCMC framework by specifying "AF_slice" as the sampler type.  

- New sampling algorithm (RW_dirichlet) added to MCMC engine, for sampling non-conjugate Dirichlet distributions.

- Added checking at time of model definition to catch errors where implied dimensions of a variable are different in different BUGS declarations.

- nimbleLists are a new list-style data structure that can now be created and used in nimbleFunctions.  nimbleList definitions can be created in R's global environment, or in setup code.  Instances of nimbleLists can be created in R's global environment, in setup code, or in run code.  nimbleLists can contain other nimbleLists.  nimbleLists can be used as arguments to nimbleFunction run and other methods, and returned from such nimbleFunctions.  

- eigen() and svd() now work to conduct eigendecompositions and singular value decompositions on matrices.  These are synonyms for nimEigen and nimSvd, respectively.  If called from R, nimEigen() and nimSvd() will execute precompiled C++ code that conducts these decompositions using the Eigen library.  Both decompositions return nimbleList objects.  nimEigen() and nimSvd() can also be used in BUGS code.

- Filtering algorithms now have "initModel" control list option.  If initModel = TRUE, the model will be initialized at the start of the filtering algorithm.  Defaults to TRUE. 

- User-defined distributions may now be defined based solely on a density ('d') function, without a simulation ('r') function. Algorithms that need to use the simulate function will fail if it is not provided.

- Checking of conjugacy is now much quicker in cases where there are recursive dependencies.

- Some refinement to the error checking of DSL code that occurs when a nimbleFunction is defined. 

## DEVELOPER LEVEL CHANGES

- Reduced copying with values and modelValues.

- Templates compatible with Eigen library to implement R's recycling rule can be used for other functions when needed.

## BUG FIXES

- Fixed bug preventing use of round() in nimbleFunction and BUGS code.

- Fixed bugs in dmnorm_chol(), dmvt_chol(), dwish_chol(), and rwish_chol() that were causing input arguments to be overwritten when used in DSL run code (but not in nodeFunctions).

- Fixed bug in multiplication and division of matrix by scalar in run code when run in R.

- Fixed bug in RW_PF sampler where dependencies were not being updated after proposal value was inserted into target node.

- Fixed bugs in RW_PF and RW_PF_block samplers where stored LP0 was not being updated with potentially new log probability for target node at each iteration.

- Fixed bug in Wishart conjugacy calculations when the Wishart-distributed matrix was multiplied by another value in the dependent node precision matrix. Conjugacy for the Wishart now is handled only when the Wishart matrix is the precision of the dependent node without any scaling.


#                           CHANGES IN VERSION 0.6-3 (December 2016)

## USER LEVEL CHANGES

- Minor bug fixes: sd() and var() were not working in BUGS code, fix to MCEM code

## DEVELOPER LEVEL CHANGES

- Removed use of compileNimble in rankSample to avoid Solaris issues
   

#                           CHANGES IN VERSION 0.6-2 (November 2016)

## USER LEVEL CHANGES

- Added ability to add conjugate (Gibbs) samplers to particular nodes in an MCMC configuration object using: conf$addSampler(nodeName, 'conjugate').

- dmulti(), rmulti(), dcat() and rcat() now allow 'probs' argument with arbitrary non-negative values and internally normalize these so they sum to one, to be consistent with R's multinom behavior.

- Added optional (turned on by default) checking of nimbleFunction run code to warn users if they are calling functions (in particular, R functions) that are not part of the NIMBLE DSL.

- Added getBound() functionality that provides dynamic access to lower and upper bounds on a node (from R or in the DSL), based on the underlying distribution and any user-defined truncation. Functionality is analogous to getParam(). Now used in the reflection sampler and the MCEM algorithm.

- cleaned up and added functionality to get information about a distribution based on its name or to get distributional information about model nodes or variables. This functionality is expected to be useful in querying nodes in nimbleFunction setup code when writing algorithms. See help on getDistributionInfo for getting information based on distribution name and help on modelBaseClass for getting distribution information about nodes and variables.

- user-defined distributions can be used in BUGS code without registration, which will be done automatically behind the scenes. Parameters and dimensions are determined automatically from the nimbleFunction density function. However, in cases where users want to provide alternative parameterizations, the range of the distribution, or indication the distribution is discrete, registerDistributions should be used in its full functionality.

## DEVELOPER LEVEL CHANGES

- Various bug fixes including: cases where compileNimble would not return a correct interface for a nimbleFunction used in multiple particular ways; access to one nimbleFunction's member data from another.

#                           CHANGES IN VERSION 0.6-1 (October 2016)

## DEVELOPER LEVEL CHANGES

- Changes to configure.ac to address Solaris installation issues and removal of some lingering std::cout calls.

- Rearranged some parts of package and compileNimble compilation.  nimble.so no longer includes redundant compilation involving classes used only at the compileNimble stage.  A sessions-specific DLL is created the first time compileNimble is used.  It contains R function registration and a finalizer management system and links to CppCode/libnimble.[a|dylib|dll].

- Added a layer to manage finalizers registered with R for nimble C++ objects.  This allows nimble to manually finalize objects and then safely handle R's calls to finalizers.

- Added an internal clearCompiled system to manually finalize C++ objects and safely manage any external pointers to them in R.  This is not yet safe on Windows (can cause crash later upon exiting R) and is intended for developer use at this stage.

- Fixed bug in making an internal function accessible to Nimble generated code.

- Fixed bug in use of nimSwitch in getParam.
                           
#			   CHANGES IN VERSION 0.6-1

## USER LEVEL CHANGES

- Added getBound() functionality that provides dynamic access to lower and upper bounds on a node (from R or in the DSL), based on the underlying distribution and any user-defined truncation. Functionality is analogous to getParam(). Now used in the reflection sampler.

#                           CHANGES IN VERSION 0.6 (August 2016)

## USER LEVEL CHANGES

- Building models is faster and uses less memory for the building steps themselves.

- The new model implementation uses less memory and seems to have faster compiled execution in many cases.

- Compiling nimbleFunctions is generally faster.

- compileNimble can handle model operation functions that use a vector of node names with a scalar index, e.g., model$calculate(nodes[i]).  Index ranges are not supported.

- The computation time spent in each sampler of an MCMC can be obtained by mcmc$run(niter, time = TRUE) and then mcmc$getTimes(). 

- New RW_multinomial sampler incorporated into MCMC engine, for sampling nodes following a dmulti() distribution.

- Progress bar added to MCMC.

- New function runMCMC() available for easily running multiple MCMC chains.

- model checking modified so that only size/dimension and NA checking is done by default. Extended checking of ability to calculate all model nodes is now invoked via model$check() only at user request.

- MCMC and modelValues configuration syntax now uses 'conf' rather than 'spec' in its naming throughout NIMBLE. 

## DEVELOPER LEVEL CHANGES

- New C++ nimbleGraph system for some steps of processing graphs.

- New design of models, specifically nodeFunctions.  This is a major implementation change.

- New C++ memory allocation system using new double[] instead of vector<double>.  This is notably more efficient since it can avoid unnecessary initialization upon allocation of local variables in otherwise fast, heavily used nimbleFunctions.

- New installation system for Windows and for binary builds on OS X. For non-Linux we now by default build a static library, libnimble.a, that is linked into compiled code for models or nimbleFunctions. 

- We now export many fewer functions, though some non-user facing functions are still exported.

- Test suite now works in Windows.

- Reduced Windows compiler warnings.

- Various modifications to packaging to align with CRAN rules.

#                           CHANGES IN VERSION 0.5-1 (May 2016)

## USER LEVEL CHANGES

- Sequential Monte Carlo algorithms added: Bootstrap filter, Auxiliary particle filter, Ensemble Kalman filter, and Liu and West filter can be accessed by calls to buildBootstrapFilter(), buildAuxiliaryFilter, buildEnsembleKF(), and buildLiuWestFilter() respectively. buildPF() is no longer available, but that name will be re-used in a future version as a wrapper to some of the filter functions listed above.

- Particle MCMC (Particle Marginal Metropolis Hastings) samplers can be used as samplers within an MCMC.  'RW_PF' can be used to sample scalar parameters, and 'RW_PF_block' can be used for multivariate parameters or multiple parameters.

- An ascent-based version of the MCEM algorithm has been added to buildMCEM(), providing automated convergence assessment and stopping.

- A block sampler that uses user-specified likelihoods has been implemented and can be used as a sampler within an MCMC, by using 'RW_llFunction_block'.

- Added functions numeric(), integer(), matrix(), and array() to the NIMBLE DSL, for creating non-scalar variables in nimbleFunctions. These are intended to replace usage of declare(), with declare() being deprecated but still currently functional. setSize() is still available to change the size of an existing object prior to assigning into specific elements, but another call to one of the new functions can also be used for that purpose.

- Added multivariate-t distribution for use in BUGS code and the DSL. 

- Added binary sampler for Gibbs sampling of discrete 0/1 nodes. This is the default sampler used for Bernoulli and Binomial(p, size=1) nodes.

- New options for Metropolis-Hastings RW sampler: log=TRUE to sample on a log scale, and reflective=TRUE to reflect normal proposal distribution to stay within the range of the target distribution.

- Improved printing in printSamplers() and getSamplers() methods of MCMC configuration objects.

- Some steps of building and compiling models and algorithms are generally faster.

- forwardsolve(), backsolve() and solve() now work more flexibly than before and can be used in BUGS code.

- A variety of small bug fixes.

## DEVELOPER LEVEL CHANGES

- Testing has been added for SMC algorithms.

#                           CHANGES IN VERSION 0.5 (March 2016)

## USER LEVEL CHANGES

- Added support for solve(), forwardsolve(), and backsolve() to the NIMBLE DSL, with the same functionality as in R.

- Improved performance of conjugate multivariate-normal - multivariate-normal sampling.

- API changes for MCMCspec objects: new method printSamplers() prints the current samplers (replacing the old getSamplers() method); new functionality for getSamplers(), which now returns a list of samplerSpec objects; new overloaded functionality of setSamplers(), which can accept a list of samplerSpec objects, replacing the current set of samplers; new method getSamplerDefinition(), which returns the nimbleFunction definition of a sampler.

- API created for modifying samplerSpec objects: setName(), setSamplerFunction(), setTarget(), setControl().

- MCMC conjugate sampling efficiency improved by ~33%.

- new system for automating the comparison of different MCMCs and generating html results pages: the main function is compareMCMCs.

- enhanced checking of model to check for size and dimension mismatch in BUGS code.

- some steps of building models and determining dependencies (e.g. getDependencies method for model classes) are faster and use less memory.

- models have a new method, getDependenciesList, that provides a list of neighbor relationships in the model graph.

- There is a new function, getParam(model, node, parameter), that returns the value of a parameter of a stochastic node.  It works in R and in the NIMBLE DSL (compilable run code).

- Models now have member functions for calculate, simulate, getLogProb, calculateDiff, and getParam.  model$calculate(nodes) is equivalent to calculate(model, nodes), etc.

- The NIMBLE DSL now allows drop=FALSE as the last argument for indices.  e.g. x[1,,drop=FALSE].  This mimics R.

- added ability to use any of the methods for R's optim function in MCEM optimization step.

- There are more flexible ways to provide arguments to declare and setSize.

- The random walk block sampler is a bit more efficient.

- There is better error trapping during some compilation steps.

## DEVELOPER LEVEL CHANGES

- various modifications to conform to CRAN rules 

## BUG FIXES

- Compilation of BUGS function inprod is fixed

- some cases of handling multivariate variables of length one (e.g. a 1x1 matrix) have been fixed.

- added warning when 'size' in multinomial is not sum of values

- added missing terms (involving only the data values) in 'ddirchmulti' in test-user.R and in manual

#                        CHANGES IN VERSION 0.4-1 (Oct. 3, 2015)

## USER LEVEL CHANGES

- added support for OpenBUGS to MCMCsuite.

- added additional option for MCMCsuite: 'calculateEfficiency' to calculate ESS and ESS/time.

- added elliptical slice sampler 'ess' to MCMC engine.

- added new MCMC option for MCMCsuite that omits conjugate samplers: 'nimble_noConj'

## BUG FIXES

- fixed bug preventing use of nimbleFunctions in packages depending on NIMBLE.

- fixed bug preventing use of nimStop in R version of nimbleFunctions.

- reduced generation of C++ compiler warnings on Windows during compileNimble.

#                        CHANGES IN VERSION 0.4 (Aug. 2, 2015)

## USER LEVEL CHANGES

- almost everything is faster in almost all cases (building models, configuring MCMCs, compiling models and nimbleFunctions), sometimes very much faster, especially during R processing.

- added DSL functions stop("Error message") and checkInterrupt() to check for a user interrupt (via R_checkUserInterrupt() in C++).

- added support for scalar and vector character strings in the compiler.  The ony meaningful use currently is as argument(s) to print or stop.

- added calculateDiff as a fourth fundamental method (after calculate, simulate, getLogProb) in model node functions.

- nearly all calls to DSL functions or other nimbleFunctions handle R-style named or ordered arguments.

- dots (".") are allowed in nimbleFunction argument names.

- most Windows compiler warnings about comparing unsigned and signed ints should be gone.  Harmless warnings remain and some objects possibly being used before initialization.

- more informative error messages added in many cases.

- nimble version of some common R functions is now consistently prefixed with "nim" instead of "nimble".  E.g. nimRound, nimStop, nimPrint, nimCopy.  These can still used without the prefix in the DSL, e.g. round, stop, print, copy.

- added nimbleOptions useMultiInterfaceForNestedNimbleFunctions (default TRUE) and (experimental) clearNimbleFunctionsAfterCompiling (default FALSE) to reduce memory use. See User Manual section 9.6.

- setSize now consistently works in both compiled and uncompiled uses as setSize(X, size1, size2, etc.), not setSize(X, c(size1, size2, etc.)).

- added support for user-defined functions in BUGS code.

- added support for use of user-defined distributions in BUGS code, with the distribution density, simulate, distribution and quantile functions coded as nimbleFunctions and registered via registerDistributions().

- added support for truncation (T(,) and I(,) syntax in BUGS) plus dinterval distribution for censoring, both following the JAGS functionality.

- added ability to impose constraints via the dconstraint syntax; these constraints act like data and therefore only take effect a posteriori and are not imposed when simulating from the model.

- added alternative distribution names (aliases) usable in BUGS code (via an extensible system for adding additional aliases as desired): dbinom for dbin, dmultinom for dmulti, dnbinom for dnegbin, dchisqr for dchisq, dweibull for dweib, ddirich for ddirch, dwishart for dwish.

- added additional parameterizations for dlnorm and dt.

- improved handling of distribution functions in nimbleFunction run code (i.e., the DSL), allowing argument matching by name and use of default values. Distribution functions can also be used as deterministic functions (e.g., pnorm to do probit calculations) in BUGS code. 

- added an optional check when building model that alerts user to presence of nodes without values and log probability calculations that return NA.

- model$checkConjugacy(nodes) now accepts a character vector argument of node/variable names, and returns a named list identifying conjugate relationships.

- incorporated automated blocking into MCMC engine: configureMCMC(model, autoBLock=TRUE) returns an MCMCspec defined by automated blocking procedure, and buildMCMC(model, autoBlock=TRUE) returns the corresponding MCMC algorithm.

- new syntax for MCMCspec$addSampler(type, target, control).  mandatory 'type' argument specifies sampler type, which may be a character string or nimbleFunction object.  mandatory 'target' argument specifices a character vector of target model nodes.  optional 'control' list argument overrides default sampler control parameters.

- default control parameters for MCMC sampling algorithms are now a NIMBLE system level option: MCMCcontrolDefaultList.

- MCMCspec$getSamplers(ind), MCMCspec$setSamplers(ind), and MCMCspec$removeSamplers(ind) can optionally also accept a character vector argument, and will act on all samplers which sample the nodes specified.

- cleaned up formatting of output when querying and setting MCMC samplers.

- NIMBLE will no longer "hang" when used in RStudio.

- reordered argument names to nimbleModel to put in more logical order of importance. first four arguments are: code, constants, data, inits.

- added documentation of NIMBLE's built-in MCMC samplers via help(samplers).

- added more informative error messages for common errors encountered when calling nimbleModel() and compileNimble().


## DEVELOPER LEVEL CHANGES

- re-wrote model definition processing and R nodeFunction instantiation to do more work at the level of a BUGS declaration rather than for each individual node.

- re-wrote checkConjugacy to use the new model definition content.

- re-wrote C++ nimCopy system and corresponding R partial evaluation to do more processing in R and set up C++ for faster processing.

- added test-copy for nimCopy, values() and values()<- .

- added nimbleProjectClass method for adding groups of nimbleFunctions from same generator at once.

- modified core samplers like RW to use the new calculateDiff nodeFunction method.

- added CmultiNimbleFunctionInterface to interface compiled nimbleFunctions contained within other nimbleFunctions.

- removed some old code deemed to be fully defunct.

- added 'range' variable to all distributions; when truncation is specified this modifies the range and also sets the 'truncated' flag to TRUE; range is not at the moment used but could be used by developers in algorithms.

- added a variety of error checking for various distribution functions and now generally use doubles in distribution functions in C++.

- modified API for nimble options to mimic R's options(), including nimbleOptions() and getOptions() and renaming of the underlying nimble options object to be .nimbleOptions.

- added an environment, nimbleUserNamespace, in the package namespace, that allows user-supplied information to control the behavior of NIMBLE. At the moment, only used for user-supplied distributions list.

- testing of user-defined distributions and functions in BUGS code (tests/test-user.R).

- full suite of testing for truncation/censoring/constraints added (tests/test-trunc.R).

- names of custom MCMC sampler nimbleFunctions need not begin with "sampler_".

- added initializeModel() nimbleFunction for use in algorithms; it performs sensible model initialization at the onset of an algorithm.

## BUG FIXES

- setSize fixed.

- generation of "1/(scale)" type bits in keyword processing for exp and gamma fixed.

- fixed bug where quitting R after dyn.unload()ing in OS X / Linux caused a segfault.  Now we rely on a finalizer that uses the base class virtual destructor, which is in the package dll, not the generated one(s).  We still have segfaults on Windows after dyn.unload()ing.

- made row subset assignment work in the compiler, e.g. x[i,] <- foo(a).

- updated package to use new igraph 1.x.x API and now explicitly link in LAPACK_LIBS and BLAS_LIBS when building libnimble.so as dpotrf not being found when using igraph 1.x.x.

- fixed bug where dexp was passing the wrong parameter (rate instead of scale) to C++; we now use our own dexp_nimble, which calls Rmath's dexp correctly.

- set default for deparse to width.cutoff=500L to avoid splitting of lines when going from expressions to string names.

#                        CHANGES IN VERSION 0.3-1 (Mar. 7, 2015)

## USER LEVEL CHANGES

- fully allow data to be provided as part of constants argument in nimbleModel for compatibility with JAGS and BUGS, which mix data and constants

- added Dirichlet-multinomial conjugacy handling

- Added 'oldSpec' argument to configureMCMC. Allows rebuilding of MCMCs without the need to recompile code if there are no new types of samplers.

- Now allow for up to 4D arrays

- Added handling of raising a vector of values to a scalar power in nimble functions (i.e., the NIMBLE DSL)

## DEVELOPER LEVEL CHANGES

- fixed naming issue in checking posterior summaries against known values in test_mcmc 

## BUG FIXES 

- fixed bug that allowed inits to overwrite data values in variables that are mixtures of data and parameters. Also now detect elements of initial values that are not variables in the model.

- fixed bug in naming of nodes in as.matrix() when there is a variable with a single indexed node (e.g. x[1], but no other nodes in the variable 'x')

- fixed bug regarding logProbs with gaps (e.g., if logProb_x[1] and logProb_x[3] are defined but not logProb_x[2])

- fixed bug in handling of unary minus for non-scalars

#		  	CHANGES IN VERSION 0.3 (Dec. 31, 2014)

## USER LEVEL CHANGES

- IMPORTANT SYNTAX CHANGE: nimbleFunctions are run via myNimbleFunction$run() instead of myNimbleFunction(). This means code written in v0.2 and earlier will not run without adding $run. 

- IMPORTANT SYNTAX CHANGES: 
  -- writing code for nimble models is now done with nimbleCode (previously modelCode)
  -- to customize MCMC use myMCMCspec <- configureMCMC(myModel) (previously MCMCspec)
  -- to build an MCMC algorithm, either myGenericMCMC <- buildMCMC(myModel) for a generic build or myCustomizedMCMC <- buildMCMC(myMCMCspec) for a customized MCMC algorithm can be used (previously one always had to build an MCMCspec, even for generic build)

- Variables and methods of a nimbleFunction can be accessed by myNimbleFunction$myVariable or myNimbleFunction$myMethod() rather than nfVar(myNimbleFunction, ‘myVariable’) or nfMethod(myNimbleFunction, ‘myMethod’)(). Similarly, names of objects/methods of a nimble function can be queried by ls(myNimbleFunction)

- faster compiling, most notably when compiling MCMC algorithms
 
## DEVELOPER LEVEL CHANGES

- Addition of graphIDs for logProbs for models (separate set than graphIDs for nodes) and graphIDs for modelValue variables

- Addition of $expandNodeNames() for modelValues

- $expandNodeNames() now evals nodeNames text array to quickly retrieve graphIDs.  expandNodeNames can also accept graphIDs and return names. Also can return multivariate nodes as either ‘x[1:10]’ or ‘x[1]’,’x[2]’, etc.

- Addition of numberedObjects (both at R and C++ level) currently used for fast construction of C++ nodeFunctionVectors, modelVariableAccessors and modelValuesAccessors

- Reduction of nimbleFunctions required for initialization functions of MCMCs

- nimCopy now only requires that the number of nodes are equal, not that the sets of consecutive nodes be equal (i.e. previous nimCopy(…, nodes = c(‘x[1]’, ‘x[2]’), nodesTo = c(‘y[1,1]’, ‘y[2,2]’)) was not allowed). 

- Addition of keywordProcessing system for more organized compiling of nimbleFunctions

Bug fixes

- Fixed bug regarding calculating log probabilities for multinomial distribution with probabilities equal to 0


#                     CHANGES IN VERSION 0.2 (Oct. 12, 2014)

- Internal changes to decrease time to build models and nimbleFunctions.

- Multivariate conjugate updaters are now included in the default MCMC implementation, as well as block updating on scalar and/or multivariate nodes.

- Vector treatment of is.na and is.nan, so is.na(values(model, nodes)) will work, equivalently to any(is.na(...)) in R.

- More extensive suite of tests.

- Extension of node maps, which can be accessed via graphIDs, to increase speed of determining node types, checking if nodes in model, etc. 

- model$getNodeNames, expandNodeNames and getDependencies all return the full node function name, i.e. if 'x[1:2]' is multivariate, will return 'x[1:2]' rather than 'x[1]', 'x[2]', unless returnScalarComponents = TRUE

- Added class of nodeVectors, which save the graphIDs of a set of nodes. Use of these is minimal at the moment. 

## BUG FIXES

- Fixed issue with specifying multivariate nodes as rows/columns of matrices

- Fixed bug in the conjugacy system

- Fixed bug in compiled version of adaptive block updater

- Fixed issue with building the NIMBLE package for Windows, as well as cleaning up causes of compiler warnings.

- Fixed issue that values(model, nodes) sorted the outcome. Remaining bug to be fixed: nimCopy(nodesFrom, nodesTo) will NOT work if nodesFrom and nodesTo have different number/different lengths of ordered contingent blocks of memory, i.e. nimCopy(nodesFrom = c('x[1]', 'x[2]'), nodesTo = c('y[2]', 'y[1]')) will fail because nodesFrom is one block but nodesTo is two blocks

- Fixed issue that chol(A+B) was being translated to chol(A) + B


#                     CHANGES IN VERSION 0.1-1 (Aug. 28, 2014)


## BUG FIXES

- handling of spaces in installed directory name on Windows

- ensured that various Makevars and Makevars.in and Makevars.win files are
  provided in the package for proper installation on all platforms

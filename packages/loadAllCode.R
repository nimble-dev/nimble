library(igraph)
mysource <- function(p) source(file.path('nimble/R', p))
mysource('all_utils.R')
mysource('options.R')
nimbleOptions$notUsingPackage <- TRUE

mysource('distributions_inputList.R')
mysource('distributions_processInputList.R')
mysource('distributions_implementations.R')
mysource('BUGS_BUGSdecl.R')
mysource('BUGS_nodeInfo.R')
mysource('BUGS_contexts.R')
mysource('BUGS_modelDef.R')
mysource('BUGS_model.R')
mysource('BUGS_graphNodeMaps.R')
mysource('BUGS_readBUGS.R')
mysource('BUGS_testBUGS.R')
mysource('BUGS_getDependencies.R')
mysource('BUGS_utils.R')
mysource('BUGS_mathCompatibility.R')

mysource('genCpp_exprClass.R')
mysource('genCpp_operatorLists.R')
mysource('genCpp_RparseTree2exprClasses.R')
mysource('genCpp_initSizes.R')
mysource('genCpp_buildIntermediates.R')
mysource('genCpp_processSpecificCalls.R')
mysource('genCpp_sizeProcessing.R')
mysource('genCpp_insertAssertions.R')
mysource('genCpp_maps.R')
mysource('genCpp_liftMaps.R')
mysource('genCpp_eigenization.R')
mysource('genCpp_addDebugMarks.R')
mysource('genCpp_generateCpp.R')

mysource('RCfunction_core.R')
mysource('RCfunction_compile.R')

mysource('nimbleFunction_util.R')
mysource('nimbleFunction_core.R')
mysource('nimbleFunction_nodeFunction.R')
mysource('nimbleFunction_Rexecution.R')
mysource('nimbleFunction_compile.R')

mysource('types_util.R')
mysource('types_symbolTable.R')
mysource('types_modelValues.R')
mysource('types_modelValuesAccessor.R')
mysource('types_modelVariableAccessor.R')
mysource('types_nimbleFunctionList.R')
mysource('types_nodeFxnVector.R')
mysource('types_numericLists.R')

mysource('cppDefs_utils.R')
mysource('cppDefs_variables.R')
mysource('cppDefs_core.R')
mysource('cppDefs_namedObjects.R')
mysource('cppDefs_BUGSmodel.R')
mysource('cppDefs_RCfunction.R')
mysource('cppDefs_nimbleFunction.R')
mysource('cppDefs_modelValues.R')
mysource('cppDefs_cppProject.R')
mysource('cppDefs_outputCppFromRparseTree.R')

mysource('cppInterfaces_utils.R')
mysource('cppInterfaces_models.R')
mysource('cppInterfaces_modelValues.R')
mysource('cppInterfaces_nimbleFunctions.R')
mysource('cppInterfaces_otherTypes.R')

mysource('nimbleProject.R')

mysource('MCEM_build.R')


mysource('MCMC_utils.R')
mysource('MCMC_spec.R')
mysource('MCMC_build.R')
mysource('MCMC_samplers.R')
mysource('MCMC_conjugacy.R')
mysource('MCMC_suite.R')
mysource('NF_utils.R')


mysource('makevars.R')
.NimbleUseRegistration = FALSE
mysource('registration.R')
mysource('zzz.R')

NeedMakevarsFile = TRUE

IncludeCodeDir <- "nimble/inst/include/nimble"
NimbleCodeDir <- "nimble/inst/CppCode"
options(nimble.Makevars.file = if(.Platform$OS.type == "windows" && file.exists("Makevars.win")) "Makevars.win" else "Makevars")

if(Sys.getenv("NIMBLE_PKG_SRC_DIR") == "") {
    path = normalizePath("nimble/inst/CppCode")
    if(.Platform$OS.type == "windows") { # check for cygwin???
         # gsub("C:", "/cygdrive/c", path)  ?
       path =  gsub("\\\\", "/", shortPathName(path))
         # You need to adjust MyMakevars to have the full path
         # to the local directory .../nimble/packages
         # Copy MyMakevars_template to MyMakevars and then edit
         # You also need to make a copy of the directory Eigen_local to Eigen
         # in nimble/packages/nimble/
         # Do not put Eigen or MyMakevars under version control
       options(nimble.Makevars.file = "MyMakevars")
    }
    Sys.setenv("NIMBLE_PKG_SRC_DIR" = path)
    NimbleCodeDir = path
}

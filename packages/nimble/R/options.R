

#' options used for NIMBLE package
#'
#' @details These options are for development use at this point.
nimbleOptions <- list(
    convertSingleVectorsToScalarsInSetupArgs = TRUE,
    messagesWhenBuildingOrFinalizingCppObjects = FALSE,
    indexDrop = TRUE,
    debugRCfunProcessing = FALSE,
    debugNFProcessing = FALSE,
    debugCppLineByLine = FALSE,
    compileAltParamFunctions = TRUE,
    verifyConjugatePosteriors = TRUE,        ## verifies the correct posterior is created for any conjugate samplers, at run-time
    includeCPPdists = TRUE,    ## includes dists.cpp and nimDists.cpp in the compilation.  Momentarily we have a problem on Windows.
    processBackwardsModelIndexRanges = FALSE,    ## if FALSE (default), for(i in 9:7) in model code becomes for(i in 7).  if TRUE, becomes for(i in c(9, 8, 7))
    prioritizeColonLikeBUGS = TRUE ## if FALSE, 1:2 + 1 evaluates to 2:3, consistent with R.  If TRUE, it evalutes to 1:3, consistent with BUGS 
)

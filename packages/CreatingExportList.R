rm(list = ls())

individualExportClassNames = c( 
 'RModelBaseClass',
 'singleModelValuesAccessClass',
 'singleVarAccessClass', 
 'CmodelBaseClass',
 'CnimbleFunctionBase', 
 'modelBaseClass',
 'modelValuesBaseClass', 
 'codeBlockClass',
 'nimbleFunctionBase')

individualExportNames = c(
 'checkInterrupt',
 'codeProcessIfThenElse',
 'parseEvalNumeric',
 'makeMV_GID_Map',
 'nimbleCode', 
 'BUGScode', 
 'nimbleModel', 
 'readBUGSmodel', 
 'testBUGSmodel',
 'populateNodeFxnVec',
 'simulate',
 'calculate',
 'calculateDiff',
 'getLogProb',
 'nfVar',
 'BUGScontextClass',
 'BUGSdeclClass',
 'BUGSsingleContextClass',
 'cppBUGSmodelClass',
 'cppClassDef',
 'dependentClass',
 'makeBUGSclassFields',
 'MakeCustomModelClass',
 'makeCustomModelValuesClass',
 'modelDefClass',
 'modelDefInfoClass',
 'nimbleFunction',
 'nimbleFunctionVirtual',
 'nimCopy',
  'nimPrint',
  'nimStop',
 '`nfVar<-`',
 'OptimReadyFunction',
 'getValues', 
 'setValues',
 'values',
 '`values<-`',
 'rankSample',
 'setSize',
    'nfMethod',
    'declare',
    'asRow',
    'getBUGSexampleDir',
    'asCol',
    'getNimbleOption',
    'nimbleOptions',
    'registerDistributions',
    'deregisterDistributions',
    'getDistribution',
    'calc_dmnormAltParams',
    'getDependencyPaths'
)

individualMaskedFunctions = c('mysource', 'individualMaskedFunctions', 'maskFileVector', 'AllFiles', 'maskSource', 'maskedFuns', 'individualExportNames', 'nfVar<-', 'getModelValuesMemberElement', 'newModelValues', 'testRows', 'individualExportClassNames',

'C_dcat', 
'C_ddirch', 
'C_dmnorm_chol',
'C_dmulti',
'C_dt_nonstandard',
'C_dwish_chol',
'C_dexp_nimble',
'C_rcat',
'C_rdirch',
'C_rmnorm_chol',
'C_rmulti',
'C_rt_nonstandard',
'C_rwish_chol',
'C_rexp_nimble',
'C_pt_nonstandard',
'C_pexp_nimble',
'C_qt_nonstandard',
'C_qexp_nimble',

#'node_determ',
'node_stoch_dbern',
'node_stoch_dbeta',
'node_stoch_dbin',
'node_stoch_dcat',
'node_stoch_dchisq',
'node_stoch_dconstraint',
'node_stoch_ddirch',
'node_stoch_dexp',
'node_stoch_dgamma',
'node_stoch_dinterval',
'node_stoch_dlnorm',
'node_stoch_dlogis',
'node_stoch_dmnorm',
'node_stoch_dmulti',
'node_stoch_dnegbin',
'node_stoch_dnorm',
'node_stoch_dpois',
'node_stoch_dt',
'node_stoch_dunif',
'node_stoch_dweib',
'node_stoch_dwish',


'makeBUGSactiveBindingDef',
'makeBUGSclassFields',
'makeCSingleModelValuesAccessor',
'makeCSingleVariableAccessor',
'MakeCustomModelClass',
'makeCustomModelValuesClass',
'makeEnvName',
'makeLogProbEnvName',
'makeLogProbName',
'makeMapAccessExpr',
'makeMapSetupCodeExprs',
'makeMapSetupCodeNames',
'makeModelBindingFields',
'makeModelCppCopyTypes',
'makeModelValuesClassFields',
'makeNameFromMatchedCall',
'makeNameName',
'makeNFBindingFields',
'makeNimbleFxnCppCopyTypes',
'makeNimbleFxnInterfaceCallMethodCode',
'makeNumericList',
'makeRowName',
'makeSingleIndexAccessCodeNames',
'makeSingleIndexAccessExpr',
'makeSingleIndexAccessorSetupCode',
'makeVecName',

'ndf_addModelDollarSignsToMethods',
'ndf_createContains',
'ndf_createMethodList',
'ndf_createSetupFunction',
'ndf_createSingleMethod',
'ndf_createStochCalculate',
'ndf_createStochSimulate',
'ndf_createVirtualNodeFunctionDefinition',
'ndf_createVirtualNodeFunctionDefinitionsList',
'ndf_generateGetParamFunction',

'replaceAllCodeSuccessfully',
'replaceWhatWeCan',
'RHSonlyInit',
'RHSonlyInit_virtual',
'activeBindingTemplate',
'addModelDollarSign',

'cc_checkLinearity',
'cc_combineExprsAddition',
'cc_combineExprsDivision',
'cc_combineExprsMultiplication',
'cc_combineExprsSubtraction',
'cc_expandDetermNodesInExpr',
'cc_getNodeDistributionText',
'cc_getNodeParamExpr',
'cc_getNodesInExpr',
'cc_getNodeValueExpr',
'cc_linkCheck',
'cc_makeConjugateSamplerName',
'cc_makeDDistributionName',
'cc_makeRDistributionName',
'cc_makeSamplerTypeName',
'cc_negateExpr',
'cc_nodeInExpr',
'cc_otherParamsCheck',

'createMakevars',
'expandMVNames',
'exprAsListDrop2',
'file',
'hasBracket',
'mcmcNodeInit',
'mcmcNodeInit_virtual',
'mcmc_findControlListNamesInCode',
'mcmc_listContentsToStr',
'modelValues2Matrix',
'modelValuesElement2Matrix',
'nameMashupFromExpr',
'NeedMakevarsFile',
'parseTreeSubstitute',
'projectNameCreator',
'UseLibraryMakevars',
'replaceDistributionAliasesNameOnly'
)

maskFileVector = c(
					'options.R',
					'distributions_inputList.R',
					 'distributions_processInputList.R',
		#			 'distributions_implementations.R',

		#			 'BUGS_BUGSdecl.R',
		#			 'BUGS_nodeInfo.R',
		#			 'BUGS_contexts.R',
					 'BUGS_modelDef.R',
		#			 'BUGS_model.R',
					 'BUGS_graphNodeMaps.R',
					 'BUGS_readBUGS.R',
					 'BUGS_testBUGS.R',
		#			 'BUGS_getDependencies.R',		Needed
		#			 'BUGS_utils.R',				Needed
		
		#			 'BUGS_mathCompatibility.R',




					'genCpp_exprClass.R',
					'genCpp_operatorLists.R',
					'genCpp_RparseTree2exprClasses.R',
					'genCpp_initSizes.R',
					'genCpp_buildIntermediates.R',
					'genCpp_processSpecificCalls.R',
					'genCpp_sizeProcessing.R',
					'genCpp_insertAssertions.R',
					'genCpp_maps.R',
					'genCpp_liftMaps.R',
					'genCpp_eigenization.R',
					'genCpp_addDebugMarks.R',
					'genCpp_generateCpp.R',
					
					
					'RCfunction_core.R',
					'RCfunction_compile.R',

		
#					 'nimbleFunction_util.R',				Needed
					 'nimbleFunction_core.R',
# 					 'nimbleFunction_nodeFunction.R',		Needed
					 'nimbleFunction_Rexecution.R',			
#					 'nimbleFunction_compile.R',			Needed

					 'types_symbolTable.R',


					 'types_numericLists.R',


#					'cppInterfaces_otherTypes.R',			Needed

					
					'cppDefs_utils.R',
					'cppDefs_variables.R',
					'cppDefs_core.R',
 					'cppDefs_namedObjects.R',
 					'cppDefs_BUGSmodel.R',
 					'cppDefs_RCfunction.R',
 					'cppDefs_nimbleFunction.R',
 					'cppDefs_modelValues.R',
 					'cppDefs_cppProject.R',
 					'cppDefs_outputCppFromRparseTree.R'
)


##	NOTE: ALL FILES VECTOR MUST BE IN CORRECT LOADING ORDER!!!!!!
  AllFiles <- readLines(file.path('nimble', 'DESCRIPTION'))
  collateStart <- grep("^Collate:", AllFiles)
  AllFiles <- AllFiles[collateStart:length(AllFiles)]
  AllFiles[1] <- gsub("Collate:", "", AllFiles[1])
  AllFiles <- gsub("\\s+", "", AllFiles)
  rm(collateStart)

maskSource <- function(fileName, maskFileVector, individualMaskedFunctions)
	{
	toMask <- fileName %in% maskFileVector
#	if(toMask)
	lsBefore <- ls(envir = .GlobalEnv)
	source(file.path('nimble', 'R', fileName))
	if(!toMask){
		lsAfter = ls(envir = .GlobalEnv)
		cat('Potential gain from masking', fileName, 'is ', length(setdiff(lsAfter, union(lsBefore, individualMaskedFunctions) ) ), '\n' ) 
		return(character(0) )
		}
	lsAfter <- ls(envir = .GlobalEnv)
	masked = setdiff(lsAfter, lsBefore)
 	return(masked)
	}

maskedFuns = character(0)
for(file in AllFiles)
	maskedFuns = c(maskedFuns, maskSource(file, maskFileVector, individualMaskedFunctions))
	
allNames <- ls()

removedFuns <- union(maskedFuns, individualMaskedFunctions)
exportNames = setdiff(allNames,  removedFuns )
exportNames = setdiff(exportNames, individualExportClassNames)
exportNames = union(exportNames, individualExportNames)

 exportClasses = individualExportClassNames

exportNames = sort(exportNames)
exportClasses = sort(exportClasses)

exportText = c(paste('export(', exportNames, ')', sep = ''), paste('exportClass(', exportClasses, ')', sep = ''))

exportText = paste(exportText, collapse = '\n')
# cat(exportText, sep = '\n')

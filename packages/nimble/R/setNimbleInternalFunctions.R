nimbleInternalFunctions <- new.env()

internalFuns <- c("CmodelValues","buildNeededObjects","checkNimbleFunctionListCpp","copyFromRobjectViaActiveBindings","dimOrLength","getBoolValue","getCharacterValue","getCharacterVectorValue","getDoubleValue","getIntValue","getMVName","getMVptr","getNimValues","getVarAndIndices","newObjElementPtr","setBoolValue","setCharacterValue","setCharacterVectorValue","setDoublePtrFromSinglePtr","setSmartPtrFromSinglePtr","setSmartPtrFromDoublePtr","setDoubleValue","setIntValue","setNimValues","setOnePtrVectorOfPtrs","setPtrVectorOfPtrs", "nimbleFinalize", "clearNeededObjects","getSetCharacterScalar","getSetCharacterVector","getSetDoubleScalar","getSetIntegerScalar","getSetLogicalScalar","getSetNimbleList","getSetNumericVector", "makeNewNimListSEXPRESSIONFromC")

for(fun in internalFuns) 
    assign(fun, get(fun), envir = nimbleInternalFunctions)

# remove these functions from the nimble namespace too?

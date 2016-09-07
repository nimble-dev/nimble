nimbleInternalFunctions <- new.env()

internalFuns <- c("CmodelValues","buildNeededObjects","checkNimbleFunctionListCpp","copyFromRobjectViaActiveBindings","dimOrLength","getBoolValue","getCharacterValue","getCharacterVectorValue","getDoubleValue","getIntValue","getMVName","getMVptr","getNimValues","getVarAndIndices","newObjElementPtr","setBoolValue","setCharacterValue","setCharacterVectorValue","setDoublePtrFromSinglePtr","setDoubleValue","setIntValue","setNimValues","setOnePtrVectorOfPtrs","setPtrVectorOfPtrs")

for(fun in internalFuns) 
    assign(fun, get(fun), envir = nimbleInternalFunctions)

# remove these functions from the nimble namespace too?



initializeModel <- nimbleFunction(
    
    setup = function(model) {
        RHSonlyNodes <- model$getMaps('nodeNamesRHSonly')
        hasRHSonlyNodes <- length(RHSonlyNodes) > 0

        allDetermNodes <- model$getNodeNames(determOnly = TRUE)
        determNodesNodeFxnVector <- nodeFunctionVector(model = model, nodeNames = allDetermNodes)
        
        stochNonDataNodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
        
        initFunctionList <- nimbleFunctionList(mcmcNodeInit_virtual)

        iter <- 1
        if(hasRHSonlyNodes) {
            initFunctionList[[iter]] <- mcmcCheckRHS_Init(model = model, node = RHSonlyNodes)
            iter <- iter + 1
        }
        
        for(i in seq_along(stochNonDataNodes))
            initFunctionList[[iter + i - 1]] <- mcmcStochNode_Init(model, stochNonDataNodes[i])
    },
    
    run = function() {
        
        for(i in seq_along(initFunctionList))  {   
            calculate(nodeFxnVector = determNodesNodeFxnVector)
            initFunctionList[[i]]$run()
        }
        calculate(model)

    },  where = getLoadingNamespace()
)



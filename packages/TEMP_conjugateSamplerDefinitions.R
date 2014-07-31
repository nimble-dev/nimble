
sampler_conjugate_dbeta <- nimbleFunction(contains = sampler_BASE, setup = function(model, mvSaved, control) {
    targetNode <- control$targetNode
    calcNodes <- model$getDependencies(targetNode)
    calcNodesDeterm <- model$getDependencies(targetNode, determOnly = TRUE)
    targetNode_nodeFunctionList <- nimbleFunctionList(node_stoch_dbeta)
    targetNode_nodeFunctionList[[1]] <- model$nodeFunctions[[targetNode]]
    dependents_dbern_nodeNames <- control$dependents_dbern
    dependents_dbern_nodeFunctions <- nimbleFunctionList(node_stoch_dbern)
    for (i in seq_along(dependents_dbern_nodeNames)) {
        dependents_dbern_nodeFunctions[[i]] <- model$nodeFunctions[[dependents_dbern_nodeNames[i]]]
    }
    dependents_dbin_nodeNames <- control$dependents_dbin
    dependents_dbin_nodeFunctions <- nimbleFunctionList(node_stoch_dbin)
    for (i in seq_along(dependents_dbin_nodeNames)) {
        dependents_dbin_nodeFunctions[[i]] <- model$nodeFunctions[[dependents_dbin_nodeNames[i]]]
    }
    dependents_dnegbin_nodeNames <- control$dependents_dnegbin
    dependents_dnegbin_nodeFunctions <- nimbleFunctionList(node_stoch_dnegbin)
    for (i in seq_along(dependents_dnegbin_nodeNames)) {
        dependents_dnegbin_nodeFunctions[[i]] <- model$nodeFunctions[[dependents_dnegbin_nodeNames[i]]]
    }
}, run = function() {
    modelLogProb0 <- getLogProb(model, calcNodes)
    origValue <- model[[targetNode]]
    prior_shape1 <- nfMethod(targetNode_nodeFunctionList[[1]], 'get_shape1')()
    prior_shape2 <- nfMethod(targetNode_nodeFunctionList[[1]], 'get_shape2')()
    declare(dependents_dbern_values, double(1, length(dependents_dbern_nodeFunctions)))
    declare(dependents_dbern_prob, double(1, length(dependents_dbern_nodeFunctions)))
    for (i in seq_along(dependents_dbern_nodeFunctions)) {
        dependents_dbern_values[i] <- mfMethod(dependents_dbern_nodeFunctions[[i]], 'get_value')()
        dependents_dbern_prob[i] <- nfMethod(dependents_dbern_nodeFunctions[[i]], 'get_prob')()
    }
    declare(dependents_dbin_values, double(1, length(dependents_dbin_nodeFunctions)))
    declare(dependents_dbin_prob, double(1, length(dependents_dbin_nodeFunctions)))
    declare(dependents_dbin_size, double(1, length(dependents_dbin_nodeFunctions)))
    for (i in seq_along(dependents_dbin_nodeFunctions)) {
        dependents_dbin_values[i] <- mfMethod(dependents_dbin_nodeFunctions[[i]], 'get_value')()
        dependents_dbin_prob[i] <- nfMethod(dependents_dbin_nodeFunctions[[i]], 'get_prob')()
        dependents_dbin_size[i] <- nfMethod(dependents_dbin_nodeFunctions[[i]], 'get_size')()
    }
    declare(dependents_dnegbin_values, double(1, length(dependents_dnegbin_nodeFunctions)))
    declare(dependents_dnegbin_prob, double(1, length(dependents_dnegbin_nodeFunctions)))
    declare(dependents_dnegbin_size, double(1, length(dependents_dnegbin_nodeFunctions)))
    for (i in seq_along(dependents_dnegbin_nodeFunctions)) {
        dependents_dnegbin_values[i] <- mfMethod(dependents_dnegbin_nodeFunctions[[i]], 'get_value')()
        dependents_dnegbin_prob[i] <- nfMethod(dependents_dnegbin_nodeFunctions[[i]], 'get_prob')()
        dependents_dnegbin_size[i] <- nfMethod(dependents_dnegbin_nodeFunctions[[i]], 'get_size')()
    }
    counter <- 1
    for (i in seq_along(dependents_dbern_nodeFunctions)) {
        if (counter == 1) {
            contribution_shape1 <- dependents_dbern_values[i]
            contribution_shape2 <- 1 - dependents_dbern_values[i]
            counter <- counter + 1
        }
        else {
            contribution_shape1 <- contribution_shape1 + dependents_dbern_values[i]
            contribution_shape2 <- contribution_shape2 + (1 - dependents_dbern_values[i])
        }
    }
    for (i in seq_along(dependents_dbin_nodeFunctions)) {
        if (counter == 1) {
            contribution_shape1 <- dependents_dbin_values[i]
            contribution_shape2 <- dependents_dbin_size[i] - dependents_dbin_values[i]
            counter <- counter + 1
        }
        else {
            contribution_shape1 <- contribution_shape1 + dependents_dbin_values[i]
            contribution_shape2 <- contribution_shape2 + (dependents_dbin_size[i] - dependents_dbin_values[i])
        }
    }
    for (i in seq_along(dependents_dnegbin_nodeFunctions)) {
        if (counter == 1) {
            contribution_shape1 <- dependents_dnegbin_size[i]
            contribution_shape2 <- dependents_dnegbin_values[i]
            counter <- counter + 1
        }
        else {
            contribution_shape1 <- contribution_shape1 + dependents_dnegbin_size[i]
            contribution_shape2 <- contribution_shape2 + dependents_dnegbin_values[i]
        }
    }
    newValue <- rbeta(1, shape1 = prior_shape1 + contribution_shape1, shape2 = prior_shape2 + contribution_shape2)
    model[[targetNode]] <<- newValue
    calculate(model, calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    modelLogProb1 <- getLogProb(model, calcNodes)
    posteriorLogDensity0 <- dbeta(origValue, shape1 = prior_shape1 + contribution_shape1, shape2 = prior_shape2 + contribution_shape2, log = 1)
    posteriorLogDensity1 <- dbeta(newValue, shape1 = prior_shape1 + contribution_shape1, shape2 = prior_shape2 + contribution_shape2, log = 1)
    posteriorVerification <- modelLogProb0 - posteriorLogDensity0 - modelLogProb1 + posteriorLogDensity1
    if (abs(posteriorVerification) > 1e-10) {
        nimPrint('conjugate posterior density appears to be wrong')
    }
}, methods = list(getPosteriorLogDensity = function() {
    prior_shape1 <- nfMethod(targetNode_nodeFunctionList[[1]], 'get_shape1')()
    prior_shape2 <- nfMethod(targetNode_nodeFunctionList[[1]], 'get_shape2')()
    declare(dependents_dbern_values, double(1, length(dependents_dbern_nodeFunctions)))
    declare(dependents_dbern_prob, double(1, length(dependents_dbern_nodeFunctions)))
    for (i in seq_along(dependents_dbern_nodeFunctions)) {
        dependents_dbern_values[i] <- mfMethod(dependents_dbern_nodeFunctions[[i]], 'get_value')()
        dependents_dbern_prob[i] <- nfMethod(dependents_dbern_nodeFunctions[[i]], 'get_prob')()
    }
    declare(dependents_dbin_values, double(1, length(dependents_dbin_nodeFunctions)))
    declare(dependents_dbin_prob, double(1, length(dependents_dbin_nodeFunctions)))
    declare(dependents_dbin_size, double(1, length(dependents_dbin_nodeFunctions)))
    for (i in seq_along(dependents_dbin_nodeFunctions)) {
        dependents_dbin_values[i] <- mfMethod(dependents_dbin_nodeFunctions[[i]], 'get_value')()
        dependents_dbin_prob[i] <- nfMethod(dependents_dbin_nodeFunctions[[i]], 'get_prob')()
        dependents_dbin_size[i] <- nfMethod(dependents_dbin_nodeFunctions[[i]], 'get_size')()
    }
    declare(dependents_dnegbin_values, double(1, length(dependents_dnegbin_nodeFunctions)))
    declare(dependents_dnegbin_prob, double(1, length(dependents_dnegbin_nodeFunctions)))
    declare(dependents_dnegbin_size, double(1, length(dependents_dnegbin_nodeFunctions)))
    for (i in seq_along(dependents_dnegbin_nodeFunctions)) {
        dependents_dnegbin_values[i] <- mfMethod(dependents_dnegbin_nodeFunctions[[i]], 'get_value')()
        dependents_dnegbin_prob[i] <- nfMethod(dependents_dnegbin_nodeFunctions[[i]], 'get_prob')()
        dependents_dnegbin_size[i] <- nfMethod(dependents_dnegbin_nodeFunctions[[i]], 'get_size')()
    }
    counter <- 1
    for (i in seq_along(dependents_dbern_nodeFunctions)) {
        if (counter == 1) {
            contribution_shape1 <- dependents_dbern_values[i]
            contribution_shape2 <- 1 - dependents_dbern_values[i]
            counter <- counter + 1
        } else {
            contribution_shape1 <- contribution_shape1 + dependents_dbern_values[i]
            contribution_shape2 <- contribution_shape2 + (1 - dependents_dbern_values[i])
        }
    }
    for (i in seq_along(dependents_dbin_nodeFunctions)) {
        if (counter == 1) {
            contribution_shape1 <- dependents_dbin_values[i]
            contribution_shape2 <- dependents_dbin_size[i] - dependents_dbin_values[i]
            counter <- counter + 1
        } else {
            contribution_shape1 <- contribution_shape1 + dependents_dbin_values[i]
            contribution_shape2 <- contribution_shape2 + (dependents_dbin_size[i] - dependents_dbin_values[i])
        }
    }
    for (i in seq_along(dependents_dnegbin_nodeFunctions)) {
        if (counter == 1) {
            contribution_shape1 <- dependents_dnegbin_size[i]
            contribution_shape2 <- dependents_dnegbin_values[i]
            counter <- counter + 1
        } else {
            contribution_shape1 <- contribution_shape1 + dependents_dnegbin_size[i]
            contribution_shape2 <- contribution_shape2 + dependents_dnegbin_values[i]
        }
    }
    targetValue <- model[[targetNode]]
    posteriorLogDensity <- dbeta(targetValue, shape1 = prior_shape1 + contribution_shape1, shape2 = prior_shape2 + contribution_shape2, log = 1)
    returnType(double())
    return(posteriorLogDensity)
}, reset = function() {
}), where = getLoadingNamespace())




sampler_conjugate_dgamma <- nimbleFunction(contains = sampler_BASE, setup = function(model, mvSaved, control) {
    targetNode <- control$targetNode
    calcNodes <- model$getDependencies(targetNode)
    calcNodesDeterm <- model$getDependencies(targetNode, determOnly = TRUE)
    targetNode_nodeFunctionList <- nimbleFunctionList(node_stoch_dgamma)
    targetNode_nodeFunctionList[[1]] <- model$nodeFunctions[[targetNode]]
    dependents_dpois_nodeNames <- control$dependents_dpois
    dependents_dpois_nodeFunctions <- nimbleFunctionList(node_stoch_dpois)
    for (i in seq_along(dependents_dpois_nodeNames)) {
        dependents_dpois_nodeFunctions[[i]] <- model$nodeFunctions[[dependents_dpois_nodeNames[i]]]
    }
    dependents_dnorm_nodeNames <- control$dependents_dnorm
    dependents_dnorm_nodeFunctions <- nimbleFunctionList(node_stoch_dnorm)
    for (i in seq_along(dependents_dnorm_nodeNames)) {
        dependents_dnorm_nodeFunctions[[i]] <- model$nodeFunctions[[dependents_dnorm_nodeNames[i]]]
    }
    dependents_dlnorm_nodeNames <- control$dependents_dlnorm
    dependents_dlnorm_nodeFunctions <- nimbleFunctionList(node_stoch_dlnorm)
    for (i in seq_along(dependents_dlnorm_nodeNames)) {
        dependents_dlnorm_nodeFunctions[[i]] <- model$nodeFunctions[[dependents_dlnorm_nodeNames[i]]]
    }
    dependents_dgamma_nodeNames <- control$dependents_dgamma
    dependents_dgamma_nodeFunctions <- nimbleFunctionList(node_stoch_dgamma)
    for (i in seq_along(dependents_dgamma_nodeNames)) {
        dependents_dgamma_nodeFunctions[[i]] <- model$nodeFunctions[[dependents_dgamma_nodeNames[i]]]
    }
    dependents_dexp_nodeNames <- control$dependents_dexp
    dependents_dexp_nodeFunctions <- nimbleFunctionList(node_stoch_dexp)
    for (i in seq_along(dependents_dexp_nodeNames)) {
        dependents_dexp_nodeFunctions[[i]] <- model$nodeFunctions[[dependents_dexp_nodeNames[i]]]
    }
}, run = function() {
    modelLogProb0 <- getLogProb(model, calcNodes)
    origValue <- model[[targetNode]]
    prior_shape <- nfMethod(targetNode_nodeFunctionList[[1]], 'get_shape')()
    prior_rate <- nfMethod(targetNode_nodeFunctionList[[1]], 'get_rate')()
    declare(dependents_dpois_values, double(1, length(dependents_dpois_nodeFunctions)))
    declare(dependents_dpois_lambda, double(1, length(dependents_dpois_nodeFunctions)))
    for (i in seq_along(dependents_dpois_nodeFunctions)) {
        dependents_dpois_values[i] <- mfMethod(dependents_dpois_nodeFunctions[[i]], 'get_value')()
        dependents_dpois_lambda[i] <- nfMethod(dependents_dpois_nodeFunctions[[i]], 'get_lambda')()
    }
    declare(dependents_dnorm_values, double(1, length(dependents_dnorm_nodeFunctions)))
    declare(dependents_dnorm_tau, double(1, length(dependents_dnorm_nodeFunctions)))
    declare(dependents_dnorm_mean, double(1, length(dependents_dnorm_nodeFunctions)))
    for (i in seq_along(dependents_dnorm_nodeFunctions)) {
        dependents_dnorm_values[i] <- mfMethod(dependents_dnorm_nodeFunctions[[i]], 'get_value')()
        dependents_dnorm_tau[i] <- nfMethod(dependents_dnorm_nodeFunctions[[i]], 'get_tau')()
        dependents_dnorm_mean[i] <- nfMethod(dependents_dnorm_nodeFunctions[[i]], 'get_mean')()
    }
    declare(dependents_dlnorm_values, double(1, length(dependents_dlnorm_nodeFunctions)))
    declare(dependents_dlnorm_tau, double(1, length(dependents_dlnorm_nodeFunctions)))
    declare(dependents_dlnorm_meanlog, double(1, length(dependents_dlnorm_nodeFunctions)))
    for (i in seq_along(dependents_dlnorm_nodeFunctions)) {
        dependents_dlnorm_values[i] <- mfMethod(dependents_dlnorm_nodeFunctions[[i]], 'get_value')()
        dependents_dlnorm_tau[i] <- nfMethod(dependents_dlnorm_nodeFunctions[[i]], 'get_tau')()
        dependents_dlnorm_meanlog[i] <- nfMethod(dependents_dlnorm_nodeFunctions[[i]], 'get_meanlog')()
    }
    declare(dependents_dgamma_values, double(1, length(dependents_dgamma_nodeFunctions)))
    declare(dependents_dgamma_rate, double(1, length(dependents_dgamma_nodeFunctions)))
    declare(dependents_dgamma_shape, double(1, length(dependents_dgamma_nodeFunctions)))
    for (i in seq_along(dependents_dgamma_nodeFunctions)) {
        dependents_dgamma_values[i] <- mfMethod(dependents_dgamma_nodeFunctions[[i]], 'get_value')()
        dependents_dgamma_rate[i] <- nfMethod(dependents_dgamma_nodeFunctions[[i]], 'get_rate')()
        dependents_dgamma_shape[i] <- nfMethod(dependents_dgamma_nodeFunctions[[i]], 'get_shape')()
    }
    declare(dependents_dexp_values, double(1, length(dependents_dexp_nodeFunctions)))
    declare(dependents_dexp_rate, double(1, length(dependents_dexp_nodeFunctions)))
    for (i in seq_along(dependents_dexp_nodeFunctions)) {
        dependents_dexp_values[i] <- mfMethod(dependents_dexp_nodeFunctions[[i]], 'get_value')()
        dependents_dexp_rate[i] <- nfMethod(dependents_dexp_nodeFunctions[[i]], 'get_rate')()
    }
    declare(dependents_dpois_offset, double(1, length(dependents_dpois_nodeFunctions)))
    declare(dependents_dpois_coeff, double(1, length(dependents_dpois_nodeFunctions)))
    declare(dependents_dnorm_offset, double(1, length(dependents_dnorm_nodeFunctions)))
    declare(dependents_dnorm_coeff, double(1, length(dependents_dnorm_nodeFunctions)))
    declare(dependents_dlnorm_offset, double(1, length(dependents_dlnorm_nodeFunctions)))
    declare(dependents_dlnorm_coeff, double(1, length(dependents_dlnorm_nodeFunctions)))
    declare(dependents_dgamma_offset, double(1, length(dependents_dgamma_nodeFunctions)))
    declare(dependents_dgamma_coeff, double(1, length(dependents_dgamma_nodeFunctions)))
    declare(dependents_dexp_offset, double(1, length(dependents_dexp_nodeFunctions)))
    declare(dependents_dexp_coeff, double(1, length(dependents_dexp_nodeFunctions)))
    x1 <- model[[targetNode]]
    model[[targetNode]] <<- 0
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dependents_dpois_nodeFunctions)) {
        dependents_dpois_offset[i] <- nfMethod(dependents_dpois_nodeFunctions[[i]], 'get_lambda')()
        dependents_dpois_coeff[i] <- (dependents_dpois_lambda[i] - dependents_dpois_offset[i])/x1
    }
    for (i in seq_along(dependents_dnorm_nodeFunctions)) {
        dependents_dnorm_offset[i] <- nfMethod(dependents_dnorm_nodeFunctions[[i]], 'get_tau')()
        dependents_dnorm_coeff[i] <- (dependents_dnorm_tau[i] - dependents_dnorm_offset[i])/x1
    }
    for (i in seq_along(dependents_dlnorm_nodeFunctions)) {
        dependents_dlnorm_offset[i] <- nfMethod(dependents_dlnorm_nodeFunctions[[i]], 'get_tau')()
        dependents_dlnorm_coeff[i] <- (dependents_dlnorm_tau[i] - dependents_dlnorm_offset[i])/x1
    }
    for (i in seq_along(dependents_dgamma_nodeFunctions)) {
        dependents_dgamma_offset[i] <- nfMethod(dependents_dgamma_nodeFunctions[[i]], 'get_rate')()
        dependents_dgamma_coeff[i] <- (dependents_dgamma_rate[i] - dependents_dgamma_offset[i])/x1
    }
    for (i in seq_along(dependents_dexp_nodeFunctions)) {
        dependents_dexp_offset[i] <- nfMethod(dependents_dexp_nodeFunctions[[i]], 'get_rate')()
        dependents_dexp_coeff[i] <- (dependents_dexp_rate[i] - dependents_dexp_offset[i])/x1
    }
    counter <- 1
    for (i in seq_along(dependents_dpois_nodeFunctions)) {
        if (counter == 1) {
            contribution_shape <- dependents_dpois_values[i]
            contribution_rate <- dependents_dpois_coeff[i]
            counter <- counter + 1
        }
        else {
            contribution_shape <- contribution_shape + dependents_dpois_values[i]
            contribution_rate <- contribution_rate + dependents_dpois_coeff[i]
        }
    }
    for (i in seq_along(dependents_dnorm_nodeFunctions)) {
        if (counter == 1) {
            contribution_shape <- 1/2
            contribution_rate <- dependents_dnorm_coeff[i]/2 * (dependents_dnorm_values[i] - dependents_dnorm_mean[i])^2
            counter <- counter + 1
        }
        else {
            contribution_shape <- contribution_shape + 1/2
            contribution_rate <- contribution_rate + dependents_dnorm_coeff[i]/2 * (dependents_dnorm_values[i] - dependents_dnorm_mean[i])^2
        }
    }
    for (i in seq_along(dependents_dlnorm_nodeFunctions)) {
        if (counter == 1) {
            contribution_shape <- 1/2
            contribution_rate <- dependents_dlnorm_coeff[i]/2 * (log(dependents_dlnorm_values[i]) - dependents_dlnorm_meanlog[i])^2
            counter <- counter + 1
        }
        else {
            contribution_shape <- contribution_shape + 1/2
            contribution_rate <- contribution_rate + dependents_dlnorm_coeff[i]/2 * (log(dependents_dlnorm_values[i]) - dependents_dlnorm_meanlog[i])^2
        }
    }
    for (i in seq_along(dependents_dgamma_nodeFunctions)) {
        if (counter == 1) {
            contribution_shape <- dependents_dgamma_shape[i]
            contribution_rate <- dependents_dgamma_coeff[i] * dependents_dgamma_values[i]
            counter <- counter + 1
        }
        else {
            contribution_shape <- contribution_shape + dependents_dgamma_shape[i]
            contribution_rate <- contribution_rate + dependents_dgamma_coeff[i] * dependents_dgamma_values[i]
        }
    }
    for (i in seq_along(dependents_dexp_nodeFunctions)) {
        if (counter == 1) {
            contribution_shape <- 1
            contribution_rate <- dependents_dexp_coeff[i] * dependents_dexp_values[i]
            counter <- counter + 1
        }
        else {
            contribution_shape <- contribution_shape + 1
            contribution_rate <- contribution_rate + dependents_dexp_coeff[i] * dependents_dexp_values[i]
        }
    }
    newValue <- rgamma(1, shape = prior_shape + contribution_shape, scale = 1/(prior_rate + contribution_rate))
    model[[targetNode]] <<- newValue
    calculate(model, calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    modelLogProb1 <- getLogProb(model, calcNodes)
    posteriorLogDensity0 <- dgamma(origValue, shape = prior_shape + contribution_shape, scale = 1/(prior_rate + contribution_rate), log = 1)
    posteriorLogDensity1 <- dgamma(newValue, shape = prior_shape + contribution_shape, scale = 1/(prior_rate + contribution_rate), log = 1)
    posteriorVerification <- modelLogProb0 - posteriorLogDensity0 - modelLogProb1 + posteriorLogDensity1
    if (abs(posteriorVerification) > 1e-10) {
        nimPrint('conjugate posterior density appears to be wrong')
    }
}, methods = list(getPosteriorLogDensity = function() {
    prior_shape <- nfMethod(targetNode_nodeFunctionList[[1]], 'get_shape')()
    prior_rate <- nfMethod(targetNode_nodeFunctionList[[1]], 'get_rate')()
    declare(dependents_dpois_values, double(1, length(dependents_dpois_nodeFunctions)))
    declare(dependents_dpois_lambda, double(1, length(dependents_dpois_nodeFunctions)))
    for (i in seq_along(dependents_dpois_nodeFunctions)) {
        dependents_dpois_values[i] <- mfMethod(dependents_dpois_nodeFunctions[[i]], 'get_value')()
        dependents_dpois_lambda[i] <- nfMethod(dependents_dpois_nodeFunctions[[i]], 'get_lambda')()
    }
    declare(dependents_dnorm_values, double(1, length(dependents_dnorm_nodeFunctions)))
    declare(dependents_dnorm_tau, double(1, length(dependents_dnorm_nodeFunctions)))
    declare(dependents_dnorm_mean, double(1, length(dependents_dnorm_nodeFunctions)))
    for (i in seq_along(dependents_dnorm_nodeFunctions)) {
        dependents_dnorm_values[i] <- mfMethod(dependents_dnorm_nodeFunctions[[i]], 'get_value')()
        dependents_dnorm_tau[i] <- nfMethod(dependents_dnorm_nodeFunctions[[i]], 'get_tau')()
        dependents_dnorm_mean[i] <- nfMethod(dependents_dnorm_nodeFunctions[[i]], 'get_mean')()
    }
    declare(dependents_dlnorm_values, double(1, length(dependents_dlnorm_nodeFunctions)))
    declare(dependents_dlnorm_tau, double(1, length(dependents_dlnorm_nodeFunctions)))
    declare(dependents_dlnorm_meanlog, double(1, length(dependents_dlnorm_nodeFunctions)))
    for (i in seq_along(dependents_dlnorm_nodeFunctions)) {
        dependents_dlnorm_values[i] <- mfMethod(dependents_dlnorm_nodeFunctions[[i]], 'get_value')()
        dependents_dlnorm_tau[i] <- nfMethod(dependents_dlnorm_nodeFunctions[[i]], 'get_tau')()
        dependents_dlnorm_meanlog[i] <- nfMethod(dependents_dlnorm_nodeFunctions[[i]], 'get_meanlog')()
    }
    declare(dependents_dgamma_values, double(1, length(dependents_dgamma_nodeFunctions)))
    declare(dependents_dgamma_rate, double(1, length(dependents_dgamma_nodeFunctions)))
    declare(dependents_dgamma_shape, double(1, length(dependents_dgamma_nodeFunctions)))
    for (i in seq_along(dependents_dgamma_nodeFunctions)) {
        dependents_dgamma_values[i] <- mfMethod(dependents_dgamma_nodeFunctions[[i]], 'get_value')()
        dependents_dgamma_rate[i] <- nfMethod(dependents_dgamma_nodeFunctions[[i]], 'get_rate')()
        dependents_dgamma_shape[i] <- nfMethod(dependents_dgamma_nodeFunctions[[i]], 'get_shape')()
    }
    declare(dependents_dexp_values, double(1, length(dependents_dexp_nodeFunctions)))
    declare(dependents_dexp_rate, double(1, length(dependents_dexp_nodeFunctions)))
    for (i in seq_along(dependents_dexp_nodeFunctions)) {
        dependents_dexp_values[i] <- mfMethod(dependents_dexp_nodeFunctions[[i]], 'get_value')()
        dependents_dexp_rate[i] <- nfMethod(dependents_dexp_nodeFunctions[[i]], 'get_rate')()
    }
    declare(dependents_dpois_offset, double(1, length(dependents_dpois_nodeFunctions)))
    declare(dependents_dpois_coeff, double(1, length(dependents_dpois_nodeFunctions)))
    declare(dependents_dnorm_offset, double(1, length(dependents_dnorm_nodeFunctions)))
    declare(dependents_dnorm_coeff, double(1, length(dependents_dnorm_nodeFunctions)))
    declare(dependents_dlnorm_offset, double(1, length(dependents_dlnorm_nodeFunctions)))
    declare(dependents_dlnorm_coeff, double(1, length(dependents_dlnorm_nodeFunctions)))
    declare(dependents_dgamma_offset, double(1, length(dependents_dgamma_nodeFunctions)))
    declare(dependents_dgamma_coeff, double(1, length(dependents_dgamma_nodeFunctions)))
    declare(dependents_dexp_offset, double(1, length(dependents_dexp_nodeFunctions)))
    declare(dependents_dexp_coeff, double(1, length(dependents_dexp_nodeFunctions)))
    x1 <- model[[targetNode]]
    model[[targetNode]] <<- 0
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dependents_dpois_nodeFunctions)) {
        dependents_dpois_offset[i] <- nfMethod(dependents_dpois_nodeFunctions[[i]], 'get_lambda')()
        dependents_dpois_coeff[i] <- (dependents_dpois_lambda[i] - dependents_dpois_offset[i])/x1
    }
    for (i in seq_along(dependents_dnorm_nodeFunctions)) {
        dependents_dnorm_offset[i] <- nfMethod(dependents_dnorm_nodeFunctions[[i]], 'get_tau')()
        dependents_dnorm_coeff[i] <- (dependents_dnorm_tau[i] - dependents_dnorm_offset[i])/x1
    }
    for (i in seq_along(dependents_dlnorm_nodeFunctions)) {
        dependents_dlnorm_offset[i] <- nfMethod(dependents_dlnorm_nodeFunctions[[i]], 'get_tau')()
        dependents_dlnorm_coeff[i] <- (dependents_dlnorm_tau[i] - dependents_dlnorm_offset[i])/x1
    }
    for (i in seq_along(dependents_dgamma_nodeFunctions)) {
        dependents_dgamma_offset[i] <- nfMethod(dependents_dgamma_nodeFunctions[[i]], 'get_rate')()
        dependents_dgamma_coeff[i] <- (dependents_dgamma_rate[i] - dependents_dgamma_offset[i])/x1
    }
    for (i in seq_along(dependents_dexp_nodeFunctions)) {
        dependents_dexp_offset[i] <- nfMethod(dependents_dexp_nodeFunctions[[i]], 'get_rate')()
        dependents_dexp_coeff[i] <- (dependents_dexp_rate[i] - dependents_dexp_offset[i])/x1
    }
    counter <- 1
    for (i in seq_along(dependents_dpois_nodeFunctions)) {
        if (counter == 1) {
            contribution_shape <- dependents_dpois_values[i]
            contribution_rate <- dependents_dpois_coeff[i]
            counter <- counter + 1
        } else {
            contribution_shape <- contribution_shape + dependents_dpois_values[i]
            contribution_rate <- contribution_rate + dependents_dpois_coeff[i]
        }
    }
    for (i in seq_along(dependents_dnorm_nodeFunctions)) {
        if (counter == 1) {
            contribution_shape <- 1/2
            contribution_rate <- dependents_dnorm_coeff[i]/2 * (dependents_dnorm_values[i] - dependents_dnorm_mean[i])^2
            counter <- counter + 1
        } else {
            contribution_shape <- contribution_shape + 1/2
            contribution_rate <- contribution_rate + dependents_dnorm_coeff[i]/2 * (dependents_dnorm_values[i] - dependents_dnorm_mean[i])^2
        }
    }
    for (i in seq_along(dependents_dlnorm_nodeFunctions)) {
        if (counter == 1) {
            contribution_shape <- 1/2
            contribution_rate <- dependents_dlnorm_coeff[i]/2 * (log(dependents_dlnorm_values[i]) - dependents_dlnorm_meanlog[i])^2
            counter <- counter + 1
        } else {
            contribution_shape <- contribution_shape + 1/2
            contribution_rate <- contribution_rate + dependents_dlnorm_coeff[i]/2 * (log(dependents_dlnorm_values[i]) - dependents_dlnorm_meanlog[i])^2
        }
    }
    for (i in seq_along(dependents_dgamma_nodeFunctions)) {
        if (counter == 1) {
            contribution_shape <- dependents_dgamma_shape[i]
            contribution_rate <- dependents_dgamma_coeff[i] * dependents_dgamma_values[i]
            counter <- counter + 1
        } else {
            contribution_shape <- contribution_shape + dependents_dgamma_shape[i]
            contribution_rate <- contribution_rate + dependents_dgamma_coeff[i] * dependents_dgamma_values[i]
        }
    }
    for (i in seq_along(dependents_dexp_nodeFunctions)) {
        if (counter == 1) {
            contribution_shape <- 1
            contribution_rate <- dependents_dexp_coeff[i] * dependents_dexp_values[i]
            counter <- counter + 1
        } else {
            contribution_shape <- contribution_shape + 1
            contribution_rate <- contribution_rate + dependents_dexp_coeff[i] * dependents_dexp_values[i]
        }
    }
    targetValue <- model[[targetNode]]
    posteriorLogDensity <- dgamma(targetValue, shape = prior_shape + contribution_shape, scale = 1/(prior_rate + contribution_rate), log = 1)
    returnType(double())
    return(posteriorLogDensity)
}, reset = function() {
}), where = getLoadingNamespace())




sampler_conjugate_dnorm <- nimbleFunction(contains = sampler_BASE, setup = function(model, mvSaved, control) {
    targetNode <- control$targetNode
    calcNodes <- model$getDependencies(targetNode)
    calcNodesDeterm <- model$getDependencies(targetNode, determOnly = TRUE)
    targetNode_nodeFunctionList <- nimbleFunctionList(node_stoch_dnorm)
    targetNode_nodeFunctionList[[1]] <- model$nodeFunctions[[targetNode]]
    dependents_dnorm_nodeNames <- control$dependents_dnorm
    dependents_dnorm_nodeFunctions <- nimbleFunctionList(node_stoch_dnorm)
    for (i in seq_along(dependents_dnorm_nodeNames)) {
        dependents_dnorm_nodeFunctions[[i]] <- model$nodeFunctions[[dependents_dnorm_nodeNames[i]]]
    }
    dependents_dlnorm_nodeNames <- control$dependents_dlnorm
    dependents_dlnorm_nodeFunctions <- nimbleFunctionList(node_stoch_dlnorm)
    for (i in seq_along(dependents_dlnorm_nodeNames)) {
        dependents_dlnorm_nodeFunctions[[i]] <- model$nodeFunctions[[dependents_dlnorm_nodeNames[i]]]
    }
}, run = function() {
    modelLogProb0 <- getLogProb(model, calcNodes)
    origValue <- model[[targetNode]]
    prior_mean <- nfMethod(targetNode_nodeFunctionList[[1]], 'get_mean')()
    prior_tau <- nfMethod(targetNode_nodeFunctionList[[1]], 'get_tau')()
    declare(dependents_dnorm_values, double(1, length(dependents_dnorm_nodeFunctions)))
    declare(dependents_dnorm_mean, double(1, length(dependents_dnorm_nodeFunctions)))
    declare(dependents_dnorm_tau, double(1, length(dependents_dnorm_nodeFunctions)))
    for (i in seq_along(dependents_dnorm_nodeFunctions)) {
        dependents_dnorm_values[i] <- mfMethod(dependents_dnorm_nodeFunctions[[i]], 'get_value')()
        dependents_dnorm_mean[i] <- nfMethod(dependents_dnorm_nodeFunctions[[i]], 'get_mean')()
        dependents_dnorm_tau[i] <- nfMethod(dependents_dnorm_nodeFunctions[[i]], 'get_tau')()
    }
    declare(dependents_dlnorm_values, double(1, length(dependents_dlnorm_nodeFunctions)))
    declare(dependents_dlnorm_meanlog, double(1, length(dependents_dlnorm_nodeFunctions)))
    declare(dependents_dlnorm_tau, double(1, length(dependents_dlnorm_nodeFunctions)))
    for (i in seq_along(dependents_dlnorm_nodeFunctions)) {
        dependents_dlnorm_values[i] <- mfMethod(dependents_dlnorm_nodeFunctions[[i]], 'get_value')()
        dependents_dlnorm_meanlog[i] <- nfMethod(dependents_dlnorm_nodeFunctions[[i]], 'get_meanlog')()
        dependents_dlnorm_tau[i] <- nfMethod(dependents_dlnorm_nodeFunctions[[i]], 'get_tau')()
    }
    declare(dependents_dnorm_offset, double(1, length(dependents_dnorm_nodeFunctions)))
    declare(dependents_dnorm_coeff, double(1, length(dependents_dnorm_nodeFunctions)))
    declare(dependents_dlnorm_offset, double(1, length(dependents_dlnorm_nodeFunctions)))
    declare(dependents_dlnorm_coeff, double(1, length(dependents_dlnorm_nodeFunctions)))
    x1 <- model[[targetNode]]
    model[[targetNode]] <<- 0
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dependents_dnorm_nodeFunctions)) {
        dependents_dnorm_offset[i] <- nfMethod(dependents_dnorm_nodeFunctions[[i]], 'get_mean')()
        dependents_dnorm_coeff[i] <- (dependents_dnorm_mean[i] - dependents_dnorm_offset[i])/x1
    }
    for (i in seq_along(dependents_dlnorm_nodeFunctions)) {
        dependents_dlnorm_offset[i] <- nfMethod(dependents_dlnorm_nodeFunctions[[i]], 'get_meanlog')()
        dependents_dlnorm_coeff[i] <- (dependents_dlnorm_meanlog[i] - dependents_dlnorm_offset[i])/x1
    }
    counter <- 1
    for (i in seq_along(dependents_dnorm_nodeFunctions)) {
        if (counter == 1) {
            contribution_mean <- dependents_dnorm_coeff[i] * (dependents_dnorm_values[i] - dependents_dnorm_offset[i]) * dependents_dnorm_tau[i]
            contribution_tau <- dependents_dnorm_coeff[i]^2 * dependents_dnorm_tau[i]
            counter <- counter + 1
        }
        else {
            contribution_mean <- contribution_mean + dependents_dnorm_coeff[i] * (dependents_dnorm_values[i] - dependents_dnorm_offset[i]) * dependents_dnorm_tau[i]
            contribution_tau <- contribution_tau + dependents_dnorm_coeff[i]^2 * dependents_dnorm_tau[i]
        }
    }
    for (i in seq_along(dependents_dlnorm_nodeFunctions)) {
        if (counter == 1) {
            contribution_mean <- dependents_dlnorm_coeff[i] * (log(dependents_dlnorm_values[i]) - dependents_dlnorm_offset[i]) * dependents_dlnorm_tau[i]
            contribution_tau <- dependents_dlnorm_coeff[i]^2 * dependents_dlnorm_tau[i]
            counter <- counter + 1
        }
        else {
            contribution_mean <- contribution_mean + dependents_dlnorm_coeff[i] * (log(dependents_dlnorm_values[i]) - dependents_dlnorm_offset[i]) * dependents_dlnorm_tau[i]
            contribution_tau <- contribution_tau + dependents_dlnorm_coeff[i]^2 * dependents_dlnorm_tau[i]
        }
    }
    newValue <- rnorm(1, mean = (prior_mean * prior_tau + contribution_mean)/(prior_tau + contribution_tau), sd = (prior_tau + contribution_tau)^(-0.5))
    model[[targetNode]] <<- newValue
    calculate(model, calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    modelLogProb1 <- getLogProb(model, calcNodes)
    posteriorLogDensity0 <- dnorm(origValue, mean = (prior_mean * prior_tau + contribution_mean)/(prior_tau + contribution_tau), sd = (prior_tau + contribution_tau)^(-0.5), log = 1)
    posteriorLogDensity1 <- dnorm(newValue, mean = (prior_mean * prior_tau + contribution_mean)/(prior_tau + contribution_tau), sd = (prior_tau + contribution_tau)^(-0.5), log = 1)
    posteriorVerification <- modelLogProb0 - posteriorLogDensity0 - modelLogProb1 + posteriorLogDensity1
    if (abs(posteriorVerification) > 1e-10) {
        nimPrint('conjugate posterior density appears to be wrong')
    }
}, methods = list(getPosteriorLogDensity = function() {
    prior_mean <- nfMethod(targetNode_nodeFunctionList[[1]], 'get_mean')()
    prior_tau <- nfMethod(targetNode_nodeFunctionList[[1]], 'get_tau')()
    declare(dependents_dnorm_values, double(1, length(dependents_dnorm_nodeFunctions)))
    declare(dependents_dnorm_mean, double(1, length(dependents_dnorm_nodeFunctions)))
    declare(dependents_dnorm_tau, double(1, length(dependents_dnorm_nodeFunctions)))
    for (i in seq_along(dependents_dnorm_nodeFunctions)) {
        dependents_dnorm_values[i] <- mfMethod(dependents_dnorm_nodeFunctions[[i]], 'get_value')()
        dependents_dnorm_mean[i] <- nfMethod(dependents_dnorm_nodeFunctions[[i]], 'get_mean')()
        dependents_dnorm_tau[i] <- nfMethod(dependents_dnorm_nodeFunctions[[i]], 'get_tau')()
    }
    declare(dependents_dlnorm_values, double(1, length(dependents_dlnorm_nodeFunctions)))
    declare(dependents_dlnorm_meanlog, double(1, length(dependents_dlnorm_nodeFunctions)))
    declare(dependents_dlnorm_tau, double(1, length(dependents_dlnorm_nodeFunctions)))
    for (i in seq_along(dependents_dlnorm_nodeFunctions)) {
        dependents_dlnorm_values[i] <- mfMethod(dependents_dlnorm_nodeFunctions[[i]], 'get_value')()
        dependents_dlnorm_meanlog[i] <- nfMethod(dependents_dlnorm_nodeFunctions[[i]], 'get_meanlog')()
        dependents_dlnorm_tau[i] <- nfMethod(dependents_dlnorm_nodeFunctions[[i]], 'get_tau')()
    }
    declare(dependents_dnorm_offset, double(1, length(dependents_dnorm_nodeFunctions)))
    declare(dependents_dnorm_coeff, double(1, length(dependents_dnorm_nodeFunctions)))
    declare(dependents_dlnorm_offset, double(1, length(dependents_dlnorm_nodeFunctions)))
    declare(dependents_dlnorm_coeff, double(1, length(dependents_dlnorm_nodeFunctions)))
    x1 <- model[[targetNode]]
    model[[targetNode]] <<- 0
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dependents_dnorm_nodeFunctions)) {
        dependents_dnorm_offset[i] <- nfMethod(dependents_dnorm_nodeFunctions[[i]], 'get_mean')()
        dependents_dnorm_coeff[i] <- (dependents_dnorm_mean[i] - dependents_dnorm_offset[i])/x1
    }
    for (i in seq_along(dependents_dlnorm_nodeFunctions)) {
        dependents_dlnorm_offset[i] <- nfMethod(dependents_dlnorm_nodeFunctions[[i]], 'get_meanlog')()
        dependents_dlnorm_coeff[i] <- (dependents_dlnorm_meanlog[i] - dependents_dlnorm_offset[i])/x1
    }
    counter <- 1
    for (i in seq_along(dependents_dnorm_nodeFunctions)) {
        if (counter == 1) {
            contribution_mean <- dependents_dnorm_coeff[i] * (dependents_dnorm_values[i] - dependents_dnorm_offset[i]) * dependents_dnorm_tau[i]
            contribution_tau <- dependents_dnorm_coeff[i]^2 * dependents_dnorm_tau[i]
            counter <- counter + 1
        } else {
            contribution_mean <- contribution_mean + dependents_dnorm_coeff[i] * (dependents_dnorm_values[i] - dependents_dnorm_offset[i]) * dependents_dnorm_tau[i]
            contribution_tau <- contribution_tau + dependents_dnorm_coeff[i]^2 * dependents_dnorm_tau[i]
        }
    }
    for (i in seq_along(dependents_dlnorm_nodeFunctions)) {
        if (counter == 1) {
            contribution_mean <- dependents_dlnorm_coeff[i] * (log(dependents_dlnorm_values[i]) - dependents_dlnorm_offset[i]) * dependents_dlnorm_tau[i]
            contribution_tau <- dependents_dlnorm_coeff[i]^2 * dependents_dlnorm_tau[i]
            counter <- counter + 1
        } else {
            contribution_mean <- contribution_mean + dependents_dlnorm_coeff[i] * (log(dependents_dlnorm_values[i]) - dependents_dlnorm_offset[i]) * dependents_dlnorm_tau[i]
            contribution_tau <- contribution_tau + dependents_dlnorm_coeff[i]^2 * dependents_dlnorm_tau[i]
        }
    }
    targetValue <- model[[targetNode]]
    posteriorLogDensity <- dnorm(targetValue, mean = (prior_mean * prior_tau + contribution_mean)/(prior_tau + contribution_tau), sd = (prior_tau + contribution_tau)^(-0.5), log = 1)
    returnType(double())
    return(posteriorLogDensity)
}, reset = function() {
}), where = getLoadingNamespace())




sampler_conjugate_dmnorm <- nimbleFunction(contains = sampler_BASE, setup = function(model, mvSaved, control) {
    targetNode <- control$targetNode
    calcNodes <- model$getDependencies(targetNode)
    calcNodesDeterm <- model$getDependencies(targetNode, determOnly = TRUE)
    targetNode_nodeFunctionList <- nimbleFunctionList(node_stoch_dmnorm)
    targetNode_nodeFunctionList[[1]] <- model$nodeFunctions[[targetNode]]
    dependents_dmnorm_nodeNames <- control$dependents_dmnorm
    dependents_dmnorm_nodeFunctions <- nimbleFunctionList(node_stoch_dmnorm)
    for (i in seq_along(dependents_dmnorm_nodeNames)) {
        dependents_dmnorm_nodeFunctions[[i]] <- model$nodeFunctions[[dependents_dmnorm_nodeNames[i]]]
    }
    d <- model$getNodeInfo()[[targetNode]]$targetNodeIndexSizes[1]
}, run = function() {
    modelLogProb0 <- getLogProb(model, calcNodes)
    origValue <- model[[targetNode]]
    prior_prec <- nfMethod(targetNode_nodeFunctionList[[1]], 'get_prec')()
    prior_mean <- nfMethod(targetNode_nodeFunctionList[[1]], 'get_mean')()
    declare(dependents_dmnorm_values, double(2, c(length(dependents_dmnorm_nodeFunctions), d)))
    declare(dependents_dmnorm_mean, double(2, c(length(dependents_dmnorm_nodeFunctions), d)))
    declare(dependents_dmnorm_prec, double(3, c(length(dependents_dmnorm_nodeFunctions), d, d)))
    for (i in seq_along(dependents_dmnorm_nodeFunctions)) {
        dependents_dmnorm_values[i, 1:d] <- mfMethod(dependents_dmnorm_nodeFunctions[[i]], 'get_value')()
        dependents_dmnorm_mean[i, 1:d] <- nfMethod(dependents_dmnorm_nodeFunctions[[i]], 'get_mean')()
        dependents_dmnorm_prec[i, 1:d, 1:d] <- nfMethod(dependents_dmnorm_nodeFunctions[[i]], 'get_prec')()
    }
    declare(dependents_dmnorm_offset, double(2, c(length(dependents_dmnorm_nodeFunctions), d)))
    declare(dependents_dmnorm_coeff, double(3, c(length(dependents_dmnorm_nodeFunctions), d, d)))
    model[[targetNode]] <<- model[[targetNode]] * 0
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dependents_dmnorm_nodeFunctions)) {
        dependents_dmnorm_offset[i, 1:d] <- nfMethod(dependents_dmnorm_nodeFunctions[[i]], 'get_mean')()
    }
    for (sizeIndex in 1:d) {
        unitVector <- model[[targetNode]] * 0
        unitVector[sizeIndex] <- 1
        model[[targetNode]] <<- unitVector
        calculate(model, calcNodesDeterm)
        for (i in seq_along(dependents_dmnorm_nodeFunctions)) {
            dependents_dmnorm_coeff[i, 1:d, sizeIndex] <- nfMethod(dependents_dmnorm_nodeFunctions[[i]], 'get_mean')() - dependents_dmnorm_offset[i, 1:d]
        }
    }
    counter <- 1
    for (i in seq_along(dependents_dmnorm_nodeFunctions)) {
        if (counter == 1) {
            contribution_prec <- t(dependents_dmnorm_coeff[i, 1:d, 1:d]) %*% dependents_dmnorm_prec[i, 1:d, 1:d] %*% dependents_dmnorm_coeff[i, 1:d, 1:d]
            contribution_mean <- t(dependents_dmnorm_coeff[i, 1:d, 1:d]) %*% dependents_dmnorm_prec[i, 1:d, 1:d] %*% asCol(dependents_dmnorm_values[i, 1:d] - dependents_dmnorm_offset[i, 1:d])
            counter <- counter + 1
        }
        else {
            contribution_prec <- contribution_prec + t(dependents_dmnorm_coeff[i, 1:d, 1:d]) %*% dependents_dmnorm_prec[i, 1:d, 1:d] %*% dependents_dmnorm_coeff[i, 1:d, 1:d]
            contribution_mean <- contribution_mean + t(dependents_dmnorm_coeff[i, 1:d, 1:d]) %*% dependents_dmnorm_prec[i, 1:d, 1:d] %*% asCol(dependents_dmnorm_values[i, 1:d] - dependents_dmnorm_offset[i, 1:d])
        }
    }
    newValue <- rmnorm(1, mean = inverse(prior_prec + contribution_prec) %*% (prior_prec %*% asCol(prior_mean) + contribution_mean), chol = chol(prior_prec + contribution_prec), prec_param = TRUE)
    model[[targetNode]] <<- newValue
    calculate(model, calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    modelLogProb1 <- getLogProb(model, calcNodes)
    posteriorLogDensity0 <- dmnorm(origValue, mean = inverse(prior_prec + contribution_prec) %*% (prior_prec %*% asCol(prior_mean) + contribution_mean), chol = chol(prior_prec + contribution_prec), prec_param = TRUE, log = 1)
    posteriorLogDensity1 <- dmnorm(newValue, mean = inverse(prior_prec + contribution_prec) %*% (prior_prec %*% asCol(prior_mean) + contribution_mean), chol = chol(prior_prec + contribution_prec), prec_param = TRUE, log = 1)
    posteriorVerification <- modelLogProb0 - posteriorLogDensity0 - modelLogProb1 + posteriorLogDensity1
    if (abs(posteriorVerification) > 1e-10) {
        nimPrint('conjugate posterior density appears to be wrong')
    }
}, methods = list(getPosteriorLogDensity = function() {
    prior_prec <- nfMethod(targetNode_nodeFunctionList[[1]], 'get_prec')()
    prior_mean <- nfMethod(targetNode_nodeFunctionList[[1]], 'get_mean')()
    declare(dependents_dmnorm_values, double(2, c(length(dependents_dmnorm_nodeFunctions), d)))
    declare(dependents_dmnorm_mean, double(2, c(length(dependents_dmnorm_nodeFunctions), d)))
    declare(dependents_dmnorm_prec, double(3, c(length(dependents_dmnorm_nodeFunctions), d, d)))
    for (i in seq_along(dependents_dmnorm_nodeFunctions)) {
        dependents_dmnorm_values[i, 1:d] <- mfMethod(dependents_dmnorm_nodeFunctions[[i]], 'get_value')()
        dependents_dmnorm_mean[i, 1:d] <- nfMethod(dependents_dmnorm_nodeFunctions[[i]], 'get_mean')()
        dependents_dmnorm_prec[i, 1:d, 1:d] <- nfMethod(dependents_dmnorm_nodeFunctions[[i]], 'get_prec')()
    }
    declare(dependents_dmnorm_offset, double(2, c(length(dependents_dmnorm_nodeFunctions), d)))
    declare(dependents_dmnorm_coeff, double(3, c(length(dependents_dmnorm_nodeFunctions), d, d)))
    model[[targetNode]] <<- model[[targetNode]] * 0
    calculate(model, calcNodesDeterm)
    for (i in seq_along(dependents_dmnorm_nodeFunctions)) {
        dependents_dmnorm_offset[i, 1:d] <- nfMethod(dependents_dmnorm_nodeFunctions[[i]], 'get_mean')()
    }
    for (sizeIndex in 1:d) {
        unitVector <- model[[targetNode]] * 0
        unitVector[sizeIndex] <- 1
        model[[targetNode]] <<- unitVector
        calculate(model, calcNodesDeterm)
        for (i in seq_along(dependents_dmnorm_nodeFunctions)) {
            dependents_dmnorm_coeff[i, 1:d, sizeIndex] <- nfMethod(dependents_dmnorm_nodeFunctions[[i]], 'get_mean')() - dependents_dmnorm_offset[i, 1:d]
        }
    }
    counter <- 1
    for (i in seq_along(dependents_dmnorm_nodeFunctions)) {
        if (counter == 1) {
            contribution_prec <- t(dependents_dmnorm_coeff[i, 1:d, 1:d]) %*% dependents_dmnorm_prec[i, 1:d, 1:d] %*% dependents_dmnorm_coeff[i, 1:d, 1:d]
            contribution_mean <- t(dependents_dmnorm_coeff[i, 1:d, 1:d]) %*% dependents_dmnorm_prec[i, 1:d, 1:d] %*% asCol(dependents_dmnorm_values[i, 1:d] - dependents_dmnorm_offset[i, 1:d])
            counter <- counter + 1
        } else {
            contribution_prec <- contribution_prec + t(dependents_dmnorm_coeff[i, 1:d, 1:d]) %*% dependents_dmnorm_prec[i, 1:d, 1:d] %*% dependents_dmnorm_coeff[i, 1:d, 1:d]
            contribution_mean <- contribution_mean + t(dependents_dmnorm_coeff[i, 1:d, 1:d]) %*% dependents_dmnorm_prec[i, 1:d, 1:d] %*% asCol(dependents_dmnorm_values[i, 1:d] - dependents_dmnorm_offset[i, 1:d])
        }
    }
    targetValue <- model[[targetNode]]
    posteriorLogDensity <- dmnorm(targetValue, mean = inverse(prior_prec + contribution_prec) %*% (prior_prec %*% asCol(prior_mean) + contribution_mean), chol = chol(prior_prec + contribution_prec), prec_param = TRUE, log = 1)
    returnType(double())
    return(posteriorLogDensity)
}, reset = function() {
}), where = getLoadingNamespace())





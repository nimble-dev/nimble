# used in altParams for dmnorm
# this needs to be sourced after nimbleFunction() is defined, so can't be done in distributions_inputList.R
calc_dmnormAltParams <- nimbleFunction(
    run = function(cholesky = double(2), prec_param = double(0), return_prec = double(0)) {
        ## original (inefficient) implementation:
        ## tmp <- t(cholesky) %*% cholesky
        ## if(prec_param == return_prec) {
        ##     return(tmp)
        ## } else {
        ##     return(inverse(tmp))
        ## }
        if(prec_param == return_prec) {
            ans <- t(cholesky) %*% cholesky
        } else {
            I <- identityMatrix(dim(cholesky)[1])
            ans <- backsolve(cholesky, forwardsolve(t(cholesky), I))
        }
        returnType(double(2))
        return(ans)
    })



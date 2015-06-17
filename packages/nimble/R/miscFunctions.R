# used in altParams for dmnorm
# this needs to be sourced after nimbleFunction() is defined, so can't be done in distributions_inputList.R
calc_dmnormAltParams <- nimbleFunction(
    run = function(cholesky = double(2), prec_param = double(0), return_prec = double(0)) {
        returnType(double(2))
        tmp <- t(cholesky) %*% cholesky
        if(prec_param == return_prec) {
            return(tmp)
        } else {
            return(inverse(tmp))
        }
    })
        


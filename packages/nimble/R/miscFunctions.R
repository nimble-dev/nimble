# used in altParams for dmnorm
# this needs to be sourced after nimbleFunction() is defined, so can't be done in distributions_inputList.R
calc_dmnormAltParams <- nimbleFunction(
    run = function(cholesky = double(2), prec_param = double(), return_prec = double()) {
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
    }
)




## used in conjugacy definition for dmnorm, to calculate 'contribution' terms;
## avoids unnecessary matrix multiplications, when 'coeff' is identity matrix
calc_dmnormConjugacyContributions <- nimbleFunction(
    run = function(coeff = double(2), prec = double(2), order = double()) {
        d <- dim(coeff)[1]
        total <- 0
        for(i in 1:d) {
            for(j in 1:d) {
                if(i == j) { total <- total + abs(coeff[i,j] - 1)
                } else     { total <- total + abs(coeff[i,j])     }
            }
        }
        if(total < d * 1E-15) return(prec)
        if(order == 1) ans <- t(coeff) %*% prec
        if(order == 2) ans <- t(coeff) %*% prec %*% coeff
        return(ans)
        returnType(double(2))
    }
)



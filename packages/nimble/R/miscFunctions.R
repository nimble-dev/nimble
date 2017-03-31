# used in altParams for dmnorm
# this needs to be sourced after nimbleFunction() is defined, so can't be done in distributions_inputList.R
calc_dmnormAltParams <- nimbleFunction(
    run = function(cholesky = double(2), prec_param = double(), return_prec = double()) {
        if(prec_param == return_prec) {
            ans <- t(cholesky) %*% cholesky
        } else {
            I <- identityMatrix(dim(cholesky)[1])
            ans <- backsolve(cholesky, forwardsolve(t(cholesky), I))
            ## Chris suggests:
            ## tmp <- forwardsolve(L, I)
            ## ans <- crossprod(tmp)
        }
        returnType(double(2))
        return(ans)
    }
)




## used in conjugacy definition for dmnorm, to calculate 'contribution' terms;
## avoids unnecessary matrix multiplications, when 'coeff' is identity matrix
calc_dmnormConjugacyContributions <- nimbleFunction(
    run = function(coeff = double(2), prec = double(2), vec = double(1), order = double()) {
        if(dim(coeff)[1] == dim(coeff)[2]) {
            d <- dim(coeff)[1]
            total <- 0
            for(i in 1:d) {
                for(j in 1:d) {
                    if(i == j) { total <- total + abs(coeff[i,j] - 1)
                             } else     { total <- total + abs(coeff[i,j])     }
                }
            }
            if(total < d * 1E-15) return(prec)
        }
        if(order == 1) ans <- t(coeff) %*% (prec %*% asCol(vec))
        if(order == 2) ans <- t(coeff) %*% prec %*% coeff
        return(ans)
        returnType(double(2))
    }
)



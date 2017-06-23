# used in altParams for dmnorm
# this needs to be sourced after nimbleFunction() is defined, so can't be done in distributions_inputList.R
calc_dmnormAltParams <- nimbleFunction(
    name = 'calc_dmnormAltParams',
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
    name = 'calc_dmnormConjugacyContributions',
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
            if(total < d * 1E-15) {
                if(order == 1) ans <- prec %*% asCol(vec)
                if(order == 2) ans <- prec
                return(ans)
            }
        }
        if(order == 1) ans <- t(coeff) %*% (prec %*% asCol(vec))
        if(order == 2) ans <- t(coeff) %*% prec %*% coeff
        return(ans)
        returnType(double(2))
    }
)

# used in altParams for dwish and dinvwish
# this needs to be sourced after nimbleFunction() is defined, so can't be done in distributions_inputList.R
calc_dwishAltParams <- nimbleFunction(
    name = 'calc_dwishAltParams',
    run = function(cholesky = double(2), scale_param = double(), return_scale = double()) {
        if(scale_param == return_scale) {
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



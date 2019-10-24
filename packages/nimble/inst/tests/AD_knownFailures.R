## try to keep in alphabetical order
AD_knownFailures <- list( ## not complete:
  "NUMDERIV_NUMERICAL_ERRORS" = list(
    ## There are many failures not documented here due to known issues with
    ## `numDeriv` numerical accuracy for hessians. See
    ## https://github.com/perrydv/nimbleCoreTeam/issues/149
  ),
  ##
  ## name the elements by operator name
  ##
  "%%" = list(
    '*' = list( ## * means all inputs fail
      compilation = TRUE
    )
  ),
  abs = list(
    "*" = list(
      input = list(
        cpp = 'negative values',
        R = NULL ## if the R version had errors too, we could say so here
      )
    )
  ),
  ceil = list(
    "*" = list(
      compilation = TRUE
    )
  ),
  ## FIXED BY MORE ACCURATE NUMERIC HESSIANS:
  ## exp = list(
  ##   ##
  ##   ## arg type when failure occured (copied from test param name)
  ##   ##
  ##   "arg1 = double(2, c(3, 4))" = list(
  ##     ##
  ##     ## first part of output that led to failure (value, jacobian, or hessian)
  ##     ##
  ##     ## numerical precision error
  ##     hessian = 'method1' ## give name of first method that led to failure
  ##   )
  ## ),
  floor = list(
    "*" = list(
      compilation = TRUE
    )
  ),
  ftrunc = list(
    "*" = list(
      compilation = TRUE
    )
  ),
  ## 'R implements lgamma as "the natural logarithm of _the absolute value of_
  ## the gamma function", and since we use exp(lgamma(x)) for gammafn(x) we get
  ## the sign of the output flipped when gamma(x) < 0.'
  factorial = list(
    "*" = list(
      input = list( ## currently doesn't do anything, just a reminder
        cpp = 'any input x where gamma(x) < 0, e.g. in (-1, 0)'
      )
    )
  ),
  gammafn = list(
    "*" = list(
      ## See entry for factorial.
    )
  ),
  inprod = list(
    "arg1 = double(1, 4) arg2 = double(1, 4)" = list(
      compilation = TRUE ## specify a compilation failure
    )
  ),
  nimRound = list(
    "*" = list(
      compilation = TRUE
    )
  ),
  var = list(
    'arg1 = double(2, c(3, 4))' = list(
      ## nimble doesn't support var of matrices, see sizeUnaryReduction
      compilation = TRUE
    )
  ),
  trace = list(
    "*" = list(
      compilation = TRUE
    )
  ),
  ### distribution tests
  dbeta = list(
    'x = double(1, 5) shape1 = double(0) shape2 = double(0) log = double(0)' = list(
      segfault = TRUE
    ),
    'x = double(0) shape1 = double(1, 3) shape2 = double(0) log = double(0)' = list(
      segfault = TRUE
    ),
    'x = double(1, 5) shape1 = double(1, 3) shape2 = double(0) log = double(0)' = list(
      segfault = TRUE
    ),
    'x = double(0) shape1 = double(0) shape2 = double(1, 7) log = double(0)' = list(
      segfault = TRUE
    ),
    'x = double(1, 5) shape1 = double(0) shape2 = double(1, 7) log = double(0)' = list(
      segfault = TRUE
    ),
    'x = double(0) shape1 = double(1, 3) shape2 = double(1, 7) log = double(0)' = list(
      segfault = TRUE
    ),
    'x = double(1, 5) shape1 = double(1, 3) shape2 = double(1, 7) log = double(0)' = list(
      segfault = TRUE
    )
  ),
  dcat = list(
    "*" = list(
      compilation = TRUE
    )
  ),
  ddexp = list(
    'x = double(0) location = double(0) scale = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) location = double(0) scale = double(0) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) location = double(1, 3) scale = double(0)' = list(
      ## compiled and uncompiled output lengths differ
      value = 'method1'
    ),
    'x = double(0) location = double(0) scale = double(1, 7)' = list(
      ## compiled and uncompiled output lengths differ
      value = 'method1'
    ),
    'x = double(1, 5) location = double(0) scale = double(1, 7)' = list(
      ## compiled and uncompiled output lengths differ
      value = 'method1'
    ),
    'x = double(0) location = double(1, 3) scale = double(1, 7)' = list(
      ## compiled and uncompiled output lengths differ
      value = 'method1'
    ),
    'x = double(1, 5) location = double(1, 3) scale = double(1, 7)' = list(
      ## compiled and uncompiled output lengths differ
      value = 'method1'
    ),
    'x = double(0) location = double(0) rate = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) location = double(0) rate = double(0) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) location = double(1, 3) rate = double(0)' = list(
      ## compiled and uncompiled output lengths differ
      value = 'method1'
    ),
    'x = double(0) location = double(0) rate = double(1, 7)' = list(
      ## Error : Number of dimensions 0 of the return() argument does not match number 1 given in the returnType() statement.
      ## This occurred for: return(out)
      ## This was part of the call:  {
      ##   out <- ddexp(x=ARG1_x_,location=ARG2_location_,scale=1,log=1)
      ##   return(out)
      ## }
      compilation = TRUE
    ),
    'x = double(0) location = double(0) rate = double(1, 7) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) location = double(0) rate = double(1, 7)' = list(
      ## compiled and uncompiled output lengths differ
      value = 'method1'
    ),
    'x = double(0) location = double(1, 3) rate = double(1, 7)' = list(
      ## compiled and uncompiled output lengths differ
      value = 'method1'
    ),
    'x = double(1, 5) location = double(1, 3) rate = double(1, 7)' = list(
      ## compiled and uncompiled output lengths differ
      value = 'method1'
    )
  ),
  ddirch = list(
    'x = double(1, 5) alpha = double(1, 5)' = list(
      ## Error : Number of dimensions 0 of the return() argument does not match number 1 given in the returnType() statement.
      ## This occurred for: return(out)
      ## This was part of the call:  {
      ##   out <- nimArr_ddirch(x=ARG1_x_,alpha=ARG2_alpha_,log=1)
      ##   return(out)
      ## }
      compilation = TRUE
    ),
    'x = double(1, 5) alpha = double(1, 5) log = double(0)' = list(
      compilation = TRUE
    )
  ),
  dexp = list(
    'x = double(0) rate = double(1, 3)' = list(
      compilation = TRUE
    )
  ),
  dexp_nimble = list(
    'x = double(0) rate = double(1, 3)' = list(
      ## compiled and uncompiled output lengths differ
      value = 'method1'
    ),
    'x = double(0) scale = double(1, 3)' = list(
      ## compiled and uncompiled output lengths differ
      value = 'method1'
    )
  ),
  dinvgamma = list(
    'x = double(0) scale = double(0) shape = double(0) log = double(0)' = list(
      segfault = TRUE
    ),
    'x = double(0) scale = double(1, 3) shape = double(0)' = list(
      ## compiled and uncompiled output lengths differ
      value = 'method1'
    ),
    'x = double(0) scale = double(1, 3) shape = double(0) log = double(0)' = list(
      segfault = TRUE
    ),
    'x = double(0) scale = double(0) shape = double(1, 7)' = list(
      ## compiled and uncompiled output lengths differ
      value = 'method1'
    ),
    'x = double(0) scale = double(0) shape = double(1, 7) log = double(0)' = list(
      segfault = TRUE
    ),
    'x = double(0) scale = double(1, 3) shape = double(1, 7)' = list(
      ## compiled and uncompiled output lengths differ
      value = 'method1'
    ),
    'x = double(0) scale = double(1, 3) shape = double(1, 7) log = double(0)' = list(
      segfault = TRUE
    ),
    'x = double(0) rate = double(0) shape = double(0) log = double(0)' = list(
      segfault = TRUE
    ),
    'x = double(0) rate = double(1, 3) shape = double(0)' = list(
      ## compiled and uncompiled output lengths differ
      value = 'method1'
    ),
    'x = double(0) rate = double(1, 3) shape = double(0) log = double(0)' = list(
      segfault = TRUE
    ),
    'x = double(0) rate = double(0) shape = double(1, 7)' = list(
      ## compiled and uncompiled output lengths differ
      value = 'method1'
    ),
    'x = double(0) rate = double(0) shape = double(1, 7) log = double(0)' = list(
      segfault = TRUE
    ),
    'x = double(0) rate = double(1, 3) shape = double(1, 7)' = list(
      ## compiled and uncompiled output lengths differ
      value = 'method1'
    ),
    'x = double(0) rate = double(1, 3) shape = double(1, 7) log = double(0)' = list(
      segfault = TRUE
    )
  ),
  dlnorm = list(
    ## NaN's produced by uncompiled output vs non-sense from compiled output
    'x = double(1, 5) meanlog = double(1, 3) sdlog = double(0)' = list(
      jacobian = 'method2' ## wrt 'x'
    ),
    'x = double(1, 5) meanlog = double(1, 3) sdlog = double(1, 7)' = list(
      jacobian = 'method2' ## wrt 'x'
    )
  ),
  dlogis = list(
    'x = double(1, 5) location = double(0) scale = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) location = double(0) scale = double(0) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) location = double(1, 3) scale = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) location = double(1, 3) scale = double(0) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) location = double(1, 3) scale = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) location = double(1, 3) scale = double(0) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) location = double(0) scale = double(1, 7)' = list(
      compilation = TRUE
    ),
    'x = double(0) location = double(0) scale = double(1, 7) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) location = double(0) scale = double(1, 7)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) location = double(0) scale = double(1, 7) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) location = double(1, 3) scale = double(1, 7)' = list(
      compilation = TRUE
    ),
    'x = double(0) location = double(1, 3) scale = double(1, 7) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) location = double(1, 3) scale = double(1, 7)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) location = double(1, 3) scale = double(1, 7) log = double(0)' = list(
      compilation = TRUE
    )
  ),
  dmnorm_chol = list(
    ## Error : Number of dimensions 0 of the return() argument does not match number 1 given in the returnType() statement.
    ## This occurred for: return(out)
    ## This was part of the call:  {
    ##   out <- nimArr_dmnorm_chol(x=ARG1_x_,mean=ARG3_mean_,cholesky=ARG2_cholesky_,prec_param=TRUE,log=1,overwrite_inputs=0)
    ##   return(out)
    ## }
    'x = double(1, 5) cholesky = double(2, c(5, 5)) mean = double(1, 5)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) cholesky = double(2, c(5, 5)) mean = double(1, 5) log = double(0)' = list(
      compilation = TRUE
    )
  ),
  dmulti = list(
    '*' = list(
      compilation = TRUE
    )
  ),
  dmvt_chol = list(
    ## Error : Number of dimensions 0 of the return() argument does not match number 1 given in the returnType() statement.
    ## This occurred for: return(out)
    ## This was part of the call:  {
    ##   out <- nimArr_dmvt_chol(x=ARG1_x_,mu=ARG3_mu_,cholesky=ARG2_cholesky_,df=ARG4_df_,prec_param=TRUE,log=1,overwrite_inputs=0)
    ##   return(out)
    ## }
    'x = double(1, 5) cholesky = double(2, c(5, 5)) mu = double(1, 5) df = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) cholesky = double(2, c(5, 5)) mu = double(1, 5) df = double(0) log = double(0)' = list(
      compilation = TRUE
    )
  ),
  dwish_chol = list(
    ## Error : Number of dimensions 0 of the return() argument does not match number 1 given in the returnType() statement.
    ## This occurred for: return(out)
    ## This was part of the call:  {
    ##   out <- nimArr_dwish_chol(x=ARG1_x_,cholesky=ARG2_cholesky_,df=ARG3_df_,scale_param=TRUE,log=1,overwrite_inputs=0)
    ##   return(out)
    ## }
    'x = double(1, 5) cholesky = double(2, c(5, 5)) df = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) cholesky = double(2, c(5, 5)) df = double(0) log = double(0)' = list(
      compilation = TRUE
    )
  ),
  dsqrtinvgamma = list(
    "*" = list(
      compilation = TRUE
    )
  ),
  ## dt fails on all inputs 'x = double(0) df = double(0)' because of missing
  ## ncp argument
  dt = list(
    'x = double(1, 5) df = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(0) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) df = double(1, 3)' = list(
      compilation = TRUE
    ),
    'x = double(0) df = double(1, 3) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(1, 3)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(1, 3) log = double(0)' = list(
      compilation = TRUE
    )
  ),
  ## dt_nonstandard fails on all inputs but 'x = double(0) df = double(0) mu =
  ## double(0) sigma = double(0)' because it has 3 parameters, but there is no
  ## MAKE_RECYCLING_RULE_CLASS4_1scalar macro yet
  dt_nonstandard = list(
    'x = double(1, 5) df = double(0) mu = double(0) sigma = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(0) mu = double(0) sigma = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) df = double(1, 3) mu = double(0) sigma = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(1, 3) mu = double(0) sigma = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) df = double(0) mu = double(1, 3) sigma = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(0) mu = double(1, 3) sigma = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) df = double(1, 3) mu = double(1, 3) sigma = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(1, 3) mu = double(1, 3) sigma = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) df = double(0) mu = double(0) sigma = double(1, 7)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(0) mu = double(0) sigma = double(1, 7)' = list(
      compilation = TRUE
    ),
    'x = double(0) df = double(1, 3) mu = double(0) sigma = double(1, 7)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(1, 3) mu = double(0) sigma = double(1, 7)' = list(
      compilation = TRUE
    ),
    'x = double(0) df = double(0) mu = double(1, 3) sigma = double(1, 7)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(0) mu = double(1, 3) sigma = double(1, 7)' = list(
      compilation = TRUE
    ),
    'x = double(0) df = double(1, 3) mu = double(1, 3) sigma = double(1, 7)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(1, 3) mu = double(1, 3) sigma = double(1, 7)' = list(
      compilation = TRUE
    ),
    ## log
    'x = double(1, 5) df = double(0) mu = double(0) sigma = double(0) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(0) mu = double(0) sigma = double(0) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) df = double(1, 3) mu = double(0) sigma = double(0) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(1, 3) mu = double(0) sigma = double(0) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) df = double(0) mu = double(1, 3) sigma = double(0) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(0) mu = double(1, 3) sigma = double(0) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) df = double(1, 3) mu = double(1, 3) sigma = double(0) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(1, 3) mu = double(1, 3) sigma = double(0) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) df = double(0) mu = double(0) sigma = double(1, 7) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(0) mu = double(0) sigma = double(1, 7) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) df = double(1, 3) mu = double(0) sigma = double(1, 7) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(1, 3) mu = double(0) sigma = double(1, 7) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) df = double(0) mu = double(1, 3) sigma = double(1, 7) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(0) mu = double(1, 3) sigma = double(1, 7) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(0) df = double(1, 3) mu = double(1, 3) sigma = double(1, 7) log = double(0)' = list(
      compilation = TRUE
    ),
    'x = double(1, 5) df = double(1, 3) mu = double(1, 3) sigma = double(1, 7) log = double(0)' = list(
      compilation = TRUE
    )
  ),
  dunif = list(
    'x = double(1, 5) min = double(1, 3) max = double(1, 7)' = list(
      value = 'method1'
    )
  ),
  dweibull = list(
    'x = double(0) scale = double(1, 3) shape = double(1, 7)' = list(
      value = 'method1'
    )
  )
)

## This test file also serves as a demo for experimental features: calling externally compiled code, and calling R functions from compiled nimbleFunctions.
## It covers rough draft functionality for calling external compiled functions and calling arbitrary R functions from within compiled nimbleFunction (including BUGS) code.

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

## Calling external C:

## Say we want a C function that adds 1.5 to a vector of values:
## We can only give non-scalars results by pointer argument.
## Any NIMBLE non-scalar can become a double* of contiguously allocated memory.
## Like in LAPACK or similar low-level libraries, you need a separate argument to say how long the allocated memory is.
sink('add1p5.h')
cat('
extern "C" {
  void my_internal_function(double *p, double*ans, int n);
}
')
sink()

sink('add1p5.cpp') 
cat('
#include <cstdio>
#include "add1p5.h"

void my_internal_function(double *p, double *ans, int n) {
  printf("In my_internal_function\\n"); /* cat reduces the double slash to single slash */ 
  for(int i = 0; i < n; i++) 
    ans[i] = p[i] + 1.5;
}
')
sink()

## A major limitation is that we cannot have compiled code return anything except a scalar.  So returned non-scalars will have to be via arguments.  That meansto use this in BUGS code we'll have to wrap it in another nimbleFunction.  It will make sense below.

## compile that to .o
## I am not an expert on all compilation twists, but on my system I need to do it with g++ instead of gcc (even though it is pure C) in order for the linker to be later happy linking it to compiled C++ from NIMBLE.  
system('g++ add1p5.cpp -c -o add1p5.o')

## now to create a nimbleFunction that interfaces to add1p5:

Radd1p5 <- nimbleExternalCall(function(x = double(1), ans = double(1), n = integer()){}, Cfun = 'my_internal_function', headerFile = 'add1p5.h', oFile = 'add1p5.o')
## The first argument uses the format of a nimbleFunction with argument type declarations, but the body of the function can be empty ('{}')
## You can choose any names for the arguments in R. They don't have to match the C code.
## Ignore the warning here and later warnings
##Radd1p5 ## this is a nimbleFunction with some internal special sauce, but from here on out it should behave like a nimbleFunction
## We'll wait to use it in BUGS code

temporarilyAssignInGlobalEnv(Radd1p5)

## Radd1p5 doesn't return anything and would only be able to return a scalar, so we can wrap it like this:
wrappedRadd1p5 <- nimbleFunction(
    run = function(x = double(1)) {
        ans <- numeric(length(x))
        Radd1p5(x, ans, length(x))
        return(ans)
        returnType(double(1))
    })

temporarilyAssignInGlobalEnv(wrappedRadd1p5)


## calling an R function
## Say we want an R function that adds 2.5 to every value in a vector
add2p5 <- function(x) {
    x + 2.5 ## This can be pure R
}
temporarilyAssignInGlobalEnv(add2p5)

## now create a nimbleFunction that interfaces to add2p5, not necessary when running uncompiled but necessary when running compiled
Radd2p5 <- nimbleRcall(function(x = double(1)){}, Rfun = 'add2p5', returnType = double(1), envir = .GlobalEnv)
## Similar to above.  The function prototype and the returnType represent a promise that add2p5 will always take and return these types.  Also an environment is needed and it defaults to R's global environment but I'm showing it explicitly.
## Again, ignore the more extensive warnings
##Radd2p5 ## this is also a nimbleFunction with more special sauce
temporarilyAssignInGlobalEnv(Radd2p5)

## Now let's use these in a model

demoCode <- nimbleCode({
    for(i in 1:4) {x[i] ~ dnorm(0,1)} ## just to get a vector
    y[1:4] <- wrappedRadd1p5(x[1:4])
    z[1:4] <- Radd2p5(x[1:4])
    })

demoModel <- nimbleModel(demoCode, inits = list(x = rnorm(4)), check = FALSE, calculate = FALSE)
## Again ignore the error during checking.  We'll have to trap and handle that, but right now I'm focused on core functionality.
## The model will not work uncompiled!

CdemoModel <- compileNimble(demoModel, dirName = '.') ## last arg puts the C++ code in your working directory so you can look at it if you like

CdemoModel$x
CdemoModel$calculate()
expect_equal(CdemoModel$y - CdemoModel$x, rep(1.5, 4), info = 'external C call works in compiled BUGS code')
expect_equal(CdemoModel$z - CdemoModel$x, rep(2.5, 4), info = 'external R call works in compiled BUGS code')

## next test alternate modes of providing arguments

sink('add1p5b.h')
cat('
extern "C" {
  double my_internal_function(double *p, double*ans, int n);
}
')
sink()

sink('add1p5b.cpp') 
cat('
#include <cstdio>
#include "add1p5b.h"

double my_internal_function(double *p, double *ans, int n) {
  printf("In my_internal_function\\n"); /* cat reduces the double slash to single slash */ 
  for(int i = 0; i < n; i++) 
    ans[i] = p[i] + 1.5;
  return(1.23);
}
')
sink()

system('g++ add1p5b.cpp -c -o add1p5b.o')
options(error = recover)
Radd1p5b <- nimbleExternalCall(list(nimbleType(name = 'x', type = 'double', dim = 1),
                                    nimbleType(name = 'ans', type = 'double', dim = 1),
                                    nimbleType(name = 'n', type = 'integer', dim = 0)),
                               Cfun = 'my_internal_function', headerFile = 'add1p5b.h', oFile = 'add1p5b.o',
                               returnType = nimbleType(name = 'NOTUSED', type = 'double', dim = 0))
temporarilyAssignInGlobalEnv(Radd1p5b)
Cadd1p5b <- compileNimble(Radd1p5b, dirName = '.')

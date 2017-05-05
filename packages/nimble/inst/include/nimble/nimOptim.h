#ifndef __NIMBLE_NIMOPTIM_H
#define __NIMBLE_NIMOPTIM_H

#include <nimble/NimArr.h>
#include <nimble/optimTypes.h>
#include <nimble/smartPtrs.h>

typedef double NimObjectiveFn(NimArr<1, double>& par);

// This behaves like R's native optim() function.
nimSmartPtr<OptimResultNimbleList> nimOptim(NimArr<1, double> &par,
                                            NimObjectiveFn fn);

#endif // __NIMBLE_NIMOPTIM_H

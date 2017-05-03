#ifndef __NIMBLE_NIMOPTIM_H
#define __NIMBLE_NIMOPTIM_H

#include <nimble/NimArr.h>
#include <nimble/optimTypes.h>
#include <nimble/smartPtrs.h>

typedef double NimObjectiveFn(const NimArr<1, double>& par);

nimSmartPtr<OptimResultNimbleList> nimFakeOptim(const NimArr<1, double> &par,
                                                void *fn);

#endif // __NIMBLE_NIMOPTIM_H

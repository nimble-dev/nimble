#include "nimble/Values.h"
#include "nimble/RcppUtils.h"
#include "nimble/Utils.h"

// ValuesFactory valuesFactory;

//int Values::getsize() { return SampleSize; }

// Values *ValuesFactory::makeNew(string &typelabel) {
//   Values *ans = (*ValuesGeneratorMap[typelabel])();
//   ValuesCollection.insert(ans);
//   return(ans);
// }

// void ValuesFactory::registerGenerator(string typelabel, ValuesGenerator vg) {
//   ValuesGeneratorMap[typelabel] = vg;
// }

// void ValuesFactory::removeValues(Values *empty) {
//   ValuesCollection.erase(empty);
// }

// ValuesFactory::~ValuesFactory() {
//   tr1::unordered_set<Values *>::iterator iVC;
//   for(iVC = ValuesCollection.begin(); iVC != ValuesCollection.end(); ++iVC) {
//     delete *iVC;
//   }
// }

// void valuesFinalizer(SEXP Sv) {
//   PRINTF("Calling valuesFinalizer\n");
//   Values *v = static_cast<Values *>(R_ExternalPtrAddr(Sv));
//   delete v;
//   valuesFactory.removeValues(v);
// }

// This uses the valuesFactory method for building a new object
// Current interface functions use getMVBuildname function
// SEXP newModelValues(SEXP Stypelabel) {
//   if(!checkString(Stypelabel,1)) RBREAK("Error, called newModelValues with an argumen that is not a string of length >= 1");
//   string typelabel = STRSEXP_2_string(Stypelabel, 0);
//   Values *ans = valuesFactory.makeNew(typelabel);
//   SEXP Sans = R_MakeExternalPtr(ans, R_NilValue, R_NilValue);
//   R_RegisterCFinalizerEx(Sans, &valuesFinalizer, TRUE); // TRUE is from enum of Rboolean
//   return(Sans);
// }


#ifndef __MODEL
#define __MODEL

#include "NamedObjects.h"
/* class GlobalObjects : public NamedObjects { */
/* public: */
/*   GlobalObjects(); */
/* }; */

//extern GlobalObjects globalObjects;

class Model : public NamedObjects{
public:
  int foo; // so it is not empty :). We might need stuff here later.
};

#endif

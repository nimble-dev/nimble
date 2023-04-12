#include <nimble/nimIntegrate.h>
#include <R_ext/Applic.h>
#include <nimble/nimIntegrate.h>
#include <string.h>
#include <algorithm>
#include <cmath>
#include <limits>

// typedef taken from "Writing R Extensions"
typedef void integr_fn(double *x, int n, void *ex);

// some_actual_integrator is a toy to be replaced by
// a call to Rdqags or Rdqagi
double some_actual_integrator(integr_fn my_fn, // my_fn will be NimIntegrateProblem::fn
                              void *ex,
                              double lower,
                              double upper) {
  int n = 1;
  std::vector<double> vals(n);
  vals[0] = (upper - lower)/2.;
  my_fn( &vals[0], n, ex);
  double sum = 0.;
  for(int i = 0; i < n; ++i) sum += vals[i];
  return sum;
}

void NimIntegrateProblem::fn(double *x, int n, void *ex) {
  NimIntegrateProblem* problem = static_cast<NimIntegrateProblem*>(ex);
  problem->par_.setSize(n, false, false);
  std::copy(x, x + n, problem->par_.getPtr());
  problem->return_vals_.setSize(n, false, false);
  problem->return_vals_ = problem->function(); // problem->function calls the actual fn provided by user
  std::copy(problem->return_vals_.getPtr(),
            problem->return_vals_.getPtr() + n,
            x);
}

double NimIntegrateProblem::integrate() {
  void* ex = this;
  double result = some_actual_integrator(NimIntegrateProblem::fn,
                                         ex,
                                         lower_, upper_);
  return result;
}

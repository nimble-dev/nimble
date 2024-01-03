/*
 * NIMBLE: an R package for programming with BUGS models.
 * Copyright (C) 2014-2017 Perry de Valpine, Christopher Paciorek,
 * Daniel Turek, Clifford Anderson-Bergman, Nick Michaud, Fritz Obermeyer,
 * Duncan Temple Lang.
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, a copy is available at
 * https://www.R-project.org/Licenses/
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <float.h>
#include <Rmath.h> /* for fmax2, fmin2, imin2 */

#include <nimble/nimIntegrate.h>
#include <string.h>
#include <algorithm>
#include <cmath>
#include <limits>

// typedef taken from "Writing R Extensions"
typedef void integr_fn(double *x, int n, void *ex);


void NimIntegrateProblem::fn(double *x, int n, void *ex) {
  NimIntegrateProblem* problem = static_cast<NimIntegrateProblem*>(ex);
  problem->x_.setSize(n, false, false);
  std::copy(x, x + n, problem->x_.getPtr());
  problem->return_vals_.setSize(n, false, false);
  problem->return_vals_ = problem->function(); // problem->function calls the actual fn provided by user
  std::copy(problem->return_vals_.getPtr(),
            problem->return_vals_.getPtr() + n,
            x);
}

double NimIntegrateProblem::integrate() {
  void* ex = this;
  result = 0.0;

  // Not sure if the allocation should be here or in constructor.
  //int* iwork = new int[subdivisions_];
  //double* work = new double[lenw];

  if (std::isfinite(lower_) && std::isfinite(upper_)) {
    // Note that `abs_tol` is before `rel_tol` in `Rqdag{s,i}` signature.
    Rdqags(NimIntegrateProblem::fn, ex, &lower_, &upper_,
           &abs_tol_, &rel_tol_, &result, &abserr, &neval, &ierr,
           &subdivisions_, &lenw, &last, iwork, work); 
  } else {
    if (std::isfinite(lower_)) {
      inf = 1;
    } else if (std::isfinite(upper_)) {
      inf = -1;
      lower_ = upper_;   // `bound` in R's `integrate`.
    } else {
      inf = 2;
      lower_ = 0.0;
    }
    Rdqagi(NimIntegrateProblem::fn, ex, &lower_, &inf,
           &abs_tol_, &rel_tol_, &result, &abserr, &neval, &ierr,
           &subdivisions_, &lenw, &last, iwork, work); 
  }
 
  //delete [] iwork;
  //delete [] work;
  return result;
}


NimIntegrateProblem::~NimIntegrateProblem() {
  delete [] work;
  delete [] iwork;
}

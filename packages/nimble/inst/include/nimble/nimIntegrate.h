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

#ifndef NIMINTEGRATE_H_
#define NIMINTEGRATE_H_

#include "nimOptim.h" // has NimBound stuff and other headers needed
// CJP note: shouldn't we just #include the specific files we need,
// namely NimArr.h and perhaps nothing else?

class NimIntegrateProblem {
public:
  NimIntegrateProblem(NimArr<1, double>& param, double lower, double upper, int subdivisions,
                      double rel_tol, double abs_tol, bool stop_on_error) :
    param_(param),
    lower_(lower),
        upper_(upper),
        subdivisions_(subdivisions),
        rel_tol_(rel_tol),
        abs_tol_(abs_tol),
        stop_on_error_(stop_on_error) {};

    double integrate();
 private:
    static void fn(double *x, int n, void *ex);

 protected:
    virtual NimArr<1, double> function()=0;

    double lower_;
    double upper_;
    int subdivisions_;
    double rel_tol_;
    double abs_tol_;
    bool stop_on_error_;
    NimArr<1, double>& param_;
    NimArr<1, double> x_;
    NimArr<1, double> return_vals_;
};

template<class Fn>
class NimIntegrateProblem_Fun : public NimIntegrateProblem {
  public:
  NimIntegrateProblem_Fun(Fn fn, NimArr<1, double>& param, double lower, double upper, int subdivisions, double rel_tol, double abs_tol, bool stop_on_error)
    : NimIntegrateProblem(param, lower, upper, subdivisions, rel_tol, abs_tol, stop_on_error), fn_(fn) {}

  protected:
  virtual NimArr<1, double> function() {return fn_(x_, param_);}

  private:
  Fn fn_;
};

template <class Fn>
inline double nimIntegrate(
    Fn fn,
    NimArr<1, double>& param,
    double lower,
    double upper,
    int subdivisions,
    double rel_tol,
    double abs_tol,
  bool stop_on_error) {

  return NimIntegrateProblem_Fun<Fn>(fn, param, lower, upper, subdivisions, rel_tol, abs_tol, stop_on_error).integrate();
}

#endif // NIMINTEGRATE_H_

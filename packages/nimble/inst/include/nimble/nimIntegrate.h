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

class NimIntegrateProblem {
public:
    NimIntegrateProblem(double lower, double upper) :
        lower_(lower),
        upper_(upper) {};

    double integrate();
 private:
    static void fn(double *x, int n, void *ex);

 protected:
    virtual NimArr<1, double> function()=0;

    double lower_;
    double upper_;
    NimArr<1, double> par_;
    NimArr<1, double> return_vals_;
};

template<class Fn>
class NimIntegrateProblem_Fun : public NimIntegrateProblem {
  public:
  NimIntegrateProblem_Fun(Fn fn, double lower, double upper)
    : NimIntegrateProblem(lower, upper), fn_(fn) {}

  protected:
  virtual NimArr<1, double> function() {return fn_(par_);}

  private:
  Fn fn_;
};

template <class Fn>
inline double nimIntegrate(
    Fn fn,
    double lower,
    double upper) {
    return NimIntegrateProblem_Fun<Fn>(fn, lower, upper).integrate();
}

#endif // NIMINTEGRATE_H_

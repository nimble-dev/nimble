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

#ifndef __EIGENUSINGCLASSES
#define __EIGENUSINGCLASSES
#include "NimArr.h"
#include "NamedObjects.h"
#include "smartPtrs.h"


class EIGEN_EIGENCLASS_R : public pointedToBase {
public:
  NimArr<1, double> values; //these will be defined in nimble.so, and nimble.a and again in on-the-fly compilation
  NimArr<2, double> vectors;
  NimArr<1, double> &getValues() {return(values);}
  NimArr<2, double> &getVectors() {return(vectors);}
  SEXP RObjectPointer;
  
  virtual SEXP copyToSEXP (   );
  void  createNewSEXP (  );
  void  copyFromSEXP ( SEXP S_nimList_ );
  EIGEN_EIGENCLASS_R(){	
    RObjectPointer = NULL;
  };
};


class EIGEN_SVDCLASS_R : public pointedToBase {
 public:
  NimArr<1, double> d;
  NimArr<2, double> u;
  NimArr<2, double> v;
  NimArr<1, double> &getD() {return(d);}
  NimArr<2, double> &getU() {return(u);}
  NimArr<2, double> &getV() {return(v);}
  SEXP RObjectPointer;
  
  virtual SEXP  copyToSEXP (   );
  void  createNewSEXP (  );
  void  copyFromSEXP ( SEXP S_nimList_ );
  EIGEN_SVDCLASS_R (  ) {
    RObjectPointer = NULL;
  };
};

extern "C" {
SEXP C_nimEigen(SEXP S_x, SEXP S_symmetric, SEXP S_valuesOnly, SEXP returnList);
SEXP C_nimSvd(SEXP S_x, SEXP S_vectors, SEXP returnList);
}


#endif

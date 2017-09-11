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

#ifndef __MODELCLASSUTILS
#define __MODELCLASSUTILS

#include "NamedObjects.h"
#include "RcppNimbleUtils.h"
#include "Values.h"

class ModelBase : public NamedObjects{
	public:
	Values* modelValues_;
	virtual Values* getModelValuesPtr(){ return modelValues_; }
	};

extern "C" {
  SEXP getModelValuesPtrFromModel (SEXP rPtr); // gets pointer to modelvalues object in model
  SEXP getModelElementPtr(SEXP Sextptr, SEXP Sname); // Gets the ptr to an element of name Sname from
  // the ModelValues object pointed to by Sextptr
  
  SEXP getMVBuildName(SEXP rPtr);			// gets character string to feed to .Call 
  // that builds a new ModelValue of the same type as
  // is pointed to by rPtr
  
  SEXP derefPtr(SEXP SmultiPtr);
 }


NimArrType** cGetModelElementPtr(SEXP Sextptr, SEXP Sname);	// Gets the ptr to an element of name Sname from
// the ModelValues object pointed to by Sextptr

#endif

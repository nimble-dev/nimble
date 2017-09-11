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

#include <iostream>
#include <nimble/smartPtrs.h>
#include <nimble/dllFinalizer.h>

//#define _DEBUG_SMARTPTR_FINALIZER

void pointedToBase_Finalizer(SEXP Snp){
  // std::cout<< "In pointedToBase_Finalizer\n";
  pointedToBase* np = static_cast<pointedToBase*>(R_ExternalPtrAddr(Snp));
  if(np) {
    np->removeWatcher(); /* object will naturally self-destruct if watcher count goes to 0*/
  }
  R_ClearExternalPtr(Snp);
}

SEXP register_pointedToBase_Finalizer(SEXP Snp, SEXP Dll, SEXP Slabel) {
  // std::cout<< "In register_pointedToBase_Finalizer\n";
  RegisterNimbleFinalizer(Snp, Dll, &pointedToBase_Finalizer, Slabel);
  return(Snp);
}

void smartPtrBase_Finalizer(SEXP Snp){
#ifdef _DEBUG_SMARTPTR_FINALIZER
  std::cout<< "In smartPtrBase_Finalizer\n";
#endif
  nimSmartPtrBase* np = static_cast<nimSmartPtrBase*>(R_ExternalPtrAddr(Snp));
  if(np) {
    delete np; /* pointed to object will naturally self-destruct if this decrements its watcher count goes to 0*/
  }
  R_ClearExternalPtr(Snp);
}

SEXP register_smartPtrBase_Finalizer(SEXP Snp, SEXP Dll, SEXP Slabel) {
#ifdef _DEBUG_SMARTPTR_FINALIZER
  std::cout<< "In register_smartPtrBase_Finalizer\n";
#endif
  RegisterNimbleFinalizer(Snp, Dll, &smartPtrBase_Finalizer, Slabel);
  return(Snp);
}

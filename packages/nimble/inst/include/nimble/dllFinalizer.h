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

#ifndef __DLLFINALIZER
#define __DLLFINALIZER

// modeled after DTL ideas in dll_auto_unload branch

#include "Utils.h"
#include <Rinternals.h>

void RegisterNimblePointer(SEXP ptr, SEXP Dll, R_CFinalizer_t finalizer, SEXP Slabel);
extern "C" {
  SEXP RNimble_Ptr_ManualFinalizer(SEXP obj);
  SEXP RNimble_Ptr_CheckAndRunAllDllFinalizers(SEXP Dll, SEXP Sforce);
  SEXP CountDllObjects();
}

#define RegisterNimbleFinalizer RegisterNimblePointer

#endif

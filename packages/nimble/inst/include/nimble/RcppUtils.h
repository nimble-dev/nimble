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

#ifndef __RCPPUTILS
#define __RCPPUTILS


#include <string>
#include <vector>
#include<iostream>
#include<sstream>
#include "R.h"
#include "Utils.h"
#include "Rdefines.h"


#include <Rinternals.h>

#include <R_ext/Applic.h>	/* this is required for optim */
#include <stdarg.h> 		/* this is required for variable number of arguments */

extern std::ostringstream _nimble_global_output;

using namespace std;
#define Inf R_PosInf
#define NA 0

class varAndIndicesClass {
public:
  std::string varName;
  std::vector< std::vector<int> > indices;
};

#define SET_NEXT_LANG_ARG(SLANG, ARGNAME) SETCAR(SLANG, ARGNAME); \
  SLANG = CDR(SLANG);

void nimble_print_to_R(std::ostringstream &input);


void multivarTestCall(double *x, int n);

vector<int> getSEXPdims(SEXP Sx);

string STRSEXP_2_string(SEXP Ss, int i = 0);
SEXP   string_2_STRSEXP(string v);
void   STRSEXP_2_vectorString(SEXP Ss, vector<string> &ans);
SEXP   vectorString_2_STRSEXP(const vector<string> &v);

vector<double> SEXP_2_vectorDouble( SEXP Sn ); /* Sn can be numeric or integer from R*/
double SEXP_2_double(SEXP Sn, int i = 0); /* Ditto */
SEXP double_2_SEXP(double v);
SEXP vectorDouble_2_SEXP(const vector<double> &v);
SEXP vectorInt_2_SEXP(const vector<int> &v);
SEXP vectorInt_2_SEXP(const vector<int> &v, int offset);

vector<int> SEXP_2_vectorInt(SEXP Sn, int offset = 0); /* Sn can be numeric or integer from R */ 
/* If Sn is numeric but not integer, a warning is issued if it contains non-integers */
int SEXP_2_int(SEXP Sn, int i = 0);
SEXP int_2_SEXP(int i);
bool SEXP_2_bool(SEXP Sn, int i = 0);
SEXP bool_2_SEXP(bool ind);

void populate_SEXP_2_double_for_copyFromRobject(void *vPtr, SEXP rScalar);
void populate_SEXP_2_int_for_copyFromRobject(void *vPtr, SEXP rScalar);
void populate_SEXP_2_bool_for_copyFromRobject(void *vPtr, SEXP rScalar);

void parseVarAndInds(const std::string &input, varAndIndicesClass &output);
extern "C" {
  SEXP populate_SEXP_2_double(SEXP rPtr, SEXP refNum, SEXP rScalar);
  SEXP extract_double_2_SEXP(SEXP rPtr, SEXP refNum);
  SEXP populate_SEXP_2_bool(SEXP rPtr, SEXP refNum, SEXP rScalar);
  SEXP extract_bool_2_SEXP(SEXP rPtr, SEXP refNum);
  SEXP populate_SEXP_2_int(SEXP rPtr, SEXP refNum, SEXP rScalar);
  SEXP extract_int_2_SEXP(SEXP rPtr, SEXP refNum);

  SEXP populate_SEXP_2_string(SEXP rPtr, SEXP rString);
  SEXP populate_SEXP_2_stringVector(SEXP rPtr, SEXP rStringVector);
  SEXP extract_string_2_SEXP(SEXP rPtr);
  SEXP extract_stringVector_2_SEXP(SEXP rPtr);

  SEXP fastMatrixInsert(SEXP matrixInto, SEXP matrix, SEXP rowStart, SEXP colStart);
  SEXP matrix2ListDouble(SEXP matrix, SEXP list, SEXP listStartIndex, SEXP RnRows,  SEXP dims);
  SEXP matrix2ListInt(SEXP matrix, SEXP list, SEXP listStartIndex, SEXP RnRows,  SEXP dims);

  SEXP C_rankSample(SEXP p, SEXP n, SEXP not_used, SEXP s);

  SEXP parseVar(SEXP Sinput);
  SEXP makeParsedVarList(SEXP Sx);
}

void rawSample(double* p, int c_samps, int N, int* ans, bool unsort, bool silent);

SEXP makeNewNimbleList(SEXP S_listName);

#endif



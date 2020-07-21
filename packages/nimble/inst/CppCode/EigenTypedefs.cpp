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

#include <nimble/EigenTypedefs.h>

SEXP  EIGEN_EIGENCLASS::copyToSEXP (  )  {
  if(!RCopiedFlag) {
    EIGEN_EIGENCLASS_R::copyToSEXP();
    RCopiedFlag = true;
  }
  return(RObjectPointer);
}

void EIGEN_EIGENCLASS::resetFlags () {
  RCopiedFlag = false;	
}

EIGEN_EIGENCLASS::EIGEN_EIGENCLASS(){
  //std::cout<<"Constructing EIGEN_EIGENCLASS\n";
  namedObjects["values"]=&values;
  namedObjects["vectors"]=&vectors;
  RCopiedFlag = false;
}

SEXP  EIGEN_SVDCLASS::copyToSEXP (  )  {
  if(!RCopiedFlag) {
    EIGEN_SVDCLASS_R::copyToSEXP();
    RCopiedFlag = true;
  }
  return(RObjectPointer);
}

EIGEN_SVDCLASS::EIGEN_SVDCLASS(){
  //std::cout<<"Constructing EIGEN_SVDCLASS\n";
  namedObjects["d"]=&d;
  namedObjects["u"]=&u;
  namedObjects["v"]=&v;
  RCopiedFlag = false;
}

void EIGEN_SVDCLASS::resetFlags () {
  RCopiedFlag = false;	
}

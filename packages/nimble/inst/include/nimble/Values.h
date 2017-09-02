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

#ifndef __VALUES
#define __VALUES

//#include<vector>
#include<string>
using std::string;
#include "NamedObjects.h"
#include "Utils.h"

class Values : public NamedObjects {
public:
  int numRows;
  int getsize(){return(numRows);};
  virtual void resize(int nrow)=0;
  
  string buildName;
  string getMVBuildName(){ return buildName;}  ;
  Values() { buildName = "missing"; numRows = 0;};
};

#endif

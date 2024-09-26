/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2012  Bauke Conijn <email>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef SOLVER_H
#define SOLVER_H

#include <fstream>
#include "structures.h"
#include <string>
#include <iostream>
#include <sstream>

extern int rr;
class BaseTask;
namespace Solver {
	typedef int (*solver_function)(BaseTask * t);
	struct type {
		const char * name;
		solver_function function;
	};
	extern type solvers[];
	static bool log=false;


};

#endif // SOLVER_H

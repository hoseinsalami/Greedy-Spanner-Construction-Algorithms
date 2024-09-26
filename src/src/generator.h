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


#ifndef DATAGENERATOR_H
#define DATAGENERATOR_H
#include "structures.h"

namespace Generator {
	typedef pointset (*generator_function)(int seed, int &count, int max_coordinate, double shape, double &t, int source);
	struct type {
		const char * name;
		generator_function function;
	};
	extern type generators[];
};

#endif // DATAGENERATOR_H

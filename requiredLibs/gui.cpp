/*
    Surface viewer - A 3D Mathematical surface inspection & renderer tool
    Copyright (C) 2011  B.J. Conijn <b.j.conijn@student.tue.nl>

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

#include <QtGui/QApplication>
#include "window.h"
#include "commandline.h"

int main(int argc, char **argv) {
	if(argc>=2 && (strcmp(argv[1],"-c")==0 || strcmp(argv[1],"--command")==0)){
		return CommandMain(argc,argv);
	}
	QApplication app(argc, argv);
	app.setApplicationName("Simulation in Computer Graphics - Project 1");
	app.setOrganizationName("TU/e");

	// Create main window
	Window window;

	// Show all
	window.show();

	// Enter mainloop
	return app.exec();
}

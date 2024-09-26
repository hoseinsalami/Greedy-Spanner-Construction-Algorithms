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


#ifndef SURFACEVIEWER_H
#define SURFACEVIEWER_H

#include <QtGui/QMainWindow>
#include "structures.h"

//Basically fixes the broken Qt macro, which interferes with C++11 user-defined literals.
#undef QLOCATION
#define QLOCATION "\0" __FILE__ ":" QTOSTRING(__LINE__)

class Task;
class QThread;
class OptionItem;
class Config;
class GLWidget;

class Window : public QMainWindow {
        Q_OBJECT
public:
	Window();
    ~Window();
	virtual void closeEvent(QCloseEvent *e);

private:
	Config *createConfig(QString name);
	void createMenu();
	void createConfigPanels();
	void createMenuOptions();
	void createContext();
	void loadSettings();

	struct {
		QMenu *file;
		QMenu *window;
		QMenu *help;
	} menu;
	GLWidget *gl;
	QThread *workThread;

private slots:
	void aboutQt();
    void generate();
    void solve();
	void createWorker();
signals:
    void solved(int);
    void newTask(Task*);
    void quitWorkerThread();
    void stopWorkerTask();
};

#endif // SURFACEVIEWER_H

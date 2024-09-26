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

#include <QtCore/QSettings>
#include <QtCore/QDebug>
#include <QtCore/QTimer>
#include <QtCore/QThread>
#include <QtGui/QLabel>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QAction>
#include <QtGui/QDockWidget>
#include <QtGui/QFileDialog>
#include <QtGui/QMessageBox>
#include <QtGui/QPushButton>

#include "window.h"
#include "glwidget.h"
#include "config.h"
#include "generator.h"
#include "solver.h"
#include "worker.h"

namespace {
int gen_alg=0, solv_alg=2;
double t=2.0, var=1.1;
int landmarks=0;
int seed=100,count=1000,maxc=1000;
int source=1;
double shape=2.0;
pointset p;
};

Window::Window() {
	setDockOptions(AllowNestedDocks | AllowTabbedDocks);
    qRegisterMetaType<pointset>("pointset");
    qRegisterMetaType<edgelist>("edgelist");
    qRegisterMetaType<edge>("edge");
    createWorker();
	createMenu();
	createContext();
	createConfigPanels();
	createMenuOptions();
	loadSettings();
}

Window::~Window() {
    // Let the worker exit its eventloop and wait for it to terminate
    emit quitWorkerThread();
    workThread->wait();
}

void Window::closeEvent(QCloseEvent *e) {
	QSettings s;
	s.setValue("geometry", saveGeometry());
	s.setValue("windowState", saveState());
	QWidget::closeEvent(e);
}

void Window::loadSettings() {
	resize(640, 480);
	QSettings s;
	restoreGeometry(s.value("geometry").toByteArray());
	restoreState(s.value("windowState").toByteArray());
}

void Window::createWorker() {
    // Create a thread for our worker.
    workThread = new QThread(this);
    workThread->start();
    // Create the worker.
    Worker * worker = new Worker();
    // Push it to his thread (Note that the thread object itself lives on our thread)
    worker->moveToThread(workThread);
    // Create the message queues
    connect(this, SIGNAL(newTask(Task*)), worker, SLOT(solve(Task*)));
    connect(worker, SIGNAL(solved(int)), this, SIGNAL(solved(int)));
    // Setup normal thread quiting sequence.
    connect(this, SIGNAL(quitWorkerThread()), worker, SLOT(deleteLater()));
    connect(worker, SIGNAL(destroyed(QObject*)), workThread, SLOT(quit()));
    qDebug("Created new Worker Thread.");
}

void Window::generate() {
    p = Generator::generators[gen_alg].function(seed,count,maxc,shape,t,source);
    gl->setPointset(p);
}

void Window::solve() {
    Task * T = new Task(p, solv_alg, t, var, landmarks);
    T->moveToThread(workThread);
    gl->resetEdges();
    gl->setT(t);
    connect(T,SIGNAL(add_edge(edge)),gl,SLOT(addEdge(edge)),Qt::QueuedConnection);
    connect(this,SIGNAL(stopWorkerTask()),T,SLOT(do_quit()),Qt::DirectConnection); 
    emit newTask(T);
}

void Window::createConfigPanels() {
    { 
        Config * c = createConfig("Generator");
        OptionItem * list = new OptionItem(c, "Generator", &gen_alg);
        for (int i=0; Generator::generators[i].name; i++) {
            list->addItem(Generator::generators[i].name);
        }
        IntegerItem * s = new IntegerItem(c, "Seed", &seed);
        s->setRange(0,10000);
        IntegerItem * m = new IntegerItem(c, "Vertices", &count);
        m->setRange(0,1000000);
        m->setSingleStep(25);
        IntegerItem * r = new IntegerItem(c, "Max coordinate", &maxc);
        r->setRange(0,1000000);
        r->setSingleStep(250);
        DoubleItem * sh = new DoubleItem(c, "Shape", &shape);
        sh->setRange(0.01,100);
        sh->setDecimals(2);
        sh->setSingleStep(0.1);
        IntegerItem * sour = new IntegerItem(c, "Source", &source);
        sour->setRange(1,6);
        sour->setSingleStep(1);
        ButtonItem * go = new ButtonItem(c, "Generate");
        connect(go,SIGNAL(clicked(bool)),this,SLOT(generate()));
        c->done();
    }
    { 
        Config * c = createConfig("Solver");
        OptionItem * list = new OptionItem(c, "Solver", &solv_alg);
        for (int i=0; Solver::solvers[i].name; i++) {
            list->addItem(Solver::solvers[i].name);
        }
        DoubleItem * s = new DoubleItem(c, "t", &t);
        s->setRange(1,1000);
        s->setDecimals(4);
        s->setSingleStep(0.1);
        DoubleItem * v = new DoubleItem(c, "var", &var);
        v->setRange(0,2000);
        v->setDecimals(3);
        v->setSingleStep(0.1);
        IntegerItem * l = new IntegerItem(c, "Landmarks", &landmarks);
        l->setRange(0,10000);
        ButtonItem * go = new ButtonItem(c, "Solve");
        connect(go,SIGNAL(clicked(bool)),this,SLOT(solve()));
        connect(go,SIGNAL(clicked(bool)),go,SLOT(disable()));
        connect(this,SIGNAL(solved(int)),go,SLOT(enable()));
        connect(this,SIGNAL(stopWorkerTask()),go,SLOT(enable()));
        LabelItem * edges = new LabelItem(c, "Edge count");
        connect(this, SIGNAL(solved(int)), edges, SLOT(setValue(int)));
        ButtonItem * verify = new ButtonItem(c, "Verify");
        connect(verify,SIGNAL(clicked(bool)),gl,SLOT(verify()));
        c->done();
        
        ButtonItem * kill = new ButtonItem(c, "Stop current computation");
        connect(kill,SIGNAL(clicked(bool)),this,SIGNAL(stopWorkerTask()));
        
    }
}

void Window::createMenu() {
	// Add menu options
	QMenuBar *bar = menuBar();
	menu.file   = bar->addMenu("File");
	menu.window = bar->addMenu("Window");
	menu.help   = bar->addMenu("Help");
}

void Window::createMenuOptions() {
	QAction *a;

	menu.file->addSeparator();
	a = new QAction("Quit", this);
	connect(a, SIGNAL(triggered()), SLOT(close()));
	menu.file->addAction(a);

	a = new QAction("About Qt4", this);
	connect(a, SIGNAL(triggered()), SLOT(aboutQt()));
	menu.help->addAction(a);
}

void Window::aboutQt() {
	QMessageBox::aboutQt(this);
}


Config *Window::createConfig(QString name) {
	// Create configuration toolbar
	QDockWidget *dock = new QDockWidget(name, this);
	dock->setObjectName(name);
	addDockWidget(Qt::LeftDockWidgetArea, dock);
	menu.window->addAction(dock->toggleViewAction());
	Config *c = new Config(dock);
	dock->setWidget(c);
	return c;
}

void Window::createContext() {
	// Create OpenGL context
	gl = new GLWidget(this);
	setCentralWidget(gl);
}

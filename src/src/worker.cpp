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


#include <QThread>
#include "worker.h"
#include "solver.h"
#include "timing.h"


Task::Task(pointset p, int algo, double t, double var, int landmarks)
    : QObject(), BaseTask(p, t, var, landmarks), edge_count(0), quit(false), algorithm(algo), required_edge_multiplicity(1) {}

uint qHash(const edge &e) {
    return e.x*(e.x+31)+e.y*(e.y+31);
}

/**
 * Creates a new edge in this tasks spanner. Returns true if the algorithm should
 * terminate.
 */
bool Task::create_edge(edge e) {
    edges[e]++;
    if (edges[e]==required_edge_multiplicity) {
#ifndef COMMANDLINE
        emit add_edge(e);
#endif
        edge_count++;
    } else if (edges[e]>required_edge_multiplicity) {
        //qDebug("Edge %d-%d repeated %d times.",e.x,e.y,edges[e]);
	}
    return quit;
}

Worker::Worker(): QObject() {
    qDebug("Worker [%p] Created", this);
}

Worker::~Worker() {
    qDebug("Worker [%p] Deleted", this);
}

void Worker::solve(Task *T) {
    qDebug("Starting to create %lf-spanner for %d points using %s",T->t,T->p.size(),Solver::solvers[T->algorithm].name);
#ifndef COMMANDLINE
    emit T->working(true);
#endif
    Timing::reset();
    Solver::solvers[T->algorithm].function(T);
    double elapsed = Timing::elapsed();
    qDebug("Created with %d edges",T->edge_count);
    qDebug("Time elapsed: %.15lfms",elapsed);
#ifndef COMMANDLINE
    emit T->working(false);
    emit solved(T->edge_count);
#endif
    T->deleteLater();
}

void Task::do_quit() {
    quit=true;
}

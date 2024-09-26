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

/*
#ifndef WORKER_H
#define WORKER_H
#include <QObject>
#include <QHash>
#include "structures.h"

uint qHash(const edge &e);

class Task : public QObject, public BaseTask {
    Q_OBJECT
    int edge_count;
public:
    QHash<edge,int> edges;
    volatile bool quit;
    int algorithm;
    int required_edge_multiplicity;
    Task(pointset p, int algo, double t, double var, int landmarks);
    bool create_edge(edge e);
public slots:
    void do_quit();
signals:
    void working(bool);
    void add_edge(edge e);
    void new_pass();
    friend class Worker;
};

class Worker : public QObject {
    Q_OBJECT
public:
	explicit Worker();
	virtual ~Worker();
public slots:
    void solve(Task *T);
signals:
    void solved(int);
};

#endif // WORKER_H
*/
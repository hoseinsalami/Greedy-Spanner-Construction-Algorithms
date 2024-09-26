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


#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>
#include "structures.h"

//Basically fixes the broken Qt macro, which interferes with C++11 user-defined literals.
#undef QLOCATION
#define QLOCATION "\0" __FILE__ ":" QTOSTRING(__LINE__)

class QTimer;

class GLWidget : public QGLWidget {
		Q_OBJECT
public:
	explicit GLWidget(QWidget *parent = 0);
	virtual ~GLWidget();
protected:
	virtual void initializeGL();
	virtual void resizeGL(int w, int h);
	virtual void paintGL();
	virtual void mousePressEvent(QMouseEvent *e);
	virtual void mouseReleaseEvent(QMouseEvent *e);
	virtual void mouseMoveEvent(QMouseEvent *e);
	virtual void keyPressEvent(QKeyEvent *e);
	virtual void keyReleaseEvent(QKeyEvent *e);
	virtual void wheelEvent(QWheelEvent *e);
private:
	QTimer * frame_timer; // Periodically updates the flow state.
	bool pressed[7]; // keys that are pressed.
	bool setKey(int keycode, bool state);
	QPoint drag;
	QPointF offset;
	float scale;

	pointset points;
	edgelist edges;

    pointset marks;
    pointset annotation_lines;

	int width;
	int height;
	int text_img;
	double t;

	void viewport();
private slots:
	void update_frame();

public slots:
	void scheduleRepaint();
	void reset();
    void setPointset(pointset);
    void resetEdges();
    void addEdge(edge);
    void setT(double v);
    void verify();

signals:
	void changeRunning(bool);
};

#endif // GLWIDGET_H

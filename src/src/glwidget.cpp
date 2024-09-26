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

#include <cmath>
#include <cassert>

#include <QTimer>
#include <QDebug>
#include <QSettings>
#include <QMouseEvent>

#if defined(_WIN32) || defined(_WIN64)
# define GL_MAX_COLOR_ATTACHMENTS          0x8CDF
#else
#include <QtOpenGL>
#endif

#include "glwidget.h"
#include "text.h"
#include "geom.h"
#include "verify.h"
#include <algorithm>


enum {
	KEY_LEFT,
	KEY_RIGHT,
	KEY_FWD,
	KEY_BACK,
	KEY_UP,
	KEY_DOWN,
	KEY_SHIFT
};

GLWidget::GLWidget(QWidget *parent) :
	QGLWidget(parent),
	frame_timer(new QTimer(this)),
	offset(),
	scale(1),
	width(-1),
	height(-1)
{
	makeCurrent(); // Ensure that openGL context is initialized.

	// Set a timer for updating the view periodically
	connect(frame_timer, SIGNAL(timeout()), SLOT(update_frame()));
	frame_timer->setSingleShot(true);
	scheduleRepaint();

	// Make the GL area focusable
	setFocusPolicy(Qt::StrongFocus);

	// Init initial key states
	for (size_t i = 0; i < sizeof(pressed) / sizeof(pressed[0]); i++)
		pressed[i] = false;

	QPixmap font("../normal.png");
	//if(font.isNull()) qDebug("Failed to load normal.png");
	text_img = bindTexture(font, GL_TEXTURE_2D);
	//if(text_img==0) qDebug("Failed to bind normal.png");
}

GLWidget::~GLWidget() {
}

void GLWidget::initializeGL() {
	QGLWidget::initializeGL();
#ifdef SPAMDEBUGMESSAGES
	qDebug("Version: %s", glGetString(GL_VERSION));
	qDebug("Vendor: %s", glGetString(GL_VENDOR));
	qDebug("Renderer: %s", glGetString(GL_RENDERER));
	qDebug("Extensions: %s", glGetString(GL_EXTENSIONS));
	GLint maxbuffers;
	glGetIntegerv(GL_MAX_COLOR_ATTACHMENTS, &maxbuffers);
	qDebug("Max Color Attachments: %d", maxbuffers);
	Q_ASSERT(maxbuffers>=2);
#endif

	glDisable(GL_DEPTH_TEST);

	glClearColor(1, 1, 1, 0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void GLWidget::viewport()
{
	glLoadIdentity();
	glScalef(scale,scale,1);
	glTranslatef(offset.rx(),-offset.ry(),0);
	scheduleRepaint();
}

void GLWidget::resizeGL(int w, int h) {
	QGLWidget::resizeGL(w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, w, h); // update context viewport size
	glOrtho(0, w, 0, h, -1, 1); // set origin to left bottom corner.
	glMatrixMode(GL_MODELVIEW);
	viewport();
	width = w;
	height = h;
}

void GLWidget::reset()
{
	offset=QPointF();
	scale=1;
	viewport();
}

/**
 * Paints the contents of the OpenGL context.
 */
void GLWidget::paintGL() {
	QGLWidget::paintGL();
	glClear(GL_COLOR_BUFFER_BIT);

	// drawing example:
	if (points.size() > 0) {
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(2, GL_DOUBLE, sizeof(vertex), &(points[0].x));

        // Draw points
        glColor3f(1,0,0);
        glPointSize(5);
        glDrawArrays(GL_POINTS,0,points.size());

        if (edges.size() > 0)
        {
            // Draw lines
            glColor4f(0,0,0,0.3);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glEnable(GL_BLEND);
            glDrawElements(GL_LINES, edges.size()*2, GL_UNSIGNED_INT, &edges[0]);
            glDisable(GL_BLEND);
        }
    }

    if (marks.size() > 0) {
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(2, GL_DOUBLE, sizeof(vertex), &(marks[0].x));

        // Draw marks
        glColor3f(0,0,1);
        glPointSize(5);
        glDrawArrays(GL_POINTS,0,marks.size());
    }
    if (annotation_lines.size() > 0) {
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(2, GL_DOUBLE, sizeof(vertex), &(annotation_lines[0].x));

        // Draw annotation lines
        glColor3f(0,1,0);
        glLineWidth(2);
        glDrawArrays(GL_LINES,0,annotation_lines.size());
        glLineWidth(1);
    }

}

void GLWidget::scheduleRepaint()
{
	if (!frame_timer->isActive())
		frame_timer->start(30);
}

struct circsort {
    point c;
    circsort(point c) : c(c) {}
    bool operator()(point a, point b) {
        return (a-c).arg() < (b-c).arg();
    }
};

void GLWidget::mousePressEvent(QMouseEvent *event) {
	QWidget::mousePressEvent(event);
	if (event->button() == Qt::LeftButton) {
		drag = event->pos();
		setCursor(QCursor(Qt::ClosedHandCursor));
	}
	if (event->buttons() == Qt::RightButton) {
        QPointF p = event->pos();
        p.setY(p.y()-height);
        p = p/scale-offset;
        double r=10/scale;
        point v(p.x(),-p.y());
        int selected=-1;
        for (int i=0; i<points.size(); i++) {
            double d=(v-points[i]).length();
            if (d<r) {
                d=r;
                selected=i;
            }
        }
        if (selected>=0) {
            marks.clear();
            annotation_lines.clear();
            QVector<point> neighbours;
            QHash<point,double> dist;
            point core = points[selected];
            marks.push_back(core);
            for (int i=0; i<edges.size(); i++) {
                point o;
                if (edges[i].x==selected) {
                    o = points[edges[i].y];
                } else if (edges[i].y==selected) {
                    o = points[edges[i].x];
                } else continue;
                neighbours.push_back(o);
            }
            std::sort(neighbours.begin(),neighbours.end(),circsort(core));
            QVector<conic> cs;
            for (int i=0; i<neighbours.size(); i++) {
                double dist=(core-neighbours[i]).length();
                conic c = conflict(core,neighbours[i],dist/t);
                cs.push_back(c);
                point r[2];
                int n=c.intersect(line(core,neighbours[i]),r);
                //assert(n==2); NOTE: commented these two lines because they gave compile errors
                point q=((core-r[0]).length() < (core-r[1]).length()) ? r[1] : r[0];
                //for (int k=0; k<2; k++) {
                //    point q=r[k];
                    point s;
                    for (int sign=-1; sign<=1; sign+=2) {
                        s=q;
                        for (int j=0; j<250; j++) {
                            annotation_lines.push_back(s);
                            point dc;
                            dc=point(2*c.a*s.x+c.b*s.y+c.d,c.b*s.x+2*c.c*s.y+c.e);
                            s=s+point(dc.y,-dc.x).normalized()*dist*.03*j*sign;
                            dc=point(2*c.a*s.x+c.b*s.y+c.d,c.b*s.x+2*c.c*s.y+c.e);
                            double v=-c(s.x,s.y)/dc.dot(dc);
                            s=s+dc*v;
                            annotation_lines.push_back(s);
                        }
                    }
                //}
                //marks.push_back(q);
                //annotation_lines.push_back(neighbours[i]);
                //annotation_lines.push_back(neighbours[(i+1)%neighbours.size()]);
            }
            for (int i=0; i<cs.size(); i++) {
                point r[4];
                int ii=(i+1)%cs.size();
                int n=cs[i].intersect(cs[ii], r);
                for (int j=0; j<n; j++) {
                    double dc=(core-r[j]).length();
                    if (dc<=(r[j]-neighbours[i ]).length()) continue;
                    if (dc<=(r[j]-neighbours[ii]).length()) continue;
                    marks.append(r[j]);
                }
            }
        }
        scheduleRepaint();
    }
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event) {
	QWidget::mouseReleaseEvent(event);
	if (event->buttons() == 0) {
		setCursor(QCursor(Qt::ArrowCursor));
	}
}

/** Uses mouse input to do (...).
 * The following actions are assigned to mouse buttons:
 *  - Left: drag viewport
 *  - Middle: -
 *  - Right: select vertex
 * The following modifiers are used:
 *  - Shift: -
 *  - Control: -
 */
void GLWidget::mouseMoveEvent(QMouseEvent *event) {
	QWidget::mouseMoveEvent(event);

	QPointF d = (event->pos() - drag); // Compute moved distance

	if (event->modifiers() & Qt::ShiftModifier) {}
	if (event->modifiers() & Qt::ControlModifier) {}

	if (event->buttons() == Qt::LeftButton) {
		offset+=d/scale;
		viewport();
		drag = event->pos();
	}
	//else if (event->buttons() == Qt::MiddleButton) {} NOTE:commented this line because it doesnt do anything but gives compile errors
	else if (event->buttons() == Qt::RightButton) {
	}
}

/** Uses mouse wheel scrolling for zooming in and out.
 */
void GLWidget::wheelEvent(QWheelEvent *event) {
	QWidget::wheelEvent(event);
	QPointF p = event->pos();
	p.ry()-=height;
	offset-=p/scale;
	scale *= pow(2,event->delta()/720.);
	offset+=p/scale;
	viewport();
}

/** Handles the key event to keep track of button press state.
 * Returns true if pressing the key requires periodic updates of the
 * OpenGL context.
 */
bool GLWidget::setKey(int keycode, bool state) {
	switch(keycode) {
		case Qt::Key_R:
			reset();
			return true;
		case Qt::Key_Left:
			pressed[KEY_LEFT] = state;
			return true;
		case Qt::Key_Right:
			pressed[KEY_RIGHT] = state;
			return true;
		case Qt::Key_A:
			pressed[KEY_FWD] = state;
			return true;
		case Qt::Key_Z:
			pressed[KEY_BACK] = state;
			return true;
		case Qt::Key_Down:
			pressed[KEY_DOWN] = state;
			return true;
		case Qt::Key_Up:
			pressed[KEY_UP] = state;
			return true;
		case Qt::Key_Shift:
			pressed[KEY_SHIFT] = state;
			return false;
		default:
			return false;
	}
}

void GLWidget::keyPressEvent(QKeyEvent *e) {
	if (!e->isAutoRepeat() && setKey(e->key(), true)) {
		scheduleRepaint();
	} else {
		QWidget::keyPressEvent(e);
	}
}

void GLWidget::keyReleaseEvent(QKeyEvent *e) {
	if (!e->isAutoRepeat() && setKey(e->key(), false)) {
	} else {
		QWidget::keyReleaseEvent(e);
	}
}

void GLWidget::update_frame() {
	bool v=false;
	float s=25.f/scale;
	offset.rx()-=width/2.f/scale;
	offset.ry()+=height/2.f/scale;
	if (pressed[KEY_FWD]) {
		scale*=1.1;
		v=true;
	}
	if (pressed[KEY_BACK]) {
		scale/=1.1;
		v=true;
	}
	offset.rx()+=width/2.f/scale;
	offset.ry()-=height/2.f/scale;
	if (pressed[KEY_LEFT]) {
		offset.rx()+=s;
		v=true;
	}
	if (pressed[KEY_RIGHT]) {
		offset.rx()-=s;
		v=true;
	}
	if (pressed[KEY_DOWN]) {
		offset.ry()-=s;
		v=true;
	}
	if (pressed[KEY_UP]) {
		offset.ry()+=s;
		v=true;
	}
	if (v) viewport();
	updateGL(); // only repaint call.
}

void GLWidget::resetEdges() {
    edges.clear();
    marks.clear();
    annotation_lines.clear();
    scheduleRepaint();
}

void GLWidget::addEdge(edge e) {
    edges.append(e);
    scheduleRepaint();
}

void GLWidget::setPointset(pointset p) {
    points = p;
    edges.clear();
    marks.clear();
    annotation_lines.clear();
    scheduleRepaint();
}

void GLWidget::setT(double v) {
    t=v;
}

void GLWidget::verify()
{
    verification_info vi = ::verify(points, edges);
	qDebug("t=%lf max average t=%lf weight=%lf max degree=%d max length=%lf", vi.t, vi.maxavgt, vi.weight, vi.maxdeg, vi.maxlength);
}

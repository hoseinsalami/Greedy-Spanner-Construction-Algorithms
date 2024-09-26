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


#ifndef GEOM_H
#define GEOM_H
#include "structures.h"

struct point {
    double x,y;
    point() {}
    point(vertex v) : x(v.x), y(v.y) {}
    point(double nx, double ny) : x(nx), y(ny) {}
    point conj() const;
    point operator+(const point &p) const;
    point operator-(const point &p) const;
    point operator*(const double c) const;
    point operator*(const point &p) const;
    point operator/(const double c) const;
    point operator/(const point &p) const;
    bool operator==(const point &p) const;
    operator vertex() const {return vertex(x,y);}
    double length() const;
    double arg() const;
    double dot(const point &p) const;
    double cross(const point &p) const;
	point normalized();
};
point polar(const double l, const double h);

struct line {
    double a, b, c;
    line() {}
    line(double na, double nb, double nc) : a(na), b(nb), c(nc) {}
    explicit line(const point &p) : a(-p.y), b(p.x), c(0) {}
    line(const point &p, const point &q) : a(q.y - p.y), b(p.x - q.x), c(a*p.x + b*p.y) {}
    line shift(const point &p) const;
    point perp() const;
    double dist(const point &p) const;
    point project(const point &p) const;
    point intersect(const line &m) const;
    bool parallel(const line &m) const;
    bool operator==(const line &k) const;
};
void bisect(line a, line b, line r[2]);

struct conic { 
    double a,b,c,d,e,f; 
    conic() {}
    conic(double a, double b, double c, double d, double e, double f) : a(a), b(b), c(c), d(d), e(e), f(f) {}
    //conic(const circle &c) : a(1), b(0), c(1), d(-2 * c.c.x), e(-2 * c.c.y), f(c.c.dot(c.c) - c.r * c.r) {}
    double operator()(double x, double y) const;
    conic operator*(double s) const;
    conic operator+(conic z) const;
    double det() const;
    double foo(const conic &z) const;
    int intersect(const line &l, point r[2]);
    int dissect(line r[2]);
    int intersect(const conic &z, point r[4]);
};
conic conflict(point p, point q, double f);
conic conflict(line l, point p);
conic conflict(line l, point p, double r);
conic conflict(line l, line m);

struct segment {
    point a;
    point b;
    segment(const point &na, const point &nb): a(na), b(nb) {}
    double length() const;
    point project(const point &p) const;
    double dist(const point &p) const;
};

#endif // GEOM_H

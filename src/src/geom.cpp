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

#include <cmath>
#include "geom.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#if defined(_WIN32) || defined(_WIN64)
#define copysign _copysign
#endif

using namespace std;

point point::conj() const {
    return point(x, -y);
}

point point::operator+(const point &p) const {
    return point(x + p.x, y + p.y);
}

point point::operator-(const point &p) const {
    return point(x - p.x, y - p.y);
}

point point::operator*(const double c) const {
    return point(x * c, y * c);
}

point point::operator*(const point &p) const {
    return point((x * p.x) - (y * p.y), (x * p.y) + (y * p.x));
}

point point::operator/(const double c) const {
    return point(x / c, y / c);
}

point point::operator/(const point &p) const {
    return (*this * p.conj()) / (p.x * p.x + p.y * p.y);
}

bool point::operator==(const point &p) const {
    return (x == p.x) && (y == p.y);
}

double point::length() const {
    return hypot(x, y);
}

double point::arg() const {
    return atan2(y, x);
}

double point::dot(const point &p) const {
    return (x * p.x) + (y * p.y);
}

double point::cross(const point &p) const {
    return (x * p.y) - (p.x * y);
}

point point::normalized() {
    return operator/(length());
}

point polar(const double l, const double h) {
    return point((l * cos(h)), (l * sin(h)));
}


line line::shift(const point &p) const {
    return line(a, b, c + (a * p.x) + (b * p.y));
}

point line::perp() const {
    return point(a,b);
}

double line::dist(const point &p) const {
    point q = perp();
    return abs(p.dot(q) - c) / q.length();
}

point line::project(const point &p) const {
    point q = perp();
    return p - q * ((p.dot(q) - c) / q.dot(q));
}

point line::intersect(const line &m) const {
    double d = (m.a * b) - (a * m.b);
    if (abs(d) < 1e-9) {} // lines are parallel
    return point( (b*m.c - m.b*c)/d , (m.a*c - a*m.c)/d );
}

bool line::parallel(const line &m) const {
    return abs((m.a * b) - (a * m.b)) < 1e-6;
}

bool line::operator==(const line &k) const {
    return (abs(a*k.c - c*k.a) + abs(b*k.c - c*k.b) + abs(a*k.b - b*k.a)) < 1e-6;
}

void bisect(line a, line b, line r[2]) {
    point c = a.intersect(b);
    point d(a.a, a.b);
    point e(b.a, b.b);
    d = d/d.length();
    e = e/e.length();
    r[0] = line(d + e).shift(c);
    r[1] = line(d - e).shift(c);
}

static int square(double a,double b,double c, double ret[]) {
    double d = b*b - 4*a*c;
    if (abs(d)<1e-6) {
        ret[0]=-b/2/a;
        return 1;
    } else if (d<0) {
        return 0;
    } else {
        d=sqrt(d); 
        ret[0]=(-b-d)/2/a;
        ret[1]=(-b+d)/2/a;
        return 2;
    }
}

static double cuberoot(double a) {
    if (a<0) return -pow(-a,1.0/3.0);
    return pow(a,1.0/3.0);
}

static int cubic(double a,double b,double c, double ret[3]) {
    double p,q,d,w1,w2,n,h;
    p = b - a*a/3;
    q = c + (2*a*a-9*b)*a/27;
    d = q*q/4 + p*p*p/27;
    if (d<0) {
        w1 = -q/2;
        w2 = sqrt(-d);
        n = 2*pow(w1*w1+w2*w2, 1.0/6.0);
        h = atan2(w2,w1)/3;
        ret[0] = n * cos(h) - a/3;
        ret[1] = n * cos(h + M_PI*2/3) - a/3;
        ret[2] = n * cos(h - M_PI*2/3) - a/3;
        return 3;
    } else {
        w1 = cuberoot(-q/2 + sqrt(d));
        w2 = cuberoot(-q/2 - sqrt(d));
        ret[0]=w1 + w2 - a/3;
        return 1;
    }
}

#define bar(a,b,c) (z.a*b*c + a*z.b*c + a*b*z.c)
double conic::operator()(double x, double y) const {
	return a * x * x + b * x * y + c * y * y + d * x + e * y + f;
}

conic conic::operator*(double s) const {
	return conic(s * a, s * b, s * c, s * d, s * e, s * f);
}

conic conic::operator+(conic z) const {
	return conic(a + z.a, b + z.b, c + z.c, d + z.d, e + z.e, f + z.f);
}

double conic::det() const {
	return 4 * a * c * f + b * d * e - a * e * e - c * d * d - f * b * b;
}

double conic::foo(const conic &z) const {
	return 4 * bar(a, c, f) + bar(b, d, e) - bar(a, e, e) - bar(c, d, d) - bar(f, b, b);
}

int conic::intersect(const line &l, point r[2]) {
	double g = l.a, h = l.b, j = l.c / (g * g + h * h), t[2];
	double A = a * h * h - b * g * h + c * g * g;
	double B = j * (2 * g * h * (c - a) + b * (g * g - h * h)) - d * h + e * g;
	double C = (a * g * g + b * g * h + c * h * h) * j * j + (d * g + e * h) * j + f;
	if (abs(A) < 1e-6) {
		if (abs(B) < 1e-6) return 0;
		r[0] = point(g * j + h * C / B, h * j - g * C / B);
		return 1;
	} else {
		int n = square(A, B, C, t);
		for(int i = 0; i < n; i++) r[i] = point(g * j - t[i] * h, h * j + t[i] * g);
		return n;
	}
}

int conic::dissect(line r[2]) {
	double u = b * b - 4 * a * c;
	if (abs(u) < 1e-6) {
		double m[2], g, h, v, w;
		if (abs(a) < 1e-9 && abs(c) < 1e-9) {
			r[0] = line(d, e, -f);
			return 1;
		}
		if (abs(a) > abs(c)) {
			g = -b / 2 / a;
			h = -1;
			v = d / -a;
			w = f / a;
		} else {
			g = 1;
			h = b / 2 / c;
			v = e / c;
			w = f / c;
		}
		int n = square(1, v, w, m);
		for(int i = 0; i < n; i++) r[i] = line(h, g, m[i]);
		return n;
	} else if (u < 0) {
		return 0;
	} else {
		point h = point(2 * c * d - b * e, 2 * a * e - d * b) / u;
		if (abs(a) > abs(c)) {
			r[0] = line(2 * a, b + sqrt(u), 0).shift(h);
			r[1] = line(2 * a, b - sqrt(u), 0).shift(h);
		} else {
			r[0] = line(b + sqrt(u), 2 * c, 0).shift(h);
			r[1] = line(b - sqrt(u), 2 * c, 0).shift(h);
		}
		return 2;
	}
}

int conic::intersect(const conic &z, point r[4]) {
	double m[3], m3 = det(), m2 = foo(z), m1 = z.foo(*this), m0 = z.det();
	conic p, q;
	if (abs(m3) > abs(m0)) {
		if (abs(m3) < 1e-9) m[0] = 0;
		else cubic(m2 / m3, m1 / m3, m0 / m3, m);
		p = *this * m[0] + z;
		q = *this + z * -m[0];
	} else {
		cubic(m1 / m0, m2 / m0, m3 / m0, m);
		p = *this + z * m[0];
		q = *this * -m[0] + z;
	}
	line u[2];
	int n = p.dissect(u), t = 0;
	for(int i = 0; i < n; i++) t += q.intersect(u[i], r + t);
	for(int i = 0; i < t; i++) {
        for(int j = 0; j < 10; j++) {
			double x = r[i].x, y = r[i].y, v[2] = { operator()(x, y), z(x, y) };
			if (!(abs(v[0]) > 1e-9 || abs(v[1]) > 1e-9)) break;
			double h[4] = {z.b *x + 2 * z.c *y + z.e, -b *x - 2 * c *y - e, -2 * z.a *x - z.b *y - z.d, 2 * a *x + b *y + d };
			double k = (h[0] * h[3] - h[1] * h[2]);
			if (abs(k) < 1e-3) break;
			r[i].x = x - (h[0] * v[0] + h[1] * v[1]) / k;
			r[i].y = y - (h[2] * v[0] + h[3] * v[1]) / k;
		}
    }
	return t;
}

conic conflict(point p, point q, double f) {
	double a = p.x - q.x, b = p.y - q.y, c = p.dot(p) - q.dot(q), d = p.dot(p) + q.dot(q);
	f *= f;
	return conic(f - a * a, -2 * a * b, f - b * b, a * c - (p.x + q.x) * f, b * c - (p.y + q.y) * f, f * d / 2 - (c * c + f * f) / 4);
}

conic conflict(line l, point p) {
	double s = 1.0 / (l.a * l.a + l.b * l.b);
	return conic(s * l.a * l.a - 1, s * l.a * l.b * 2, s * l.b * l.b - 1,
	             2 * (p.x - l.a * l.c * s), 2 * (p.y - l.b * l.c * s), s * l.c * l.c - p.dot(p));
}

conic conflict(line l, point p, double r) {
	point t = l.perp();
	line u = l.shift(t * copysign(r / t.length(), -t.dot(p)));
	return conflict(u, p);
}

conic conflict(line l, line m) {
	double s = 1.0 / (l.a * l.a + l.b * l.b), t = 1.0 / (m.a * m.a + m.b * m.b);
	return conic(s * l.a * l.a - t * m.a * m.a, s * l.a * l.b * 2 - t * m.a * m.b * 2, s * l.b * l.b - t * m.b * m.b,
	             2 * (t * m.a * m.c - s * l.a * l.c), 2 * (t * m.b * m.c - s * l.b * l.c), s * l.c * l.c - t * m.c * m.c);
}


double segment::length() const {
    return (a-b).length();
}

point segment::project(const point &p) const {
    double l = (p-a).dot(b-a)/(b-a).dot(b-a);
    if (l<0) l=0;
    if (l>1) l=1;
    return a+(b-a)*l;
}

double segment::dist(const point &p) const {
    return (p-project(p)).length();
}


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
#include <iostream>
#include <fstream>
#include <sstream>

#include "random.h"
#include "generator.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

pointset random(int seed, int& count, int max_coordinate, double shape, double& t, int source) {
    random_init(seed);
    pointset r;
    for (int i = 0; i < count; i++) {
        r.push_back(vertex(
            random_range_d(0, max_coordinate),
            random_range_d(0, max_coordinate)
        ));
    }

    return r;
}

pointset normal(int seed, int& count, int max_coordinate, double shape, double& t, int source) {
    random_init(seed);
    pointset r;
    max_coordinate/=2;
    for (int i = 0; i < count; i++) {
        double U = random_range_d(0,1);
        double V = random_range_d(0,2*M_PI);
        U = sqrt(-2*log(U)) * max_coordinate;
        r.push_back(vertex(
            U*sin(V) + max_coordinate,
            U*cos(V) + max_coordinate
        ));
    }

    return r;
}

pointset normal_spots(int seed, int& count, int max_coordinate, double shape, double& t, int source) {
    random_init(seed);
    int sn = sqrt((double)count);
    pointset r;
    int x=0,y=0;
    int smc = sqrt((double)max_coordinate)/2;
    max_coordinate/=2;
    for (int i = 0; i < count; i++) {
        if (i%sn == 0) {
            double U = random_range_d(0,1);
            double V = random_range_d(0,2*M_PI);
            U = sqrt(-2*log(U)) * max_coordinate;
            x=U*sin(V) + max_coordinate;
            y=U*cos(V) + max_coordinate;
        }
        double U = random_range_d(0,1);
        double V = random_range_d(0,2*M_PI);
        U = sqrt(-2*log(U)) * smc;
        r.push_back(vertex(
            U*sin(V) + x,
            U*cos(V) + y
        ));
    }

    return r;
}

pointset squares(int seed, int& count, int max_coordinate, double shape, double& t, int source) {
    random_init(seed);
    int sn = sqrt((double)count);
    pointset r;
    int x=0,y=0;
    double smc = sqrt((double)max_coordinate);
    for (int i = 0; i < count; i++) {
        if (i%sn == 0) {
            x=random_range_d(0, max_coordinate-smc);
            y=random_range_d(0, max_coordinate-smc);
        }
        r.push_back(vertex(
            x+random_range_d(0, smc),
            y+random_range_d(0, smc)
        ));
    }
    return r;
}

pointset circles(int seed, int& count, int max_coordinate, double shape, double& t, int source) {
    random_init(seed);
    int sn = count/(sqrt(sqrt((double)count))*shape);
    pointset res;
    double x=0,y=0,r=0;
    for (int i = 0; i < count; i++) {
        if (i%sn == 0) {
            r = random_range_d(1,random_range_d(1,max_coordinate/2));
            x=random_range_d(r, max_coordinate-r);
            y=random_range_d(r, max_coordinate-r);
        }
        double phi = random_range_d(0,2*M_PI);
        res.push_back(vertex(
            x+r*cos(phi),
            y+r*sin(phi)
        ));
    }
    return res;
}

pointset M(int seed, int& count, int max_coordinate, double shape, double& t, int source) {
    random_init(seed);
    pointset r;
    for (int i = 0; i < count / 4; i++) {
        r.push_back(vertex(random_range_d(0, max_coordinate),max_coordinate));
        r.push_back(vertex(random_range_d(0, max_coordinate),0));
    }
    for (int i = 0; i < count / 8; i++) {
        r.push_back(vertex(max_coordinate,random_range_d(max_coordinate * 3 / 4, max_coordinate)));
        r.push_back(vertex(max_coordinate,random_range_d(0, max_coordinate / 4)));
        r.push_back(vertex(max_coordinate / 2,random_range_d(max_coordinate / 4, max_coordinate * 3 / 4)));
        r.push_back(vertex(random_range_d(max_coordinate / 2, max_coordinate),max_coordinate / 4));
        r.push_back(vertex(random_range_d(max_coordinate / 2, max_coordinate),max_coordinate * 3 / 4));
    }
    return r;
}

pointset parabola(int seed, int& count, int max_coordinate, double shape, double& t, int source) {
    random_init(seed);
    pointset r;
	const double exponent = shape;
	const double s = pow((double)max_coordinate, (exponent-1) / exponent);
    for (int i = 0; i < count / 2; i++) {
		const double x = random_range_d(0, max_coordinate);
        r.push_back(vertex(x,(unsigned int)(pow((double)x, 1.0 / exponent) * s / 2.0) + max_coordinate / 2));
        r.push_back(vertex(x,max_coordinate / 2 - (unsigned int)(pow((double)x, 1.0 / exponent) * s / 2.0)));
    }
    return r;
}

pointset circles2(int seed, int& count, int max_coordinate, double shape, double& t, int source) {
	using namespace std;
    random_init(seed);
    pointset res;
	ofstream myfile;
	myfile.open("challenge.txt");
	myfile << count << endl;
	myfile << 11 << " " << 10 << endl;
	for (int i = 0; i < count / 8; i++) {
		double x = ((double)max_coordinate) + ((double)max_coordinate) / 2.0 * (1+cos(2.0*M_PI*((double)i)/(((double)count)/8.0))*1.3);
		double y = ((double)max_coordinate) + ((double)max_coordinate) / 2.0 * (1+sin(2.0*M_PI*((double)i)/(((double)count)/8.0)));
        res.push_back(vertex(
            x,
            y
        ));
		myfile << x << " " << y << endl;
    }
	for (int i = 0; i < count * 7 / 8; i++) {
		double x = ((double)max_coordinate) + ((double)max_coordinate) / 2.0 * (1.0 / shape + cos(2.0*M_PI*((double)i)/(((double)count)*7.0/8.0)) / shape / 2.0);
		double y = ((double)max_coordinate) + ((double)max_coordinate) / 2.0 * (1.0 + sin(2.0*M_PI*((double)i)/(((double)count)*7.0/8.0)) / shape / 2.0);
        res.push_back(vertex(
            x,
            y
        ));
		myfile << x << " " << y << endl;
    }
	myfile.close();
	/*
	for (int i = 0; i < count * 7 / 32; i++) {
        double phi = random_range_d(0,2*M_PI);
        res.push_back(vertex(
            (int)(max_coordinate / 2.0 * (2 - 1 / shape - cos(phi) / shape / 2.0)),
            (int)(max_coordinate / 2.0 * (1 + sin(phi) / shape / 2.0))
        ));
    }
	for (int i = 0; i < count * 7 / 32; i++) {
        double phi = random_range_d(0,2*M_PI);
        res.push_back(vertex(
            (int)(max_coordinate / 2.0 * (1 + sin(phi) / shape / 2.0)),
            (int)(max_coordinate / 2.0 * (1 / shape + cos(phi) / shape / 2.0))
        ));
    }
	for (int i = 0; i < count * 7 / 32; i++) {
        double phi = random_range_d(0,2*M_PI);
        res.push_back(vertex(
            (int)(max_coordinate / 2.0 * (1 + sin(phi) / shape / 2.0)),
            (int)(max_coordinate / 2.0 * (2 - 1 / shape - cos(phi) / shape / 2.0))
        ));
    }
	*/
    return res;
}

pointset martin(int seed, int& count, int max_coordinate, double shape, double& t, int source) {
    random_init(seed);
    pointset res;
	unsigned int squares = 8;
	for (unsigned int k = 0; k < squares; k++) {
		double dk = k;
		double dmax_coordinate = max_coordinate;
		double dcount = count;
		double dsquares = squares;
		for (unsigned int i = 0; i < count / squares / 4; i++) {
			double di = i;
			// dmax_coordinate * (1 - dk / 8.0)
			// dmax_coordinate * (dk / 8.0 / 2.0)
			res.push_back(vertex(
				dmax_coordinate * (dk / 8.0 / 2.0) + di / (dcount / dsquares / 4) * dmax_coordinate * (1 - dk / 8.0),
				dmax_coordinate * (dk / 8.0 / 2.0)
			));
		}
		for (unsigned int i = 0; i < count / squares / 4; i++) {
			double di = i;
			res.push_back(vertex(
				dmax_coordinate * (dk / 8.0 / 2.0),
				dmax_coordinate * (dk / 8.0 / 2.0) + di / (dcount / dsquares / 4) * dmax_coordinate * (1 - dk / 8.0)
			));
		}
		for (unsigned int i = 0; i < count / squares / 4; i++) {
			double di = i;
			res.push_back(vertex(
				dmax_coordinate * (1 - dk / 8.0 / 2.0) - di / (dcount / dsquares / 4) * dmax_coordinate * (1 - dk / 8.0),
				dmax_coordinate * (1 - dk / 8.0 / 2.0)
			));
		}
		for (unsigned int i = 0; i < count / squares / 4; i++) {
			double di = i;
			res.push_back(vertex(
				dmax_coordinate * (1 - dk / 8.0 / 2.0),
				dmax_coordinate * (1 - dk / 8.0 / 2.0) - di / (dcount / dsquares / 4) * dmax_coordinate * (1 - dk / 8.0)
			));
		}
	}
    return res;
}

double generateErlang(int k, double mu) {
	double result = 0;
	for (int i = 0; i < k; i++) {
		result -= log(random_range_d(0, 1)) * mu;
	}
	return result;
}

const double E = 2.7182818284590452353602874;

double generateGammaPartial(double delta) {
    double v0 = E / (E - delta);
    while (true)
    {
        double v1 = random_range_d(0, 1);
        double v2 = random_range_d(0, 1);
        double v3 = random_range_d(0, 1);
        double xi, eta;
        if (v1 < v0)
        {
            xi = pow(v2, 1 / delta);
            eta = v3 * pow(xi, delta - 1);
        }
        else
        {
            xi = 1 - log(v2);
            eta = v3 * pow(E, -xi);
        }
        if (eta <= pow(xi, delta - 1) * pow(E, -xi))
        {
            return xi;
        }
    }
}

double generateGamma(double alpha, double beta) {
	int alphaFloor = (int) floor(alpha);
	double alphaRest = alpha - alphaFloor;
	return beta * generateGammaPartial(alphaRest) - generateErlang(alphaFloor, 1);
}

pointset gamma(int seed, int &count, int max_coordinate, double shape, double& t, int source) {
    random_init(seed);
    pointset r;

    for (int i = 0; i < count; i++) {
        r.push_back(vertex(
            generateGamma(shape, max_coordinate),
            generateGamma(shape, max_coordinate)
        ));
    }

    return r;
}

pointset file(int seed, int &count, int max_coordinate, double shape, double& t, int source) {
	using namespace std;
    random_init(seed);
    const int bits = random_count_bits(0, max_coordinate);
    pointset r;
        char filename[80];
	sprintf(filename,"input.txt");
	ifstream myfile(filename);
	myfile >> count;
	double x, y;
    for (int i = 0; i < count; i++) {
		myfile >> x;
		myfile >> y;
        r.push_back(vertex(
            x,
            y
        ));
    }
	myfile.close();

    return r;
}

pointset circles3(int seed, int& count, int max_coordinate, double shape, double& t, int source) {
    random_init(seed);
    pointset res;
	res.push_back(vertex(max_coordinate/2, max_coordinate/2));
	for (unsigned int i = 0; i < count; i++) {
		int x = max_coordinate/2 * (1.0 + cos(2.0 * M_PI * ((double)i) / ((double)count)));
		int y = max_coordinate/2 * (1.0 + sin(2.0 * M_PI * ((double)i) / ((double)count)));
        res.push_back(vertex(
            x,
            y
        ));
    }
	return res;
}

pointset counter(int seed, int& count, int max_coordinate, double shape, double& t, int source) {
    pointset res;

	for (int i = 0; i < count / 2; i++)
	{
		res.push_back(vertex(i, 2+0.001*i));
		res.push_back(vertex(i, 1-0.001*i));
	}

	return res;
}

Generator::type Generator::generators[] = {
	{"Uniform random", random}, // 0
    {"Squares", squares}, // 1
    {"Normal", normal}, // 2
    {"Normal spots", normal_spots}, // 3
    {"Circles", circles}, // 4
    {"M", M}, // 5
    {"Parabola", parabola}, // 6
    {"Circles2", circles2}, // 7
    {"Martin", martin}, // 8
    {"Gamma", gamma}, // 9
    {"From file", file}, // 10
    {"Circles3", circles3}, // 11
    {"Counter", counter}, // 12
		{nullptr,nullptr}, // Last element of list.
};

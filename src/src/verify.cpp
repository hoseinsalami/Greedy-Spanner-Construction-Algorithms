/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2012  <copyright holder> <email>

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

#include <vector>

#include "verify.h"
#include "geom.h"
#include "binaryheap.h"

#define D(x,y) distances[(x)+(y)*N]

inline bool is_infinite(const double x) {
    double y = x - x;
    return x == x && y != y;
}

void updateDistances(pointset& vertices, std::vector< std::pair<unsigned int, double> >* myEdges, unsigned int from, std::vector<double>& distances) {
    unsigned int N = vertices.size();
    BinHeap<double> heap(N, 0);
    for (unsigned int i = 0; i < N; i++) {
        heap.insert(i, std::numeric_limits<double>::infinity());
    }
    heap.decreaseKey(from, 0);
    while (heap.getCount() > 0) {
        std::pair<unsigned int, double> pair = heap.getMin();
        if (is_infinite(pair.second))
            break;
        D(from, pair.first) = pair.second;
        D(pair.first, from) = pair.second;
        heap.extractMin();
        for (unsigned int i = 0; i < myEdges[pair.first].size(); i++) {
            std::pair<unsigned int, double> edge = myEdges[pair.first][i];
            double alt = pair.second + edge.second;
            if (alt < heap.getValue(edge.first)) {
                heap.decreaseKey(edge.first, alt);
            }
        }
    }
}

verification_info verify(pointset p, edgelist e)
{
    verification_info r;
    r.maxdeg=0;
    r.t=0;
    r.maxavgt=0;
    r.weight=0;
    r.maxlength=0;
    
    int N = p.size();
    std::vector<double> degrees(N);
	for (size_t i = 0; i < N; i++)
	{
		degrees[i] = 0;
	}

    //degrees..fill(0,N);
    for (unsigned int i = 0; i < e.size(); i++) {
		edge m = e[i];
        double d = (point(p[m.x])-p[m.y]).length();
        r.weight+=d;
        degrees[m.x]++;
        degrees[m.y]++;
		if (d > r.maxlength)
			r.maxlength = d;
    }       
    if (p.size()<=dilationVerificationBound) {
        std::vector<double> distances(N*N);
		for (size_t i = 0; i < N; i++)
		{
			distances[i] = std::numeric_limits<double>::infinity();
		}
        //distances.fill(std::numeric_limits<double>::infinity(),N*N);
		std::vector< std::pair<unsigned int, double> >* myEdges = new std::vector< std::pair<unsigned int, double> > [N];
        for (int i=0; i<N; i++) 
            distances[N*i+i] = 0;
        for (unsigned int i = 0; i < e.size(); i++) {
			edge m = e[i];
            double d = (point(p[m.x])-p[m.y]).length();
            D(m.x,m.y) = d;
            D(m.y,m.x) = d;
			myEdges[m.x].push_back(std::pair<unsigned int, double>(m.y, d));
			myEdges[m.y].push_back(std::pair<unsigned int, double>(m.x, d));
        }       
        
        for (int i=0; i<N; i++) {
			updateDistances(p, myEdges, i, distances);
        }
		delete[] myEdges;
        for (int j=0; j<N; j++) {
			double avgt = 0;
			int countervoordezekerheid = 0;
            for (int k=0; k<N; k++) {
                double d = (point(p[j])-p[k]).length();
                if (d<1e-9) continue;
                if (D(j,k)>1e49) r.t=10000;
                else {
                    double at = D(j,k)/d;
					avgt += at;
					countervoordezekerheid++;
                    if (r.t < at)
                        r.t = at;
                }
            }
			avgt /= countervoordezekerheid;
			if (r.maxavgt < avgt)
				r.maxavgt = avgt;
        }
    }
    /*
    for (int j=0; j<N; j++) {
        for (int i=0; i<N; i++) {
            if (D(i,j)<10000000) 
                printf("%8.2lf ",D(i,j));
            else
                printf("  ****** ");
        }
        printf("\n");
    }
    */
    
    for (int j=0; j<N; j++) {
        if (r.maxdeg < degrees[j])
            r.maxdeg = degrees[j];
    }
    return r;
}

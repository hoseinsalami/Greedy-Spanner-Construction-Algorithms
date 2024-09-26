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

#define NOMINMAX

#include <limits>
#include <algorithm>
#include <cmath>
#include <stack>
#include <iostream>
#include <unordered_set>
#include <fstream>
#include <set>

#include "solver.h"
#include "binaryheap.h"
#include "rangetree.h"
#include "mystack.h"
#include "random.h"
#include "timing.h"

#include "treap.h"
#include "scapegoat.h"

//#include "del_interface.h"

//#include "memory.h"
#include "logstuff.h"

#include "windows.h"
#include "psapi.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define	FirstCandidate	1
#define	BestCandidate	2
#define	NoneVertex	0
#define	FirstVertex	1
#define	SecondVertex	2
#define	BothVertex	3

//candidatePolicy:
//1:FirstCandidate
//2:BestCandidate

using std::cout;
using std::endl;

extern int	qcount = 0;
int rr = 0;

inline edge order(edge v) {
	if (v.x <= v.y) return v;
	else return edge(v.y, v.x);
}

struct EdgeInfo : public edge
{
	EdgeInfo(unsigned int x, unsigned int y, double distance) : edge(x, y), distance(distance) { }
    double distance;
    inline bool operator < (const EdgeInfo &o) const {
		if (distance != o.distance)//(abs(distance - o.distance) >= 1e-6)
			return distance < o.distance;
		if (x != o.x)
			return x < o.x;
		return y < o.y;
	}
	inline EdgeInfo Invert() const {
		EdgeInfo result(y, x, distance);
		return result;
	}
};

inline double distance(const vertex &a, const double x, const double y)
{
    return hypot(a.x-x,a.y-y);
}

inline double distance(const vertex &a, const vertex &b)
{
    return distance(a, b.x, b.y);
}

inline double NotSquaredDistance(const vertex& a, const vertex& b)
{
	double dx = (a.x - b.x) * (a.x - b.x);
	double dy = (a.y - b.y) * (a.y - b.y);
	return dx + dy;
}

inline bool is_infinite(const double x) {
    double y = x - x;
    return x == x && y != y;
}

template<class Heap>
struct Dijkstra {
static double computeDistance(pointset vertices, std::vector< std::pair<unsigned int, double> >* myEdges, unsigned int from, unsigned int to)
{
    unsigned int N = vertices.size();
    Heap heap(N, 0);
    for (unsigned int i = 0; i < N; i++) {
        heap.insert(i, std::numeric_limits<double>::infinity());
    }
    heap.decreaseKey(from, 0);
    double result = std::numeric_limits<double>::infinity();
    while (heap.getCount() > 0) {
        std::pair<unsigned int, double> pair = heap.getMin();
        if (pair.first == to)
        {
            result = pair.second;
            break;
        }
        heap.extractMin();
        for (unsigned int i = 0; i < myEdges[pair.first].size(); i++) {
            std::pair<unsigned int, double> edge = myEdges[pair.first][i];
            double alt = pair.second + edge.second;
            if (alt < heap.getValue(edge.first)) {
                heap.decreaseKey(edge.first, alt);
            }
        }
    }
    return result;
}
};

template<class Heap>
struct AStar {
static double computeDistance(pointset vertices, std::vector< std::pair<unsigned int, double> >* myEdges, unsigned int from, unsigned int to)
{
    unsigned int N = vertices.size();
    double* realDistance = new double[N];
    Heap heap(N, 0);
    for (unsigned int i = 0; i < N; i++) {
        realDistance[i] = std::numeric_limits<double>::infinity();
        heap.insert(i, std::numeric_limits<double>::infinity());
    }
    realDistance[from] = 0;
    heap.decreaseKey(from, 0);
    double result = std::numeric_limits<double>::infinity();
    while (heap.getCount() > 0) {
        std::pair<unsigned int, double> pair = heap.getMin();
        if (pair.first == to)
        {
            result = realDistance[pair.first];
            break;
        }
        heap.extractMin();
        for (unsigned int i = 0; i < myEdges[pair.first].size(); i++) {
            std::pair<unsigned int, double> edge = myEdges[pair.first][i];
            double alt = realDistance[pair.first] + edge.second;
            if (alt < realDistance[edge.first]) {
                realDistance[edge.first] = alt;
                heap.decreaseKey(edge.first, alt + distance(vertices[edge.first], vertices[to]));
            }
        }
    }
    delete[] realDistance;
    return result;
}
};

template<class DistAlgo>
int GreedySpanner(BaseTask * T)
{
	STARTLOG
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
    std::vector<EdgeInfo> edgeList;
    std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];

    edgeList.reserve((N-1) * N / 2);
    for (unsigned int y = 0; y < N; y++) {
        for (unsigned int x = y+1; x < N; x++) {
            EdgeInfo e(x, y, distance(vertices[x],vertices[y]));
            edgeList.push_back(e);
        }
    }
    std::sort(edgeList.begin(), edgeList.end());

    int edgeCount = 0;
    for (unsigned int k = 0; k < (N-1) * N / 2; k++) {
        double computedDist = DistAlgo::computeDistance(vertices, myEdges, edgeList[k].x, edgeList[k].y);
        if (edgeList[k].distance * T->t < computedDist) {
            std::pair<unsigned int, double> pair(edgeList[k].x, edgeList[k].distance);
            myEdges[edgeList[k].y].push_back(pair);
            pair.first = edgeList[k].y;
            myEdges[edgeList[k].x].push_back(pair);
            if (T->create_edge(edgeList[k])) return 0;
            edgeCount++;
        }
    }
	ENDLOG

    delete[] myEdges;
    return edgeCount;
}

template<class Heap>
void updateDistances(pointset vertices, std::vector< std::pair<unsigned int, double> >* myEdges, unsigned int from, double* matrix) {
	qcount++;
    unsigned int N = vertices.size();
    Heap heap(N, 0);
    for (unsigned int i = 0; i < N; i++) {
        heap.insert(i, std::numeric_limits<double>::infinity());
    }
    heap.decreaseKey(from, 0);
    while (heap.getCount() > 0) {
        std::pair<unsigned int, double> pair = heap.getMin();
        if (is_infinite(pair.second))
            break;
        matrix[from + N * pair.first] = pair.second;
        matrix[pair.first + N * from] = pair.second;
        heap.extractMin();
        for (unsigned int i = 0; i < myEdges[pair.first].size(); i++) {
            std::pair<unsigned int, double> edge = myEdges[pair.first][i];
            double alt = pair.second + edge.second;
            if (alt < heap.getValue(edge.first)) {
                heap.decreaseKey(edge.first, alt);
            }
        }
    }
}//

 //********************************************************************************************************
template<class Heap>
void SSSP(pointset vertices, std::vector< std::pair<unsigned int, double> >* myEdges, unsigned int from, double* distVector, double maxDist) {
	qcount++;
	unsigned int N = vertices.size();
	double maxDistPow2 = maxDist * maxDist;
	Heap heap(N, 0);
	for (unsigned int i = 0; i < N; i++) {
		heap.insert(i, std::numeric_limits<double>::infinity());
		distVector[i] = std::numeric_limits<double>::infinity();
	}
	distVector[from] = 0;
	heap.decreaseKey(from, 0);
	while (heap.getCount() > 0) {
		std::pair<unsigned int, double> pair = heap.getMin();
		if (is_infinite(pair.second)) //|| maxDistPow2 < NotSquaredDistance(vertices[pair.first], vertices[from]))
			break;
		distVector[pair.first] = pair.second;
		
		heap.extractMin();
		for (unsigned int i = 0; i < myEdges[pair.first].size(); i++) {
			std::pair<unsigned int, double> edge = myEdges[pair.first][i];
			double alt = pair.second + edge.second;
			if (alt < heap.getValue(edge.first)) {
				heap.decreaseKey(edge.first, alt);
			}
		}
	}
}//
 //*********************************************************************************************************
template<class Heap>
void ModifiedUpdateDistances(pointset vertices, std::vector< std::pair<unsigned int, double> >* myEdges, unsigned int from, double* matrix, double t) {
	qcount++;
	unsigned int N = vertices.size();
	Heap heap(N, 0);
	for (unsigned int i = 0; i < N; i++) {
		heap.insert(i, std::numeric_limits<double>::infinity());
	}
	heap.decreaseKey(from, 0);
	while (heap.getCount() > 0) {
		std::pair<unsigned int, double> pair = heap.getMin();
		if (is_infinite(pair.second))
			break;
		matrix[from + N * pair.first] = pair.second;
		matrix[pair.first + N * from] = pair.second;
		heap.extractMin();
		for (unsigned int i = 0; i < myEdges[pair.first].size(); i++) {
			std::pair<unsigned int, double> edge = myEdges[pair.first][i];
			double alt = pair.second + edge.second;
			if (alt < heap.getValue(edge.first) && alt <= distance(vertices[from], vertices[edge.first]) * t) {
				heap.decreaseKey(edge.first, alt);
			}
		}
	}
}
///******************************************************
template<class Heap>
static double computeDistance(const pointset &vertices, std::vector< std::pair<unsigned int, double> >* myEdges, unsigned int from, unsigned int to)
{
	qcount++;

	unsigned int N = vertices.size();
	Heap heap(N, 0);
	for (unsigned int i = 0; i < N; i++) {
		heap.insert(i, std::numeric_limits<double>::infinity());
	}
	heap.decreaseKey(from, 0);
	double result = std::numeric_limits<double>::infinity();
	while (heap.getCount() > 0) {
		std::pair<unsigned int, double> pair = heap.getMin();
		if (pair.first == to)
		{
			result = pair.second;
			break;
		}
		heap.extractMin();
		for (unsigned int i = 0; i < myEdges[pair.first].size(); i++) {
			std::pair<unsigned int, double> edge = myEdges[pair.first][i];
			double alt = pair.second + edge.second;
			if (alt < heap.getValue(edge.first)) {
				heap.decreaseKey(edge.first, alt);
			}
		}
	}
	return result;
}

//********************************************************************************************************

void ComputeGreedySpannerED(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>& edgeList) {
	STARTLOG

	//std::ofstream myfile;
	//myfile.open ("output.txt");
	//myfile << "edges" << endl;

    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
    double* matrix = new double[N * N];

    for (unsigned int x = 0; x < N; x++)
		for (unsigned int y = 0; y < N; y++)
            matrix[y + N * x] = std::numeric_limits<double>::infinity();
    for (unsigned int y = 0; y < N; y++)
        matrix[y + N * y] = 0;

    std::sort(edgeList.begin(), edgeList.end());

    std::vector<unsigned int> indices;
    indices.reserve(N-1);
    std::vector<unsigned int> indices2;
    indices2.reserve(N-1);

    for (unsigned int k = 0; k < edgeList.size(); k++) {
        double edgeDist = edgeList[k].distance;
		unsigned int ex = edgeList[k].x;
		unsigned int ey = edgeList[k].y;
        if (edgeDist * T->t < matrix[ey + N * ex]) {
			if (T->create_edge(edgeList[k])) { edgeCount = 0; goto cleanup; }
            edgeCount++;
			//myfile << std::fixed << std::setprecision(15) << edgeDist << endl;
			//int edgeDisti = (int)floor(edgeDist);
			//myfile << edgeDisti << endl;

            indices.clear();
			indices2.clear();
			for (unsigned int y = 0; y < N; y++) {
				unsigned int x = y;
				double oldDistxa = matrix[x + N * ex];
				double oldDistby = matrix[ey + N * y];
				if (oldDistxa + edgeDist < oldDistby) {
					indices.push_back(x);
				}
				if (edgeDist + oldDistby < oldDistxa) {
					indices2.push_back(y);
				}
			}
			for (unsigned int i = 0; i < indices2.size(); i++) {
				unsigned int y = indices2[i];
				double oldDistby = matrix[ey + N * y];
				for (unsigned int j = 0; j < indices.size(); j++) {
					unsigned int x = indices[j];
					double oldDistxy = matrix[x + N * y];
					double oldDistxa = matrix[x + N * ex];
					double newDist = oldDistxa + edgeDist + oldDistby;
					if (oldDistxy > newDist) {
						matrix[y + N * x] = newDist;
						matrix[x + N * y] = newDist;
					}
				}
			}
        }
    }
cleanup:
	ENDLOG
    delete[] matrix;
	//myfile.close();
}

void ComputeGreedySpannerApprox(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>& edgeList) {
	STARTLOG
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
    double* matrix = new double[N * N];

    for (unsigned int x = 0; x < N; x++)
		for (unsigned int y = 0; y < N; y++)
            matrix[y + N * x] = std::numeric_limits<double>::infinity();
    for (unsigned int y = 0; y < N; y++)
        matrix[y + N * y] = 0;

    std::sort(edgeList.begin(), edgeList.end());
    std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];

    for (unsigned int k = 0; k < edgeList.size(); k++) {
        double edgeDist = edgeList[k].distance;
		double approxDist = std::numeric_limits<double>::infinity();
		for (unsigned int i = 0; i < myEdges[edgeList[k].x].size(); i++) {
			unsigned int y = myEdges[edgeList[k].x][i].first;
			double distby = matrix[edgeList[k].x + N * y];
			for (unsigned int j = 0; j < myEdges[edgeList[k].y].size(); j++) {
				unsigned int x = myEdges[edgeList[k].y][j].first;
				double approxDistxy = matrix[x + N * y];
				double distxa = matrix[x + N * edgeList[k].y];
				double newDist = distxa + approxDistxy + distby;
				if (approxDist > newDist) {
					approxDist = newDist;
				}
			}
		}
        if (edgeDist * T->t < approxDist) {
            std::pair<unsigned int, double> pair(edgeList[k].x, edgeList[k].distance);
            myEdges[edgeList[k].y].push_back(pair);
            pair.first = edgeList[k].y;
            myEdges[edgeList[k].x].push_back(pair);
			if (T->create_edge(edgeList[k])) { edgeCount = 0; goto cleanup; }
            edgeCount++;
			matrix[edgeList[k].y + N * edgeList[k].x] = edgeDist;
			matrix[edgeList[k].x + N * edgeList[k].y] = edgeDist;
        }
		else
		{
			matrix[edgeList[k].y + N * edgeList[k].x] = approxDist;
			matrix[edgeList[k].x + N * edgeList[k].y] = approxDist;
		}
    }
cleanup:
	ENDLOG
    delete[] matrix;
}

template<class Heap, bool doDouble>
void ComputeGreedySpannerFG(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>& edgeList) {
	STARTLOG
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
    double* matrix = new double[N * N];

    for (unsigned int x = 0; x < N; x++)
		for (unsigned int y = 0; y < N; y++)
            matrix[y + N * x] = std::numeric_limits<double>::infinity();
    for (unsigned int y = 0; y < N; y++)
        matrix[y + N * y] = 0;

    std::sort(edgeList.begin(), edgeList.end());
    std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
    for (unsigned int k = 0; k < edgeList.size(); k++) {
        if (edgeList[k].distance * T->t < matrix[edgeList[k].y + N * edgeList[k].x]) {
            updateDistances<Heap>(vertices, myEdges, edgeList[k].x, matrix);
			if (doDouble) {
				updateDistances<Heap>(vertices, myEdges, edgeList[k].y, matrix);
			}
            if (edgeList[k].distance * T->t < matrix[edgeList[k].y + N * edgeList[k].x]) {
                std::pair<unsigned int, double> pair(edgeList[k].x, edgeList[k].distance);
                myEdges[edgeList[k].y].push_back(pair);
                pair.first = edgeList[k].y;
                myEdges[edgeList[k].x].push_back(pair);
				if (T->create_edge(edgeList[k])) { edgeCount = 0; goto cleanup; }
                edgeCount++;
            }
        }
    }
cleanup:
	ENDLOG
	double vm, rss;
    delete[] matrix;
    delete[] myEdges;
}

struct EDGreedySpanner {
	static void ComputeGreedySpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>& edgeList) {
		ComputeGreedySpannerED(T, edgeCount, edgeList);
	}
};

struct ApproxGreedySpanner {
	static void ComputeGreedySpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>& edgeList) {
		ComputeGreedySpannerApprox(T, edgeCount, edgeList);
	}
};

template<class Heap>
struct FGGreedySpanner {
	static void ComputeGreedySpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>& edgeList) {
		ComputeGreedySpannerFG<Heap, false>(T, edgeCount, edgeList);
	}
};

template<class Heap>
struct FGDGreedySpanner {
	static void ComputeGreedySpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>& edgeList) {
		ComputeGreedySpannerFG<Heap, true>(T, edgeCount, edgeList);
	}
};

template<class Algorithm>
int GreedySpannerOnCompleteGraph(BaseTask * T)
{
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
    std::vector<EdgeInfo> edgeList;
    edgeList.reserve((N-1) * N / 2);
    for (unsigned int y = 0; y < N; y++) {
        for (unsigned int x = y+1; x < N; x++) {
            EdgeInfo e(x, y, distance(vertices[x],vertices[y]));
            edgeList.push_back(e);
        }
    }
    unsigned int edgeCount = 0;

	Algorithm::ComputeGreedySpanner(T, edgeCount, edgeList);

    return edgeCount;
}

struct NeighborInfo
{
	int indexInLargeBucket;
	int indexInBucket;
	int bucketIndex;
};

template<class Heap>
int ComputeGreedySpannerBuckets(BaseTask * T,
								unsigned int& edgeCount,
								std::vector<std::pair<unsigned int, double> >* myEdges,
								double t,
								std::vector<int>* buckets,
								double estimatedLongestEdgeLength,
								double estimatedLongestTPath,
								int* bucketLookupX,
								int* bucketLookupY,
								int* bucketLookupZ,
								int bucketeers
								) {
    const pointset &vertices = T->p;
    int N = vertices.size();
	std::vector<int>* bucketsLarge = new std::vector<int>[bucketeers * bucketeers];
	double** matrices = new double*[bucketeers * bucketeers];
	NeighborInfo (*neighbors)[9] = new NeighborInfo[N][9];
	for (int i = 0; i < N; i++) {
		int bucketX = bucketLookupX[i];
		int bucketY = bucketLookupY[i];
		int bucketCoord = bucketX + bucketeers * bucketY;
		for (int a = -1; a <= 1; a++) {
			for (int b = -1; b <= 1; b++) {
				NeighborInfo ninfo;
				int bucketLargeIndex = bucketX + a + bucketeers * (bucketY + b);
				if (bucketX + a >= 0 && bucketX + a < bucketeers
					&&
					bucketY + b >= 0 && bucketY + b < bucketeers) {
					ninfo.indexInLargeBucket = bucketsLarge[bucketLargeIndex].size();
					ninfo.indexInBucket = bucketLookupZ[i];
					ninfo.bucketIndex = bucketLargeIndex;
					bucketsLarge[bucketLargeIndex].push_back(i);
				}
				else
				{
					ninfo.indexInLargeBucket = -1;
					ninfo.indexInBucket = -1;
					ninfo.bucketIndex = -1;
				}
				neighbors[i][(a+1) + 3 * (b+1)] = ninfo;
			}
		}
	}

	unsigned int edgeCount2 = 0;
    for (int bucketLargeCoord = 0; bucketLargeCoord < bucketeers * bucketeers; bucketLargeCoord++) {
		int bucketLargeSize = bucketsLarge[bucketLargeCoord].size();
		int bucketSize = buckets[bucketLargeCoord].size();
		for (int i = 0; i < bucketSize; i++) {
			int x = buckets[bucketLargeCoord][i];
			for (int j = 0; j < bucketLargeSize; j++) {
				int y = bucketsLarge[bucketLargeCoord][j];
				if (x < y && distance(vertices[x],vertices[y]) <= estimatedLongestEdgeLength) {
					edgeCount2++;
				}
			}
		}
	}

    std::vector<EdgeInfo> edgeList;
	edgeList.reserve(edgeCount2);

    for (int bucketLargeCoord = 0; bucketLargeCoord < bucketeers * bucketeers; bucketLargeCoord++) {
		int bucketLargeSize = bucketsLarge[bucketLargeCoord].size();
		int bucketSize = buckets[bucketLargeCoord].size();
		double* matrix = matrices[bucketLargeCoord] = new double[bucketSize * bucketLargeSize];
		for (int i = 0; i < bucketSize; i++) {
			int x = buckets[bucketLargeCoord][i];
			for (int j = 0; j < bucketLargeSize; j++) {
				int y = bucketsLarge[bucketLargeCoord][j];
				matrix[i + bucketSize * j] = std::numeric_limits<double>::infinity();
				if (x < y && distance(vertices[x],vertices[y]) <= estimatedLongestEdgeLength) {
					EdgeInfo e(x, y, distance(vertices[x],vertices[y]));
					edgeList.push_back(e);
				}
				else if (x == y) {
					matrix[i + bucketSize * j] = 0;
				}
			}
		}
    }

    std::sort(edgeList.begin(), edgeList.end());
    std::vector<unsigned int> indices;
    std::vector<unsigned int> indices2;

    for (unsigned int k = 0; k < edgeList.size(); k++) {
        double edgeDist = edgeList[k].distance;
		int bucketIndexC = neighbors[edgeList[k].x][4].bucketIndex;
		int CBucketSize = buckets[bucketIndexC].size();
		int indexCInBucketC = neighbors[edgeList[k].x][4].indexInBucket;
		for (int ni0 = 0; ni0 < 9; ni0++) {
			if (neighbors[edgeList[k].y][ni0].bucketIndex == bucketIndexC) {
				if (edgeDist * t < matrices[bucketIndexC][indexCInBucketC + CBucketSize * neighbors[edgeList[k].y][ni0].indexInLargeBucket]) {
					if (T->create_edge(edgeList[k])) return 0;
					myEdges[edgeList[k].x].push_back(std::pair<unsigned int, double>(edgeList[k].y, edgeList[k].distance));
					myEdges[edgeList[k].y].push_back(std::pair<unsigned int, double>(edgeList[k].x, edgeList[k].distance));
					edgeCount++;

					// Just added edge (c, d)
					int bucketIndexD = neighbors[edgeList[k].y][4].bucketIndex;
					int indexDInBucketD = neighbors[edgeList[k].y][4].indexInBucket;
					double* matrixC = matrices[bucketIndexC];
					double* matrixD = matrices[bucketIndexD];
					int DBucketSize = buckets[bucketIndexD].size();
					int CBucketLargeSize = bucketsLarge[bucketIndexC].size();
					int DBucketLargeSize = bucketsLarge[bucketIndexD].size();
					indices.clear();
					indices2.clear();
					for (unsigned int a = 0; a < CBucketLargeSize; a++) {
						double oldDistac = matrixC[indexCInBucketC + CBucketSize * a];
						int Ai = bucketsLarge[bucketIndexC][a];
						for (int ni = 0; ni < 9; ni++) {
							if (neighbors[Ai][ni].bucketIndex == bucketIndexD) {
								double oldDistad = matrixD[indexDInBucketD + DBucketSize * neighbors[Ai][ni].indexInLargeBucket];
								if (oldDistac + edgeDist < oldDistad) {
									indices.push_back(a);
								}
								break;
							}
						}
					}
					for (unsigned int b = 0; b < DBucketLargeSize; b++) {
						double oldDistdb = matrixD[indexDInBucketD + DBucketSize * b];
						int Bi = bucketsLarge[bucketIndexD][b];
						for (int ni = 0; ni < 9; ni++) {
							if (neighbors[Bi][ni].bucketIndex == bucketIndexC) {
								double oldDistcb = matrixC[indexCInBucketC + CBucketSize * neighbors[Bi][ni].indexInLargeBucket];
								if (edgeDist + oldDistdb < oldDistcb) {
									indices2.push_back(b);
								}
								break;
							}
						}
					}

					for (unsigned int i = 0; i < indices2.size(); i++) {
						unsigned int b = indices2[i];
						int Bi = bucketsLarge[bucketIndexD][b];
						double oldDistdb = matrixD[indexDInBucketD + DBucketSize * b];
						int bucketIndexB = neighbors[Bi][4].bucketIndex;
						double* matrixB = matrices[bucketIndexB];
						int indexBInBucketB = neighbors[Bi][4].indexInBucket;
						int BBucketSize = buckets[bucketIndexB].size();
						for (unsigned int j = 0; j < indices.size(); j++) {
							unsigned int a = indices[j];
							int Ai = bucketsLarge[bucketIndexC][a];
							int bucketIndexA = neighbors[Ai][4].bucketIndex;
							for (int ni = 0; ni < 9; ni++) {
								if (neighbors[Bi][ni].bucketIndex == bucketIndexA) {
									double oldDistac = matrixC[indexCInBucketC + CBucketSize * a];
									double* matrixA = matrices[bucketIndexA];
									int indexAInBucketA = neighbors[Ai][4].indexInBucket;
									int ABucketSize = buckets[bucketIndexA].size();
									double oldDistab = matrixA[indexAInBucketA + ABucketSize * neighbors[Bi][ni].indexInLargeBucket];
									double newDist = oldDistac + edgeDist + oldDistdb;
									if (oldDistab > newDist) {
										matrixA[indexAInBucketA + ABucketSize * neighbors[Bi][ni].indexInLargeBucket] = newDist;
										for (int ni2 = 0; ni2 < 9; ni2++) {
											if (neighbors[Ai][ni2].bucketIndex == bucketIndexB) {
												matrixB[indexBInBucketB + BBucketSize * neighbors[Ai][ni2].indexInLargeBucket] = newDist;
												break;
											}
										}
									}
									break;
								}
							}
						}
					}
				}
				break;
			}
		}
    }

	delete[] bucketsLarge;
	delete[] neighbors;
	for (int i = 0; i < bucketeers * bucketeers; i++)
		delete[] matrices[i];
	delete[] matrices;
    return edgeCount;
}

inline double minimum(const double a, const double b) {
	return a < b ? a : b;
}

inline double maximum(const double a, const double b) {
	return a > b ? a : b;
}

inline unsigned int maximum(const unsigned int a, const unsigned int b) {
	return a > b ? a : b;
}

inline bool hypComp(const std::pair<unsigned int, double>& a, const std::pair<unsigned int, double>& b) { return a.second < b.second; }

int HyperbolaSpanner(BaseTask * T) {
	STARTLOG
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
	unsigned int edgeCount = 0;

	std::vector<std::pair<unsigned int, double> > points;
	std::vector<std::pair<unsigned int, double> > edges;
	for (unsigned int i = 0; i < N; i++) {
		points.clear();
		edges.clear();
		for (unsigned int j = 0; j < N; j++) {
			if (distance(vertices[i], vertices[j]) > 1e-10)
				points.push_back(std::pair<unsigned int, double>(j, distance(vertices[i], vertices[j])));
		}
		std::sort(points.begin(), points.end(), hypComp);
		for (unsigned int k = 0; k < points.size(); k++) {
			std::pair<unsigned int, double> min = points[k];
			unsigned int v = min.first;

			bool addEdge = true;
			for (unsigned int j = 0; j < edges.size(); j++) {
				unsigned int u = edges[j].first;
				double iu = edges[j].second;
				double uv = distance(vertices[u], vertices[v]);
				double iv = min.second;
				if (iu + uv * T->t <= iv * T->t) {
					addEdge = false;
					break;
				}
			}

			if (addEdge) {
				edge e;
				e.x = i;
				e.y = min.first;
				if (T->create_edge(e)) return 0;
				edgeCount++;
				edges.push_back(min);
			}
		}
	}
	ENDLOG
    return edgeCount;
}

int HyperbolaSpannerFast(BaseTask * T) {
	STARTLOG
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
	unsigned int edgeCount = 0;

	std::vector<std::pair<unsigned int, double> > edges;
	for (unsigned int i = 0; i < N; i++) {
		edges.clear();
		bool done = false;
		while (!done) {
			done = true;
			unsigned int closest = 0;
			double closestDist = std::numeric_limits<double>::infinity();
			for (unsigned int v = 0; v < N; v++) {
				if (v != i) {
					bool addEdge = true;
					double iv = distance(vertices[i], vertices[v]);
					for (unsigned int j = 0; j < edges.size(); j++) {
						unsigned int u = edges[j].first;
						double iu = edges[j].second;
						double uv = distance(vertices[u], vertices[v]);
						if (iu + uv * T->t <= iv * T->t) {
							addEdge = false;
							break;
						}
					}

					if (addEdge && closestDist > iv) {
						closestDist = iv;
						closest = v;
						done = false;
					}
				}
			}
			if (!done) {
				edge e;
				e.x = i;
				e.y = closest;
				if (T->create_edge(e)) return 0;
				edgeCount++;
				edges.push_back(std::pair<unsigned int, double>(closest, closestDist));
			}
		}
	}
	ENDLOG
    return edgeCount;
}

bool liesInNoHyperbola(EdgeInfo e, double t, std::vector<std::pair<unsigned int, double> >* myEdges, const pointset &vertices) {
	unsigned int v = e.x;
	unsigned int w = e.y;
	for (unsigned int j = 0; j < myEdges[v].size(); j++) {
		unsigned int u = myEdges[v][j].first;
		double vu = myEdges[v][j].second;
		double uw = distance(vertices[u], vertices[w]);
		double vw = e.distance;
		if (vu + uw * t <= vw * t)
			return false;
	}
	return true;
}

void ComputeSymmetricHyperbolaSpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>& edgeList) {
	STARTLOG
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
    std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
    std::sort(edgeList.begin(), edgeList.end());
    for (unsigned int k = 0; k < edgeList.size(); k++) {
		EdgeInfo e = edgeList[k];
		unsigned int v = e.x;
		unsigned int w = e.y;
		if (v != w
			&&
			liesInNoHyperbola(e, T->t, myEdges, vertices)
			&&
			liesInNoHyperbola(e.Invert(), T->t, myEdges, vertices)) {
			if (T->create_edge(e)) {
				edgeCount = 0;
				goto cleanup;
			}
			edgeCount++;
			myEdges[v].push_back(std::pair<unsigned int, double>(w, e.distance));
			myEdges[w].push_back(std::pair<unsigned int, double>(v, e.distance));
		}
	}
	cleanup:
	ENDLOG
	delete[] myEdges;
}

int SymmetricHyperbolaSpanner(BaseTask * T) {
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
    std::vector<EdgeInfo> edgeList;

    //edgeList.reserve((N-1) * N / 2);
    for (unsigned int y = 0; y < N; y++) {
        for (unsigned int x = y+1; x < N; x++) {
            EdgeInfo e(x, y, distance(vertices[x],vertices[y]));
			edgeList.push_back(e);
        }
    }

    unsigned int edgeCount = 0;
	ComputeSymmetricHyperbolaSpanner(T, edgeCount, edgeList);
    return edgeCount;
}

template<class Heap, bool Insert>
void FindClosestUncoveredPoint(unsigned int i, double t, std::vector<int>& QContent, Heap& Q, std::vector<std::pair<unsigned int, double> >* myEdges, const pointset &vertices) {
    unsigned int N = vertices.size();
	bool foundNew = false;
	unsigned int closest = 0;
	double closestDist = std::numeric_limits<double>::infinity();
	for (unsigned int v = 0; v < N; v++) {
		if (v != i) {
			double iv = distance(vertices[i], vertices[v]);
			EdgeInfo e(i, v, iv);
			if (closestDist > iv
				&&
				liesInNoHyperbola(e, t, myEdges, vertices)
				&&
				liesInNoHyperbola(e.Invert(), t, myEdges, vertices)) {
				closestDist = iv;
				closest = v;
				foundNew = true;
			}
		}
	}
	if (foundNew) {
		QContent[i] = closest;
		if (Insert)
			Q.insert(i, closestDist);
		else
			Q.increaseKey(i, closestDist);
	}
	else
		Q.extractMin();
}

// MinHeap, T->t, extractMin

template<class Heap>
int SymmetricHyperbolaSpannerFast(BaseTask * T) {
	STARTLOG
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
	double t = T->t;

	std::vector<int> QContent(N, -1);
	Heap Q(N, -std::numeric_limits<double>::infinity());
	unsigned int edgeCount = 0;
    std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	
	for (unsigned int i = 0; i < N; i++)
		FindClosestUncoveredPoint<Heap, true>(i, t, QContent, Q, myEdges, vertices);
	
	while (Q.getCount() > 0) {
		bool foundOne = false;
		std::pair<unsigned int, double> minQ;
		while (!foundOne)
		{
			minQ = Q.getMin();
			if (liesInNoHyperbola(EdgeInfo(QContent[minQ.first], minQ.first, minQ.second), t, myEdges, vertices))
				foundOne = true;
			else
				FindClosestUncoveredPoint<Heap, false>(minQ.first, t, QContent, Q, myEdges, vertices);
			if (Q.getCount() == 0)
				goto cleanup;
		}

		int minRepFrom = minQ.first;
		int minRepTo = QContent[minQ.first];
		double minRepDist = minQ.second;

		std::pair<unsigned int, double> pair(minRepFrom, minRepDist);
        myEdges[minRepTo].push_back(pair);
        pair.first = minRepTo;
        myEdges[minRepFrom].push_back(pair);
        if (T->create_edge(edge(minRepFrom, minRepTo))) return 0;
        edgeCount++;
		FindClosestUncoveredPoint<Heap, false>(minQ.first, t, QContent, Q, myEdges, vertices);
	}
	cleanup:
	ENDLOG
	delete[] myEdges;
    return edgeCount;
}

struct MyData {
	double intersectionX2;
	vertex point;
	unsigned int pointIndex;
	vertex minimum;
	unsigned int minimumIndex;
};

class MyAugmentation {
public:
	static void process(ScapegoatNode<double, MyData>* current) {
		vertex minV = current->value.point;
		unsigned int minI = current->value.pointIndex;
		if (current->left != nullptr && current->left->value.minimum.x < minV.x) {
			minV = current->left->value.minimum;
			minI = current->left->value.minimumIndex;
		}
		if (current->right != nullptr && current->right->value.minimum.x < minV.x) {
			minV = current->right->value.minimum;
			minI = current->right->value.minimumIndex;
		}
		current->value.minimum = minV;
		current->value.minimumIndex = minI;
	}
};

class MySearchResult {
public:
	MySearchResult() : result(std::numeric_limits<double>::infinity()), foundSomething(false) {
	}
	double result;
	unsigned int resultI;
	bool foundSomething;
	void update(ScapegoatNode<double, MyData>* x, bool goingLeft) {
		if (goingLeft) {
			if (result >= x->value.point.x) {
				result = x->value.point.x;
				resultI = x->value.pointIndex;
				foundSomething = true;
			}
			if (x->right != nullptr && result >= x->right->value.minimum.x) {
				result = x->right->value.minimum.x;
				resultI = x->right->value.minimumIndex;
				foundSomething = true;
			}
		}
	}
};

bool myOrdering(const ScapegoatNode<double, MyData>* a, const ScapegoatNode<double, MyData>* b) {
	return a->value.intersectionX2 < b->value.intersectionX2;
}

inline vertex rotate(vertex v, double theta) {
	return vertex(v.x * cos(theta) - v.y * sin(theta), v.x * sin(theta) + v.y * cos(theta));
}

void ComputeThetaSpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>* edges) {
	STARTLOG
	const pointset &vertices = T->p;
    unsigned int N = vertices.size();
	unsigned int k = 8;
	double theta;
	if (T->t > 0) {
		if (T->t >= 8 * (29 + 23 * sqrt(2.0))) {
			k = 4;
		}
		else if (T->t >= 2) {
			k = 6;
		}
		else {
			double t = std::numeric_limits<double>::infinity();
			theta = std::numeric_limits<double>::infinity();
			while (t > T->t) {
				k++;
				theta =  (2.0 * M_PI) / ((double)k);
				t = 1.0 / (cos(theta) - sin(theta));
			}
		}
	}
	else {
		k = (int)-T->t;
	}
	theta = (2.0 * M_PI) / ((double)k);
	std::vector<vertex> verticesCopy(N);
	for (unsigned int i = 0; i < N; i++)
		verticesCopy[i] = vertices[i];

	double eps = 1.01;
	double lineA = -tan(theta / 2.0 * eps);
	double lineAR = tan(theta / 2.0 * eps);
	for (unsigned int c = 0; c < k; c++) {
		for (unsigned int i = 0; i < N; i++)
			verticesCopy[i] = rotate(verticesCopy[i], theta);

		ScapegoatHeap<double, MyData, MyAugmentation> heap;
		std::vector<ScapegoatNode<double, MyData>*> nodes;

		for (unsigned int i = 0; i < N; i++) {
			vertex v = verticesCopy[i];
			double intersectionX = v.x - v.y / lineA;
			MyData data;
			data.point = v;
			data.pointIndex = i;
			data.intersectionX2 = v.x - v.y / lineAR;
			nodes.push_back(heap.insert(intersectionX, data));
		}

		std::sort(nodes.begin(), nodes.end(), myOrdering);

		for (unsigned int i = 0; i < nodes.size(); i++) {
			ScapegoatNode<double, MyData>* node = nodes[i];
			heap.remove(node);
			MySearchResult result;
			heap.augmentedSearch(node->key, result);

			if (result.foundSomething) {
				if (edges == nullptr) {
					edge e;
					e.x = node->value.pointIndex;
					e.y = result.resultI;
					if (T->create_edge(e)) { edgeCount = 0; return; }
				}
				else {
					int ex = node->value.pointIndex;
					int ey = result.resultI;
					EdgeInfo e(ex, ey, distance(vertices[ex], vertices[ey]));
					edges->push_back(e);
				}
				edgeCount++;
			}
			node->reset();
			delete node;
		}
	}
	ENDLOG
}

int ThetaSpanner(BaseTask * T) {
    unsigned int edgeCount = 0;
	ComputeThetaSpanner(T, edgeCount, nullptr);
    return edgeCount;
}

struct BoundingBox {
	BoundingBox(const pointset &vertices, std::vector<unsigned int>* indices, unsigned int begin, unsigned int end) {
		left = std::numeric_limits<double>::infinity();
		right = -std::numeric_limits<double>::infinity();
		bottom = std::numeric_limits<double>::infinity();
		top = -std::numeric_limits<double>::infinity();
		maxDist = 0;
		for (unsigned int i = begin; i != end; i++) {
			vertex v = vertices[(*indices)[i]];
			if (left > v.x)
				left = v.x;
			if (right < v.x)
				right = v.x;
			if (bottom > v.y)
				bottom = v.y;
			if (top < v.y)
				top = v.y;
			if (distance(v, center()) > maxDist)
				maxDist = distance(v, center());
		}
		longest = maximum(top - bottom, right - left);
	}
	double left, right, bottom, top, longest, maxDist;
	inline double lrMean() const {
		return (right + left) / 2.0;
	}
	inline double btMean() const {
		return (top + bottom) / 2.0;
	}
	inline bool contains(vertex v) const {
		return v.x >= left && v.x <= right && v.y >= bottom && v.y <= top;
	}
	inline bool leftRightSplit() const {
		return right - left >= top - bottom;
	}
	inline vertex center() const {
		vertex result(lrMean(), btMean());
		return result;
	}
	inline double size() const {
		return maxDist;
	}
	inline double radius() const {
		vertex corner(left, bottom);
		return distance(corner, center());
	}
};

double distance(BoundingBox* leftBoundingBox, BoundingBox* rightBoundingBox) {
	return distance(leftBoundingBox->center(), rightBoundingBox->center());
}

double maxDistance(BoundingBox* leftBoundingBox, BoundingBox* rightBoundingBox) {
	return distance(leftBoundingBox, rightBoundingBox) + leftBoundingBox->radius() + rightBoundingBox->radius();
}

double minDistance(BoundingBox* leftBoundingBox, BoundingBox* rightBoundingBox) {
	return distance(leftBoundingBox, rightBoundingBox) - leftBoundingBox->radius() - rightBoundingBox->radius();
}

bool WellSeparated(BoundingBox* leftBoundingBox, BoundingBox* rightBoundingBox, const double s) {
	double R = maximum(leftBoundingBox->radius(), rightBoundingBox->radius());
	return distance(leftBoundingBox, rightBoundingBox) >= (s + 2) * R;
}
/*
double NotSquaredDistance(BoundingBox* leftBoundingBox, BoundingBox* rightBoundingBox) {
	return NotSquaredDistance(leftBoundingBox->center(), rightBoundingBox->center());
}

double NotSquaredMaxDistance(BoundingBox* leftBoundingBox, BoundingBox* rightBoundingBox) {
	return NotSquaredDistance(leftBoundingBox, rightBoundingBox) + leftBoundingBox->radius() + rightBoundingBox->radius();
}

double NotSquaredMinDistance(BoundingBox* leftBoundingBox, BoundingBox* rightBoundingBox) {
	return NotSquaredDistance(leftBoundingBox, rightBoundingBox) - leftBoundingBox->radius() - rightBoundingBox->radius();
}*/

template<bool IsImplicit>
class SplitTree;

template<bool IsImplicit>
struct WspdEdge {
	SplitTree<IsImplicit>* first;
	SplitTree<IsImplicit>* second;
	double dist;
    inline bool operator < (const WspdEdge &o) const {
        return dist < o.dist;
	}
};

template<bool IsImplicit>
struct WellSeparatedPair {
	WellSeparatedPair(SplitTree<IsImplicit>* first, SplitTree<IsImplicit>* second, int firstSize) : first(first), second(second)
	{
		interestingCache = new bool[firstSize];
		for (int i = 0; i < firstSize; i++) {
			interestingCache[i] = true;
		}
		shortestPathSoFar = std::numeric_limits<double>::infinity();
		shortestPathSoFarIndex = -1;
	}
	SplitTree<IsImplicit>* first;
	SplitTree<IsImplicit>* second;
	bool* interestingCache;
	double shortestPathSoFar;
	int shortestPathSoFarIndex;
	int shortestPathSoFarOtherIndex;
	double length() const;
	double minlength() const;
	double maxlength() const;
    inline bool operator < (const WellSeparatedPair &o) const {
        return minlength() < o.minlength();
	}
	~WellSeparatedPair() {
		delete[] interestingCache;
	}
};

template<bool IsImplicit>
struct WellSeparatedPairPointerComparer {
	inline bool operator()(WellSeparatedPair<IsImplicit> const * i, WellSeparatedPair<IsImplicit> const * j) const {
		return (*i) < (*j);
	}
};

// Fix: parent + count
template<bool IsImplicit>
class SplitTree {
private:
	unsigned int representative;
	SplitTree<IsImplicit>* left;
	SplitTree<IsImplicit>* right;
	unsigned int starti;
	unsigned int endi;
	SplitTree() {
	}
	SplitTree(unsigned int representative, BoundingBox* box, SplitTree<IsImplicit>* left, SplitTree<IsImplicit>* right, std::vector<unsigned int>* indices, unsigned int start, unsigned int end)
	: representative(representative), box(box), left(left), right(right), indices(indices), starti(start), endi(end), isTop(false) { }
	std::vector<unsigned int>* indices;
	bool isTop;
public:
	unsigned int count() {
		return endi - starti;
	}
	class Iterator;
	Iterator begin();
	Iterator end();
	static SplitTree<IsImplicit>* createSplitTreeRec(const pointset &vertices, std::vector<unsigned int>* indices, unsigned int start, unsigned int end, unsigned int& rightmost);
	static SplitTree<IsImplicit>* createSplitTree(const pointset &vertices) {
		std::vector<unsigned int>* indices = new std::vector<unsigned int>();
		for (unsigned int i = 0; i < vertices.size(); i++)
			indices->push_back(i);
		unsigned int rightmost;
		SplitTree<IsImplicit>* result = createSplitTreeRec(vertices, indices, 0, indices->size(), rightmost);
		result->isTop = true;
		return result;
	}
	~SplitTree() {
		delete box;
		if (left != nullptr) {
			delete left;
		}
		if (right != nullptr){
			delete right;
		}
		//if (isTop)
			delete indices;
	}
	BoundingBox* box;
	void AddPairs(SplitTree<IsImplicit>* other, BaseTask * T, double s, unsigned int &edgeCount, std::vector<EdgeInfo>* edges, const pointset &vertices) {
		if (WellSeparated(box, other->box, s)) {
			if (edges == nullptr) {
				edge e;
				e.x = representative;
				e.y = other->representative;
				T->create_edge(e);
			}
			else {
				int ex = representative;
				int ey = other->representative;
				EdgeInfo e(ex, ey, distance(vertices[ex], vertices[ey]));
				edges->push_back(e);
			}
			edgeCount++;
		}
		else {
			if (box->longest > other->box->longest) {
				left->AddPairs(other, T, s, edgeCount, edges, vertices);
				right->AddPairs(other, T, s, edgeCount, edges, vertices);
			}
			else {
				other->left->AddPairs(this, T, s, edgeCount, edges, vertices);
				other->right->AddPairs(this, T, s, edgeCount, edges, vertices);
			}
		}
	}
	void ComputeWspdSpanner(BaseTask * T, double s, unsigned int &edgeCount, std::vector<EdgeInfo>* edges, const pointset &vertices) {
		if (left != nullptr) {
			left->AddPairs(right, T, s, edgeCount, edges, vertices);
			left->ComputeWspdSpanner(T, s, edgeCount, edges, vertices);
			right->ComputeWspdSpanner(T, s, edgeCount, edges, vertices);
		}
	}
	std::vector<std::pair<SplitTree<IsImplicit>*, double> > pairs;
	void AddHypPairs(SplitTree<IsImplicit>* other, double s, std::vector<WspdEdge<IsImplicit>>* wspdEdges) {
		if (WellSeparated(box, other->box, s)) {
			double dist = distance(box->center(), other->box->center());
			if (wspdEdges == nullptr) {
				pairs.push_back(std::pair<SplitTree<IsImplicit>*, double>(other, dist));
				other->pairs.push_back(std::pair<SplitTree<IsImplicit>*, double>(this, dist));
			}
			else {
				WspdEdge<IsImplicit> e;
				e.first = this;
				e.second = other;
				e.dist = dist;
				wspdEdges->push_back(e);
			}
		}
		else {
			if (box->longest > other->box->longest) {
				if (left != nullptr) {
					left->AddHypPairs(other, s, wspdEdges);
					right->AddHypPairs(other, s, wspdEdges);
				}
			}
			else {
				if (other->left != nullptr) {
					other->left->AddHypPairs(this, s, wspdEdges);
					other->right->AddHypPairs(this, s, wspdEdges);
				}
			}
		}
	}
	void ComputeDecomposition(double s, std::vector<WspdEdge<IsImplicit>>* wspdEdges = nullptr) {
		if (left != nullptr) {
			left->AddHypPairs(right, s, wspdEdges);
			left->ComputeDecomposition(s, wspdEdges);
			right->ComputeDecomposition(s, wspdEdges);
		}
	}
	std::vector<std::pair<SplitTree<IsImplicit>*, double> > edges;
	static inline bool wspdComp(const std::pair<SplitTree<IsImplicit>*, double>& a, const std::pair<SplitTree<IsImplicit>*, double>& b) { return a.second < b.second; }
	void AddWspdHypEdges(BaseTask * T, double s, unsigned int& edgeCount, bool symmetric) {
		std::sort(pairs.begin(), pairs.end(), wspdComp);
		for (unsigned int k = 0; k < pairs.size(); k++) {
			SplitTree<IsImplicit>* min = pairs[k].first;
			double dist = pairs[k].second;

			bool addEdge = true;
			for (unsigned int j = 0; j < edges.size(); j++) {
				SplitTree<IsImplicit>* u = edges[j].first;
				double tu = edges[j].second;
				double um = distance(u->box->center(), min->box->center());
				double tm = dist;
				if (tu * (1 + 2 * T->t / s) + um * T->t * (1 + 2 / s) <= tm * T->t) {
					addEdge = false;
					break;
				}
			}
			if (addEdge) {
				if (symmetric) {
					for (unsigned int j = 0; j < min->edges.size(); j++) {
						SplitTree<IsImplicit>* u = min->edges[j].first;
						double tu = min->edges[j].second;
						double um = distance(u->box->center(), min->box->center());
						double tm = dist;
						if (tu * (1 + 2 * T->t / s) + um * T->t * (1 + 2 / s) <= tm * T->t) {
							addEdge = false;
							break;
						}
					}
				}
				if (addEdge) {
					edge e;
					e.x = representative;
					e.y = min->representative;
					T->create_edge(e);
					edgeCount++;
					edges.push_back(std::pair<SplitTree<IsImplicit>*, double>(min, dist));
				}
			}
		}

		if (left != nullptr) {
			left->AddWspdHypEdges(T, s, edgeCount, symmetric);
			right->AddWspdHypEdges(T, s, edgeCount, symmetric);
		}
	}
	static void AddWspdSymHypEdges(BaseTask * T, double s, unsigned int& edgeCount, std::vector<WspdEdge<IsImplicit>>& wspdEdges) {
		std::sort(wspdEdges.begin(), wspdEdges.end());
		for (unsigned int k = 0; k < wspdEdges.size(); k++) {
			SplitTree<IsImplicit>* min1 = wspdEdges[k].first;
			SplitTree<IsImplicit>* min2 = wspdEdges[k].second;
			double dist = wspdEdges[k].dist;

			bool addEdge = true;
			for (unsigned int j = 0; j < min1->edges.size(); j++) {
				SplitTree<IsImplicit>* u = min1->edges[j].first;
				double tu = min1->edges[j].second;
				double um = distance(u->box->center(), min2->box->center());
				double tm = dist;
				if (tu * (1 + 2 * T->t / s) + um * T->t * (1 + 2 / s) <= tm * T->t) {
					addEdge = false;
					break;
				}
			}
			if (addEdge) {
				for (unsigned int j = 0; j < min2->edges.size(); j++) {
					SplitTree<IsImplicit>* u = min2->edges[j].first;
					double tu = min2->edges[j].second;
					double um = distance(u->box->center(), min1->box->center());
					double tm = dist;
					if (tu * (1 + 2 * T->t / s) + um * T->t * (1 + 2 / s) <= tm * T->t) {
						addEdge = false;
						break;
					}
				}
				if (addEdge) {
					edge e;
					e.x = min1->representative;
					e.y = min2->representative;
					T->create_edge(e);
					edgeCount++;
					min1->edges.push_back(std::pair<SplitTree<IsImplicit>*, double>(min2, dist));
				}
			}
		}
	}
	void ComputeWspdHypSpanner(BaseTask * T, double s, unsigned int& edgeCount, bool symmetric = false) {
		ComputeDecomposition(s);
		AddWspdHypEdges(T, s, edgeCount, symmetric);
	}
	void ComputeWspdSymHypSpanner(BaseTask * T, double s, unsigned int& edgeCount) {
		std::vector<WspdEdge<IsImplicit> > wspdEdges;
		ComputeDecomposition(s, &wspdEdges);
		AddWspdSymHypEdges(T, s, edgeCount, wspdEdges);
	}
	void GetPairs(SplitTree<IsImplicit>* other, double s, std::vector<WellSeparatedPair<IsImplicit>*>& pairs, double sizeLowerBound) {
		if (maxDistance(box, other->box) >= sizeLowerBound) {
			if (WellSeparated(box, other->box, s)) {
				WellSeparatedPair<IsImplicit>* p;
				if (count() <= other->count())
				{
					p = new WellSeparatedPair<IsImplicit>(this, other, count());
				}
				else
				{
					p = new WellSeparatedPair<IsImplicit>(other, this, other->count());
				}
				pairs.push_back(p);
			}
			else {
				if (box->longest > other->box->longest) {
					if (left != nullptr) {
						left->GetPairs(other, s, pairs, sizeLowerBound);
						right->GetPairs(other, s, pairs, sizeLowerBound);
					}
				}
				else {
					if (other->left != nullptr) {
						other->left->GetPairs(this, s, pairs, sizeLowerBound);
						other->right->GetPairs(this, s, pairs, sizeLowerBound);
					}
				}
			}
		}
	}
	void ComputeWspdPairs(double s, std::vector<WellSeparatedPair<IsImplicit>*>& pairs, double sizeLowerBound = 0.0) {
		if (left != nullptr) {
			left->GetPairs(right, s, pairs, sizeLowerBound);
			left->ComputeWspdPairs(s, pairs, sizeLowerBound);
			right->ComputeWspdPairs(s, pairs, sizeLowerBound);
		}
	}
};

template<>
class SplitTree<false>::Iterator {
private:
	unsigned int i;
	SplitTree<false>* owner;
public:
	Iterator(SplitTree<false>* owner, unsigned int i) : owner(owner), i(i) {
	}
	bool operator==(const Iterator& rhs) {
		return owner == rhs.owner && i == rhs.i;
	}
	bool operator!=(const Iterator& rhs) {
		return owner != rhs.owner || i != rhs.i;
	}
	void operator++() {
		i++;
	}
	unsigned int operator*()
	{
		return (*owner->indices)[i];
	}
};

template<>
SplitTree<false>::Iterator SplitTree<false>::begin() {
	return SplitTree<false>::Iterator(this, 0);
}
template<>
SplitTree<false>::Iterator SplitTree<false>::end() {
	return SplitTree<false>::Iterator(this, this->indices->size());
}

template<>
class SplitTree<true>::Iterator {
private:
	unsigned int i;
	SplitTree<true>* owner;
public:
	Iterator(SplitTree<true>* owner, unsigned int i) : owner(owner), i(i) {
	}
	bool operator==(const Iterator& rhs) {
		return owner == rhs.owner && i == rhs.i;
	}
	bool operator!=(const Iterator& rhs) {
		return owner != rhs.owner || i != rhs.i;
	}
	void operator++() {
		i++;
	}
	unsigned int operator*()
	{
		return (*owner->indices)[i];
	}
};

template<>
SplitTree<true>::Iterator SplitTree<true>::begin() {
	return Iterator(this, starti);
}
template<>
SplitTree<true>::Iterator SplitTree<true>::end() {
	return Iterator(this, endi);
}

template<>
SplitTree<false>* SplitTree<false>::createSplitTreeRec(const pointset &vertices, std::vector<unsigned int>* indices, unsigned int start, unsigned int end, unsigned int& rightmost) {
	unsigned int representative;
	vertex representativeV;
	BoundingBox* box = new BoundingBox(vertices, indices, 0, indices->size());
	SplitTree<false>* left = nullptr;
	SplitTree<false>* right = nullptr;
	if (indices->size() > 1) {
		std::vector<unsigned int>* pointsLeft = new std::vector<unsigned int>();
		std::vector<unsigned int>* pointsRight = new std::vector<unsigned int>();
		for (unsigned int i = 0; i < indices->size(); i++) {
			vertex v = vertices[(*indices)[i]];
			if (box->leftRightSplit()) {
				if (v.x < box->lrMean())
					pointsLeft->push_back((*indices)[i]);
				else
					pointsRight->push_back((*indices)[i]);
			}
			else {
				if (v.y < box->btMean())
					pointsLeft->push_back((*indices)[i]);
				else
					pointsRight->push_back((*indices)[i]);
			}
		}
		left = createSplitTreeRec(vertices, pointsLeft, 0, pointsLeft->size(), rightmost);
		representative = rightmost;
		right = createSplitTreeRec(vertices, pointsRight, 0, pointsRight->size(), rightmost);

		//delete pointsLeft;
		//delete pointsRight;
	}
	else if (indices->size() == 1) {
		left = right = nullptr;
		rightmost = representative = (*indices)[0];
	}
	else
		std::cout << "FAAL" << endl;
	return new SplitTree(representative, box, left, right, indices, 0, indices->size());
}

struct SortIndicesByX {
	const pointset &vertices;
	SortIndicesByX(const pointset &vertices) : vertices(vertices) {
	}
    bool operator()(const unsigned int& a, const unsigned int& b) {
		return vertices[a].x < vertices[b].x;
    }
};

struct SortIndicesByY {
	const pointset &vertices;
	SortIndicesByY(const pointset &vertices) : vertices(vertices) {
	}
    bool operator()(const unsigned int& a, const unsigned int& b) {
		return vertices[a].y < vertices[b].y;
    }
};

template<>
SplitTree<true>* SplitTree<true>::createSplitTreeRec(const pointset &vertices, std::vector<unsigned int>* indices, unsigned int start, unsigned int end, unsigned int& rightmost) {
	unsigned int representative;
	vertex representativeV;
	BoundingBox* box = new BoundingBox(vertices, indices, start, end);
	SplitTree<true>* left = nullptr;
	SplitTree<true>* right = nullptr;
	if (end - start > 1) {
		if (box->leftRightSplit())
			sort(indices->begin() + start, indices->begin() + end, SortIndicesByX(vertices));
		else
			sort(indices->begin() + start, indices->begin() + end, SortIndicesByY(vertices));
		unsigned int index = start;
		for (; index < end; index++) {
			vertex v = vertices[(*indices)[index]];
			if (box->leftRightSplit()) {
				if (v.x > box->lrMean() || (v.x == box->lrMean() && index != start))
					break;
			}
			else {
				if (v.y > box->btMean() || (v.y == box->btMean() && index != start))
					break;
			}
		}
		left = createSplitTreeRec(vertices, indices, start, index, rightmost);
		representative = rightmost;
		right = createSplitTreeRec(vertices, indices, index, end, rightmost);
		return new SplitTree<true>(representative, box, left, right, indices, start, end);
	}
	else if (end - start == 1) {
		left = right = nullptr;
		rightmost = representative = (*indices)[start];
	}
	else
		std::cout << "FAAL" << endl;
	return new SplitTree<true>(representative, box, left, right, indices, start, end);
}

template<bool IsImplicit>
double dist(SplitTree<IsImplicit>* first, SplitTree<IsImplicit>* second) {
	return distance(first->box, second->box);
}

template<bool IsImplicit>
double WellSeparatedPair<IsImplicit>::minlength() const {
	return minDistance(first->box, second->box);
}

template<bool IsImplicit>
double WellSeparatedPair<IsImplicit>::maxlength() const {
	return maxDistance(first->box, second->box);
}

template<bool IsImplicit>
double WellSeparatedPair<IsImplicit>::length() const {
	return dist(first, second);
}

template<bool IsImplicit>
void ComputeWspdSpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>* edges) {
	STARTLOG
    const pointset &vertices = T->p;
	double s = 4.0 * (T->t + 1.0) / (T->t - 1.0);
	SplitTree<IsImplicit>* splitTree = SplitTree<IsImplicit>::createSplitTree(vertices);
	splitTree->ComputeWspdSpanner(T, s, edgeCount, edges, vertices);
	ENDLOG
	delete splitTree;
}

template<bool IsImplicit>
int WspdSpanner(BaseTask * T) {
    unsigned int edgeCount = 0;
	ComputeWspdSpanner<IsImplicit>(T, edgeCount, nullptr);
    return edgeCount;
}

struct WellSeparatedPairRepresentative {
	WellSeparatedPairRepresentative(int u, int v, double dist, int pairIndex) : u(u), v(v), dist(dist), pairIndex(pairIndex) { }
	int u, v;
	double dist;
	int pairIndex;
};

template<class Heap, bool IsImplicit>
int closestDilationHigherT(const pointset& vertices,
						   std::vector< std::pair<unsigned int, double> >* myEdges,
						   unsigned int from,
						   double t,
						   SplitTree<IsImplicit>* set,
						   Heap& myHeap,
						   WellSeparatedPair<IsImplicit>* wspair,
						   double* realDistance,
						   bool* mySet
						   ) {
	qcount++;
	unsigned int N = vertices.size();
	unsigned int foundCount = 0;
	myHeap.clear();
	for (unsigned int i = 0; i < N; i++) {
		realDistance[i] = std::numeric_limits<double>::infinity();
		myHeap.insert(i, std::numeric_limits<double>::infinity());
	}
	realDistance[from] = 0;
    myHeap.decreaseKey(from, 0);
	int result = -1;
	double minDist = std::numeric_limits<double>::infinity();
    while (myHeap.getCount() > 0) {
        std::pair<unsigned int, double> pair = myHeap.getMin();
        if (is_infinite(pair.second))
            break;
        myHeap.extractMin();

		if (mySet[pair.first]) {
			double dist = distance(vertices[pair.first], vertices[from]);
			double dilation = realDistance[pair.first] / dist;
			if (realDistance[pair.first] < wspair->shortestPathSoFar) {
				wspair->shortestPathSoFar = realDistance[pair.first];
				wspair->shortestPathSoFarIndex = from;
				wspair->shortestPathSoFarOtherIndex = pair.first;
			}
			if (dilation > t
				&&
				dist < minDist) {
				minDist = dist;
				result = pair.first;
			}
			foundCount++;
			if (foundCount == set->count()) {
				return result;
			}
		}

		for (unsigned int i = 0; i < myEdges[pair.first].size(); i++) {
            std::pair<unsigned int, double> edge = myEdges[pair.first][i];
            double alt = realDistance[pair.first] + edge.second;
            if (alt < realDistance[edge.first]) {
                realDistance[edge.first] = alt;
				double estimate = distance(vertices[edge.first], wspair->second->box->center()) - wspair->second->box->radius();
                myHeap.decreaseKey(edge.first, alt + (estimate >= 0 ? estimate : 0));
            }
        }
    }
	if (foundCount < set->count()) {
		for (typename SplitTree<IsImplicit>::Iterator it = set->begin(); it != set->end(); ++it) {
			int index = *it;
			double dist = distance(vertices[index], vertices[from]);
			if (is_infinite(realDistance[index])
				&&
				dist < minDist)
			{
				minDist = dist;
				result = index;
			}
		}
	}
	return result;
}

template<class Heap, bool IsImplicit>
WellSeparatedPairRepresentative* ClosestPair(unsigned int pairIndex,
											 WellSeparatedPair<IsImplicit>* pair,
											 const pointset &vertices,
											 std::vector< std::pair<unsigned int, double> >* myEdges,
											 double t,
											 Heap& myHeap,
											 bool* mySet,
											 double* realDistance
											 ) {
	unsigned int N = vertices.size();

	for (unsigned int i = 0; i < N; i++)
		mySet[i] = false;
	for (typename SplitTree<IsImplicit>::Iterator it = pair->second->begin(); it != pair->second->end(); ++it) {
		unsigned int index = *it;
		mySet[index] = true;
	}
	int best = -1;
	int index = -1;
	double bestDist = std::numeric_limits<double>::infinity();
	unsigned int j = 0;
	for (typename SplitTree<IsImplicit>::Iterator it = pair->first->begin(); it != pair->first->end(); ++it, j++) {
		if (pair->interestingCache[j]) {
			unsigned int indexj = *it;
			if (pair->shortestPathSoFarIndex == -1
				||
				t * distance(vertices[indexj], vertices[pair->shortestPathSoFarIndex])
				+ pair->shortestPathSoFar
				+ t * (pair->second->box->size()
				       + distance(vertices[pair->shortestPathSoFarOtherIndex], pair->second->box->center()))
				+ 1e-6 >
				t * (distance(vertices[indexj], pair->second->box->center())
				       - pair->second->box->size())) {
				int result = closestDilationHigherT<Heap, IsImplicit>(vertices, myEdges, indexj, t, pair->second, myHeap, pair, realDistance, mySet);
				if (result != -1) {
					double dist = distance(vertices[indexj], vertices[result]);
					if (dist < bestDist) {
						bestDist = dist;
						best = result;
						index = indexj;
					}
				}
				else
					pair->interestingCache[j] = false;
			}
			else
				pair->interestingCache[j] = false;
		}
	}

	if (best == -1)
		return nullptr;
	else
		return new WellSeparatedPairRepresentative(index, best, bestDist, pairIndex);
}

template<class Heap, bool IsImplicit>
void FillQueue(Heap& Q,
			   std::vector<WellSeparatedPairRepresentative*>& QContent,
			   unsigned int& i,
			   std::vector<WellSeparatedPair<IsImplicit>*>& pairs,
			   const pointset &vertices,
			   std::vector< std::pair<unsigned int, double> >* myEdges,
			   double t,
			   double s,
			   Heap& myHeap,
			   int edgeCount,
			   bool* mySet,
			   double* realDistance
			   ) {
	while (i < pairs.size() && (Q.getCount() == 0 || pairs[i]->minlength() <= Q.getMin().second)) {
		WellSeparatedPairRepresentative* p = ClosestPair<Heap, IsImplicit>(i, pairs[i], vertices, myEdges, t, myHeap, mySet, realDistance);
		if (i % 1000 == 0)
			std::cout << i << "/" << pairs.size() << " " << Q.getCount() << " " << edgeCount << endl;
		if (p != nullptr) {
			Q.insert(i, p->dist);
			QContent[i] = p;
		}
		i++;
	}
}

template<class Heap, bool IsImplicit>
int GreedyLinspace(BaseTask * T) {
	STARTLOG
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
	unsigned int edgeCount = 0;
	unsigned int savedDijkstraCounter = 0;
	unsigned int discardedPairsCounter = 0;
	double t = T->t < 2.0 ? T->t : 2.0;
	double s = 4.0 * t / (t - 1.0);
	Heap myHeap(N, 0);
	bool* mySet = new bool[N];
	double* realDistance = new double[N];

	SplitTree<IsImplicit>* splitTree = SplitTree<IsImplicit>::createSplitTree(vertices);
	std::vector<WellSeparatedPair<IsImplicit>*> pairs;
	splitTree->ComputeWspdPairs(s, pairs);
	sort(pairs.begin(), pairs.end(), WellSeparatedPairPointerComparer<IsImplicit>());

    std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	std::vector<WellSeparatedPairRepresentative*> QContent(pairs.size(), nullptr);
	Heap Q(pairs.size(), -std::numeric_limits<double>::infinity());

	unsigned int index = 0;
	FillQueue<Heap>(Q, QContent, index, pairs, vertices, myEdges, t, s, myHeap, edgeCount, mySet, realDistance);

	while (Q.getCount() > 0) {
        std::pair<unsigned int, double> minQ = Q.getMin();

        Q.extractMin();
		WellSeparatedPairRepresentative* minRepresentative = QContent[minQ.first];
		QContent[minQ.first] = nullptr;

		std::pair<unsigned int, double> pair(minRepresentative->u, minRepresentative->dist);
        myEdges[minRepresentative->v].push_back(pair);
        pair.first = minRepresentative->v;
        myEdges[minRepresentative->u].push_back(pair);
        if (T->create_edge(edge(minRepresentative->u, minRepresentative->v))) return 0;
        edgeCount++;

		WellSeparatedPair<IsImplicit>* first = pairs[minRepresentative->pairIndex];
		double length = first->length();

		for (int j = 0; j < pairs.size(); j++) {
			if (QContent[j] != nullptr) {
				WellSeparatedPair<IsImplicit>* second = pairs[j];
				double distffsf = minDistance(first->first->box, second->first->box);
				double distffss = minDistance(first->first->box, second->second->box);
				double distfssf = minDistance(first->second->box, second->first->box);
				double distfsss = minDistance(first->second->box, second->second->box);
				double lengthf = length;
				double lengths = maxDistance(second->first->box, second->second->box) * t + 1e-6;
				if (distffsf + lengthf + distfsss <= lengths
					||
					distffss + lengthf + distfssf <= lengths) {
					WellSeparatedPairRepresentative* p = ClosestPair<Heap>(j, second, vertices, myEdges, t, myHeap, mySet, realDistance);
					if (p == nullptr) {
						Q.remove(j);
						QContent[j] = nullptr;
					}
					else {
						Q.increaseKey(j, p->dist);
						QContent[j] = p;
					}
				}
			}
		}

		FillQueue<Heap>(Q, QContent, index, pairs, vertices, myEdges, t, s, myHeap, edgeCount, mySet, realDistance);
	}
	ENDLOG

	delete[] mySet;
	delete[] realDistance;
	delete splitTree;
    delete[] myEdges;
    return edgeCount;
}

template<class Heap, bool IsImplicit>
void UpdatePair(unsigned int j,
				Heap& Q,
				std::vector<WellSeparatedPairRepresentative*>& QContent,
				WellSeparatedPair<IsImplicit>* pair,
				const pointset &vertices,
				std::vector< std::pair<unsigned int, double> >* myEdges,
				double t,
				Heap& myHeap,
				bool* mySet,
				bool* dirtyBits,
				double* realDistance
				) {
	WellSeparatedPairRepresentative* p = ClosestPair<Heap, IsImplicit>(j, pair, vertices, myEdges, t, myHeap, mySet, realDistance);
	if (p == nullptr) {
		Q.remove(j);
		QContent[j] = nullptr;
	}
	else {
		Q.increaseKey(j, p->dist);
		QContent[j] = p;
	}
	dirtyBits[j] = false;
}

template<class Heap, bool IsImplicit>
void ReinsertPair(unsigned int j,
				  Heap& Q,
				  std::vector<WellSeparatedPairRepresentative*>& QContent,
				  WellSeparatedPair<IsImplicit>* pair,
				  const pointset &vertices,
				  std::vector< std::pair<unsigned int, double> >* myEdges,
				  double t,
				  Heap& myHeap,
				  bool* mySet,
				  bool* dirtyBits,
				  double* realDistance
				  ) {
	WellSeparatedPairRepresentative* p = ClosestPair<Heap, IsImplicit>(j, pair, vertices, myEdges, t, myHeap, mySet, realDistance);
	if (p != nullptr) {
		Q.insert(j, p->dist);
		QContent[j] = p;
	}
	dirtyBits[j] = false;
}

template<class Heap, bool IsImplicit>
void FillQueue2(Heap& Q,
				std::vector<WellSeparatedPairRepresentative*>& QContent,
				std::vector<WellSeparatedPair<IsImplicit>*>& pairs,
				const pointset &vertices
				) {
	for (unsigned int i = 0; i < pairs.size(); i++) {
		double bestDist = std::numeric_limits<double>::infinity();
		unsigned int bestFirst = 0;
		unsigned int bestSecond = 0;
		for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->first->begin(); it != pairs[i]->first->end(); ++it) {
			unsigned int firstI = *it;
			for (typename SplitTree<IsImplicit>::Iterator jt = pairs[i]->second->begin(); jt != pairs[i]->second->end(); ++jt) {
				unsigned int secondI = *jt;
				double nextDist = distance(vertices[firstI], vertices[secondI]);
				if (bestDist > nextDist) {
					bestDist = nextDist;
					bestFirst = firstI;
					bestSecond = secondI;
				}
			}
		}
		Q.insert(i, bestDist);
		QContent[i] = new WellSeparatedPairRepresentative(bestFirst, bestSecond, bestDist, i);
	}
}

template<class Heap, bool IsImplicit>
int GreedyLinspace2(BaseTask * T) {
	STARTLOG
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
	unsigned int edgeCount = 0;
	unsigned int savedDijkstraCounter = 0;
	unsigned int discardedPairsCounter = 0;
	double t = T->t < 2.0 ? T->t : 2.0;
	double s = 5.0;//T->var;
	Heap myHeap(N, 0);
	bool* mySet = new bool[N];
	double* realDistance = new double[N];

	SplitTree<IsImplicit>* splitTree = SplitTree<IsImplicit>::createSplitTree(vertices);
	std::vector<WellSeparatedPair<IsImplicit>*> pairs;
	splitTree->ComputeWspdPairs(s, pairs);
	unsigned int m = pairs.size();

    std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	std::vector<WellSeparatedPairRepresentative*> QContent(m, nullptr);
	Heap Q(m, -std::numeric_limits<double>::infinity());
	bool* dirtyBits = new bool[m];
	for (unsigned int i = 0; i < m; i++)
		dirtyBits[i] = false;

	FillQueue2<Heap>(Q, QContent, pairs, vertices);

	while (Q.getCount() > 0) {
		bool foundOne = false;
		std::pair<unsigned int, double> minQ;
		while (!foundOne) {
			minQ = Q.getMin();

			if (dirtyBits[minQ.first])
				UpdatePair(minQ.first, Q, QContent, pairs[minQ.first], vertices, myEdges, t, myHeap, mySet, dirtyBits, realDistance);
			else
				foundOne = true;
			if (Q.getCount() == 0)
				goto end;
		}
        Q.extractMin();
		WellSeparatedPairRepresentative* minRepresentative = QContent[minQ.first];

		std::pair<unsigned int, double> pair(minRepresentative->u, minRepresentative->dist);
        myEdges[minRepresentative->v].push_back(pair);
        pair.first = minRepresentative->v;
        myEdges[minRepresentative->u].push_back(pair);
        if (T->create_edge(edge(minRepresentative->u, minRepresentative->v))) return 0;
        edgeCount++;

		for (int j = 0; j < m; j++)
			if (QContent[j] != nullptr)
				dirtyBits[j] = true;

		ReinsertPair(minQ.first, Q, QContent, pairs[minQ.first], vertices, myEdges, t, myHeap, mySet, dirtyBits, realDistance);
	}
end:
	ENDLOG

	delete[] mySet;
	delete[] realDistance;
	delete splitTree;
    delete[] myEdges;
	delete[] dirtyBits;
    return edgeCount;
}

template<class Heap>
void DoDijkstra(unsigned int j,
				const pointset& vertices,
				std::vector< std::pair<unsigned int, double> >* myEdges,
				Heap &myHeap,
				double t,
				std::vector<int>& QContent,
				Heap& Q,
				bool* dirtyBits
				) {
	qcount++;
    unsigned int N = vertices.size();
	QContent[j] = -1;
	dirtyBits[j] = false;
	myHeap.clear();
    for (unsigned int i = 0; i < N; i++)
        myHeap.insert(i, std::numeric_limits<double>::infinity());
    myHeap.decreaseKey(j, 0);
    while (myHeap.getCount() > 0) {
        std::pair<unsigned int, double> pair = myHeap.getMin();
		double directDist = distance(vertices[pair.first], vertices[j]);
        if (pair.second > t * directDist)
        {
			if (vertices[pair.first].x >= vertices[j].x) {
				if (QContent[j] == -1 || directDist < distance(vertices[j], vertices[QContent[j]]))
					QContent[j] = pair.first;
			}
        }
        myHeap.extractMin();
        for (unsigned int i = 0; i < myEdges[pair.first].size(); i++) {
            std::pair<unsigned int, double> edge = myEdges[pair.first][i];
            double alt = pair.second + edge.second;
            if (alt < myHeap.getValue(edge.first))
                myHeap.decreaseKey(edge.first, alt);
        }
    }
	if (QContent[j] != -1) {
		if (!Q.contains(j))
			Q.insert(j, distance(vertices[QContent[j]], vertices[j]));
		else {
			if (Q.getValue(j) < distance(vertices[QContent[j]], vertices[j]))
				Q.increaseKey(j, distance(vertices[QContent[j]], vertices[j]));
			else
				Q.decreaseKey(j, distance(vertices[QContent[j]], vertices[j]));
		}
	}
	else if (Q.contains(j))
		Q.remove(j);
}

template<class Heap>
int GreedyLinspace3(BaseTask * T) {
	STARTLOG
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
	unsigned int edgeCount = 0;
	unsigned int savedDijkstraCounter = 0;
	unsigned int discardedPairsCounter = 0;
	double t = T->t < 2.0 ? T->t : 2.0;
	Heap myHeap(N, 0);

    std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	std::vector<int> QContent(N, -1);
	Heap Q(N, -std::numeric_limits<double>::infinity());
	bool* dirtyBits = new bool[N];
	for (unsigned int j = 0; j < N; j++)
		DoDijkstra<Heap>(j, vertices, myEdges, myHeap, t, QContent, Q, dirtyBits);

	while (Q.getCount() > 0) {
		bool foundOne = false;
		std::pair<unsigned int, double> minQ;
		while (!foundOne)
		{
			minQ = Q.getMin();
			if (dirtyBits[minQ.first])
				DoDijkstra<Heap>(minQ.first, vertices, myEdges, myHeap, t, QContent, Q, dirtyBits);
			else
				foundOne = true;
			if (Q.getCount() == 0)
				goto end;
		}
		int minRepFrom = minQ.first;
		int minRepTo = QContent[minQ.first];
		double minRepDist = minQ.second;//distance(vertices[minRepFrom], vertices[minRepTo]);

		std::pair<unsigned int, double> pair(minRepFrom, minRepDist);
        myEdges[minRepTo].push_back(pair);
        pair.first = minRepTo;
        myEdges[minRepFrom].push_back(pair);
        if (T->create_edge(edge(minRepFrom, minRepTo))) return 0;
        edgeCount++;

		//DoDijkstra<Heap>(minQ.first, vertices, myEdges, myHeap, t, QContent, Q);

		for (int j = 0; j < N; j++)
			if (QContent[j] != -1)
				dirtyBits[j] = true;
	}
end:
	ENDLOG

	delete[] dirtyBits;
    delete[] myEdges;
    return edgeCount;
}

template<class Heap>
void FindClosest(unsigned int j,
				 const pointset& vertices,
				 std::vector<int>& QContent,
				 Heap& Q,
				 int bucketeers,
				 std::vector<int>* buckets,
				 int fromBucketIndexX,
				 int fromBucketIndexY
				 ) {
	int closestPoint = -1;
	double closestPointDist = std::numeric_limits<double>::infinity();
	for (int a = -1; a <= 1; a++) {
		for (int b = -1; b <= 1; b++) {
			if (fromBucketIndexX + a >= 0 && fromBucketIndexX + a < bucketeers
				&&
				fromBucketIndexY + b >= 0 && fromBucketIndexY + b < bucketeers) {
				int bucketIndex = fromBucketIndexX + a + bucketeers * (b + fromBucketIndexY);
				int bucketSize = buckets[bucketIndex].size();
				for (int i = 0; i < bucketSize; i++) {
					int newIndex = buckets[bucketIndex][i];
					if (newIndex != j) {
						double newDist = distance(vertices[j], vertices[newIndex]);
						if (newDist < closestPointDist) {
							closestPointDist = newDist;
							closestPoint = newIndex;
						}
					}
				}
			}
		}
	}
	if (closestPoint != -1) {
		QContent[j] = closestPoint;
		Q.insert(j, closestPointDist);
	}
}


template<class Heap>
void DoBoundedDijkstra(unsigned int j,
					   const pointset& vertices,
					   std::vector< std::pair<unsigned int, double> >* const myEdges,
					   Heap &myHeap,
					   double t,
					   std::vector<int>& QContent,
					   Heap& Q,
					   double maxtdist,
					   double maxedgelength,
					   int bucketeers,
					   std::vector<int>* buckets,
					   int fromBucketIndexX,
					   int fromBucketIndexY,
					   bool* dirtyBits
					   ) {
	QContent[j] = -1;
	dirtyBits[j] = false;
	myHeap.clear(false);
	for (int a = -1; a <= 1; a++) {
		for (int b = -1; b <= 1; b++) {
			if (fromBucketIndexX + a >= 0 && fromBucketIndexX + a < bucketeers
				&&
				fromBucketIndexY + b >= 0 && fromBucketIndexY + b < bucketeers) {
				int bucketIndex = fromBucketIndexX + a + bucketeers * (b + fromBucketIndexY);
				int bucketSize = buckets[bucketIndex].size();
				for (int i = 0; i < bucketSize; i++) {
					if (distance(vertices[j], vertices[buckets[bucketIndex][i]]) <= maxtdist) {
						myHeap.insert(buckets[bucketIndex][i], std::numeric_limits<double>::infinity());
					}
				}
			}
		}
	}
    myHeap.decreaseKey(j, 0);
    while (myHeap.getCount() > 0) {
        std::pair<unsigned int, double> pair = myHeap.getMin();
		double directDist = distance(vertices[pair.first], vertices[j]);
        if (pair.second > t * directDist
			&&
			directDist <= maxedgelength) {
			if (vertices[pair.first].x >= vertices[j].x) {
				if (QContent[j] == -1 || directDist < distance(vertices[j], vertices[QContent[j]]))
					QContent[j] = pair.first;
			}
        }
        myHeap.extractMin();
        for (unsigned int i = 0; i < myEdges[pair.first].size(); i++) {
            std::pair<unsigned int, double> edge = myEdges[pair.first][i];
			if (distance(vertices[j], vertices[edge.first]) <= maxtdist) {
				double alt = pair.second + edge.second;
				if (alt < myHeap.getValue(edge.first))
					myHeap.decreaseKey(edge.first, alt);
			}
        }
    }
	if (QContent[j] != -1) {
		if (!Q.contains(j))
			Q.insert(j, distance(vertices[QContent[j]], vertices[j]));
		else {
			if (Q.getValue(j) < distance(vertices[QContent[j]], vertices[j]))
				Q.increaseKey(j, distance(vertices[QContent[j]], vertices[j]));
			else
				Q.decreaseKey(j, distance(vertices[QContent[j]], vertices[j]));
		}
	}
	else if (Q.contains(j))
		Q.remove(j);
}

template<class Heap>
int ComputeGreedySpannerBuckets2(BaseTask * T,
								 unsigned int& edgeCount,
								 std::vector<std::pair<unsigned int, double> >* myEdges,
								 double t,
								 std::vector<int>* buckets,
								 double estimatedLongestEdgeLength,
								 double estimatedLongestTPath,
								 int* bucketLookupX,
								 int* bucketLookupY,
								 int* bucketLookupZ,
								 int bucketeers
								 ) {
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
	Heap myHeap(N, 0);
	std::vector<int> QContent(N, -1);
	Heap Q(N, -std::numeric_limits<double>::infinity());
	bool* dirtyBits = new bool[N];
	for (unsigned int j = 0; j < N; j++) {
		FindClosest<Heap>(j, vertices, QContent, Q, bucketeers, buckets, bucketLookupX[j], bucketLookupY[j]);
		dirtyBits[j] = false;
	}

	while (Q.getCount() > 0) {
		bool foundOne = false;
		std::pair<unsigned int, double> minQ;
		while (!foundOne) {
			minQ = Q.getMin();
			if (dirtyBits[minQ.first])
				DoBoundedDijkstra<Heap>(minQ.first, vertices, myEdges, myHeap, t, QContent, Q, estimatedLongestTPath, estimatedLongestEdgeLength, bucketeers, buckets, bucketLookupX[minQ.first], bucketLookupY[minQ.first], dirtyBits);
			else
				foundOne = true;
			if (Q.getCount() == 0)
				goto end;
		}
		int minRepFrom = minQ.first;
		int minRepTo = QContent[minQ.first];
		double minRepDist = minQ.second;
		Q.extractMin();

		std::pair<unsigned int, double> pair(minRepFrom, minRepDist);
        myEdges[minRepTo].push_back(pair);
        pair.first = minRepTo;
        myEdges[minRepFrom].push_back(pair);
        if (T->create_edge(edge(minRepFrom, minRepTo))) return 0;
        edgeCount++;

		int fromBucketIndexX = bucketLookupX[minRepFrom];
		int fromBucketIndexY = bucketLookupY[minRepFrom];
		for (int a = -1; a <= 1; a++) {
			for (int b = -1; b <= 1; b++) {
				if (fromBucketIndexX + a >= 0 && fromBucketIndexX + a < bucketeers
					&&
					fromBucketIndexY + b >= 0 && fromBucketIndexY + b < bucketeers) {
					int bucketIndex = fromBucketIndexX + a + bucketeers * (b + fromBucketIndexY);
					int bucketSize = buckets[bucketIndex].size();
					for (int i = 0; i < bucketSize; i++) {
						dirtyBits[buckets[bucketIndex][i]] = true;
					}
				}
			}
		}
		DoBoundedDijkstra<Heap>(minRepFrom, vertices, myEdges, myHeap, t, QContent, Q, estimatedLongestTPath, estimatedLongestEdgeLength, bucketeers, buckets, bucketLookupX[minRepFrom], bucketLookupY[minRepFrom], dirtyBits);
	}
end:
	delete[] dirtyBits;

    return edgeCount;
}

template<class Heap>
struct GreedySpannerBucketsWrapper {
	static void ComputeShortEdges(BaseTask * T,
								  unsigned int& edgeCount,
								  std::vector<std::pair<unsigned int, double> >* edges,
								  double t,
								  std::vector<int>* buckets,
								  double estimatedLongestEdgeLength,
								  double estimatedLongestTPath,
								  int* bucketLookupX,
								  int* bucketLookupY,
								  int* bucketLookupZ,
								  int bucketeers
								  ) {
		ComputeGreedySpannerBuckets<Heap>(T, edgeCount, edges, t, buckets, estimatedLongestEdgeLength, estimatedLongestTPath, bucketLookupX, bucketLookupY, bucketLookupZ, bucketeers);
	}
};

template<class Heap>
struct GreedySpannerBuckets2Wrapper {
	static void ComputeShortEdges(BaseTask * T,
								  unsigned int& edgeCount,
								  std::vector<std::pair<unsigned int, double> >* edges,
								  double t,
								  std::vector<int>* buckets,
								  double estimatedLongestEdgeLength,
								  double estimatedLongestTPath,
								  int* bucketLookupX,
								  int* bucketLookupY,
								  int* bucketLookupZ,
								  int bucketeers
								  ) {
		ComputeGreedySpannerBuckets2<Heap>(T, edgeCount, edges, t, buckets, estimatedLongestEdgeLength, estimatedLongestTPath, bucketLookupX, bucketLookupY, bucketLookupZ, bucketeers);
	}
};

template<class Algorithm>
int GreedySpannerBuckets(BaseTask * T) {
	STARTLOG
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
	double t = T->t;
    std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
    unsigned int edgeCount = 0;
	
	double estimatedLongestEdgeLengthUnit = log((double)N) / sqrt((double)N) / log(log((double)N)) / sqrt(sqrt(t - 1.0)) * T->var;
	double minX = vertices[0].x, minY = vertices[0].y, maxX = vertices[0].x, maxY = vertices[0].y;
	for (int i = 1; i < N; i++) {
		if (vertices[i].x < minX) minX = vertices[i].x;
		if (vertices[i].x > maxX) maxX = vertices[i].x;
		if (vertices[i].y < minY) minY = vertices[i].y;
		if (vertices[i].y > maxY) maxY = vertices[i].y;
	}
	double xwidth = maxX - minX;
	double ywidth = maxY - minY;
	double maxwidth = xwidth > ywidth ? xwidth : ywidth;
	double estimatedLongestEdgeLength = estimatedLongestEdgeLengthUnit * maxwidth;
	double estimatedLongestTPath = estimatedLongestEdgeLength * (1 + T->t / 2.0);
    //cout << "estimatedLongestEdgeLength: " << std::fixed << estimatedLongestEdgeLength << endl;
    //cout << "estimatedLongestTPath: " << std::fixed << estimatedLongestTPath << endl;
	int bucketeers = (int)(1.0 / (estimatedLongestEdgeLengthUnit * (1 + T->t / 2.0)));
	if(bucketeers<1){bucketeers=1;}
    //cout << "bucketSide: " << bucketeers << endl;
	double bucketWidth = xwidth / ((double)bucketeers) + 1e-6;
	double bucketHeight = ywidth / ((double)bucketeers) + 1e-6;
	if (bucketWidth < 1.0) bucketWidth = 1e-6;
	if (bucketHeight < 1.0) bucketHeight = 1e-6;
    //cout << "bucketWidth: " << std::fixed << bucketWidth << endl;
    //cout << "bucketHeight: " << std::fixed << bucketHeight << endl;

	std::vector<int>* buckets = new std::vector<int>[bucketeers * bucketeers];
	int* bucketLookupX = new int[N];
	int* bucketLookupY = new int[N];
	int* bucketLookupZ = new int[N];
	for (int i = 0; i < N; i++) {
		int bucketX = (int)((vertices[i].x - minX) / bucketWidth);
		int bucketY = (int)((vertices[i].y - minY) / bucketHeight);
		int bucketCoord = bucketX + bucketeers * bucketY;
		bucketLookupZ[i] = buckets[bucketCoord].size();
		buckets[bucketCoord].push_back(i);
		bucketLookupX[i] = bucketX;
		bucketLookupY[i] = bucketY;
	}

	Algorithm::ComputeShortEdges(T, edgeCount, myEdges, t, buckets, estimatedLongestEdgeLength, estimatedLongestTPath, bucketLookupX, bucketLookupY, bucketLookupZ, bucketeers);
	ENDLOG
	delete[] myEdges;
	delete[] buckets;
	delete[] bucketLookupX;
	delete[] bucketLookupY;
	delete[] bucketLookupZ;

    return edgeCount;
}

template<class Heap, bool IsImplicit>
WellSeparatedPairRepresentative* ClosestPair2(unsigned int pairNo,
											  WellSeparatedPair<IsImplicit>* pair,
											  const pointset &vertices,
											  std::vector< std::pair<unsigned int, double> >* myEdges,
											  double t,
											  Heap& myHeap,
											  bool* mySet,
											  double* realDistance,
											  bool** discounters
											  ) {
	unsigned int N = vertices.size();

	for (unsigned int i = 0; i < N; i++)
		mySet[i] = false;
	for (typename SplitTree<IsImplicit>::Iterator it = pair->second->begin(); it != pair->second->end(); ++it) {
		unsigned int index = *it;
		mySet[index] = true;
	}
	int best = -1;
	int index = -1;
	double bestDist = std::numeric_limits<double>::infinity();
	unsigned int j = 0;
	for (typename SplitTree<IsImplicit>::Iterator it = pair->first->begin(); it != pair->first->end(); ++it, j++) {
		if (pair->interestingCache[j]
			&&
			!discounters[pairNo][j]) {
			unsigned int indexj = *it;
			if (pair->shortestPathSoFarIndex == -1
				||
				t * distance(vertices[indexj], vertices[pair->shortestPathSoFarIndex])
				+ pair->shortestPathSoFar
				+ t * (pair->second->box->size()
				       + distance(vertices[pair->shortestPathSoFarOtherIndex], pair->second->box->center()))
				+ 1e-6 >
				t * (distance(vertices[indexj], pair->second->box->center())
				       - pair->second->box->size())) {
				int result = closestDilationHigherT<Heap, IsImplicit>(vertices, myEdges, indexj, t, pair->second, myHeap, pair, realDistance, mySet);
				if (result != -1) {
					double dist = distance(vertices[indexj], vertices[result]);
					if (dist < bestDist) {
						bestDist = dist;
						best = result;
						index = indexj;
					}
				}
				else
					pair->interestingCache[j] = false;
			}
			else
				pair->interestingCache[j] = false;
		}
	}

	if (best == -1)
		return nullptr;
	else
		return new WellSeparatedPairRepresentative(index, best, bestDist, pairNo);
}

template<class Heap, bool IsImplicit>
void UpdatePair2(unsigned int j,
				 Heap& Q,
				 std::vector<WellSeparatedPairRepresentative*>& QContent,
				 WellSeparatedPair<IsImplicit>* pair,
				 const pointset &vertices,
				 std::vector< std::pair<unsigned int, double> >* myEdges,
				 double t,
				 Heap& myHeap,
				 bool* mySet,
				 bool* dirtyBits,
				 double* realDistance,
				 bool** discounters
				 ) {
	WellSeparatedPairRepresentative* p = ClosestPair2<Heap, IsImplicit>(j, pair, vertices, myEdges, t, myHeap, mySet, realDistance, discounters);
	if (p == nullptr) {
		Q.remove(j);
		QContent[j] = nullptr;
	}
	else {
		Q.increaseKey(j, p->dist);
		QContent[j] = p;
	}
	dirtyBits[j] = false;
}

template<class Heap, bool IsImplicit>
void ReinsertPair2(unsigned int j,
				   Heap& Q,
				   std::vector<WellSeparatedPairRepresentative*>& QContent,
				   WellSeparatedPair<IsImplicit>* pair,
				   const pointset &vertices,
				   std::vector< std::pair<unsigned int, double> >* myEdges,
				   double t,
				   Heap& myHeap,
				   bool* mySet,
				   bool* dirtyBits,
				   double* realDistance,
				   bool** discounters
				   ) {
	WellSeparatedPairRepresentative* p = ClosestPair2<Heap, IsImplicit>(j, pair, vertices, myEdges, t, myHeap, mySet, realDistance, discounters);
	if (p != nullptr) {
		Q.insert(j, p->dist);
		QContent[j] = p;
	}
	dirtyBits[j] = false;
}

template<class Heap, bool IsImplicit>
void FillQueue3(Heap& Q,
				std::vector<WellSeparatedPairRepresentative*>& QContent,
				std::vector<WellSeparatedPair<IsImplicit>*>& pairs,
				const pointset &vertices,
				std::vector< std::pair<unsigned int, double> >* myEdges,
				double t,
				Heap& myHeap,
				int edgeCount,
				bool* discountedPairs,
				bool* mySet,
				bool* dirtyBits,
				double* realDistance,
				bool** discounters
				) {
	for (unsigned int j = 0; j < pairs.size(); j++) {
		if (!discountedPairs[j]) {
			WellSeparatedPairRepresentative* p = ClosestPair2<Heap, IsImplicit>(j, pairs[j], vertices, myEdges, t, myHeap, mySet, realDistance, discounters);
			if (p == nullptr)
				QContent[j] = nullptr;
			else {
				Q.insert(j, p->dist);
				QContent[j] = p;
			}
		}
	}
}

double atan3(vertex point, vertex origin) {
	double x = point.x - origin.x;
	double y = point.y - origin.y;
	double z = atan2(y, x); // WTF?
	if (z < 0) z = z + 2 * M_PI;
	return z;
}

class AngleStorer {
private:
	ScapegoatHeap<double, double> heap;
public :
	vertex origin;
	void Add(double start, double end) {
		if (heap.getCount() == 0) {
			heap.insert(start, end);
			return;
		}
		ScapegoatNode<double, double>* closest = heap.searchClosest(start);
		if (closest->key > start) {
			if (closest->prev == nullptr) {
				while (closest->value < end && closest->next != nullptr) {
					closest = closest->next;
					heap.remove(closest->prev);
				}
				if (closest->value < end && closest->next == nullptr) {
					closest->key = start;
					closest->value = end;
					return;
				}
				if (closest->key > end) {
					heap.insert(start, end);
					return;
				}
				closest->key = start;
				return;
			}
			closest = closest->prev;
		}

		ScapegoatNode<double, double>* rightClosest = closest->next;
		while (rightClosest != nullptr && rightClosest->value < end) {
			ScapegoatNode<double, double>* previous = rightClosest;
			rightClosest = rightClosest->next;
			heap.remove(previous);
		}
		if (rightClosest == nullptr
			||
			rightClosest->key > end) {
			if (start <= closest->value) {
				if (closest->value < end)
					closest->value = end;
			}
			else
				heap.insert(start, end);
		}
		else {
			if (start <= closest->value) {
				closest->value = rightClosest->value;
				heap.remove(rightClosest);
			}
			else
				rightClosest->key = start;
		}
	}
	bool Covers(ScapegoatNode<double, double>* node, double start, double end) {
		return node->key <= start && node->value >= end;
	}
	bool Covers(double start, double end) {
		if (heap.getCount() == 0)
			return false;
		ScapegoatNode<double, double>* closest = heap.searchClosest(start);
		if (closest->value < start)
			return closest->next != nullptr && Covers(closest->next, start, end);
		else if (end < closest->key)
			return closest->prev != nullptr && Covers(closest->prev, start, end);
		else
			return Covers(closest, start, end);
	}
	bool Covers(BoundingBox* box) {
		double angles[4];
		angles[0] = atan3(vertex(box->left, box->bottom), origin);
		angles[1] = atan3(vertex(box->right, box->bottom), origin);
		angles[2] = atan3(vertex(box->left, box->top), origin);
		angles[3] = atan3(vertex(box->right, box->top), origin);
		std::sort(angles, angles + 4);

		if (angles[0] + M_PI > angles[3]) {
			return Covers(angles[0], angles[3]);
		}
		else {
			if (angles[0] + M_PI < angles[1]) {
				return Covers(0.0, angles[0]) && Covers(angles[1], 2 * M_PI);
			}
			else if (angles[1] + M_PI < angles[2]) {
				return Covers(0.0, angles[1]) && Covers(angles[2], 2 * M_PI);
			}
			return Covers(0.0, angles[2]) && Covers(angles[3], 2 * M_PI);
		}
	}
};

template<class Algorithm, class Heap, bool IsImplicit>
int GreedySpannerBucketsFixed(BaseTask * T) {
	//TODO:
	// Make ESA' local
	STARTLOG
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
    std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
    unsigned int edgeCount = 0;

	// Bucketing stuff
	double t = T->t;
	double estimatedLongestEdgeLengthUnit = log((double)N) / sqrt((double)N) / log(log((double)N)) / sqrt(sqrt(t - 1.0)) * T->var;
	double minX = vertices[0].x, minY = vertices[0].y, maxX = vertices[0].x, maxY = vertices[0].y;
	for (int i = 1; i < N; i++) {
		if (vertices[i].x < minX) minX = vertices[i].x;
		if (vertices[i].x > maxX) maxX = vertices[i].x;
		if (vertices[i].y < minY) minY = vertices[i].y;
		if (vertices[i].y > maxY) maxY = vertices[i].y;
	}
	double xwidth = maxX - minX;
	double ywidth = maxY - minY;
	double maxwidth = xwidth > ywidth ? xwidth : ywidth;
	double estimatedLongestEdgeLength = estimatedLongestEdgeLengthUnit * maxwidth;
	double estimatedLongestTPath = estimatedLongestEdgeLength * (1 + t / 2.0);
    //cout << "estimatedLongestEdgeLength: " << std::fixed << estimatedLongestEdgeLength << endl;
    //cout << "estimatedLongestTPath: " << std::fixed << estimatedLongestTPath << endl;
	int bucketeers = (int)(1.0 / (estimatedLongestEdgeLengthUnit * (1 + t / 2.0)));
	if(bucketeers<1){bucketeers=1;}
    //cout << "bucketSide: " << bucketeers << endl;
	double bucketWidth = xwidth / ((double)bucketeers) + 1e-6;
	double bucketHeight = ywidth / ((double)bucketeers) + 1e-6;
	if (bucketWidth < 1.0) bucketWidth = 1e-6;
	if (bucketHeight < 1.0) bucketHeight = 1e-6;
    //cout << "bucketWidth: " << std::fixed << bucketWidth << endl;
   // cout << "bucketHeight: " << std::fixed << bucketHeight << endl;

	std::vector<int>* buckets = new std::vector<int>[bucketeers * bucketeers];
	int* bucketLookupX = new int[N];
	int* bucketLookupY = new int[N];
	int* bucketLookupZ = new int[N];
	for (int i = 0; i < N; i++) {
		int bucketX = (int)((vertices[i].x - minX) / bucketWidth);
		int bucketY = (int)((vertices[i].y - minY) / bucketHeight);
		int bucketCoord = bucketX + bucketeers * bucketY;
		bucketLookupZ[i] = buckets[bucketCoord].size();
		buckets[bucketCoord].push_back(i);
		bucketLookupX[i] = bucketX;
		bucketLookupY[i] = bucketY;
	}

	// Call base algorithm
	Algorithm::ComputeShortEdges(T, edgeCount, myEdges, t, buckets, estimatedLongestEdgeLength, estimatedLongestTPath, bucketLookupX, bucketLookupY, bucketLookupZ, bucketeers);
	//cout << edgeCount << endl;
	//cout << Timing::elapsed() << endl;

	// Find bridging points
	//هدف از این بخش این است تا برای هر نقطه، به ازای هر مسیر پوششی که از آن به یکی از نقاط داخل سلولهای 
	//مجاور (و خود سلول مربوط به نقطه) 
	//وجود دارد یک مخروط ایجاد شود که اندازه آن بستگی به طول مسیر پوششی بین آن دو دارد
	AngleStorer* certificates = new AngleStorer[N];
	Heap myHeap(N, 0);
	for (int i = 0; i < N; i++) {
		certificates[i].origin = vertices[i];
		myHeap.clear(false);
		int fromBucketIndexX = bucketLookupX[i];
		int fromBucketIndexY = bucketLookupY[i];
		for (int a = -1; a <= 1; a++) {
			for (int b = -1; b <= 1; b++) {
				if (fromBucketIndexX + a >= 0 && fromBucketIndexX + a < bucketeers
					&&
					fromBucketIndexY + b >= 0 && fromBucketIndexY + b < bucketeers) {
					int bucketIndex = fromBucketIndexX + a + bucketeers * (b + fromBucketIndexY);
					int bucketSize = buckets[bucketIndex].size();
					for (int j = 0; j < bucketSize; j++) {
						if (distance(vertices[i], vertices[buckets[bucketIndex][j]]) <= estimatedLongestTPath) {
							myHeap.insert(buckets[bucketIndex][j], std::numeric_limits<double>::infinity());
						}
					}
				}
			}
		}
		myHeap.decreaseKey(i, 0);
		while (myHeap.getCount() > 0) {
			std::pair<unsigned int, double> pair = myHeap.getMin();
			double directDist = distance(vertices[pair.first], vertices[i]);
			if (pair.first != i
				&&
				directDist <= estimatedLongestEdgeLength
				&&
				t * directDist > pair.second) {
				double inner = (t * directDist - pair.second) / (2.0 * t * directDist);
				double angle = 2.0 * asin(inner);
				vertex target = vertices[pair.first];
				vertex origin = vertices[i];
				double startAngle = atan3(target, origin);
				double start = startAngle - angle;
				double end = startAngle + angle;
				if (start < 0) {
					certificates[i].Add(start + 2*M_PI, 2*M_PI);
					certificates[i].Add(0.0, end);
				}
				else if (end >= 2*M_PI) {
					certificates[i].Add(start, 2*M_PI);
					certificates[i].Add(0.0, end - 2*M_PI);
				}
				else
					certificates[i].Add(start, end);
			}

			myHeap.extractMin();
			for (unsigned int j = 0; j < myEdges[pair.first].size(); j++) {
				std::pair<unsigned int, double> edge = myEdges[pair.first][j];
				if (myHeap.contains(edge.first)) {
					double alt = pair.second + edge.second;
					if (alt < myHeap.getValue(edge.first))
						myHeap.decreaseKey(edge.first, alt);
				}
			}
		}
	}
	
	delete[] buckets;
	delete[] bucketLookupX;
	delete[] bucketLookupY;

	// Compute WSPD
	double s = 5.0;
	SplitTree<IsImplicit>* splitTree = SplitTree<IsImplicit>::createSplitTree(vertices);
	std::vector<WellSeparatedPair<IsImplicit>*> pairs;
	splitTree->ComputeWspdPairs(s, pairs, estimatedLongestEdgeLength);
	//cout << "Long pairs: " << pairs.size() << endl;

	// Discount pairs
	bool* discountedPairs = new bool[pairs.size()];
	bool** discounters = new bool*[pairs.size()];
	unsigned int discountedPairsCounter = 0;
	unsigned int discountersCounter = 0;
	unsigned int discountersTotalCounter = 0;
	for (unsigned int i = 0; i < pairs.size(); i++) {
		discountedPairs[i] = false;
		bool foundBadPoint = false;
		for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->first->begin(); it != pairs[i]->first->end(); ++it) {
			unsigned int a = *it;
			if (!certificates[a].Covers(pairs[i]->second->box)) {
				foundBadPoint = true;
				break;
			}
		}
		if (!foundBadPoint) {
			discountedPairs[i] = true;
			discountedPairsCounter++;
		}
		else {
			foundBadPoint = false;
			for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->second->begin(); it != pairs[i]->second->end(); ++it) {
				unsigned int b = *it;
				if (!certificates[b].Covers(pairs[i]->first->box)) {
					foundBadPoint = true;
					break;
				}
			}
			if (!foundBadPoint) {
				discountedPairs[i] = true;
				discountedPairsCounter++;
			}
			else {
				discounters[i] = new bool[pairs[i]->first->count()];
				discountersTotalCounter += pairs[i]->first->count();
				unsigned int index = 0;
				for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->first->begin(); it != pairs[i]->first->end(); ++it, index++) {
					unsigned int a = *it;
					bool discount = certificates[a].Covers(pairs[i]->second->box);
					discounters[i][index] = discount;
					if (discount) discountersCounter++;
				}
			}
		}
	}
	//cout << "Discounted " << discountedPairsCounter << "/" << pairs.size() << endl;
	//cout << "Discounted " << discountersCounter << "/" << discountersTotalCounter << endl;

	delete[] certificates;

	// Do ESA', ignore discounted pairs, don't do dijkstras on discounted points
	bool* mySet = new bool[N];
	double* realDistance = new double[N];
	
	unsigned int m = pairs.size();
	std::vector<WellSeparatedPairRepresentative*> QContent(m, nullptr);
	Heap Q(m, -std::numeric_limits<double>::infinity());
	bool* dirtyBits = new bool[m];
	for (unsigned int i = 0; i < m; i++)
		dirtyBits[i] = false;

	FillQueue3<Heap>(Q, QContent, pairs, vertices, myEdges, t, myHeap, edgeCount, discountedPairs, mySet, dirtyBits, realDistance, discounters);

	while (Q.getCount() > 0) {
		bool foundOne = false;
		std::pair<unsigned int, double> minQ;
		while (!foundOne) {
			minQ = Q.getMin();

			if (dirtyBits[minQ.first])
				UpdatePair2(minQ.first, Q, QContent, pairs[minQ.first], vertices, myEdges, t, myHeap, mySet, dirtyBits, realDistance, discounters);
			else
				foundOne = true;
			if (Q.getCount() == 0)
				goto end;
		}
        Q.extractMin();
		WellSeparatedPairRepresentative* minRepresentative = QContent[minQ.first];

		std::pair<unsigned int, double> pair(minRepresentative->u, minRepresentative->dist);
        myEdges[minRepresentative->v].push_back(pair);
        pair.first = minRepresentative->v;
        myEdges[minRepresentative->u].push_back(pair);
        if (T->create_edge(edge(minRepresentative->u, minRepresentative->v))) return 0;
        edgeCount++;

		for (int j = 0; j < m; j++)
			if (QContent[j] != nullptr)
				dirtyBits[j] = true;

		ReinsertPair2(minQ.first, Q, QContent, pairs[minQ.first], vertices, myEdges, t, myHeap, mySet, dirtyBits, realDistance, discounters);
	}
end:
	ENDLOG
	delete[] myEdges;
	delete splitTree;
	for (unsigned int i = 0; i < m; i++) {
		if (!discountedPairs[i]) {
			delete[] discounters[i];
		}
	}
	delete[] discounters;
	delete[] discountedPairs;
	delete[] mySet;
	delete[] realDistance;
	delete[] dirtyBits;

    return edgeCount;
}

template<class DistAlgo, class Heap>
int GreedySpannerPrelim(BaseTask * T) {
	STARTLOG
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
    std::vector<EdgeInfo> edgeList;
	const unsigned int totalEdges = (N-1) * N / 2;
    edgeList.reserve(totalEdges);
    for (unsigned int y = 0; y < N; y++)
        for (unsigned int x = y+1; x < N; x++) {
            EdgeInfo e(x, y, distance(vertices[x],vertices[y]));
            edgeList.push_back(e);
        }
    std::sort(edgeList.begin(), edgeList.end());

    double* matrix = new double[N * N];
    for (unsigned int y = 0; y < N; y++)
        for (unsigned int x = y+1; x < N; x++) {
            matrix[y + N * x] = std::numeric_limits<double>::infinity();
            matrix[x + N * y] = std::numeric_limits<double>::infinity();
		}
	for (unsigned int y = 0; y < N; y++)
        matrix[y + N * y] = 0;

	std::vector<unsigned int> Eindices;
	double Lcurrent = T->var;
	unsigned int index = 0;
	while (index < totalEdges) {
		Eindices.push_back(index);
		while (index < totalEdges
			   &&
			   edgeList[index].distance < Lcurrent)
			index++;
		Lcurrent *= 2;
	}

    unsigned int edgeCount = 0;
	unsigned int E1 = Eindices.size() > 1 ? Eindices[1] : totalEdges;
    std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
    for (unsigned int k = 0; k < E1; k++) {
        double computedDist = DistAlgo::computeDistance(vertices, myEdges, edgeList[k].x, edgeList[k].y);
        if (edgeList[k].distance * T->t < computedDist) {
            std::pair<unsigned int, double> pair(edgeList[k].x, edgeList[k].distance);
            myEdges[edgeList[k].y].push_back(pair);
            pair.first = edgeList[k].y;
            myEdges[edgeList[k].x].push_back(pair);
            if (T->create_edge(edgeList[k])) { edgeCount = 0; goto cleanup; }
            edgeCount++;
        }
    }
	
	for (unsigned int i = 1; i < Eindices.size(); i++) {
		double Li = (1 << i - 1) * T->var;
		for (unsigned int j = 0; j < N; j++)
            updateDistances<Heap>(vertices, myEdges, j, matrix);
		unsigned int next = totalEdges;
		for (unsigned int k = Eindices[i]; k < next; k++) {
			next = i < Eindices.size() - 1 ? Eindices[i+1] : totalEdges;

            if (edgeList[k].distance * T->t < matrix[edgeList[k].y + N * edgeList[k].x]) {
                std::pair<unsigned int, double> pair(edgeList[k].x, edgeList[k].distance);
                myEdges[edgeList[k].y].push_back(pair);
                pair.first = edgeList[k].y;
                myEdges[edgeList[k].x].push_back(pair);
				if (T->create_edge(edgeList[k])) { edgeCount = 0; goto cleanup; }
                edgeCount++;

				for (unsigned int p = 0; p < N; p++)
					if (distance(vertices[p], vertices[edgeList[k].x]) < (T->t-0.5) * Li
						||
						distance(vertices[p], vertices[edgeList[k].y]) < (T->t-0.5) * Li)
						updateDistances<Heap>(vertices, myEdges, p, matrix);
            }
		}
	}

cleanup:
	ENDLOG
    delete[] matrix;
    delete[] myEdges;
    return edgeCount;
}

struct DijkstraInformation {
	DijkstraInformation(bool extracted, unsigned int first, double second) : extracted(extracted), first(first), second(second) { }
	bool extracted;
	unsigned int first;
	double second;
};

double& matrixA(double* matrix, unsigned int from, unsigned int to, unsigned int N) {
	return matrix[from + N * to];
}

template<class Heap>
void DijkstraBounded(unsigned int N,
					 std::vector<std::pair<unsigned int, double> >* myEdges,
					 unsigned int from,
					 double L,
					 Heap* queue,
					 std::stack<DijkstraInformation>* undoInformation,
					 double* matrix
					 ) {
	qcount++;
    while (queue->getCount() > 0
		   &&
		   queue->getMin().second <= L) {
        std::pair<unsigned int, double> pair = queue->getMin();
		undoInformation->push(DijkstraInformation(true, pair.first, pair.second));
        queue->extractMin();
        for (unsigned int i = 0; i < myEdges[pair.first].size(); i++) {
            std::pair<unsigned int, double> edge = myEdges[pair.first][i];
			if (edge.first != from) {
				double alt = pair.second + edge.second;
				double old = matrixA(matrix, from, edge.first, N);
				if (alt < old) {
					undoInformation->push(DijkstraInformation(false, edge.first, old));
					queue->decreaseKey(edge.first, alt);
					matrixA(matrix, from, edge.first, N) = alt;
				}
			}
        }
    }
}

template<class Heap>
void DijkstraUndo(unsigned int N,
				  unsigned int from,
				  double L,
				  Heap* queue,
				  std::stack<DijkstraInformation>* undoInformation,
				  double* matrix
				  ) {
	while (!undoInformation->empty()
		   &&
		   queue->getMin().second > L) {
		DijkstraInformation c = undoInformation->top();
		undoInformation->pop();
		if (c.extracted) {
			queue->insert(c.first, c.second);
		}
		else {
			queue->increaseKey(c.first, c.second);
			matrixA(matrix, from, c.first, N) = c.second;
		}
	}
}

template<class Heap>
int GreedySpannerNew(BaseTask * T) {
	STARTLOG
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
    std::vector<EdgeInfo> edgeList;
	const unsigned int totalEdges = (N-1) * N / 2;
    edgeList.reserve(totalEdges);
    for (unsigned int y = 0; y < N; y++) {
        for (unsigned int x = y+1; x < N; x++) {
            EdgeInfo e(x, y, distance(vertices[x],vertices[y]));
            edgeList.push_back(e);
        }
    }
    std::sort(edgeList.begin(), edgeList.end());

	std::vector<unsigned int> Eindices;
	std::vector<double> Li;
	unsigned int index = 0;
	while (index < totalEdges) {
		Eindices.push_back(index);
		Li.push_back(edgeList[index].distance);
		double Lcurrent = edgeList[index].distance * 2;
		while (index < totalEdges
			   &&
			   edgeList[index].distance < Lcurrent)
			index++;
	}
	Eindices.push_back(totalEdges);

    std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	Heap** queues = new Heap*[N];
    double* matrix = new double[N * N];
	for (unsigned int i = 0; i < N; i++) {
		Heap* newHeap = new Heap(N, -std::numeric_limits<double>::infinity());
		for (int j = 0; j < N; j++) {
			if (i != j) {
				newHeap->insert(j, std::numeric_limits<double>::infinity());
				matrixA(matrix, i, j, N) = std::numeric_limits<double>::infinity();
			}
			else
				matrixA(matrix, i, j, N) = 0;
		}
		queues[i] = newHeap;
	}

	char filename[80];
	sprintf(filename, "BCFMS(%d-%f).txt", N, T->t);
	std::ofstream myfile(filename);

	std::stack<DijkstraInformation>* undoInformations = new std::stack<DijkstraInformation>[N];

    unsigned int edgeCount = 0;
	for (unsigned int i = 0; i < Eindices.size() - 1; i++) {
		for (unsigned int u = 0; u < N; u++) {
			if (i > 0) {
				DijkstraUndo<Heap>(N, u, (T->t - 0.5) * Li[i-1], queues[u], undoInformations + u, matrix);
			}
			DijkstraBounded<Heap>(N, myEdges, u, 2 * T->t * Li[i], queues[u], undoInformations + u, matrix);
		}

		for (unsigned int k = Eindices[i]; k < Eindices[i+1]; k++) {
			unsigned int u = edgeList[k].x;
			unsigned int v = edgeList[k].y;
            if (edgeList[k].distance * T->t < matrixA(matrix, v, u, N)
				&&
				edgeList[k].distance * T->t < matrixA(matrix, u, v, N)) {
                std::pair<unsigned int, double> pair(u, edgeList[k].distance);
                myEdges[v].push_back(pair);
                pair.first = v;
                myEdges[u].push_back(pair);
				if (T->create_edge(edgeList[k])) { edgeCount = 0; goto cleanup; }
                edgeCount++;
				
				//cout << edgeCount<<"\t" <<edgeList[k].x <<","<< edgeList[k].y << endl;
				myfile << "-" << edgeList[k].x << "&" << edgeList[k].y << "-" << "(" << vertices[edgeList[k].x].x << ", " << vertices[edgeList[k].x].y << ") ---> (" << vertices[edgeList[k].y].x << ", " << vertices[edgeList[k].y].y << ")" << endl;

				for (unsigned int p = 0; p < N; p++) {
					if (distance(vertices[p], vertices[u]) < (T->t - 0.5) * Li[i]
						||
						distance(vertices[p], vertices[v]) < (T->t - 0.5) * Li[i]) {
						double alt1 = matrixA(matrix, p, u, N) + edgeList[k].distance;
						double old1 = matrixA(matrix, p, v, N);
						if (alt1 < old1) {
							DijkstraUndo<Heap>(N, p, minimum((T->t - 0.5) * Li[i], Li[i]), queues[p], undoInformations + p, matrix);
							queues[p]->decreaseKey(v, alt1);
							matrixA(matrix, p, v, N) = alt1;
							DijkstraBounded<Heap>(N, myEdges, p, 2 * T->t * Li[i], queues[p], undoInformations + p, matrix);
						}
						double alt2 = matrixA(matrix, p, v, N) + edgeList[k].distance;
						double old2 = matrixA(matrix, p, u, N);
						if (alt2 < old2) {
							DijkstraUndo<Heap>(N, p, minimum((T->t - 0.5) * Li[i], Li[i]), queues[p], undoInformations + p, matrix);
							queues[p]->decreaseKey(u, alt2);
							matrixA(matrix, p, u, N) = alt2;
							DijkstraBounded<Heap>(N, myEdges, p, 2 * T->t * Li[i], queues[p], undoInformations + p, matrix);
						}
					}
				}
			}
		}
	}

cleanup:
	ENDLOG
		//std::cout << "Cleaning up..." << std::endl;
		myfile.close();
    delete[] matrix;
    delete[] myEdges;
    delete[] queues;
    delete[] undoInformations;
    return edgeCount;
}

struct vertex2 {
	vertex2() {}
	vertex2(vertex v, unsigned int i) : v(v), index(i) { }
	vertex v;
	unsigned int index;
};

struct BB {
	BB() {
		left = std::numeric_limits<double>::infinity();
		right = -std::numeric_limits<double>::infinity();
		bottom = std::numeric_limits<double>::infinity();
		top = -std::numeric_limits<double>::infinity();
	}
	BB(double left, double right, double top, double bottom)
		: left(left), right(right), top(top), bottom(bottom) {}
	double left, right, top, bottom;
	inline BB splitLeft(const vertex& v) const { return BB(left, v.x, top, bottom); }
	inline BB splitRight(const vertex& v) const { return BB(v.x, right, top, bottom); }
	inline BB splitTop(const vertex& v) const { return BB(left, right, top, v.y	); }
	inline BB splitBottom(const vertex& v) const { return BB(left, right, v.y, bottom); }
	inline bool contains(const vertex& v) const {
		return v.x >= left && v.x <= right && v.y >= bottom && v.y <= top;
	}
	inline vertex LBCorner() const { return vertex(left, bottom); }
	inline vertex RBCorner() const { return vertex(right, bottom); }
	inline vertex LTCorner() const { return vertex(left, top); }
	inline vertex RTCorner() const { return vertex(right, top); }
	inline vertex Closest(vertex v) const
	{
		if (v.x < left) {
			if (v.y < bottom)
				return LBCorner();
			else if (v.y > top)
				return LTCorner();
			else
				return vertex(left, v.y);
		}
		else if (v.x > right) {
			if (v.y < bottom)
				return RBCorner();
			else if (v.y > top)
				return RTCorner();
			else
				return vertex(right, v.y);
		}
		else if (v.y < bottom)
			return vertex(v.x, bottom);
		else if (v.y > top)
			return vertex(v.x, top);
		else
			return v;
	}
};

inline bool Xpred(const vertex2& a, const vertex2& b) { return (a.v.x < b.v.x); }
inline bool Ypred(const vertex2& a, const vertex2& b) { return (a.v.y < b.v.y); }

class KDTree {
private:
	KDTree(vertex2 location, bool sortOnX, BB box) : location(location), sortOnX(sortOnX), box(box) {
		left = right = nullptr;
	}
	bool IsNotContainedInHyperbola(vertex w, vertex v, std::vector<std::pair<unsigned int, double> >& edges, double t, const pointset &vertices) {
		bool isNotContained = true;
		for (unsigned int j = 0; j < edges.size(); j++) {
			unsigned int u = edges[j].first;
			double vu = edges[j].second;
			double uw = distance(vertices[u], w);
			double vw = distance(v, w);
			if (vu + uw * t <= vw * t) {
				isNotContained = false;
				break;
			}
		}
		return isNotContained;
	}
	bool IsFeasiblePoint(vertex w, vertex v, double& best, std::vector<std::pair<unsigned int, double> >& edges, double t, const pointset &vertices) {
		return IsNotContainedInHyperbola(w, v, edges, t, vertices)
				&&
				distance(w, v) < best;
	}
	void ContinueConstrainedNNSearch2(vertex v, unsigned int& result, double& best, std::vector<std::pair<unsigned int, double> >& edges, double t, const pointset &vertices) {
		if (box.contains(v))
			ConstrainedNNSearch2(v, result, best, edges, t, vertices);
		else if (IsFeasiblePoint(box.LBCorner(), v, best, edges, t, vertices)
				||
				IsFeasiblePoint(box.RBCorner(), v, best, edges, t, vertices)
				||
				IsFeasiblePoint(box.LTCorner(), v, best, edges, t, vertices)
				||
				IsFeasiblePoint(box.RTCorner(), v, best, edges, t, vertices)
				||
				IsFeasiblePoint(box.Closest(v), v, best, edges, t, vertices))
				ConstrainedNNSearch2(v, result, best, edges, t, vertices);
	}
	void ConstrainedNNSearch2(vertex v, unsigned int& result, double& best, std::vector<std::pair<unsigned int, double> >& edges, double t, const pointset &vertices) {
		double dist = distance(location.v, v);
		if (dist > 1e-10) {
			if (dist < best && IsNotContainedInHyperbola(location.v, v, edges, t, vertices)) {
				best = dist;
				result = location.index;
			}
		}
		if (left != nullptr) {
			if (right != nullptr) {
				if (distance(left->location.v, v) < distance(right->location.v, v)) {
					left->ContinueConstrainedNNSearch2(v, result, best, edges, t, vertices);
					right->ContinueConstrainedNNSearch2(v, result, best, edges, t, vertices);
				}
				else {
					right->ContinueConstrainedNNSearch2(v, result, best, edges, t, vertices);
					left->ContinueConstrainedNNSearch2(v, result, best, edges, t, vertices);
				}
			}
			else
				left->ContinueConstrainedNNSearch2(v, result, best, edges, t, vertices);
		}
		else if (right != nullptr)
			right->ContinueConstrainedNNSearch2(v, result, best, edges, t, vertices);
	}
	vertex2 location;
	bool sortOnX;
	KDTree* left;
	KDTree* right;
	BB box;
public:
	KDTree(std::vector<vertex2>& vertices, bool sortOnX, BB box) : box(box), sortOnX(sortOnX) {
		unsigned int medianlocation = vertices.size() / 2;
		if (sortOnX)
			std::sort(vertices.begin(), vertices.end(), Xpred);
		else
			std::sort(vertices.begin(), vertices.end(), Ypred);
		location = vertices[medianlocation];
		if (vertices.size() > 1) {
			std::vector<vertex2> leftV;
			std::vector<vertex2> rightV;
			for (unsigned int i = 0; i < medianlocation; i++)
				leftV.push_back(vertices[i]);
			for (unsigned int i = medianlocation + 1; i < vertices.size(); i++)
				rightV.push_back(vertices[i]);
			left = leftV.size() > 0 ? new KDTree(leftV, !sortOnX, sortOnX ? box.splitLeft(location.v) : box.splitBottom(location.v)) : nullptr;
			right = rightV.size() > 0 ? new KDTree(rightV, !sortOnX, sortOnX ? box.splitRight(location.v) : box.splitTop(location.v)) : nullptr;
		}
		else
			left = right = nullptr;
	}
	bool ConstrainedNNSearch(vertex v, unsigned int& result, std::vector<std::pair<unsigned int, double> >& edges, double t, const pointset &vertices) {
		double best = std::numeric_limits<double>::infinity();
		ConstrainedNNSearch2(v, result, best, edges, t, vertices);
		return !is_infinite(best);
	}
};

int HyperbolaKDSpanner(BaseTask * T) {
	STARTLOG
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
	unsigned int edgeCount = 0;

	BB box;
	std::vector<vertex2> vertices2;
	for (unsigned int i = 0; i < N; i++) {
		vertex v = vertices[i];
		vertices2.push_back(vertex2(v, i));
		if (box.left > v.x)
			box.left = v.x;
		if (box.right < v.x)
			box.right = v.x;
		if (box.bottom > v.y)
			box.bottom = v.y;
		if (box.top < v.y)
			box.top = v.y;
	}
	KDTree tree(vertices2, true, box);
	for (unsigned int i = 0; i < N; i++) {
		vertex v = vertices[i];
		std::vector<std::pair<unsigned int, double> > edges;
		unsigned int nearestNeighbor = 0;
		while (tree.ConstrainedNNSearch(v, nearestNeighbor, edges, T->t, vertices)) {
			edge e;
			e.x = i;
			e.y = nearestNeighbor;
			if (T->create_edge(e)) return 0;
			edgeCount++;
			std::pair<unsigned int, double> p;
			p.first = nearestNeighbor;
			p.second = distance(v, vertices[nearestNeighbor]);
			edges.push_back(p);
		}
	}
	ENDLOG
    return edgeCount;
}

struct ThetaEDHybrid {
	static void ComputeFirstSpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>* edges) {
		ComputeThetaSpanner(T, edgeCount, edges);
	}
	static void ComputeSecondSpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>& edgeList) {
		ComputeGreedySpannerED(T, edgeCount, edgeList);
	}
};

template<class Heap>
struct ThetaFGHybrid {
	static void ComputeFirstSpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>* edges) {
		ComputeThetaSpanner(T, edgeCount, edges);
	}
	static void ComputeSecondSpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>& edgeList) {
		ComputeGreedySpannerFG<Heap, false>(T, edgeCount, edgeList);
	}
};

struct ThetaSymHypHybrid {
	static void ComputeFirstSpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>* edges) {
		ComputeThetaSpanner(T, edgeCount, edges);
	}
	static void ComputeSecondSpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>& edgeList) {
		ComputeSymmetricHyperbolaSpanner(T, edgeCount, edgeList);
	}
};

template<bool IsImplicit>
struct WspdEDHybrid {
	static void ComputeFirstSpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>* edges) {
		ComputeWspdSpanner<IsImplicit>(T, edgeCount, edges);
	}
	static void ComputeSecondSpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>& edgeList) {
		ComputeGreedySpannerED(T, edgeCount, edgeList);
	}
};

template<class Heap, bool IsImplicit>
struct WspdFGHybrid {
	static void ComputeFirstSpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>* edges) {
		ComputeWspdSpanner<IsImplicit>(T, edgeCount, edges);
	}
	static void ComputeSecondSpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>& edgeList) {
		ComputeGreedySpannerFG<Heap, false>(T, edgeCount, edgeList);
	}
};

template<bool IsImplicit>
struct WspdSymHypHybrid {
	static void ComputeFirstSpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>* edges) {
		ComputeWspdSpanner<IsImplicit>(T, edgeCount, edges);
	}
	static void ComputeSecondSpanner(BaseTask * T, unsigned int& edgeCount, std::vector<EdgeInfo>& edgeList) {
		ComputeSymmetricHyperbolaSpanner(T, edgeCount, edgeList);
	}
};

template<class Algorithms>
int HybridSpanner(BaseTask * T) {
	unsigned int edgeCount = 0;
	std::vector<EdgeInfo> edges;
	double t = T->t;
	T->t = pow(t, T->var);
	Algorithms::ComputeFirstSpanner(T, edgeCount, &edges);
	if (edgeCount > 0) {
		T->t = pow(t, 1.0 - T->var);
		Algorithms::ComputeSecondSpanner(T, edgeCount, edges);
	}
    return edgeCount;
}

bool controlPointHyperbolaTest(unsigned int v, vertex w1, vertex w2, double t, std::vector<std::pair<unsigned int, double> >* myEdges, const pointset &vertices) {
	for (unsigned int j = 0; j < myEdges[v].size(); j++) {
		unsigned int u = myEdges[v][j].first;
		double vu = myEdges[v][j].second;
		double uw1 = distance(vertices[u], w1);
		double vw1 = distance(vertices[v], w1);
		if (vu + uw1 * t <= vw1 * t) {
			double vu = myEdges[v][j].second;
			double uw2 = distance(vertices[u], w2);
			double vw2 = distance(vertices[v], w2);
			if (vu + uw2 * t <= vw2 * t)
				return true;
		}
	}
	return false;
}

int ThetaSymHypSparsificationSpanner(BaseTask * T) {
	STARTLOG
	unsigned int edgeCount = 0;
	std::vector<EdgeInfo> edges;
	std::vector<std::pair<vertex, vertex>> controlPoints;
	const pointset &vertices = T->p;
    unsigned int N = vertices.size();
	unsigned int k = 8;
	double t = std::numeric_limits<double>::infinity();
	double theta = std::numeric_limits<double>::infinity();
	while (t > T->t) {
		k++;
		theta = (2.0 * M_PI) / ((double)k);
		t = 1.0 / (cos(theta)-sin(theta));
	}
	std::vector<vertex> verticesCopy(N);
	for (unsigned int i = 0; i < N; i++)
		verticesCopy[i] = vertices[i];

	double eps = 1.01;
	double lineA = -tan(theta / 2.0 * eps);
	double lineAR = tan(theta / 2.0 * eps);
	for (unsigned int c = 0; c < k; c++) {
		for (unsigned int i = 0; i < N; i++)
			verticesCopy[i] = rotate(verticesCopy[i], theta);

		ScapegoatHeap<double, MyData, MyAugmentation> heap;
		std::vector<ScapegoatNode<double, MyData>*> nodes;

		for (unsigned int i = 0; i < N; i++) {
			vertex v = verticesCopy[i];
			double intersectionX = v.x - v.y / lineA;
			MyData data;
			data.point = v;
			data.pointIndex = i;
			data.intersectionX2 = v.x - v.y / lineAR;
			nodes.push_back(heap.insert(intersectionX, data));
		}

		std::sort(nodes.begin(), nodes.end(), myOrdering);

		for (unsigned int i = 0; i < nodes.size(); i++) {
			ScapegoatNode<double, MyData>* node = nodes[i];
			heap.remove(node);
			MySearchResult result;
			heap.augmentedSearch(node->key, result);

			if (result.foundSomething) {
				int ex = node->value.pointIndex;
				int ey = result.resultI;
				EdgeInfo e(ex, ey, distance(vertices[ex], vertices[ey]));
				edges.push_back(e);
				edgeCount++;
				vertex xv = verticesCopy[e.x];
				vertex yv = verticesCopy[e.y];
				double dist = yv.x - xv.x;
				vertex c1(yv.x, xv.y + dist * tan(theta / 2));
				vertex c2(yv.x, xv.y - dist * tan(theta / 2));
				controlPoints.push_back(std::pair<vertex, vertex>(rotate(c1, -theta * (c+1)), rotate(c2, -theta * (c+1))));
			}
			node->reset();
			delete node;
		}
	}
	if (edgeCount > 0) {
		const pointset &vertices = T->p;
		unsigned int N = vertices.size();
		std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
		std::sort(edges.begin(), edges.end());
		for (unsigned int k = 0; k < edges.size(); k++) {
			EdgeInfo e = edges[k];
			unsigned int v = e.x;
			unsigned int w = e.y;
			if (distance(vertices[v], vertices[w]) < 1e-10) continue;
			if (liesInNoHyperbola(e, T->t, myEdges, vertices)
				&&
				liesInNoHyperbola(e.Invert(), T->t, myEdges, vertices)) {
				if (T->create_edge(e)) { edgeCount = 0; goto cleanup; }
				edgeCount++;
				myEdges[v].push_back(std::pair<unsigned int, double>(w, e.distance));
				myEdges[w].push_back(std::pair<unsigned int, double>(v, e.distance));
			}
		}
		for (unsigned int k = 0; k < edges.size(); k++) {
			EdgeInfo e = edges[k];
			unsigned int v = e.x;
			vertex c1 = controlPoints[k].first;
			vertex w = vertices[e.y];
			vertex c2 = controlPoints[k].second;
			if (!controlPointHyperbolaTest(v, c1, w, T->t, myEdges, vertices)
				||
				!controlPointHyperbolaTest(v, w, c2, T->t, myEdges, vertices)) {
				if (T->create_edge(e)) { edgeCount = 0; goto cleanup; }
			}
		}
		cleanup:
		delete[] myEdges;
	}
	ENDLOG
    return edgeCount;
}

int ThetaSymHypSparsificationSpanner2(BaseTask * T) {
	unsigned int edgeCount = 0;
	std::vector<EdgeInfo> edges;
	ComputeThetaSpanner(T, edgeCount, &edges);
	if (edgeCount > 0) {
		ComputeSymmetricHyperbolaSpanner(T, edgeCount, edges);
	}
    return edgeCount;
}

int ThetaSymHypSparsificationSpanner3(BaseTask * T) {
	unsigned int edgeCount = 0;
	std::vector<EdgeInfo> edges;
	T->t /= T->var;
	ComputeThetaSpanner(T, edgeCount, &edges);
	if (edgeCount > 0) {
		T->t = T->var;
		ComputeSymmetricHyperbolaSpanner(T, edgeCount, edges);
	}
    return edgeCount;
}

bool addEdge(BaseTask * T, unsigned int i, unsigned int j) {
	edge e;
	e.x = i;
	e.y = j;
	return T->create_edge(e);
}

template<bool IsImplicit>
int WspdHypSparsificationSpanner(BaseTask * T) {
    unsigned int edgeCount = 0;
    const pointset &vertices = T->p;
	double s = 4.0 * (T->t + 1.0) / (T->t - 1.0);
	SplitTree<IsImplicit>* splitTree = SplitTree<IsImplicit>::createSplitTree(vertices);
	splitTree->ComputeWspdHypSpanner(T, s, edgeCount);
	delete splitTree;
    return edgeCount;
}

template<bool IsImplicit>
int WspdSymHypSparsificationSpanner(BaseTask * T) {
    unsigned int edgeCount = 0;
    const pointset &vertices = T->p;
	double s = 4.0 * (T->t + 1.0) / (T->t - 1.0);
	SplitTree<IsImplicit>* splitTree = SplitTree<IsImplicit>::createSplitTree(vertices);
	splitTree->ComputeWspdHypSpanner(T, s, edgeCount, true);
	delete splitTree;
    return edgeCount;
}

template<bool IsImplicit>
int WspdSymHypSparsificationSpanner2(BaseTask * T) {
    unsigned int edgeCount = 0;
    const pointset &vertices = T->p;
	double s = 4.0 * (T->t + 1.0) / (T->t - 1.0);
	SplitTree<IsImplicit>* splitTree = SplitTree<IsImplicit>::createSplitTree(vertices);
	splitTree->ComputeWspdSymHypSpanner(T, s, edgeCount);
	delete splitTree;
    return edgeCount;
}

int YaoSpanner(BaseTask * T) {
	STARTLOG
    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
	unsigned int edgeCount = 0;
	unsigned int k = 8;
	double t = std::numeric_limits<double>::infinity();
	double theta = std::numeric_limits<double>::infinity();
	while (t > T->t) {
		k++;
		theta = (2.0 * M_PI) / ((double)k);
		t = 1.0 / (cos(theta)-sin(theta));
	}

	bool* coneBools = new bool[k];
	unsigned int* coneInts = new unsigned int[k];
	double* coneDists = new double[k];
	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < k; j++)
			coneBools[j] = false;
		for (unsigned int j = 0; j < N; j++) {
			double dist = distance(vertices[i], vertices[j]);
			if (dist > 1e-10) {
				double angle = atan2(vertices[j].y - vertices[i].y, vertices[j].x - vertices[i].x);
				if (angle < 0) angle = 2.0 * M_PI + angle;
				unsigned int bucket = floor(angle / theta);
				if (!coneBools[bucket]
					||
					coneDists[bucket] > dist) {
					coneBools[bucket] = true;
					coneDists[bucket] = dist;
					coneInts[bucket] = j;
				}
			}
		}
		for (unsigned int j = 0; j < k; j++) {
			if (coneBools[j]) {
				edge e;
				e.x = i;
				e.y = coneInts[j];
				if (T->create_edge(e)) return 0;
				edgeCount++;
			}
		}
	}
	ENDLOG
	delete[] coneBools;
	delete[] coneInts;
	delete[] coneDists;
    return edgeCount;
}

int DelaunaySpanner(BaseTask * T) {
/*    const pointset &vertices = T->p;
    unsigned int N = vertices.size();
	unsigned int edgeCount = 0;

	std::vector< Delaunay::Point > v;
	for (unsigned int i = 0; i < N; i++) {
		Delaunay::Point tempP;
		tempP[0] = vertices[i].x;
		tempP[1] = vertices[i].y;
		v.push_back(tempP);
	}
	Delaunay delobject(v);
	delobject.Triangulate();
	for (Delaunay::fIterator fit = delobject.fbegin(); fit != delobject.fend(); ++fit) {
		if (addEdge(T, delobject.Org(fit), delobject.Dest(fit))) return 0;
		if (addEdge(T, delobject.Dest(fit), delobject.Apex(fit))) return 0;
		if (addEdge(T, delobject.Apex(fit), delobject.Org(fit))) return 0;
	}

    return edgeCount;*/return 0;
}

int FromFile(BaseTask * T) {
	using namespace std;
    char filename[80];
	sprintf(filename,"graph.txt");
	ifstream myfile(filename);
	unsigned int edgeCount;
	myfile >> edgeCount;
	unsigned int x, y;
    for (int i = 0; i < edgeCount; i++) {
		myfile >> x;
		myfile.ignore();
		myfile >> y;
		T->create_edge(edge(x, y));
    }
	myfile.close();

    return edgeCount;
}

//****************************************
double shortest_T_Path( std::vector<std::pair<unsigned int, double>>* myEdges, unsigned int N, double* matrix, int source, int target)
{
	double nearestNeighborToTargetDistance = INFINITY, d;// , tpath = distance(vertices[source], vertices[target])*t;

	int ns = myEdges[source].size();

	for (size_t i = 0; i < ns; i++)
	{
		d = matrix[source * N + myEdges[source][i].first] + matrix[myEdges[source][i].first * N + target];

		
		if (d < nearestNeighborToTargetDistance)
			nearestNeighborToTargetDistance = d;
	}

	return nearestNeighborToTargetDistance;
}
/**********************************************/
bool can_have_t_Path(BaseTask * T, std::vector<std::pair<unsigned int, double> >* myEdges, int source, int target)
{
	double tPathStartToTarget = T->t * distance(T->p[source], T->p[target]);
	int v, v1;
	std::vector<std::pair<unsigned int, double>> np, nq;

	for (size_t i = 0; i < myEdges[source].size(); i++)
	{
		if (distance(T->p[source], T->p[myEdges[source][i].first]) + distance(T->p[myEdges[source][i].first], T->p[target]) <= tPathStartToTarget)
			np.push_back(myEdges[source][i]);
	}

	for (size_t i = 0; i < myEdges[target].size(); i++)
	{
		if (distance(T->p[target], T->p[myEdges[target][i].first]) + distance(T->p[myEdges[target][i].first], T->p[source]) <= tPathStartToTarget)
			nq.push_back(myEdges[target][i]);
	}

	for (size_t i = 0; i < np.size(); i++)
	{
		for (size_t j = 0; j < nq.size(); j++)
		{
			if (distance(T->p[source], T->p[np[i].first]) + distance(T->p[np[i].first], T->p[nq[j].first]) + distance(T->p[target], T->p[nq[j].first]) <= tPathStartToTarget)
				return true;
		}
	}

	return false;
}

//**************************************************************************
template<class Heap>
int ComputeGreedySpannerIFGBN1(BaseTask * T) {
	STARTLOG
		const pointset &vertices = T->p;
	unsigned int N = vertices.size();
	double* matrix = new double[N * N];
	double dist12, dist21;

	std::vector<EdgeInfo> edgeList;
	edgeList.reserve((N - 1) * N / 2);
	for (unsigned int y = 0; y < N; y++) {
		for (unsigned int x = y + 1; x < N; x++) {
			EdgeInfo e(x, y, distance(vertices[x], vertices[y]));
			edgeList.push_back(e);
		}
	}

	for (unsigned int x = 0; x < N; x++)
		for (unsigned int y = 0; y < N; y++)
			matrix[y + N * x] = std::numeric_limits<double>::infinity();
	for (unsigned int y = 0; y < N; y++)
		matrix[y + N * y] = 0;

	std::sort(edgeList.begin(), edgeList.end());
	std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	for (unsigned int k = 0; k < edgeList.size(); k++) {
		if (edgeList[k].distance * T->t < matrix[edgeList[k].y + N * edgeList[k].x]) {

			dist12 = shortest_T_Path(myEdges, N, matrix, edgeList[k].y, edgeList[k].x);
			dist21 = shortest_T_Path(myEdges, N, matrix, edgeList[k].x, edgeList[k].y);
			if (dist12 > edgeList[k].distance * T->t && dist21 > edgeList[k].distance * T->t)
			{
				updateDistances<Heap>(vertices, myEdges, edgeList[k].x, matrix);

				if (edgeList[k].distance * T->t < matrix[edgeList[k].y + N * edgeList[k].x]) {
					std::pair<unsigned int, double> pair(edgeList[k].x, edgeList[k].distance);
					myEdges[edgeList[k].y].push_back(pair);
					pair.first = edgeList[k].y;
					myEdges[edgeList[k].x].push_back(pair);
					if (T->create_edge(edgeList[k])) { goto cleanup; }
				}
			}
			else
			{
				matrix[edgeList[k].y + N * edgeList[k].x] = matrix[edgeList[k].x + N * edgeList[k].y] = dist12<dist21 ? dist12 : dist21;
			}
		}
	}
cleanup:
	ENDLOG
		double vm, rss;
	delete[] matrix;
	delete[] myEdges;
	return  0;
}

//**************************************************************************
template<class Heap>
int ComputeGreedySpannerIFGBN2(BaseTask * T) {
	STARTLOG
		const pointset &vertices = T->p;
	unsigned int N = vertices.size();
	double* matrix = new double[N * N];

	std::vector<EdgeInfo> edgeList;
	edgeList.reserve((N - 1) * N / 2);
	for (unsigned int y = 0; y < N; y++) {
		for (unsigned int x = y + 1; x < N; x++) {
			EdgeInfo e(x, y, distance(vertices[x], vertices[y]));
			edgeList.push_back(e);
		}
	}

	for (unsigned int x = 0; x < N; x++)
		for (unsigned int y = 0; y < N; y++)
			matrix[y + N * x] = std::numeric_limits<double>::infinity();
	for (unsigned int y = 0; y < N; y++)
		matrix[y + N * y] = 0;

	std::sort(edgeList.begin(), edgeList.end());
	std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	for (unsigned int k = 0; k < edgeList.size(); k++) {
		if (edgeList[k].distance * T->t < matrix[edgeList[k].y + N * edgeList[k].x]) {
			if (can_have_t_Path(T, myEdges, edgeList[k].y, edgeList[k].x))
			{
				updateDistances<Heap>(vertices, myEdges, edgeList[k].x, matrix);

				if (edgeList[k].distance * T->t < matrix[edgeList[k].y + N * edgeList[k].x]) {
					std::pair<unsigned int, double> pair(edgeList[k].x, edgeList[k].distance);
					myEdges[edgeList[k].y].push_back(pair);
					pair.first = edgeList[k].y;
					myEdges[edgeList[k].x].push_back(pair);
					if (T->create_edge(edgeList[k])) { goto cleanup; }

				}
			}
			else
			{
				std::pair<unsigned int, double> pair(edgeList[k].x, edgeList[k].distance);
				myEdges[edgeList[k].y].push_back(pair);
				pair.first = edgeList[k].y;
				myEdges[edgeList[k].x].push_back(pair);
				if (T->create_edge(edgeList[k])) {  goto cleanup; }

				matrix[edgeList[k].y + N * edgeList[k].x] = matrix[edgeList[k].x + N * edgeList[k].y] = edgeList[k].distance;
			}
		}
	}
cleanup:
	ENDLOG
		double vm, rss;
	delete[] matrix;
	delete[] myEdges;
	return 0;
}

//**************************************************************************
template<class Heap>
int ComputeGreedySpannerIFGBN(BaseTask * T) {
	STARTLOG
		const pointset &vertices = T->p;
	unsigned int N = vertices.size();
	double* matrix = new double[N * N];
	double dist12, dist21;
	std::vector<EdgeInfo> edgeList;
	edgeList.reserve((N - 1) * N / 2);
	for (unsigned int y = 0; y < N; y++) {
		for (unsigned int x = y + 1; x < N; x++) {
			EdgeInfo e(x, y, distance(vertices[x], vertices[y]));
			edgeList.push_back(e);
		}
	}

	for (unsigned int x = 0; x < N; x++)
		for (unsigned int y = 0; y < N; y++)
			matrix[y + N * x] = std::numeric_limits<double>::infinity();
	for (unsigned int y = 0; y < N; y++)
		matrix[y + N * y] = 0;

	char filename[80];
	sprintf(filename, "IFGBN(%d-%f).txt", N, T->t);
	std::ofstream myfile(filename, std::ofstream::out);

	std::sort(edgeList.begin(), edgeList.end());
	std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	for (unsigned int k = 0; k < edgeList.size(); k++) {
		if (edgeList[k].distance * T->t < matrix[edgeList[k].y + N * edgeList[k].x]) {
			
			dist12 = shortest_T_Path(myEdges, N, matrix, edgeList[k].y, edgeList[k].x);
			dist21 = shortest_T_Path(myEdges, N, matrix, edgeList[k].x, edgeList[k].y);
			if (dist12 > edgeList[k].distance * T->t && dist21 > edgeList[k].distance * T->t)
			{
				if (can_have_t_Path(T, myEdges, edgeList[k].y, edgeList[k].x))
				{
					updateDistances<Heap>(vertices, myEdges, edgeList[k].x, matrix);

					if (edgeList[k].distance * T->t < matrix[edgeList[k].y + N * edgeList[k].x]) {
						std::pair<unsigned int, double> pair(edgeList[k].x, edgeList[k].distance);
						myEdges[edgeList[k].y].push_back(pair);
						pair.first = edgeList[k].y;
						myEdges[edgeList[k].x].push_back(pair);
						if (T->create_edge(edgeList[k])) {  goto cleanup; }
						myfile << "-" << edgeList[k].x << "&" << edgeList[k].y << "-" << "(" << vertices[edgeList[k].x].x << ", " << vertices[edgeList[k].x].y << ") ---> (" << vertices[edgeList[k].y].x << ", " << vertices[edgeList[k].y].y << ")" << endl;
					}
				}
				else
				{
					std::pair<unsigned int, double> pair(edgeList[k].x, edgeList[k].distance);
					myEdges[edgeList[k].y].push_back(pair);
					pair.first = edgeList[k].y;
					myEdges[edgeList[k].x].push_back(pair);
					if (T->create_edge(edgeList[k])) { goto cleanup; }
					myfile << "-" << edgeList[k].x << "&" << edgeList[k].y << "-" << "(" << vertices[edgeList[k].x].x << ", " << vertices[edgeList[k].x].y << ") ---> (" << vertices[edgeList[k].y].x << ", " << vertices[edgeList[k].y].y << ")" << endl;

					matrix[edgeList[k].y + N * edgeList[k].x] = matrix[edgeList[k].x + N * edgeList[k].y] = edgeList[k].distance;
				}
			}
			else
				matrix[edgeList[k].y + N * edgeList[k].x] = matrix[edgeList[k].x + N * edgeList[k].y] = dist12<dist21 ? dist12 : dist21;
		}
	}
cleanup:
	ENDLOG
		myfile.close();
		double vm, rss;
	delete[] matrix;
	delete[] myEdges;
	return 0;
}

template<class Heap>
int ComputeGreedySpannerIFGBC(BaseTask * T) {
	STARTLOG
		const pointset &vertices = T->p;
	unsigned int N = vertices.size();
	double* matrix = new double[N * N];
	double dist12, dist21;

	std::vector<EdgeInfo> edgeList;
	edgeList.reserve((N - 1) * N / 2);
	for (unsigned int y = 0; y < N; y++) {
		for (unsigned int x = y + 1; x < N; x++) {
			EdgeInfo e(x, y, distance(vertices[x], vertices[y]));
			edgeList.push_back(e);
		}
	}

	for (unsigned int x = 0; x < N; x++)
		for (unsigned int y = 0; y < N; y++)
			matrix[y + N * x] = std::numeric_limits<double>::infinity();
	for (unsigned int y = 0; y < N; y++)
		matrix[y + N * y] = 0;

	std::sort(edgeList.begin(), edgeList.end());
	std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	for (unsigned int k = 0; k < edgeList.size(); k++) {
		if (edgeList[k].distance * T->t < matrix[edgeList[k].y + N * edgeList[k].x]) {
			ModifiedUpdateDistances<Heap>(vertices, myEdges, edgeList[k].x, matrix, T->t);

			if (edgeList[k].distance * T->t < matrix[edgeList[k].y + N * edgeList[k].x]) {
				std::pair<unsigned int, double> pair(edgeList[k].x, edgeList[k].distance);
				myEdges[edgeList[k].y].push_back(pair);
				pair.first = edgeList[k].y;
				myEdges[edgeList[k].x].push_back(pair);
				if (T->create_edge(edgeList[k])) { goto cleanup; }
			}
		}
	}
cleanup:
	ENDLOG
		double vm, rss;
	delete[] matrix;
	delete[] myEdges;
	return  0;
}

//************************************************************************************************
const double PI = 3.14159265358979323846264338327950288419716939937510L;
class Cone
{
public:
	vertex apex;
	vertex bisector;
	vertex start;
	vertex end;
	double	angle;

	Cone(const vertex& ap, const vertex& bs, double ang)
	{
		apex = ap;
		bisector = bs;
		angle = ang;

		double  angleToRadian = angle; //(PI / 180.0) * angle;
		vertex   b;
		b.x = bisector.x - apex.x;
		b.y = bisector.y - apex.y;

		start.x = b.x * cos(angleToRadian / 2.0) - b.y * sin(angleToRadian / 2.0);
		start.y = b.x * sin(angleToRadian / 2.0) + b.y * cos(angleToRadian / 2.0);

		end.x = b.x * cos(-angleToRadian / 2.0) - b.y * sin(-angleToRadian / 2.0);
		end.y = b.x * sin(-angleToRadian / 2.0) + b.y * cos(-angleToRadian / 2.0);

	}

	bool IsIn(const vertex* p)
	{
		vertex relPoint;
		relPoint.x = p->x - apex.x;
		relPoint.y = p->y - apex.y;

		bool k = areClockwise(start, relPoint) && !areClockwise(end, relPoint);

		return k;
	}

	bool areClockwise(vertex v1, vertex v2)
	{
		return -v1.x * v2.y + v1.y * v2.x > 0;
	}
};

class cone_list_node {
public:
	Cone *data;
	cone_list_node *next;
};

class cone_list {
public:
	cone_list_node *first;
	int size;
	cone_list() {
		first = NULL;
		size = 0;
	}
	~cone_list() {
		cone_list_node *t, *f = first;
		while (f != NULL) {
			t = f;
			f = f->next;
			delete t;
		}
	}
	void add(Cone *d) {
		cone_list_node *p = first;
		first = new cone_list_node;
		first->next = p;
		first->data = d;
		size++;
	}
};

//*****************************************
void	AddCone(cone_list* list, const pointset &vertices, int apex, int bisector, double	angle)
{
	Cone	*cone = new Cone(vertices[apex], vertices[bisector], angle);

	list[apex].add(cone);
}
//***************************************
bool IsInCone(cone_list*	Cones, const pointset &vertices, int apex, int p)
{
	cone_list_node* v_item = Cones[apex].first;
	Cone*	cone;
	int c = 0;
	while (v_item != NULL)
	{
		c++;
		cone = v_item->data;
		v_item = v_item->next;
		if (cone->IsIn(&vertices[p]))
		{
			//cout<<"True," <<"apex:"<<apex<<",point:"<<p<<",cone count:"<<c<<endl;
			return true;
		}
	}
	//cout <<"false,"<<"apex:"<<apex<<",point:"<<p<<",cone count:"<<c<<endl;
	return false;
}
//************************************************
template<class Heap>

int DeltaGreedy(BaseTask* T)
{
	STARTLOG
		pointset vertices = T->p;
	cone_list* Cones = new cone_list[T->p.size()];
	double d, g, asn, tetha, delta = T->t;
	double* dist = new double[T->p.size()];
	int c = 0, N = T->p.size();
	std::vector<EdgeInfo> edgeList;

	edgeList.reserve((N - 1) * N / 2);
	for (unsigned int y = 0; y < N; y++) {
		for (unsigned int x = y + 1; x < N; x++) {
			EdgeInfo e(x, y, distance(vertices[x], vertices[y]));
			edgeList.push_back(e);
		}
	}

	std::sort(edgeList.begin(), edgeList.end());
	std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];

	for (int i = 0; i < (N * (N - 1)) / 2; i++)
	{
		if (!IsInCone(Cones, vertices, edgeList[i].x, edgeList[i].y) && !IsInCone(Cones, vertices, edgeList[i].y, edgeList[i].x))
		{
			c++;
			double dist = computeDistance<Heap>(vertices, myEdges, edgeList[i].x, edgeList[i].y);
			d = dist / edgeList[i].distance;

			if (d > delta)
			{
				std::pair<unsigned int, double> pair(edgeList[i].x, edgeList[i].distance);
				myEdges[edgeList[i].y].push_back(pair);
				pair.first = edgeList[i].y;
				myEdges[edgeList[i].x].push_back(pair);
				if (T->create_edge(edgeList[i])) { goto cleanup; }
				d = 1;
			}
			g = d / (sqrt(2.0) * T->t);
			asn = asin(g) * (180.0 / PI);
			tetha = 45.0 - asn;

			AddCone(Cones, vertices, edgeList[i].x, edgeList[i].y, 2 * tetha);
			AddCone(Cones, vertices, edgeList[i].y, edgeList[i].x, 2 * tetha);
		}
	}
cleanup:
	ENDLOG
		cout << "cleaning ...." << endl;
	delete[]dist;
	delete[]Cones;
	delete[] myEdges;
	return c;
}


//////////////////////Gap Greedy Construction Algorithms\\\\\\\\\\\\\\\\\\\\\\\\

#pragma region Gap Greedy Construction Algorithms

class HeapIndexNode
{
public:
	unsigned int	index;
	HeapIndexNode* next;
};

class HeapIndex
{
public:
	HeapIndexNode** map;
	HeapIndexNode* first, *last;
	unsigned int size;
	unsigned int capacity;

	HeapIndex(unsigned int count)
	{
		capacity = count;
		size = 0;
		map = new HeapIndexNode* [capacity];
		for (unsigned i = 0; i < capacity; i++)
		{
			map[i] = nullptr;
		}
		first = new HeapIndexNode();
		last = first;
		first->next = nullptr;
	}

	~HeapIndex()
	{
		//cout << "clear linked list " << endl;
		Dispose();
		
	}

	void Add(unsigned int index)
	{
		//cout << "Begin adding pair (" << index << "), List Size:" << size << endl;

		HeapIndexNode* n = new HeapIndexNode();
		n->index = index;
		n->next = nullptr;

		map[index] = last;
		last->next = n;
		last = n;
		size++;

		//cout << "pair (" << index << ") is Added, List Size:" << size << endl;
	}

	void Remove(unsigned int index)
	{
		//cout << "Begin remove pair ("<< index <<"), List size:" << size << endl;

		HeapIndexNode* prev = map[index];
		HeapIndexNode* nodeToBeRemoved = prev->next;
		HeapIndexNode* nextNode = nodeToBeRemoved->next;

		if (nodeToBeRemoved == last)
			last = prev;
		else
			map[nextNode->index] = prev;

		prev->next = nextNode;

		delete nodeToBeRemoved;
		map[index] = nullptr;
		size--;

		//cout << "pair (" << index << ") Removed, List size:" << size << endl;
	}

	void Dispose()
	{
		for (size_t i = 0; i < capacity; i++)
		{
			if (map[i] != nullptr)
				delete map[i];
		}
		delete[] map;
	}

	HeapIndexNode* getNode(unsigned int index)
	{
		if (map[index] == nullptr)
			return nullptr;
		else
			return map[index]->next;
	}
};

//*******************************************************************
double getSlope(double dy, double dx)
{
	if (dx == 0)
	{
		if (dy > 0)
			return INFINITY;
		return -1 * INFINITY;
	}
	return dy / dx;
}
//*******************************************************************
bool  haveGapProperty(double theta,
	double w,
	const pointset &vertices,
	unsigned int v1,
	unsigned int v2,
	std::vector<std::pair<unsigned int, unsigned int>>& edges) {

	qcount++;
	double  angle, distss, distee, d;
	double ang1, ang2, ang1_1, ang2_1;
	/*
	ang1 = atan3(vertices[v2], vertices[v1]);

	for (int i = 0; i < edges.size(); i++) {

		ang2 = atan3(vertices[edges[i].second], vertices[edges[i].first]);
		angle = abs(ang1 - ang2);
		*/
	ang1 = atan2(vertices[v1].y - vertices[v2].y, vertices[v1].x - vertices[v2].x);
	if (ang1 < 0)
		ang1_1 = ang1 + 2 * PI;
	else
		ang1_1 = ang1;

	for (int i = 0; i < edges.size(); i++) {

		ang2 = atan2(vertices[edges[i].first].y - vertices[edges[i].second].y, vertices[edges[i].first].x - vertices[edges[i].second].x);
		if (ang2 < 0)
			ang2_1 = ang2 + 2 * PI;
		else
			ang2_1 = ang2;

		if (abs(ang1) < PI / 2.0 &&  abs(ang2) < PI / 2.0)
			angle = abs(ang1 - ang2);
		else
			angle = abs(ang1_1 - ang2_1);
		if (angle <= theta)
		{
			d = w*w* NotSquaredDistance(vertices[edges[i].first], vertices[edges[i].second]);
			distss = NotSquaredDistance(vertices[edges[i].first], vertices[v1]);
			distee = NotSquaredDistance(vertices[edges[i].second], vertices[v2]);

			if (distss <= d || distee <= d)
				return false;
		}
	}

	return true;
}

//***************************************************************
template<class Heap, bool IsImplicit>
int GapGreedyOriginal(BaseTask * T) {
	STARTLOG
		const pointset &vertices = T->p;
	unsigned int N = vertices.size();
	unsigned int edgeCount = 0;
	double w = T->var, theta = T->t;
	unsigned int TotalEdges = ((N - 1) * N) / 2;

	std::vector<EdgeInfo> edgeList;
	edgeList.reserve(TotalEdges);
	for (unsigned int x = 0; x < N; x++) {
		for (unsigned int y = x + 1; y < N; y++) {
			EdgeInfo e(x, y, distance(vertices[x], vertices[y]));
			edgeList.push_back(e);
		}
	}

	std::sort(edgeList.begin(), edgeList.end());

	std::vector<std::pair<unsigned int, unsigned int>> edges;
	
	
	for (unsigned int i = 0; i < edgeList.size(); i++)
	{
		if (haveGapProperty(theta, w, vertices, edgeList[i].x, edgeList[i].y, edges))
		{
			std::pair<unsigned int, unsigned int> e1(edgeList[i].x, edgeList[i].y);
			edges.push_back(e1);

			std::pair<unsigned int, unsigned int> e2(edgeList[i].y, edgeList[i].x);
			edges.push_back(e2);

			if (T->create_edge(edge(edgeList[i].x, edgeList[i].y))) return 0;
			if (T->create_edge(edge(edgeList[i].y, edgeList[i].x))) return 0;
			edgeCount += 2;
			
			//cout << k << "=> (" << vertices[edgeList[i].x].x << ", " << vertices[edgeList[i].x].y << ") ---> (" << vertices[edgeList[i].y].x << ", " << vertices[edgeList[i].y].y << ")" << endl;

			//myfile << "-" << edgeList[i].x << "&" << edgeList[i].y << "-" <<"(" << vertices[edgeList[i].x].x << ", " << vertices[edgeList[i].x].y << ") ---> (" << vertices[edgeList[i].y].x << ", " << vertices[edgeList[i].y].y << ")" << endl;
			
		}
	}

	
	char filename[80];
	sprintf(filename, "GapOrig(%d-%f-%f).txt", N, theta, w);
	std::ofstream myfile(filename);
	PROCESS_MEMORY_COUNTERS_EX pmc;
	GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));

	//SIZE_T physMemUsedByMe = pmc.WorkingSetSize;
	myfile << pmc.WorkingSetSize;

	myfile.close();

	edgeList.clear();
	edges.clear();

	return edgeCount;
}

//******************************************

#pragma region Bakhshesh methods

#pragma region Method 1

//method1
template<class Heap, bool IsImplicit>
WellSeparatedPairRepresentative* SelectedPair(	unsigned int pairIndex,
												WellSeparatedPair<IsImplicit>* pair,
												const pointset &vertices,
												double theta,
												double w,
												std::vector<std::pair<unsigned int, unsigned int> >& edges) {
	unsigned int N = vertices.size();
	rr++;

	unsigned int shortestPairInedx1, shortestPairInedx2;
	double bestDist = std::numeric_limits<double>::infinity();
	unsigned int 	index1, index2;

	for (typename SplitTree<IsImplicit>::Iterator it1 = pair->first->begin(); it1 != pair->first->end(); ++it1) {

		index1 = *it1;
		for (typename SplitTree<IsImplicit>::Iterator it2 = pair->second->begin(); it2 != pair->second->end(); ++it2) {
			index2 = *it2;

			bool r = haveGapProperty(theta, w, vertices, index1, index2, edges);// || haveGapProperty(theta, w, vertices, index2, index1, edges);

			if (r)
			{
				double dist = NotSquaredDistance(vertices[index1], vertices[index2]);
				if (dist < bestDist)
				{
					bestDist = dist;
					shortestPairInedx1 = index1;
					shortestPairInedx2 = index2;
				}
			}
		}
	}
	if (std::isinf(bestDist))
		return nullptr;
	else
		return new WellSeparatedPairRepresentative(shortestPairInedx1, shortestPairInedx2, sqrt(bestDist), pairIndex);
}

template<class Heap, bool IsImplicit>
void GapFillQueue1(	Heap& Q,
					std::vector<WellSeparatedPairRepresentative*>& QContent,
					unsigned int& i,
					std::vector<WellSeparatedPair<IsImplicit>*>& pairs,
					const pointset &vertices,
					double theta,
					double w,
					int edgeCount,
					std::vector<std::pair<unsigned int, unsigned int> > &edges, 
					HeapIndex& heapIndex) {
	while (i < pairs.size() && (Q.getCount() == 0 || pairs[i]->minlength() <= Q.getMin().second)) {
		WellSeparatedPairRepresentative* p = SelectedPair<Heap, IsImplicit>(i, pairs[i], vertices, theta, w, edges);
		//if (i % 40000 == 0)
			//std::cout << i << "/" << pairs.size() << " " << Q.getCount() << " " << edgeCount << endl;
		if (p != nullptr) {
			Q.insert(i, p->dist);
			QContent[i] = p;
			heapIndex.Add(i);
		}
		i++;
	}
}

template<class Heap, bool IsImplicit>
int BakhsheshMethod1(BaseTask * T) {
	STARTLOG
		rr = 0;
		const pointset &vertices = T->p;
	unsigned int N = vertices.size();
	unsigned int edgeCount = 0;

	double w = T->var;
	double theta = T->t;

	double s = maximum(4.0 / sin(theta / 2.0), 2.0 / w); //4.0 * t / (t - 1.0);

	SplitTree<IsImplicit>* splitTree = SplitTree<IsImplicit>::createSplitTree(vertices);
	std::vector<WellSeparatedPair<IsImplicit>*> pairs;
	splitTree->ComputeWspdPairs(s, pairs);
	sort(pairs.begin(), pairs.end(), WellSeparatedPairPointerComparer<IsImplicit>());

	//std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	std::vector<WellSeparatedPairRepresentative*> QContent(pairs.size(), nullptr);
	Heap Q(pairs.size(), -std::numeric_limits<double>::infinity());

	HeapIndex heapIndex(pairs.size());

	std::vector<std::pair<unsigned int, unsigned int>> edges;

	unsigned int index = 0, minIndex =0;
	GapFillQueue1<Heap>(Q, QContent, index, pairs, vertices, theta, w, edgeCount, edges, heapIndex);

	while (Q.getCount() > 0) {
		std::pair<unsigned int, double> minQ = Q.getMin();

		Q.extractMin();
		WellSeparatedPairRepresentative* minRepresentative = QContent[minQ.first];
		QContent[minQ.first] = nullptr;
		heapIndex.Remove(minQ.first);

		/*
		std::pair<unsigned int, double> pair(minRepresentative->u, minRepresentative->dist);
		myEdges[minRepresentative->v].push_back(pair);
		pair.first = minRepresentative->v;
		myEdges[minRepresentative->u].push_back(pair);*/

		std::pair<unsigned int, unsigned int> e1(minRepresentative->u, minRepresentative->v);
		edges.push_back(e1);

		std::pair<unsigned int, unsigned int> e2(minRepresentative->v, minRepresentative->u);
		edges.push_back(e2);

		if (T->create_edge(edge(minRepresentative->u, minRepresentative->v))) return 0;
		if (T->create_edge(edge(minRepresentative->v, minRepresentative->u))) return 0;
		edgeCount += 2;


		WellSeparatedPair<IsImplicit>* first = pairs[minRepresentative->pairIndex];

		//for (int j = 0; j < index; j++) {
		HeapIndexNode* node = heapIndex.first->next;
		while (node != nullptr)
		{
			//if (QContent[j] != nullptr) {
			int j = node->index;
			node = node->next;

			WellSeparatedPair<IsImplicit>* second = pairs[j];
			double distffsf = minDistance(first->first->box, second->first->box);
			double distffss = minDistance(first->first->box, second->second->box);
			double distfssf = minDistance(first->second->box, second->first->box);
			double distfsss = minDistance(first->second->box, second->second->box);
			double lengths = maxDistance(second->first->box, second->second->box) * w + 1e-6;
			if (distffsf <= lengths || distfsss <= lengths || distffss <= lengths || distfssf <= lengths) {
				WellSeparatedPairRepresentative* p = SelectedPair<Heap>(j, second, vertices, theta, w, edges);
				if (p == nullptr) {
					Q.remove(j);
					QContent[j] = nullptr;
					heapIndex.Remove(j);
				}
				else {
					Q.increaseKey(j, p->dist);
					QContent[j] = p;
				}
			}
		}

		GapFillQueue1<Heap>(Q, QContent, index, pairs, vertices, theta, w, edgeCount, edges, heapIndex);

		delete minRepresentative;
	}
	ENDLOG
		//cout << " qc = " << rr << endl;

	char filename[80];
	sprintf(filename, "BM1(%d-%f-%f).txt", N, theta, w);
	std::ofstream myfile(filename);

	PROCESS_MEMORY_COUNTERS_EX pmc;
	GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
	SIZE_T physMemUsedByMe = pmc.WorkingSetSize;

	myfile << pmc.WorkingSetSize;
	myfile << "\t" << pairs.size() << endl;
	myfile.close();

	delete splitTree;

	for (size_t i = 0; i < pairs.size(); i++)
	{
		delete pairs[i];
	}
	pairs.clear();
	edges.clear();
	QContent.clear();
	//delete[] myEdges;
	return edgeCount;
}

#pragma endregion

#pragma region Bakhshesh method 2

template<class Heap, bool IsImplicit>
WellSeparatedPairRepresentative* GapImprovedSelectedPair(	unsigned int pairIndex,
														WellSeparatedPair<IsImplicit>* pair,
														const pointset &vertices,
														double theta,
														double w,
														std::vector<std::pair<unsigned int, unsigned int>> edges,
														std::vector<std::pair<unsigned int, unsigned int>> &Ni)
{
	unsigned int N = vertices.size();
	rr++;
	unsigned int shortestPairInedx1, shortestPairInedx2;
	double bestDist = std::numeric_limits<double>::infinity();
	std::vector<unsigned int> toRemove;

	for (unsigned int i = 0; i < Ni.size(); i++) {
		bool r = haveGapProperty(theta, w, vertices, Ni[i].first, Ni[i].second, edges);

		if (r)
		{
			double dist = distance(vertices[Ni[i].first], vertices[Ni[i].second]);
			if (dist < bestDist)
			{
				bestDist = dist;
				shortestPairInedx1 = Ni[i].first;
				shortestPairInedx2 = Ni[i].second;
			}
		}
		else
		{
			toRemove.push_back(i);
		}
	}

	for (int i = 0; i < toRemove.size(); i++) {
		Ni.erase(Ni.begin() + (toRemove[i] - i));
	}

	if (std::isinf(bestDist))
		return nullptr;
	else
		return new WellSeparatedPairRepresentative(shortestPairInedx1, shortestPairInedx2, bestDist, pairIndex);
}

template<class Heap, bool IsImplicit>
void GapFillQueue2(	Heap& Q,
					std::vector<WellSeparatedPairRepresentative*>& QContent,
					unsigned int& i,
					std::vector<WellSeparatedPair<IsImplicit>*>& pairs,
					const pointset &vertices,
					double theta,
					double w,
					int edgeCount,
					std::vector<std::pair<unsigned int, unsigned int> > &edges,
					std::vector<std::pair<unsigned int, unsigned int>>* Nis)
{
	while (i < pairs.size() && (Q.getCount() == 0 || pairs[i]->minlength() <= Q.getMin().second)) {
		WellSeparatedPairRepresentative* p = GapImprovedSelectedPair<Heap, IsImplicit>(i, pairs[i], vertices, theta, w, edges, Nis[i]);
		//if (i % 40000 == 0)
			//std::cout << i << "/" << pairs.size() << " " << Q.getCount() << " " << edgeCount << endl;
		if (p != nullptr) {
			Q.insert(i, p->dist);
			QContent[i] = p;
		}
		i++;
	}
}

template<class Heap, bool IsImplicit>
int BakhsheshMethod2(BaseTask * T) {
	STARTLOG
		rr = 0;
		const pointset &vertices = T->p;
	unsigned int N = vertices.size();
	unsigned int edgeCount = 0;

	double w = T->var;
	double theta = T->t;

	//double t = T->t < 2.0 ? T->t : 2.0;
	double s = maximum(4.0 / sin(theta / 2.0), 2.0 / w); //4.0 * t / (t - 1.0);
	Heap myHeap(N, 0);

	SplitTree<IsImplicit>* splitTree = SplitTree<IsImplicit>::createSplitTree(vertices);
	std::vector<WellSeparatedPair<IsImplicit>*> pairs;
	splitTree->ComputeWspdPairs(s, pairs);
	sort(pairs.begin(), pairs.end(), WellSeparatedPairPointerComparer<IsImplicit>());

	std::vector<std::pair<unsigned int, unsigned int>> edges;

	std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	std::vector<WellSeparatedPairRepresentative*> QContent(pairs.size(), nullptr);
	Heap Q(pairs.size(), -std::numeric_limits<double>::infinity());
	std::vector<std::pair<unsigned int, unsigned int>>* Nis =
		new std::vector<std::pair<unsigned int, unsigned int>>[pairs.size()];

	char filename[80];
	sprintf(filename, "BM2-pairs(%d).txt", N);

	std::ofstream myfile(filename);

	for (int i = 0; i < pairs.size(); i++) {
		myfile << "pairs (" << i << ")" << endl;
		for (typename SplitTree<IsImplicit>::Iterator it1 = pairs[i]->first->begin(); it1 != pairs[i]->first->end(); ++it1) {
			for (typename SplitTree<IsImplicit>::Iterator it2 = pairs[i]->second->begin(); it2 != pairs[i]->second->end(); ++it2) {
				std::pair<unsigned int, unsigned int> pair(*it1, *it2);
				Nis[i].push_back(pair);
				
				myfile << "(" << vertices[*it1].x << ", " << vertices[*it1].y << ") ---> (" << vertices[*it2].x << ", " << vertices[*it2].y << ")" << endl;
			}
		}
		myfile << endl;
	}

	unsigned int index = 0;
	GapFillQueue2<Heap>(Q, QContent, index, pairs, vertices, theta, w, edgeCount, edges, Nis);

	

	while (Q.getCount() > 0) {
		std::pair<unsigned int, double> minQ = Q.getMin();

		Q.extractMin();
		WellSeparatedPairRepresentative* minRepresentative = QContent[minQ.first];
		QContent[minQ.first] = nullptr;

		std::pair<unsigned int, double> pair(minRepresentative->u, minRepresentative->dist);
		myEdges[minRepresentative->v].push_back(pair);
		pair.first = minRepresentative->v;
		myEdges[minRepresentative->u].push_back(pair);

		/*if (minRepresentative->v > minRepresentative->u)
		{
			myfile << "-" << minRepresentative->v << "&" << minRepresentative->u << "-" << "(" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ") ---> (" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ")" << endl;
			myfile << "-" << minRepresentative->u << "&" << minRepresentative->v << "-" << "(" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ") ---> (" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ")" << endl;
		}
		else
		{
			myfile << "-" << minRepresentative->u << "&" << minRepresentative->v << "-" << "(" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ") ---> (" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ")" << endl;
			myfile << "-" << minRepresentative->v << "&" << minRepresentative->u << "-" << "(" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ") ---> (" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ")" << endl;
		}*/

		std::pair<unsigned int, unsigned int> e1(minRepresentative->u, minRepresentative->v);
		edges.push_back(e1);

		if (T->create_edge(edge(minRepresentative->u, minRepresentative->v))) return 0;
		edgeCount++;

		WellSeparatedPair<IsImplicit>* first = pairs[minRepresentative->pairIndex];
		std::vector<std::pair<unsigned int, unsigned int>> newEdges;
		newEdges.push_back(e1);

		for (int j = 0; j < pairs.size(); j++) {
			if (QContent[j] != nullptr) {
				WellSeparatedPairRepresentative* minRepresentative1 = QContent[j];
				

				if (!haveGapProperty(theta, w, vertices, minRepresentative1->u, minRepresentative1->v, newEdges)) {
					WellSeparatedPair<IsImplicit>* second = pairs[j];
					WellSeparatedPairRepresentative* p = GapImprovedSelectedPair<Heap>(j, second, vertices, theta, w, newEdges, Nis[j]);
					if (p == nullptr) {
						Q.remove(j);
						QContent[j] = nullptr;
					}
					else {
						Q.increaseKey(j, p->dist);
						QContent[j] = p;
					}
				}
			}
		}

		GapFillQueue2<Heap>(Q, QContent, index, pairs, vertices, theta, w, edgeCount, newEdges, Nis);
	}
	ENDLOG
		//cout << " qc = " << rr << endl;
	
		myfile.close();

	delete splitTree;
	delete[] myEdges;
	delete[] Nis;
	return edgeCount;
}

#pragma endregion

#pragma region Improved_Bakhshesh Methods (changing their close criteria)

template<class Heap, bool IsImplicit>
int ImprovedMethod1BakhsheshWithChangingCloseCriteria(BaseTask * T) {
	unsigned int edgeCount = 0;
	/*STARTLOG
		rr = 0;
	const pointset &vertices = T->p;
	unsigned int N = vertices.size();
	

	double w = T->var;
	double theta = T->t;

	//double t = T->t < 2.0 ? T->t : 2.0;
	double s = maximum(4.0 / sin(theta / 2.0), 2.0 / w); //4.0 * t / (t - 1.0);

	SplitTree<IsImplicit>* splitTree = SplitTree<IsImplicit>::createSplitTree(vertices);
	std::vector<WellSeparatedPair<IsImplicit>*> pairs;
	splitTree->ComputeWspdPairs(s, pairs);
	sort(pairs.begin(), pairs.end(), WellSeparatedPairPointerComparer<IsImplicit>());

	std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	std::vector<WellSeparatedPairRepresentative*> QContent(pairs.size(), nullptr);
	Heap Q(pairs.size(), -std::numeric_limits<double>::infinity());

	std::vector<std::pair<unsigned int, unsigned int>> edges;

	unsigned int index = 0;
	GapFillQueue1<Heap>(Q, QContent, index, pairs, vertices, theta, w, edgeCount, edges);

	while (Q.getCount() > 0) {
		std::pair<unsigned int, double> minQ = Q.getMin();

		Q.extractMin();
		WellSeparatedPairRepresentative* minRepresentative = QContent[minQ.first];
		QContent[minQ.first] = nullptr;

		std::pair<unsigned int, double> pair(minRepresentative->u, minRepresentative->dist);
		myEdges[minRepresentative->v].push_back(pair);
		pair.first = minRepresentative->v;
		myEdges[minRepresentative->u].push_back(pair);

		std::pair<unsigned int, unsigned int> e1(minRepresentative->u, minRepresentative->v);
		edges.push_back(e1);

		if (T->create_edge(edge(minRepresentative->u, minRepresentative->v))) return 0;
		edgeCount++;

		WellSeparatedPair<IsImplicit>* first = pairs[minRepresentative->pairIndex];

		for (int j = 0; j < pairs.size(); j++) {
			if (QContent[j] != nullptr) {
				WellSeparatedPairRepresentative* minRepresentative1 = QContent[j];
				std::vector<std::pair<unsigned int, unsigned int>> newEdges;
				newEdges.push_back(e1);

				if (!haveGapProperty(theta, w, vertices, minRepresentative1->u, minRepresentative1->v, newEdges)) {
					WellSeparatedPair<IsImplicit>* second = pairs[j];
					WellSeparatedPairRepresentative* p = SelectedPair<Heap>(j, second, vertices, theta, w, edges);
					if (p == nullptr) {
						Q.remove(j);
						QContent[j] = nullptr;
					}
					else {
						Q.increaseKey(j, p->dist);
						QContent[j] = p;
					}
				}
			}
		}

		GapFillQueue1<Heap>(Q, QContent, index, pairs, vertices, theta, w, edgeCount, edges);
	}
	ENDLOG
		cout << " qc = " << rr << endl;
	delete splitTree;
	delete[] myEdges;*/
	return edgeCount;
}

#pragma endregion

#pragma endregion

//*************************//Our Methods\\****************************

#pragma region My Methods

#pragma region My Linear Space Method

//**********************Linear Space Method For Constructing Gap Greedy *******************
template<class Heap>
void getClosestGapedNeighbor(unsigned int j,
							const pointset& vertices,
							double theta,
							double w,
							std::vector<int>& QContent,
							Heap& Q,
							std::vector<std::pair<unsigned int, unsigned int>> &edges,
							bool* dirtyBits,
							double& smallesCleanCandidate)
{
	unsigned int N = vertices.size();
	std::vector< unsigned int> investNeighbors;
	double closestGapedNeighbor = std::numeric_limits<double>::infinity(),
			smallestNighborLargerThan_smallesCleanCandidate = std::numeric_limits<double>::infinity();
	
	int closestGapedNeighborIndex = -1, smallestNighborLargerThan_smallesCleanCandidateIndex = -1;
	double h = distance(vertices[j], vertices[QContent[j]]);

	for (unsigned int i = j + 1; i < N; i++)
	{
		double d = distance(vertices[j], vertices[i]);

		if (d >= h && d <= smallesCleanCandidate)
			investNeighbors.push_back(i);

		if (d > smallesCleanCandidate && d < smallestNighborLargerThan_smallesCleanCandidate)
		{
			smallestNighborLargerThan_smallesCleanCandidate = d;
			smallestNighborLargerThan_smallesCleanCandidateIndex = i;
		}
	}

	for (unsigned int i = 0; i < investNeighbors.size(); i++)
	{
		if (haveGapProperty(theta, w, vertices, j, investNeighbors[i], edges))
		{
			double d = distance(vertices[j], vertices[investNeighbors[i]]);
			if (closestGapedNeighbor > d)
			{
				closestGapedNeighbor = d;
				closestGapedNeighborIndex = investNeighbors[i];
			}
		}
	}

	if (std::isinf(closestGapedNeighbor))
	{
		QContent[j] = smallestNighborLargerThan_smallesCleanCandidateIndex;
		
		if (smallestNighborLargerThan_smallesCleanCandidateIndex == -1)
			Q.remove(j);
		else
			Q.increaseKey(j, smallestNighborLargerThan_smallesCleanCandidate);
	}
	else
	{
		if (closestGapedNeighbor < smallesCleanCandidate)
			smallesCleanCandidate = closestGapedNeighbor;

		dirtyBits[j] = false;
		Q.increaseKey(j, closestGapedNeighbor);
		QContent[j] = closestGapedNeighborIndex;
	}
}

template<class Heap, bool IsImplicit>
int MyGapGreedyLinspace(BaseTask * T) {
	STARTLOG
	const pointset &vertices = T->p;
	unsigned int N = vertices.size();
	unsigned int edgeCount = 0;
	double w = T->var, theta = T->t, smallestCleanCandidate = std::numeric_limits<double>::infinity();

	std::vector<std::pair<unsigned int, unsigned int>> edges, newEdge;
	std::vector<int> QContent(N - 1, -1);
	Heap Q(N - 1, -std::numeric_limits<double>::infinity());
	bool* dirtyBits = new bool[N - 1];
	int closestNeighborIndex; 
	for (unsigned int i = 0; i < N - 1; i++)
	{
		double closestDist = std::numeric_limits<double>::infinity();
		for (unsigned int j = i + 1; j < N; j++)
		{
			double dist = distance(vertices[i], vertices[j]);
			if (dist < closestDist)
			{
				closestDist = dist;
				closestNeighborIndex = j;
			}
		}
		dirtyBits[i] = false;
		QContent[i] = closestNeighborIndex;
		Q.insert(i, closestDist);
	}

	while (Q.getCount() > 0) {
		bool foundOne = false;
		std::pair<unsigned int, double> minQ;
		while (!foundOne)
		{
			minQ = Q.getMin();
			if (dirtyBits[minQ.first])
				getClosestGapedNeighbor(minQ.first, vertices, theta, w, QContent, Q, edges ,dirtyBits, smallestCleanCandidate);
			else
				foundOne = true;
			if (Q.getCount() == 0)
				goto end;
		}
		dirtyBits[minQ.first] = true;

		int minRepFrom = minQ.first;
		int minRepTo = QContent[minQ.first];
		double minRepDist = minQ.second;//distance(vertices[minRepFrom], vertices[minRepTo]);

		std::pair<unsigned int, unsigned int> e1(minRepFrom, minRepTo);
		edges.push_back(e1);

		newEdge.clear();
		newEdge.push_back(e1);

		if (T->create_edge(edge(minRepFrom, minRepTo))) return 0;
		edgeCount++;

		smallestCleanCandidate = std::numeric_limits<double>::infinity();
		for (int j = 0; j < N - 1; j++)
			if (QContent[j] != -1 && dirtyBits[j] == false)
			{
				if (!haveGapProperty(theta, w, vertices, j, QContent[j], newEdge))
					dirtyBits[j] = true;
				else
				{
					double d = distance(vertices[j], vertices[QContent[j]]);
					if (smallestCleanCandidate > d)
						smallestCleanCandidate = d;
				}
			}
	}
end:
	ENDLOG

	delete[] dirtyBits;
	return edgeCount;
}
#pragma endregion


#pragma region Alpha Gap Greedy Methods

#pragma region Structures And Common Functions

class NewCone
{
public:
	vertex apex;
	vertex bisector;
	double	angle;
	double  angle2;

	NewCone(const vertex& ap, const vertex& bs, double ang)
	{
		apex = ap;
		bisector = bs;
		angle = ang;
		angle2 = atan3(bisector, apex);
	}

	bool IsIn(const vertex* p)
	{
		double   ang1, ang;

		vertex p1(p->x, p->y);
		ang1 = atan3(p1, apex);

		ang = abs(ang1-angle2);
		
		if (ang <= angle)
			return true;
		else
			return false;
	}
};

//*****************************************
void	NewAddCone(std::vector<NewCone> *	Cones, const pointset &vertices, int apex, int bisector, double	angle)
{
	NewCone	cone(vertices[apex], vertices[bisector], angle);

	Cones[apex].push_back(cone);
}
//++++++++++++++++++++++++++++++++++++++++
bool NewIsInCone(std::vector<NewCone> *	Cones, const pointset &vertices, int apex, int p)
{
	for (unsigned int i = 0; i < Cones[apex].size(); i++)
	{
		if (Cones[apex][i].IsIn(&vertices[p]))
			return true;
	}
	return false;
}
//************************************************

#pragma endregion

//************************************************
bool  NewHaveGapProperty(double theta,
	double alpha,
	double w,
	const pointset &vertices,
	unsigned int v1,
	unsigned int v2,
	std::vector<std::pair<unsigned int, unsigned int>>& edges,
	int	candidatePolicy,
	double&	maxAngle,
	bool& sourceIsClose) {

	double  angle, distss, distee, d;
	maxAngle = -1;
	double ang1, ang2, ang1_1, ang2_1;

	ang1 = atan2(vertices[v1].y - vertices[v2].y, vertices[v1].x - vertices[v2].x);
	if (ang1 < 0)
		ang1_1 = ang1 + 2 * PI;
	else
		ang1_1 = ang1;

	for (int i = 0; i < edges.size(); i++) {

		ang2 = atan2(vertices[edges[i].first].y - vertices[edges[i].second].y, vertices[edges[i].first].x - vertices[edges[i].second].x);
		if (ang2 < 0)
			ang2_1 = ang2 + 2 * PI;
		else
			ang2_1 = ang2;

		if (abs(ang1) < PI / 2.0 &&  abs(ang2) < PI / 2.0)
			angle = abs(ang1 - ang2);
		else
			angle = abs(ang1_1 - ang2_1);

		if (angle <= theta - alpha)
		{
			d = w*distance(vertices[edges[i].first], vertices[edges[i].second]);
			distss = distance(vertices[edges[i].first], vertices[v1]);
			distee = distance(vertices[edges[i].second], vertices[v2]);

			if (distss <= d || distee <= d)
			{
				sourceIsClose = (distss <= d);

				if (candidatePolicy == FirstCandidate)
				{
					maxAngle = angle;
					return false;
				}
				else // Best Candidate Policy
				{
					if (maxAngle < angle)
					{
						maxAngle = angle;
					}
				}
			}
		}
	}
	if (candidatePolicy == BestCandidate && maxAngle != -1)
		return false;

	return true;
}

#pragma region Alpha Gap Greedy

template<class Heap, bool IsImplicit>
int NewAlphaGapGreedyFirstCandidate(BaseTask * T) {
	STARTLOG
		const pointset &vertices = T->p;
	std::vector<NewCone> *Cones = new std::vector<NewCone>[T->p.size()];
	unsigned int N = vertices.size();
	unsigned int edgeCount = 0;
	double w = T->var, theta = T->t, alpha = T->alpha, angle;

	int candidatePolicy = FirstCandidate;

	std::vector<EdgeInfo> edgeList;
	edgeList.reserve((N - 1) * N);
	for (unsigned int x = 0; x < N; x++) {
		for (unsigned int y = 0; y < N; y++) {
			if (x != y) {
				EdgeInfo e(x, y, distance(vertices[x], vertices[y]));
				edgeList.push_back(e);
			}
		}
	}

	std::sort(edgeList.begin(), edgeList.end());
	bool sourceIsClose;
	std::vector<std::pair<unsigned int, unsigned int>> edges;
	for (unsigned int i = 0; i < edgeList.size(); i++)
	{
		if (NewIsInCone(Cones, vertices, edgeList[i].x, edgeList[i].y))
			continue;

		if (NewHaveGapProperty(theta, alpha, w, vertices, edgeList[i].x, edgeList[i].y, edges, candidatePolicy, angle, sourceIsClose))
		{
			std::pair<unsigned int, unsigned int> e1(edgeList[i].x, edgeList[i].y);
			edges.push_back(e1);

			NewAddCone(Cones, vertices, edgeList[i].x, edgeList[i].y, theta - alpha);
			if (T->create_edge(edge(edgeList[i].x, edgeList[i].y))) return 0;
			edgeCount++;
		}
		else if (sourceIsClose)
		{
			NewAddCone(Cones, vertices, edgeList[i].x, edgeList[i].y, theta - angle);
		}
	}

	delete[] Cones;

	return edgeCount;
}

template<class Heap, bool IsImplicit>
int NewAlphaGapGreedyBestCandidate(BaseTask * T) {
	STARTLOG
		const pointset &vertices = T->p;
	std::vector<NewCone> *Cones = new std::vector<NewCone>[T->p.size()];
	unsigned int N = vertices.size();
	unsigned int edgeCount = 0;
	double w = T->var, theta = T->t, alpha = T->alpha, angle;

	int candidatePolicy = BestCandidate;

	std::vector<EdgeInfo> edgeList;
	edgeList.reserve((N - 1) * N);
	for (unsigned int x = 0; x < N; x++) {
		for (unsigned int y = 0; y < N; y++) {
			if (x != y) {
				EdgeInfo e(x, y, distance(vertices[x], vertices[y]));
				edgeList.push_back(e);
			}
		}
	}

	std::sort(edgeList.begin(), edgeList.end());
	bool sourceIsClose;
	std::vector<std::pair<unsigned int, unsigned int>> edges;
	for (unsigned int i = 0; i < edgeList.size(); i++)
	{
		if (NewIsInCone(Cones, vertices, edgeList[i].x, edgeList[i].y))
			continue;

		if (NewHaveGapProperty(theta, alpha, w, vertices, edgeList[i].x, edgeList[i].y, edges, candidatePolicy, angle, sourceIsClose))
		{
			std::pair<unsigned int, unsigned int> e1(edgeList[i].x, edgeList[i].y);
			edges.push_back(e1);

			NewAddCone(Cones, vertices, edgeList[i].x, edgeList[i].y, theta - alpha);
			if (T->create_edge(edge(edgeList[i].x, edgeList[i].y))) return 0;
			edgeCount++;
		}
		else if (sourceIsClose)
		{
			NewAddCone(Cones, vertices, edgeList[i].x, edgeList[i].y, theta - angle);
		}
	}

	delete[] Cones;

	return edgeCount;
}


#pragma endregion

#pragma region Linear Space Alpha Gap Greedy

template<class Heap>
void NewAlphaGapClosestGapedNeighbor(unsigned int j,
	const pointset& vertices,
	double theta,
	double alpha,
	double w,
	std::vector<int>& QContent,
	Heap& Q,
	Heap& myHeap,
	std::vector<std::pair<unsigned int, unsigned int>> &edges,
	std::vector<NewCone>* Cones,
	bool* dirtyBits,
	double& smallesCleanCandidate,
	int candidatePolicy)
{

	unsigned int N = vertices.size();
	double closestGapedNeighbor = std::numeric_limits<double>::infinity(),
		smallestNighborLargerThan_smallesCleanCandidate = std::numeric_limits<double>::infinity();

	int closestGapedNeighborIndex = -1, smallestNighborLargerThan_smallesCleanCandidateIndex = -1;
	double h = distance(vertices[j], vertices[QContent[j]]), largestInvestigatedNeighbor = 0, angle;
	myHeap.clear();

	for (unsigned int i = 0; i < N; i++)
	{
		if (i == j)
			continue;

		double d = distance(vertices[j], vertices[i]);

		if (d >= h && d <= smallesCleanCandidate)
		{
			myHeap.insert(i, d);

			if (d > largestInvestigatedNeighbor)
				largestInvestigatedNeighbor = d;
		}
		if (d > smallesCleanCandidate && d < smallestNighborLargerThan_smallesCleanCandidate)
		{
			smallestNighborLargerThan_smallesCleanCandidate = d;
			smallestNighborLargerThan_smallesCleanCandidateIndex = i;
		}
	}
	bool sourceIsClose;
	while (myHeap.getCount() > 0)
	{
		std::pair<unsigned int, double> pair = myHeap.getMin();
		myHeap.extractMin();

		if (NewIsInCone(Cones, vertices, j, pair.first))
			continue;

		//int whichVertex = WhichVertexRuinGapProperty(theta, alpha, w, vertices, j, pair.first, edges, candidatePolicy, angle);
		if (NewHaveGapProperty(theta, alpha, w, vertices, j, pair.first, edges, candidatePolicy, angle, sourceIsClose))
		{
			closestGapedNeighbor = pair.second;
			closestGapedNeighborIndex = pair.first;
			break;
		}
		else
		{
			if (sourceIsClose)
				NewAddCone(Cones, vertices, j, pair.first, theta - angle);
		}
	}

	if (std::isinf(closestGapedNeighbor))
	{
		QContent[j] = smallestNighborLargerThan_smallesCleanCandidateIndex;

		if (smallestNighborLargerThan_smallesCleanCandidateIndex == -1)
			Q.remove(j);
		else
			Q.increaseKey(j, smallestNighborLargerThan_smallesCleanCandidate);
	}
	else
	{
		if (closestGapedNeighbor < smallesCleanCandidate)
			smallesCleanCandidate = closestGapedNeighbor;

		dirtyBits[j] = false;
		Q.increaseKey(j, closestGapedNeighbor);
		QContent[j] = closestGapedNeighborIndex;
	}
}

template<class Heap, bool IsImplicit>
int NewLinearSpaceAlphaGapGreedyFirstCandidate(BaseTask * T) {
	STARTLOG
		const pointset &vertices = T->p;
	unsigned int N = vertices.size();
	unsigned int edgeCount = 0;
	double w = T->var, theta = T->t, alpha = T->alpha, angle, smallestCleanCandidate = std::numeric_limits<double>::infinity();

	int candidatePolicy = FirstCandidate;

	std::vector<NewCone> *Cones = new std::vector<NewCone>[N];
	std::vector<std::pair<unsigned int, unsigned int>> edges, newEdge;
	std::vector<int> QContent(N, -1);
	Heap Q(N, -std::numeric_limits<double>::infinity());
	Heap myHeap(N, 0);

	bool* dirtyBits = new bool[N];
	int closestNeighborIndex;
	for (unsigned int i = 0; i < N; i++)
	{
		double closestDist = std::numeric_limits<double>::infinity();
		for (unsigned int j = 0; j < N; j++)
		{
			if (i != j)
			{
				double dist = distance(vertices[i], vertices[j]);
				if (dist < closestDist)
				{
					closestDist = dist;
					closestNeighborIndex = j;
				}
			}
		}
		dirtyBits[i] = false;
		QContent[i] = closestNeighborIndex;
		Q.insert(i, closestDist);
	}

	while (Q.getCount() > 0) {
		bool foundOne = false;
		std::pair<unsigned int, double> minQ;
		while (!foundOne)
		{
			minQ = Q.getMin();
			if (dirtyBits[minQ.first])
				NewAlphaGapClosestGapedNeighbor(minQ.first, vertices, theta, alpha, w, QContent, Q, myHeap, edges, Cones, dirtyBits, smallestCleanCandidate, candidatePolicy);
			else
				foundOne = true;
			if (Q.getCount() == 0)
				goto end;
		}
		dirtyBits[minQ.first] = true;

		int minRepFrom = minQ.first;
		int minRepTo = QContent[minQ.first];
		double minRepDist = minQ.second;

		std::pair<unsigned int, unsigned int> e1(minRepFrom, minRepTo);
		edges.push_back(e1);

		NewAddCone(Cones, vertices, minRepFrom, minRepTo, theta - alpha);

		newEdge.clear();
		newEdge.push_back(e1);

		if (T->create_edge(edge(minRepFrom, minRepTo))) return 0;
		edgeCount++;

		smallestCleanCandidate = std::numeric_limits<double>::infinity();
		for (int j = 0; j < N; j++)
			if (QContent[j] != -1 && dirtyBits[j] == false)
			{
				bool sourceIsClose;
				if (!NewHaveGapProperty(theta, alpha, w, vertices, j, QContent[j], newEdge, candidatePolicy, angle, sourceIsClose))
					dirtyBits[j] = true;
				else
				{
					double d = distance(vertices[j], vertices[QContent[j]]);
					if (smallestCleanCandidate > d)
						smallestCleanCandidate = d;
				}
			}
	}
end:
	ENDLOG

		delete[] dirtyBits;
	delete[] Cones;
	return edgeCount;
}

template<class Heap, bool IsImplicit>
int NewLinearSpaceAlphaGapGreedyBestCandidate(BaseTask * T) {
	STARTLOG
		const pointset &vertices = T->p;
	unsigned int N = vertices.size();
	unsigned int edgeCount = 0;
	double w = T->var, theta = T->t, alpha = T->alpha, angle, smallestCleanCandidate = std::numeric_limits<double>::infinity();
	// to have a gap greedy Spanner, alpha must be 0
	//alpha = 0;
	int candidatePolicy = BestCandidate;

	std::vector<NewCone> *Cones = new std::vector<NewCone>[N];
	std::vector<std::pair<unsigned int, unsigned int>> edges, newEdge;
	std::vector<int> QContent(N, -1);
	Heap Q(N, -std::numeric_limits<double>::infinity());
	Heap myHeap(N, 0);

	bool* dirtyBits = new bool[N];
	int closestNeighborIndex;
	for (unsigned int i = 0; i < N; i++)
	{
		double closestDist = std::numeric_limits<double>::infinity();
		for (unsigned int j = 0; j < N; j++)
		{
			if (i != j)
			{
				double dist = distance(vertices[i], vertices[j]);
				if (dist < closestDist)
				{
					closestDist = dist;
					closestNeighborIndex = j;
				}
			}
		}
		dirtyBits[i] = false;
		QContent[i] = closestNeighborIndex;
		Q.insert(i, closestDist);
	}

	while (Q.getCount() > 0) {
		bool foundOne = false;
		std::pair<unsigned int, double> minQ;
		while (!foundOne)
		{
			minQ = Q.getMin();
			if (dirtyBits[minQ.first])
				NewAlphaGapClosestGapedNeighbor(minQ.first, vertices, theta, alpha, w, QContent, Q, myHeap, edges, Cones, dirtyBits, smallestCleanCandidate, candidatePolicy);
			else
				foundOne = true;
			if (Q.getCount() == 0)
				goto end;
		}
		dirtyBits[minQ.first] = true;

		int minRepFrom = minQ.first;
		int minRepTo = QContent[minQ.first];
		double minRepDist = minQ.second;

		std::pair<unsigned int, unsigned int> e1(minRepFrom, minRepTo);
		edges.push_back(e1);

		NewAddCone(Cones, vertices, minRepFrom, minRepTo, theta - alpha);

		newEdge.clear();
		newEdge.push_back(e1);

		if (T->create_edge(edge(minRepFrom, minRepTo))) return 0;
		edgeCount++;

		smallestCleanCandidate = std::numeric_limits<double>::infinity();
		for (int j = 0; j < N; j++)
			if (QContent[j] != -1 && dirtyBits[j] == false)
			{
				bool sourceIsClose;
				if (!NewHaveGapProperty(theta, alpha, w, vertices, j, QContent[j], newEdge, candidatePolicy, angle, sourceIsClose))
					dirtyBits[j] = true;
				else
				{
					double d = distance(vertices[j], vertices[QContent[j]]);
					if (smallestCleanCandidate > d)
						smallestCleanCandidate = d;
				}
			}
	}
end:
	ENDLOG

		delete[] dirtyBits;
	delete[] Cones;
	return edgeCount;
}
#pragma endregion

#pragma endregion

#pragma region MyNewMethod

#pragma region Structures And Common Functions



class CircleSector
{
public:
	static int count;
	vertex apex;
	double 	startAngle, endAngle;
	int id;

	CircleSector(const vertex& ap, double   startAngle, double   endAngle)
	{
		apex = ap;
		this->startAngle = startAngle;
		this->endAngle = endAngle;
		id = count;
		count++;
	}

	bool IsIn(const vertex* p)
	{
		double   ang = atan3(*p, apex);
		/*ang = atan2((double)(p->y - apex.y), (double)(p->x - apex.x));
		if (ang < 0)
			ang += 2.0 * PI;*/

		if (startAngle < endAngle)
		{
			return (ang >= startAngle && ang <= endAngle);
		}
		else
		{
			return (ang >= startAngle || ang <= endAngle);
		}
	}
};

class CircleSector_list_node {
public:
	CircleSector *data;
	CircleSector_list_node *next;
};

class CircleSector_list {
public:
	CircleSector_list_node *first;
	int size;
	CircleSector_list() {
		first = nullptr;
		size = 0;
	}
	~CircleSector_list() {
		CircleSector_list_node *t, *f = first;
		while (f != nullptr) {
			t = f;
			f = f->next;
			//if (t->data != nullptr)
				delete t->data;
			delete t;
		}
	}
	void add(CircleSector *d) {
		CircleSector_list_node *p = first;
		first = new CircleSector_list_node;
		first->next = p;
		first->data = d;
		size++;
	}

	void remove(int id) {
		CircleSector_list_node *p = first, *q = nullptr;
		while (p != nullptr && p->data->id != id)
		{
			q = p;
			p = p->next;
		}

		if (p != nullptr)
		{
			if (p == first)
				first = first->next;
			else
				q->next = p->next;

			delete p->data;
			delete p;
			size--;
		}
	}
};


int CircleSector::count = 1;

//*****************************************
void	AddCircleSector(CircleSector_list *	Cones, const pointset &vertices, int apex, double  startAngle, double 	endAngle)
{
	CircleSector	*cone = new CircleSector(vertices[apex], startAngle, endAngle);
	double w;
	if (endAngle > startAngle)
		w = (endAngle - startAngle) * (180.0 / PI);
	else
		w = (2.0 * PI - startAngle + endAngle) * (180.0 / PI);

	//if (apex == 8)
	//cout << apex << " New (sa:" << startAngle * (180.0 / PI) << ", ea:" << endAngle * (180.0 / PI) << ")\t "<<w << endl;
	Cones[apex].add(cone);
}
///*****************************************
void	RemoveCircleSector(CircleSector_list *	Cones,int apex, int id)
{
	Cones[apex].remove(id); 
}
//***************************************
bool IsInCircleSector(CircleSector_list *	Cones, const pointset &vertices, int apex, int p)
{
	CircleSector_list_node* v_item = Cones[apex].first;
	CircleSector*	cone;
	int c = 0;
	while (v_item != nullptr)
	{
		c++;
		cone = v_item->data;
		v_item = v_item->next;
		if (cone->IsIn(&vertices[p]))
		{
			return true;
		}
	}
	return false;
}

CircleSector* GetConePointIn(CircleSector_list *	Cones, const pointset &vertices, int apex, double   ang)
{
	CircleSector_list_node* v_item = Cones[apex].first;
	CircleSector*	cone, *findCone = nullptr;
	while (v_item != nullptr)
	{
		cone = v_item->data;
		v_item = v_item->next;
		if (cone->startAngle < cone->endAngle)
		{
			if (ang >= cone->startAngle && ang <= cone->endAngle)
				return cone;
		}
		else
		{
			if (ang >= cone->startAngle || ang <= cone->endAngle)
				return cone;
		}
	}
	return nullptr;
}

#pragma endregion

/*
vertex	rotatePoint(const pointset &vertices, int apex, double slope, double theta, double sign)
{
	vertex point;
	double l = 5.0;

	if (slope == std::numeric_limits<float>
		::max()) {
		point.x = vertices[apex].x;
		point.y = vertices[apex].y + l*sign;
	}
	else {
		float dx = (l / sqrt(1 + (slope * slope)));
		float dy = slope * dx;
		point.x = vertices[apex].x + sign*dx;
		point.y = vertices[apex].y + sign*dy;
	}
	double angle = (theta) * (PI / 180.0); // Convert to radians

	double rotatedX = cos(angle) * (point.x - vertices[apex].x) - sin(angle) * (point.y - vertices[apex].y) + vertices[apex].x;

	double rotatedY = sin(angle) * (point.x - vertices[apex].x) + cos(angle) * (point.y - vertices[apex].y) + vertices[apex].y;
	return point;
}*/


bool NewGetConePointIn(CircleSector_list *	Cones, double ang)
{
	CircleSector_list_node* v_item = Cones[0].first;
	CircleSector*	cone, *findCone = nullptr;
	while (v_item != nullptr)
	{
		cone = v_item->data;
		v_item = v_item->next;
		if (cone->startAngle < cone->endAngle)
		{
			if (ang >= cone->startAngle && ang <= cone->endAngle)
				return true;
		}
		else
		{
			if (ang >= cone->startAngle || ang <= cone->endAngle)
				return cone;
		}
	}
	return false;
}


double toAngle(double rad)
{
	return rad *(180.0 / PI);
}
int hh = 1;
int	MergeCircleSector(CircleSector_list*	Cones, const pointset &vertices, int apex, double   startAngle, double   endAngle,  CircleSector* fc,  CircleSector* sc)
{
	double sa = startAngle, ea = endAngle;
	
	if (fc == nullptr)
	{
		sc->startAngle = sa;
		return 1;
	}
	else if (sc == nullptr)
	{
		fc->endAngle = ea;
		return 2;
	}
	else
	{
		sa = fc->startAngle;
		ea = sc->endAngle;

		if (fc->id != sc->id)
		{
			RemoveCircleSector(Cones, apex, fc->id);
			RemoveCircleSector(Cones, apex, sc->id);
			AddCircleSector(Cones, vertices, apex, sa, ea);
			return 3;
		}
		else
		{
			CircleSector_list *ca = new CircleSector_list[1];
			CircleSector	*cone = new CircleSector(vertices[apex], startAngle, endAngle);
			ca[0].add(cone);

			bool a = NewGetConePointIn(ca, fc->startAngle);
			bool b = NewGetConePointIn(ca, fc->endAngle);

			delete[] ca;

			if (!(fc->endAngle == endAngle && fc->startAngle == startAngle) && a && b )
			{
				fc->startAngle = 0.0;
				fc->endAngle = 2.0 * PI;
				return 4;
			}
			else
			{
				return 5;
			}
		}
	}
}

int	AddOrMergeCones(CircleSector_list *	Cones, const pointset &vertices, int apex, double   ang, double   theta)
{
	double  	sa = ang - theta,
			ea = ang + theta;
	if (sa < 0)
		sa += 2.0 * M_PI;

	if (ea > 2.0 * M_PI)
		ea -= 2.0 * M_PI;

	CircleSector* fc = GetConePointIn(Cones, vertices, apex, sa);
	CircleSector* sc = GetConePointIn(Cones, vertices, apex, ea);

	if (fc == nullptr && sc == nullptr)
	{
		AddCircleSector(Cones, vertices, apex, sa, ea);
		return 3;
	}
	else
	{
		return MergeCircleSector(Cones, vertices, apex, sa, ea, fc, sc);
	}
}

//method1
template<class Heap, bool IsImplicit>
WellSeparatedPairRepresentative* MNM_SelectedPair(unsigned int pairIndex,
	WellSeparatedPair<IsImplicit>* pair,
	const pointset &vertices,
	double   theta,
	double   w,
	CircleSector_list *Cones) {
	
	unsigned int N = vertices.size();
	rr++;

	unsigned int shortestPairInedx1, shortestPairInedx2;
	double   bestDist = std::numeric_limits<double>::infinity();
	unsigned int 	index1, index2;

	for (typename SplitTree<IsImplicit>::Iterator it1 = pair->first->begin(); it1 != pair->first->end(); ++it1) {

		index1 = *it1;
		for (typename SplitTree<IsImplicit>::Iterator it2 = pair->second->begin(); it2 != pair->second->end(); ++it2) {
			index2 = *it2;

			if (!IsInCircleSector(Cones, vertices, index1, index2) && !IsInCircleSector(Cones, vertices, index2, index1))
			{
				double   dist = NotSquaredDistance(vertices[index1], vertices[index2]);
				if (dist < bestDist)
				{
					bestDist = dist;
					shortestPairInedx1 = index1;
					shortestPairInedx2 = index2;
				}
			}
		}
	}
	if (std::isinf(bestDist))
		return nullptr;
	else
		return new WellSeparatedPairRepresentative(shortestPairInedx1, shortestPairInedx2, sqrt(bestDist), pairIndex);
}

template<class Heap, bool IsImplicit>
void MNM_GapFillQueue1(Heap& Q,
	std::vector<WellSeparatedPairRepresentative*>& QContent,
	unsigned int& i,
	std::vector<WellSeparatedPair<IsImplicit>*>& pairs,
	const pointset &vertices,
	double   theta,
	double   w,
	int edgeCount,
	std::vector<std::pair<unsigned int, unsigned int> > &edges,
	CircleSector_list *Cones,
	HeapIndex& heapIndex) {

	while (i < pairs.size() && (Q.getCount() == 0 || pairs[i]->minlength() <= Q.getMin().second)) {
		WellSeparatedPairRepresentative* p = MNM_SelectedPair<Heap, IsImplicit>(i, pairs[i], vertices, theta, w, Cones);
		//if (i % 10000000 == 0)
		//std::cout << i << "/" << pairs.size() << " " << Q.getCount() << " " << edgeCount << endl;
		if (p != nullptr) {
			Q.insert(i, p->dist);
			QContent[i] = p;
			heapIndex.Add(i);
		}
		i++;
	}
}


bool	printError = false;

inline void PrintError(char	message[]) {

	if (printError)
		cout << message;
}


template<class Heap, bool IsImplicit>
int MNM_Method(BaseTask * T) {
	STARTLOG
		rr = 0;
	const pointset &vertices = T->p;
	unsigned int N = vertices.size();
	unsigned int edgeCount = 0;
	
	CircleSector_list *Cones = new CircleSector_list[T->p.size()];
	
	double   w = T->var;
	double   theta = T->t;

	double   s = maximum((double) (4.0 / sin(theta / 2.0)), (double)(2.0 / w)); //4.0 * t / (t - 1.0);

	//PrintError("Before CreateSplit Tree..\n");

	SplitTree<IsImplicit>* splitTree = SplitTree<IsImplicit>::createSplitTree(vertices);
	//cout << "CreateSplit Tree Done..\n";
	std::vector<WellSeparatedPair<IsImplicit>*> pairs;
	
	//cout << "Before Computing WSPD (s:" << s << ",w:" << w << ",theta:" << theta << ")..\n";
	splitTree->ComputeWspdPairs(s, pairs);
	
	//cout << "Computing WSPD Done (Count:"<< pairs.size()<<")\n";

	sort(pairs.begin(), pairs.end(), WellSeparatedPairPointerComparer<IsImplicit>());

	//cout << "Sort Done\n";

	std::vector<WellSeparatedPairRepresentative*> QContent(pairs.size(), nullptr);
	Heap Q(pairs.size(), -std::numeric_limits<double>::infinity());
	
	std::vector<std::pair<unsigned int, unsigned int>> edges;

	HeapIndex heapIndex(pairs.size());

	unsigned int index = 0, i=0;
	MNM_GapFillQueue1<Heap>(Q, QContent, index, pairs, vertices, theta, w, edgeCount, edges, Cones, heapIndex);
	
	while (Q.getCount() > 0) {
		
		//cout << "next Pair.." << endl;
		
		std::pair<unsigned int, double > minQ = Q.getMin();
		
		Q.extractMin();
		WellSeparatedPairRepresentative* minRepresentative = QContent[minQ.first];
		
		QContent[minQ.first] = nullptr;

		heapIndex.Remove(minQ.first);

		std::pair<unsigned int, unsigned int> e1(minRepresentative->u, minRepresentative->v);
		edges.push_back(e1);

		std::pair<unsigned int, unsigned int> e2(minRepresentative->v, minRepresentative->u);
		edges.push_back(e2);
		/*
		if (minRepresentative->v > minRepresentative->u)
		{
			myfile << "-" << minRepresentative->v << "&" << minRepresentative->u << "-" << "(" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ") ---> (" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ")" << endl;
			myfile << "-" << minRepresentative->u << "&" << minRepresentative->v << "-" << "(" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ") ---> (" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ")" << endl;
		}
		else
		{
			myfile << "-" << minRepresentative->u << "&" << minRepresentative->v << "-" << "(" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ") ---> (" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ")" << endl;
			myfile << "-" << minRepresentative->v << "&" << minRepresentative->u << "-" << "(" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ") ---> (" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ")" << endl;
		}*/

		
		if (T->create_edge(edge(minRepresentative->u, minRepresentative->v))) return 0;
		if (T->create_edge(edge(minRepresentative->v, minRepresentative->u))) return 0; 
		edgeCount+=2;
		
		//AddCircleSector(Cones, vertices, 0, 26.0/18.0*PI, 28.0 / 18.0 * PI);
		//IsInCircleSector(Cones, vertices, 0, 1);

		double  wd = w * minQ.second,
			ang = atan2((double)(vertices[minRepresentative->u].y - vertices[minRepresentative->v].y), 
				(double)(vertices[minRepresentative->u].x - vertices[minRepresentative->v].x));

		if (ang < 0)
			ang += 2.0 * PI;
		wd = wd * wd;

		for (int i = 0; i < vertices.size(); i++)
		{
			double   dist = NotSquaredDistance(vertices[minRepresentative->v], vertices[i]);
			if (dist <= wd)
			{
				AddOrMergeCones(Cones, vertices, i, ang, theta);
				//cout << "Cones " << i << ":" << Cones[i].size() << endl;
			}
		}

		ang += PI;
		if (ang > 2.0 * PI)
			ang -= 2.0 * PI;

		for (int i = 0; i < vertices.size(); i++)
		{
			double   dist = NotSquaredDistance(vertices[minRepresentative->u], vertices[i]);
			if (dist <= wd)
			{
				AddOrMergeCones(Cones, vertices, i, ang, theta);
				//cout << "Cones " << i << ":" << Cones[i].size()<<endl;
			}
		}
		
		WellSeparatedPair<IsImplicit>* first = pairs[minRepresentative->pairIndex];
		
		//for (int j = 0; j < index; j++) {
		HeapIndexNode* node = heapIndex.first->next;
		while (node != nullptr)
		 {
			//if (QContent[j] != nullptr) {
			int j = node->index;
			node = node->next;
			WellSeparatedPair<IsImplicit>* second = pairs[j];
			double   distffsf = minDistance(first->first->box, second->first->box);
			double   distffss = minDistance(first->first->box, second->second->box);
			double   distfssf = minDistance(first->second->box, second->first->box);
			double   distfsss = minDistance(first->second->box, second->second->box);
			double   lengths = maxDistance(second->first->box, second->second->box) * w + 1e-6;
			if (distffsf <= lengths || distfsss <= lengths || distffss <= lengths || distfssf <= lengths) {
				WellSeparatedPairRepresentative* p = MNM_SelectedPair<Heap, IsImplicit>(j, pairs[j], vertices, theta, w,Cones);
				if (p == nullptr) {
					Q.remove(j);
					QContent[j] = nullptr;
					heapIndex.Remove(j);
				}
				else {
					Q.increaseKey(j, p->dist);
					QContent[j] = p;
				}
			}
		   //}
		}

		MNM_GapFillQueue1<Heap>(Q, QContent, index, pairs, vertices, theta, w, edgeCount, edges, Cones, heapIndex);

		delete minRepresentative;
	}
	ENDLOG
	
	char filename[80];
	sprintf(filename, "MNM(%d-%f-%f).txt", N, theta, w);
	std::ofstream myfile(filename);
	
	
	
	PROCESS_MEMORY_COUNTERS_EX pmc;
	GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
	SIZE_T physMemUsedByMe = pmc.WorkingSetSize;
	
	myfile << pmc.WorkingSetSize;
	myfile << "\t" << pairs.size() << endl;
	myfile.close();
	
	delete[] Cones;
	delete splitTree;

	for (size_t i = 0; i < pairs.size(); i++)
	{
		delete pairs[i];
	}
	pairs.clear();
	edges.clear();
	QContent.clear();
	//heapIndex.Dispose();

	//GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));

	//SIZE_T physMemUsedByMe = pmc.WorkingSetSize;
//	size_t vv = pmc.WorkingSetSize;
	
	return edgeCount;
}

//*********************************************************************
template<class Heap, bool IsImplicit>
WellSeparatedPairRepresentative* GapClosestPair2(unsigned int pairNo,
	WellSeparatedPair<IsImplicit>* pair,
	const pointset &vertices,
	double t,
	double w,
	bool** discounters,
	CircleSector_list *Cones
) {
	int shortestPairInedx1, shortestPairInedx2;
	double bestDist = std::numeric_limits<double>::infinity();
	unsigned int j = 0;
	for (typename SplitTree<IsImplicit>::Iterator it = pair->first->begin(); it != pair->first->end(); ++it, j++) {
		unsigned int indexj = *it;
		if (discounters[pairNo][j] == false) {
			bool discounted = true;

			for (typename SplitTree<IsImplicit>::Iterator it1 = pair->second->begin(); it1 != pair->second->end(); ++it1) {
				unsigned int indexk = *it1;
				if (!IsInCircleSector(Cones, vertices, indexj, indexk) && !IsInCircleSector(Cones, vertices, indexk, indexj)) {
					double   dist = NotSquaredDistance(vertices[indexj], vertices[indexk]);//distance(vertices[indexj], vertices[indexk]);
					discounted = false;
					if (dist < bestDist)
					{
						bestDist = dist;
						shortestPairInedx1 = indexj;
						shortestPairInedx2 = indexk;
					}
				}
			}
			discounters[pairNo][j] = discounted;
		}
	}

	if (std::isinf(bestDist))
		return nullptr;
	else
		return new WellSeparatedPairRepresentative(shortestPairInedx1, shortestPairInedx2, sqrt(bestDist), pairNo);
}

template<class Heap, bool IsImplicit>
void GapUpdatePair2(unsigned int j,
	Heap& Q,
	std::vector<WellSeparatedPairRepresentative*>& QContent,
	WellSeparatedPair<IsImplicit>* pair,
	const pointset &vertices,
	std::vector< std::pair<unsigned int, double> >* myEdges,
	double t,
	double w,
	bool* dirtyBits,
	bool** discounters,
	AngleStorer* certificates,
	CircleSector_list *Cones
) {
	WellSeparatedPairRepresentative* p = GapClosestPair2<Heap, IsImplicit>(j, pair, vertices, t, w, discounters, Cones);
	if (p == nullptr) {
		Q.remove(j);
		QContent[j] = nullptr;
	}
	else {
		Q.increaseKey(j, p->dist);
		QContent[j] = p;
	}
	dirtyBits[j] = false;
}

template<class Heap, bool IsImplicit>
void GapReinsertPair2(unsigned int j,
	Heap& Q,
	std::vector<WellSeparatedPairRepresentative*>& QContent,
	WellSeparatedPair<IsImplicit>* pair,
	const pointset &vertices,
	std::vector< std::pair<unsigned int, double> >* myEdges,
	double t,
	double w,
	bool* dirtyBits,
	bool** discounters,
	AngleStorer* certificates,
	CircleSector_list *Cones
) {
	WellSeparatedPairRepresentative* p = GapClosestPair2<Heap, IsImplicit>(j, pair, vertices, t, w, discounters, Cones);
	if (p != nullptr) {
		Q.insert(j, p->dist);
		QContent[j] = p;
	}
	/*else
		QContent[j] = nullptr;*/
	dirtyBits[j] = false;
}

template<class Heap, bool IsImplicit>
void GapFillQueue3(Heap& Q,
	std::vector<WellSeparatedPairRepresentative*>& QContent,
	std::vector<WellSeparatedPair<IsImplicit>*>& pairs,
	const pointset &vertices,
	std::vector< std::pair<unsigned int, double> >* myEdges,
	double t,
	double w,
	int edgeCount,
	bool* discountedPairs,
	bool* dirtyBits,
	bool** discounters,
	AngleStorer* certificates,
	CircleSector_list *Cones
) {
	for (unsigned int j = 0; j < pairs.size(); j++) {
		if (!discountedPairs[j]) {
			WellSeparatedPairRepresentative* p = GapClosestPair2<Heap, IsImplicit>(j, pairs[j], vertices, t, w, discounters, Cones);
			if (p == nullptr)
				QContent[j] = nullptr;
			else {
				Q.insert(j, p->dist);
				QContent[j] = p;
			}
		}
	}
}

template<class Heap>
SIZE_T ComputeGapSpannerBuckets(BaseTask * T,
	unsigned int& edgeCount,
	std::vector<std::pair<unsigned int, double> >* myEdges,
	double t,
	std::vector<int>* buckets,
	double estimatedLongestEdgeLength,
	double estimatedLongestTPath,
	int* bucketLookupX,
	int* bucketLookupY,
	int bucketeers,
	CircleSector_list *Cones,
	std::vector<int>* bucketsLarge
) {
	const pointset &vertices = T->p;
	int N = vertices.size();
	double   w = T->var;
	double   theta = T->t;
	int k = 0;
	unsigned int edgeCount2 = 0;
	for (int bucketLargeCoord = 0; bucketLargeCoord < bucketeers * bucketeers; bucketLargeCoord++) {
		int bucketLargeSize = bucketsLarge[bucketLargeCoord].size();
		int bucketSize = buckets[bucketLargeCoord].size();
		for (int i = 0; i < bucketSize; i++) {
			int x = buckets[bucketLargeCoord][i];
			for (int j = 0; j < bucketLargeSize; j++) {
				int y = bucketsLarge[bucketLargeCoord][j];
				if (x < y && distance(vertices[x], vertices[y]) <= estimatedLongestEdgeLength) {
					edgeCount2++;
				}
					//if (distance(vertices[x], vertices[y]) > estimatedLongestEdgeLength)
						//cout << k++ <<":"<< distance(vertices[x], vertices[y]) << endl;
			}
		}
	}

	std::vector<EdgeInfo> edgeList;
	edgeList.reserve(edgeCount2);

	for (int bucketLargeCoord = 0; bucketLargeCoord < bucketeers * bucketeers; bucketLargeCoord++) {
		int bucketLargeSize = bucketsLarge[bucketLargeCoord].size();
		int bucketSize = buckets[bucketLargeCoord].size();
		for (int i = 0; i < bucketSize; i++) {
			int x = buckets[bucketLargeCoord][i];
			for (int j = 0; j < bucketLargeSize; j++) {
				int y = bucketsLarge[bucketLargeCoord][j];
				if (x < y && distance(vertices[x], vertices[y]) <= estimatedLongestEdgeLength) {
					EdgeInfo e(x, y, distance(vertices[x], vertices[y]));
					edgeList.push_back(e);
				}
			}
		}
	}

	/*
	char filename[80];
	sprintf(filename, "NMNM(%d).txt", N);
	remove(filename);
	std::ofstream myfile(filename, std::ofstream::out | std::ofstream::app);
	*/
	std::sort(edgeList.begin(), edgeList.end());

	for (unsigned int k = 0; k < edgeList.size(); k++) {
		double edgeDist = w*edgeList[k].distance;
		double  angyx = atan3(vertices[edgeList[k].x], vertices[edgeList[k].y]);
		double  angxy = atan3(vertices[edgeList[k].y], vertices[edgeList[k].x]);

		if (!IsInCircleSector(Cones, vertices, edgeList[k].x, edgeList[k].y) && !IsInCircleSector(Cones, vertices, edgeList[k].y, edgeList[k].x)) {
		//if (!certificates[edgeList[k].x].Covers(angxy, angxy) && !certificates[edgeList[k].y].Covers(angyx, angyx)) {
			if (T->create_edge(edgeList[k])) return 0;
			if (T->create_edge(edgeList[k].Invert())) return 0;
			myEdges[edgeList[k].x].push_back(std::pair<unsigned int, double>(edgeList[k].y, edgeList[k].distance));
			myEdges[edgeList[k].y].push_back(std::pair<unsigned int, double>(edgeList[k].x, edgeList[k].distance));
			edgeCount+=2;


			/*if (edgeList[k].y > edgeList[k].x)
			{
				myfile << "-" << edgeList[k].y << "&" << edgeList[k].x << "-" << "(" << vertices[edgeList[k].y].x << ", " << vertices[edgeList[k].y].y << ") ---> (" << vertices[edgeList[k].x].x << ", " << vertices[edgeList[k].x].y << ")" << endl;
				myfile << "-" << edgeList[k].x << "&" << edgeList[k].y << "-" << "(" << vertices[edgeList[k].x].x << ", " << vertices[edgeList[k].x].y << ") ---> (" << vertices[edgeList[k].y].x << ", " << vertices[edgeList[k].y].y << ")" << endl;
			}
			else
			{
				myfile << "-" << edgeList[k].x << "&" << edgeList[k].y << "-" << "(" << vertices[edgeList[k].x].x << ", " << vertices[edgeList[k].x].y << ") ---> (" << vertices[edgeList[k].y].x << ", " << vertices[edgeList[k].y].y << ")" << endl;
				myfile << "-" << edgeList[k].y << "&" << edgeList[k].x << "-" << "(" << vertices[edgeList[k].y].x << ", " << vertices[edgeList[k].y].y << ") ---> (" << vertices[edgeList[k].x].x << ", " << vertices[edgeList[k].x].y << ")" << endl;
			}*/

			

			// Just added edge (c, d)

			int bucketIndexC = bucketLookupX[edgeList[k].x] + bucketeers * bucketLookupY[edgeList[k].x];
			int bucketIndexD = bucketLookupX[edgeList[k].y] + bucketeers * bucketLookupY[edgeList[k].y];

			int CBucketLargeSize = bucketsLarge[bucketIndexC].size();
			int DBucketLargeSize = bucketsLarge[bucketIndexD].size();

			for (unsigned int a = 0; a < CBucketLargeSize; a++) {
				int Ai = bucketsLarge[bucketIndexC][a];
				if (distance(vertices[Ai], vertices[edgeList[k].x]) <= edgeDist) {
					AddOrMergeCones(Cones, vertices, Ai, angxy, theta);
				}
			}

			for (unsigned int b = 0; b < DBucketLargeSize; b++) {
				int Bi = bucketsLarge[bucketIndexD][b];
				if (distance(vertices[Bi], vertices[edgeList[k].y]) <= edgeDist) {
					AddOrMergeCones(Cones, vertices, Bi, angyx, theta);
				}
			}
		}
	}

	PROCESS_MEMORY_COUNTERS_EX pmc;
	GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
	SIZE_T physMemUsedByMe = pmc.WorkingSetSize;

	return pmc.WorkingSetSize;

	//delete[] bucketsLarge;
	//myfile.close();
	//return edgeCount;
}

template<class Heap, bool IsImplicit>
int GapSpannerBucketsFixed(BaseTask * T) {
	//TODO:
	// Make ESA' local
	STARTLOG
		const pointset &vertices = T->p;
	unsigned int N = vertices.size();
	std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	unsigned int edgeCount = 0;

	// Bucketing stuff
	//double w_ratio = T->t;
	double theta = T->t;
	double w = T->var;
	double t = 1.0 / (cos(theta) - sin(theta) - 2 * w);
	double estimatedLongestEdgeLengthUnit = (log((double)N) / sqrt((double)N) / log(log((double)N)) / sqrt(sqrt(t - 1.0)));
	double minX = vertices[0].x, minY = vertices[0].y, maxX = vertices[0].x, maxY = vertices[0].y;
	for (int i = 1; i < N; i++) {
		if (vertices[i].x < minX) minX = vertices[i].x;
		if (vertices[i].x > maxX) maxX = vertices[i].x;
		if (vertices[i].y < minY) minY = vertices[i].y;
		if (vertices[i].y > maxY) maxY = vertices[i].y;
	}
	double xwidth = maxX - minX;
	double ywidth = maxY - minY;
	double maxwidth = xwidth > ywidth ? xwidth : ywidth;
	double estimatedLongestEdgeLength = estimatedLongestEdgeLengthUnit * maxwidth;
	double estimatedLongestTPath = estimatedLongestEdgeLength * (1 + t / 2.0);

	/*cout << "log N: " << std::fixed << log((double)N)  << endl;
	cout << "Sqrt N: " << std::fixed << sqrt((double)N) << endl;
	cout << "log log N: " << std::fixed << log(log((double)N)) << endl;
	cout << "sqrt(sqrt(t - 1.0)): " << std::fixed << sqrt(sqrt(t - 1.0)) << endl;
	
	cout << "estimatedLongestEdgeLengthUnit: " << std::fixed << estimatedLongestEdgeLengthUnit << endl;
	cout << "estimatedLongestEdgeLength: " << std::fixed << estimatedLongestEdgeLength << endl;
	cout << "estimatedLongestTPath: " << std::fixed << estimatedLongestTPath << endl;*/

	int bucketeers = (int)(1.0 / (estimatedLongestEdgeLengthUnit * (1 + t / 2.0)));
	
	//cout << "xwidth: " << std::fixed << xwidth << endl;
	//cout << "ywidth: " << std::fixed << ywidth << endl;

	if (bucketeers<1) { bucketeers = 1; }
	//cout << "bucketSide: " << bucketeers << endl;
	
	double bucketWidth = xwidth / ((double)bucketeers) + 1e-6;
	double bucketHeight = ywidth / ((double)bucketeers) + 1e-6;
	if (bucketWidth < 1.0) bucketWidth = 1e-6;
	if (bucketHeight < 1.0) bucketHeight = 1e-6;
	//cout << "bucketWidth: " << std::fixed << bucketWidth << endl;
	//cout << "bucketHeight: " << std::fixed << bucketHeight << endl;

	std::vector<int>* buckets = new std::vector<int>[bucketeers * bucketeers];
	int* bucketLookupX = new int[N];
	int* bucketLookupY = new int[N];
	//int* bucketLookupZ = new int[N];
	for (int i = 0; i < N; i++) {
		int bucketX = (int)((vertices[i].x - minX) / bucketWidth);
		int bucketY = (int)((vertices[i].y - minY) / bucketHeight);
		int bucketCoord = bucketX + bucketeers * bucketY;
		//bucketLookupZ[i] = buckets[bucketCoord].size();
		buckets[bucketCoord].push_back(i);
		bucketLookupX[i] = bucketX;
		bucketLookupY[i] = bucketY;
	}

	std::vector<int>* bucketsLarge = new std::vector<int>[bucketeers * bucketeers];
	for (int i = 0; i < N; i++) {
		int bucketX = bucketLookupX[i];
		int bucketY = bucketLookupY[i];
		int bucketCoord = bucketX + bucketeers * bucketY;
		for (int a = -1; a <= 1; a++) {
			for (int b = -1; b <= 1; b++) {
				int bucketLargeIndex = bucketX + a + bucketeers * (bucketY + b);
				if (bucketX + a >= 0 && bucketX + a < bucketeers
					&&
					bucketY + b >= 0 && bucketY + b < bucketeers) {
					bucketsLarge[bucketLargeIndex].push_back(i);
				}
			}
		}
	}

	CircleSector_list *Cones = new CircleSector_list[N];

	// Call base algorithm
	SIZE_T memUsed1 = ComputeGapSpannerBuckets<Heap>(T, edgeCount, myEdges, t, buckets, estimatedLongestEdgeLength, estimatedLongestTPath, bucketLookupX, bucketLookupY, 
		 bucketeers,Cones, bucketsLarge);

	int shortEdges = edgeCount;
	//cout << "Short Edges: " << edgeCount << endl;
	//cout << Timing::elapsed() << endl;

	AngleStorer* certificates = new AngleStorer[N];
	for (int i = 0; i < N; i++) {
		certificates[i].origin = vertices[i];

		CircleSector_list_node* v_item = Cones[i].first;
		CircleSector* cone;
		while (v_item != nullptr)
		{
			cone = v_item->data;
			certificates[i].Add(cone->startAngle, cone->endAngle);
			v_item = v_item->next;
		}
	}

	

	// Compute WSPD
	//double s = maximum((double)(4.0 / sin(theta / 2.0)), (double)(2.0 / w)); //4.0 * t / (t - 1.0);//5.0;
	double s = 5.0;
	SplitTree<IsImplicit>* splitTree = SplitTree<IsImplicit>::createSplitTree(vertices);
	std::vector<WellSeparatedPair<IsImplicit>*> pairs;
	splitTree->ComputeWspdPairs(s, pairs, estimatedLongestEdgeLength);
	//cout << "Long pairs: " << pairs.size() << endl;

	std::vector<WellSeparatedPair<IsImplicit>*> pairs1;
	splitTree->ComputeWspdPairs(s, pairs1);

	// Discount pairs
	bool* discountedPairs = new bool[pairs.size()];
	bool** discounters = new bool* [pairs.size()];
	unsigned int discountedPairsCounter = 0;
	unsigned int discountersCounter = 0;
	unsigned int discountersTotalCounter = 0;
	for (unsigned int i = 0; i < pairs.size(); i++) {
		discountedPairs[i] = false;
		bool foundBadPoint = false;
		for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->first->begin(); it != pairs[i]->first->end(); ++it) {
			unsigned int a = *it;
			if (!certificates[a].Covers(pairs[i]->second->box)) {
				foundBadPoint = true;
				break;
			}
		}
		if (!foundBadPoint) {
			discountedPairs[i] = true;
			discountedPairsCounter++;
		}
		else {
			foundBadPoint = false;
			for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->second->begin(); it != pairs[i]->second->end(); ++it) {
				unsigned int b = *it;
				if (!certificates[b].Covers(pairs[i]->first->box)) {
					foundBadPoint = true;
					break;
				}
			}
			if (!foundBadPoint) {
				discountedPairs[i] = true;
				discountedPairsCounter++;
			}
			else {
				discounters[i] = new bool[pairs[i]->first->count()];
				discountersTotalCounter += pairs[i]->first->count();
				unsigned int index = 0;
				for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->first->begin(); it != pairs[i]->first->end(); ++it, index++) {
					unsigned int a = *it;
					bool discount = certificates[a].Covers(pairs[i]->second->box);
					discounters[i][index] = discount;
					if (discount) discountersCounter++;
				}
			}
		}
	}
	//cout << "Discounted " << discountedPairsCounter << "/" << pairs.size() << endl;
	//cout << "Discounted " << discountersCounter << "/" << discountersTotalCounter << endl;
	//char filename[80];
	//sprintf(filename, "NMNM(%d-%f-%f).txt", N, theta, w);
	//std::ofstream myfile(filename, std::ofstream::out);

	unsigned int m = pairs.size();
	std::vector<WellSeparatedPairRepresentative*> QContent(m, nullptr);
	Heap Q(m, -std::numeric_limits<double>::infinity());
	bool* dirtyBits = new bool[m];
	for (unsigned int i = 0; i < m; i++)
		dirtyBits[i] = false;

	GapFillQueue3<Heap>(Q, QContent, pairs, vertices, myEdges, t, w, edgeCount, discountedPairs, dirtyBits, discounters, certificates, Cones);

	while (Q.getCount() > 0) {
		bool foundOne = false;
		std::pair<unsigned int, double> minQ;
		while (!foundOne) {
			minQ = Q.getMin();

			if (dirtyBits[minQ.first])
				GapUpdatePair2(minQ.first, Q, QContent, pairs[minQ.first], vertices, myEdges, t, w, dirtyBits, discounters, certificates, Cones);
			else
				foundOne = true;
			if (Q.getCount() == 0)
				goto end;
		}
		Q.extractMin();
		WellSeparatedPairRepresentative* minRepresentative = QContent[minQ.first];

		std::pair<unsigned int, double> pair(minRepresentative->u, minRepresentative->dist);
		myEdges[minRepresentative->v].push_back(pair);
		pair.first = minRepresentative->v;
		myEdges[minRepresentative->u].push_back(pair);
		if (T->create_edge(edge(minRepresentative->u, minRepresentative->v))) return 0;
		if (T->create_edge(edge(minRepresentative->v, minRepresentative->u))) return 0;
		edgeCount+=2;

		
		/*if (minRepresentative->v > minRepresentative->u)
		{
			myfile << "-" << minRepresentative->v << "&" << minRepresentative->u << "-" << "(" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ") ---> (" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ")" << endl;
			myfile << "-" << minRepresentative->u << "&" << minRepresentative->v << "-" << "(" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ") ---> (" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ")" << endl;
		}
		else
		{
			myfile << "-" << minRepresentative->u << "&" << minRepresentative->v << "-" << "(" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ") ---> (" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ")" << endl;
			myfile << "-" << minRepresentative->v << "&" << minRepresentative->u << "-" << "(" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ") ---> (" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ")" << endl;
		}*/
		

		double edgeDist = w*distance(vertices[minRepresentative->u], vertices[minRepresentative->v]);
		double  angyx = atan3(vertices[minRepresentative->u], vertices[minRepresentative->v]);
		double  angxy = atan3(vertices[minRepresentative->v], vertices[minRepresentative->u]);

		for (unsigned int a = 0; a < N; a++) {
			if (distance(vertices[a], vertices[minRepresentative->u]) <= edgeDist) {
				AddOrMergeCones(Cones, vertices, a, angxy, theta);
			}
		}

		for (unsigned int b = 0; b < N; b++) {
			if (distance(vertices[b], vertices[minRepresentative->v]) <= edgeDist) {
				AddOrMergeCones(Cones, vertices, b, angyx, theta);
			}
		}

		for (int j = 0; j < m; j++)
			if (QContent[j] != nullptr)
				dirtyBits[j] = true;

		GapReinsertPair2(minQ.first, Q, QContent, pairs[minQ.first], vertices, myEdges, t, w, dirtyBits, discounters, certificates, Cones);
	}
end:
	
	char filename[80];
	sprintf(filename, "AM(%d-%f-%f).txt", N, theta, w);
	std::ofstream myfile(filename, std::ofstream::app);
	
	PROCESS_MEMORY_COUNTERS_EX pmc;
	GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
	SIZE_T physMemUsedByMe = pmc.WorkingSetSize;

	if (memUsed1 > physMemUsedByMe)
		physMemUsedByMe = memUsed1;

	myfile << physMemUsedByMe;
	myfile << "\t" << pairs.size() << "\t" << discountedPairsCounter << "\t" << pairs1.size() << endl;
	myfile.close();


	delete[] buckets;
	delete[] bucketLookupX;
	delete[] bucketLookupY;

	delete[] myEdges;
	delete splitTree;
	for (unsigned int i = 0; i < m; i++) {
		if (!discountedPairs[i]) {
			delete[] discounters[i];
		}
	}
	delete[] discounters;
	delete[] discountedPairs;
	delete[] dirtyBits;
	delete[] bucketsLarge;
	delete[] Cones;

	/*
	cout <<endl<< "Long Pairs%\t\tDiscounted\t\t\tShort Edges" << endl;
	cout <<  (double)pairs.size() / (double)pairs1.size()<<"\t\t"
	<< discountedPairsCounter << "/" << pairs.size() << ", %=" << (double)discountedPairsCounter / (double)pairs.size() << "\t\t"
	<< (double)shortEdges / (double)edgeCount <<endl << endl;*/
	//cout << "Discounted " << discountersCounter << "/" << discountersTotalCounter << ", %=" << (double)discountersCounter / (double)discountersTotalCounter << endl;

	

	/*myfile << shortEdges << "\t" << edgeCount <<"\t" 
		<< (double)pairs.size() / (double)pairs1.size() << "\t"
		<< discountedPairsCounter << "/" << pairs.size() << ", %=" << (double)discountedPairsCounter / (double)pairs.size() << "\t"
		<< (double)shortEdges / (double)edgeCount << endl;*/

	myfile.close();

	return edgeCount;
}

//***************************************************************
//***************************************************************
template<class Heap, bool IsImplicit>
int GapSpannerBucketsFixed1(BaseTask* T) {
	//TODO:
	// Make ESA' local
	STARTLOG
	const pointset& vertices = T->p;
	unsigned int N = vertices.size();
	std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	unsigned int edgeCount = 0;

	// Bucketing stuff
	//double w_ratio = T->t;
	double theta = T->t;
	double w = T->var;
	double t = 1.0 / (cos(theta) - sin(theta) - 2 * w);
	double estimatedLongestEdgeLengthUnit = (log((double)N) / sqrt((double)N) / log(log((double)N)) / sqrt(sqrt(t - 1.0)));
	double minX = vertices[0].x, minY = vertices[0].y, maxX = vertices[0].x, maxY = vertices[0].y;
	for (int i = 1; i < N; i++) {
		if (vertices[i].x < minX) minX = vertices[i].x;
		if (vertices[i].x > maxX) maxX = vertices[i].x;
		if (vertices[i].y < minY) minY = vertices[i].y;
		if (vertices[i].y > maxY) maxY = vertices[i].y;
	}
	double xwidth = maxX - minX;
	double ywidth = maxY - minY;
	double maxwidth = xwidth > ywidth ? xwidth : ywidth;
	double estimatedLongestEdgeLength = estimatedLongestEdgeLengthUnit * maxwidth;
	double estimatedLongestTPath = estimatedLongestEdgeLength * (1 + t / 2.0);

	/*cout << "log N: " << std::fixed << log((double)N)  << endl;
	cout << "Sqrt N: " << std::fixed << sqrt((double)N) << endl;
	cout << "log log N: " << std::fixed << log(log((double)N)) << endl;
	cout << "sqrt(sqrt(t - 1.0)): " << std::fixed << sqrt(sqrt(t - 1.0)) << endl;

	cout << "estimatedLongestEdgeLengthUnit: " << std::fixed << estimatedLongestEdgeLengthUnit << endl;
	cout << "estimatedLongestTPath: " << std::fixed << estimatedLongestTPath << endl;*/
	
	//cout << "est Long Edge Len: " << std::fixed << estimatedLongestEdgeLength << endl;

	int bucketeers = (int)(1.0 / (estimatedLongestEdgeLengthUnit * (1 + t / 2.0)));

	//cout << "xwidth: " << std::fixed << xwidth << endl;
	//cout << "ywidth: " << std::fixed << ywidth << endl;

	if (bucketeers < 1) { bucketeers = 1; }
	//cout << "bucketSide: " << bucketeers << endl;

	double bucketWidth = xwidth / ((double)bucketeers) + 1e-6;
	double bucketHeight = ywidth / ((double)bucketeers) + 1e-6;
	if (bucketWidth < 1.0) bucketWidth = 1e-6;
	if (bucketHeight < 1.0) bucketHeight = 1e-6;
	//cout << "bucketWidth: " << std::fixed << bucketWidth << endl;
	//cout << "bucketHeight: " << std::fixed << bucketHeight << endl;

	std::vector<int>* buckets = new std::vector<int>[bucketeers * bucketeers];
	int* bucketLookupX = new int[N];
	int* bucketLookupY = new int[N];
	//int* bucketLookupZ = new int[N];
	for (int i = 0; i < N; i++) {
		int bucketX = (int)((vertices[i].x - minX) / bucketWidth);
		int bucketY = (int)((vertices[i].y - minY) / bucketHeight);
		int bucketCoord = bucketX + bucketeers * bucketY;
		//bucketLookupZ[i] = buckets[bucketCoord].size();
		buckets[bucketCoord].push_back(i);
		bucketLookupX[i] = bucketX;
		bucketLookupY[i] = bucketY;
	}

	std::vector<int>* bucketsLarge = new std::vector<int>[bucketeers * bucketeers];
	for (int i = 0; i < N; i++) {
		int bucketX = bucketLookupX[i];
		int bucketY = bucketLookupY[i];
		int bucketCoord = bucketX + bucketeers * bucketY;
		for (int a = -1; a <= 1; a++) {
			for (int b = -1; b <= 1; b++) {
				int bucketLargeIndex = bucketX + a + bucketeers * (bucketY + b);
				if (bucketX + a >= 0 && bucketX + a < bucketeers
					&&
					bucketY + b >= 0 && bucketY + b < bucketeers) {
					bucketsLarge[bucketLargeIndex].push_back(i);
				}
			}
		}
	}

	CircleSector_list* Cones = new CircleSector_list[N];

	// Call base algorithm
	ComputeGapSpannerBuckets<Heap>(T, edgeCount, myEdges, t, buckets, 
		estimatedLongestEdgeLength, estimatedLongestTPath, bucketLookupX, bucketLookupY,
		bucketeers, Cones, bucketsLarge);

	int shortEdges = edgeCount;
	//cout << "Short Edges: " << edgeCount << endl;
	//cout << Timing::elapsed() << endl;

	AngleStorer* certificates = new AngleStorer[N];
	for (int i = 0; i < N; i++) {
		certificates[i].origin = vertices[i];

		CircleSector_list_node* v_item = Cones[i].first;
		CircleSector* cone;
		while (v_item != nullptr)
		{
			cone = v_item->data;
			certificates[i].Add(cone->startAngle, cone->endAngle);
			v_item = v_item->next;
		}
	}

	delete[] buckets;
	delete[] bucketLookupX;
	delete[] bucketLookupY;

	// Compute WSPD
	double s = 5.0;
	SplitTree<IsImplicit>* splitTree = SplitTree<IsImplicit>::createSplitTree(vertices);
	std::vector<WellSeparatedPair<IsImplicit>*> pairs;
	splitTree->ComputeWspdPairs(s, pairs, estimatedLongestEdgeLength);
	//cout << "Long pairs: " << pairs.size() << endl;

	std::vector<WellSeparatedPair<IsImplicit>*> pairs1;
	splitTree->ComputeWspdPairs(s, pairs1);

	// Discount pairs
	bool* discountedPairs = new bool[pairs.size()];
	bool** discounters = new bool* [pairs.size()];
	unsigned int discountedPairsCounter = 0;
	unsigned int discountersCounter = 0;
	unsigned int discountersTotalCounter = 0;
	for (unsigned int i = 0; i < pairs.size(); i++) {
		discountedPairs[i] = false;
		bool foundBadPoint = false;
		for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->first->begin(); it != pairs[i]->first->end(); ++it) {
			unsigned int a = *it;
			if (!certificates[a].Covers(pairs[i]->second->box)) {
				foundBadPoint = true;
				break;
			}
		}
		if (!foundBadPoint) {
			discountedPairs[i] = true;
			discountedPairsCounter++;
		}
		else {
			foundBadPoint = false;
			for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->second->begin(); it != pairs[i]->second->end(); ++it) {
				unsigned int b = *it;
				if (!certificates[b].Covers(pairs[i]->first->box)) {
					foundBadPoint = true;
					break;
				}
			}
			if (!foundBadPoint) {
				discountedPairs[i] = true;
				discountedPairsCounter++;
			}
			else {
				discounters[i] = new bool[pairs[i]->first->count()];
				discountersTotalCounter += pairs[i]->first->count();
				unsigned int index = 0;
				for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->first->begin(); it != pairs[i]->first->end(); ++it, index++) {
					unsigned int a = *it;
					bool discount = certificates[a].Covers(pairs[i]->second->box);
					discounters[i][index] = discount;
					if (discount) discountersCounter++;
				}
			}
		}
	}
	//cout << "Discounted " << discountedPairsCounter << "/" << pairs.size() << endl;
	//cout << "Discounted " << discountersCounter << "/" << discountersTotalCounter << endl;
	//char filename[80];
	//sprintf(filename, "NMNM(%d-%f-%f).txt", N, theta, w);
	//std::ofstream myfile(filename, std::ofstream::out);

	std::vector<int>* hostPairs = new std::vector<int>[N];

	for (unsigned int i = 0; i < pairs.size(); i++) {
		if (discountedPairs[i] == true)
			continue;

		for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->first->begin(); it != pairs[i]->first->end(); ++it) {
			unsigned int a = *it;
			hostPairs[a].push_back(i);
		}
		for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->second->begin(); it != pairs[i]->second->end(); ++it) {
			unsigned int b = *it;
			hostPairs[b].push_back(i);
		}
	}

	unsigned int m = pairs.size();
	
	std::vector<WellSeparatedPairRepresentative*> QContent(m, nullptr);
	Heap Q(m, -std::numeric_limits<double>::infinity());
	bool* dirtyBits = new bool[m];
	for (unsigned int i = 0; i < m; i++)
		dirtyBits[i] = false;

	GapFillQueue3<Heap>(Q, QContent, pairs, vertices, myEdges, t, w, edgeCount, discountedPairs, dirtyBits, discounters, certificates, Cones);

	while (Q.getCount() > 0) {
		bool foundOne = false;
		std::pair<unsigned int, double> minQ;
		while (!foundOne) {
			minQ = Q.getMin();

			if (dirtyBits[minQ.first])
				GapUpdatePair2(minQ.first, Q, QContent, pairs[minQ.first], vertices, myEdges, t, w, dirtyBits, discounters, certificates, Cones);
			else
				foundOne = true;
			if (Q.getCount() == 0)
				goto end;
		}
		Q.extractMin();
		WellSeparatedPairRepresentative* minRepresentative = QContent[minQ.first];

		std::pair<unsigned int, double> pair(minRepresentative->u, minRepresentative->dist);
		myEdges[minRepresentative->v].push_back(pair);
		pair.first = minRepresentative->v;
		myEdges[minRepresentative->u].push_back(pair);
		if (T->create_edge(edge(minRepresentative->u, minRepresentative->v))) return 0;
		if (T->create_edge(edge(minRepresentative->v, minRepresentative->u))) return 0;
		edgeCount += 2;


		/*if (minRepresentative->v > minRepresentative->u)
		{
			myfile << "-" << minRepresentative->v << "&" << minRepresentative->u << "-" << "(" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ") ---> (" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ")" << endl;
			myfile << "-" << minRepresentative->u << "&" << minRepresentative->v << "-" << "(" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ") ---> (" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ")" << endl;
		}
		else
		{
			myfile << "-" << minRepresentative->u << "&" << minRepresentative->v << "-" << "(" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ") ---> (" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ")" << endl;
			myfile << "-" << minRepresentative->v << "&" << minRepresentative->u << "-" << "(" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ") ---> (" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ")" << endl;
		}*/


		double edgeDist = w * distance(vertices[minRepresentative->u], vertices[minRepresentative->v]);
		edgeDist = edgeDist * edgeDist;
		double  angyx = atan3(vertices[minRepresentative->u], vertices[minRepresentative->v]);
		double  angxy = atan3(vertices[minRepresentative->v], vertices[minRepresentative->u]);

		for (unsigned int a = 0; a < N; a++) {
			if (NotSquaredDistance(vertices[a], vertices[minRepresentative->u]) <= edgeDist) {
				if (AddOrMergeCones(Cones, vertices, a, angxy, theta) != 5)
				{
					for (size_t s = 0; s < hostPairs[a].size(); s++)
					{
						dirtyBits[hostPairs[a][s]] = true;
					}
				}
			}
		}

		for (unsigned int b = 0; b < N; b++) {
			if (NotSquaredDistance(vertices[b], vertices[minRepresentative->v]) <= edgeDist) {
				if (AddOrMergeCones(Cones, vertices, b, angyx, theta) != 5)
				{
					for (size_t s = 0; s < hostPairs[b].size(); s++)
					{
						dirtyBits[hostPairs[b][s]] = true;
					}
				}
			}
		}

		/*for (int j = 0; j < m; j++)
			if (QContent[j] != nullptr)
				dirtyBits[j] = true;*/

		GapReinsertPair2(minQ.first, Q, QContent, pairs[minQ.first], vertices, myEdges, t, w, dirtyBits, discounters, certificates, Cones);
	}
end:


	std::vector<std::pair<std::pair<int, int>, double>> partitionEdgeCount;
	std::vector<double> edgeList;
	double currentPartitionMinEdgeLen;

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < myEdges[i].size(); j++)
		{
			edgeList.push_back(myEdges[i][j].second);
		}
	}

	std::sort(edgeList.begin(), edgeList.end());
	currentPartitionMinEdgeLen = edgeList[0];
	int curPartition = 0;
	partitionEdgeCount
		.push_back(std::pair<std::pair<int, int>, double>(std::pair<int, int>(1,0), currentPartitionMinEdgeLen));

	for (size_t i = 1; i < edgeList.size(); i++)
	{
		if (edgeList[i] > 2 * currentPartitionMinEdgeLen)
		{
			curPartition++;
			currentPartitionMinEdgeLen = edgeList[i];
			partitionEdgeCount//.push_back(std::pair<int, double>(1, currentPartitionMinEdgeLen));
				.push_back(std::pair<std::pair<int, int>, double>(std::pair<int, int>(1, 0), currentPartitionMinEdgeLen));
		}
		else
		{
			partitionEdgeCount[curPartition].first.first++;
		}
	}
	
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = i+1; j < N; j++)
		{
			double edge_size = distance(vertices[i], vertices[j]);
			int k = 0;
			while (k < partitionEdgeCount.size() && edge_size >= partitionEdgeCount[k].second)
				k++;

			partitionEdgeCount[k--].first.second++;
		}
	}
	delete[] myEdges;
	delete splitTree;
	for (unsigned int i = 0; i < m; i++) {
		if (!discountedPairs[i]) {
			delete[] discounters[i];
		}
	}
	delete[] discounters;
	delete[] discountedPairs;
	delete[] dirtyBits;
	delete[] bucketsLarge;
	delete[] Cones;
	delete[] hostPairs;

	/*
	cout <<endl<< "Long Pairs%\t\tDiscounted\t\t\tShort Edges" << endl;
	cout <<  (double)pairs.size() / (double)pairs1.size()<<"\t\t"
	<< discountedPairsCounter << "/" << pairs.size() << ", %=" << (double)discountedPairsCounter / (double)pairs.size() << "\t\t"
	<< (double)shortEdges / (double)edgeCount <<endl << endl;*/
	//cout << "Discounted " << discountersCounter << "/" << discountersTotalCounter << ", %=" << (double)discountersCounter / (double)discountersTotalCounter << endl;

	char filename[80];
	sprintf(filename, "NMNM(%d-%f-%f).txt", N, theta, w);
	std::ofstream myfile(filename, std::ofstream::app);

	/*myfile << shortEdges << "\t" << edgeCount << "\t"
		<< (double)pairs.size() / (double)pairs1.size() << "\t"
		<< discountedPairsCounter << "/" << pairs.size() << ", %=" << (double)discountedPairsCounter / (double)pairs.size() << "\t"
		<< (double)shortEdges / (double)edgeCount << endl;*/
	double total = 0, totalPairs = 0;
	int allPairs = N * (N-1)/2;
	for (size_t i = 0; i < partitionEdgeCount.size(); i++)
	{
		total += (double)partitionEdgeCount[i].first.first / edgeCount;
		totalPairs += (double)partitionEdgeCount[i].first.second / allPairs;
		myfile 
			<< std::fixed 
			<< partitionEdgeCount[i].second
			<< "(" 
			<< std::setprecision(4)
			<< totalPairs
			<<"%-"
			//<< (double) partitionEdgeCount[i].first.first / edgeCount 
			//<< "% - "
			<< total
			<< "%), ";
	}
	myfile << endl;

	myfile.close();

	return edgeCount;
}

//***************************************************************

template<class Heap, bool IsImplicit>
int GapSpannerBucketsFixed2(BaseTask* T) {
	//TODO:
	// Make ESA' local
	STARTLOG
		const pointset& vertices = T->p;
	unsigned int N = vertices.size();
	std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	unsigned int edgeCount = 0;

	// Bucketing stuff
	//double w_ratio = T->t;
	double theta = T->t;
	double w = T->var;
	double t = 1.0 / (cos(theta) - sin(theta) - 2 * w);
	double estimatedLongestEdgeLengthUnit = (log((double)N) / sqrt((double)N) / log(log((double)N)) / sqrt(sqrt(t - 1.0)));
	double minX = vertices[0].x, minY = vertices[0].y, maxX = vertices[0].x, maxY = vertices[0].y;
	for (int i = 1; i < N; i++) {
		if (vertices[i].x < minX) minX = vertices[i].x;
		if (vertices[i].x > maxX) maxX = vertices[i].x;
		if (vertices[i].y < minY) minY = vertices[i].y;
		if (vertices[i].y > maxY) maxY = vertices[i].y;
	}
	double xwidth = maxX - minX;
	double ywidth = maxY - minY;
	double maxwidth = xwidth > ywidth ? xwidth : ywidth;
	double estimatedLongestEdgeLength = estimatedLongestEdgeLengthUnit * maxwidth;
	double estimatedLongestTPath = estimatedLongestEdgeLength * (1 + t / 2.0);



	int bucketeers = (int)(1.0 / (estimatedLongestEdgeLengthUnit * (1 + t / 2.0)));

	//cout << "xwidth: " << std::fixed << xwidth << endl;
	//cout << "ywidth: " << std::fixed << ywidth << endl;

	if (bucketeers < 1) { bucketeers = 1; }
	//cout << "bucketSide: " << bucketeers << endl;

	double bucketWidth = xwidth / ((double)bucketeers) + 1e-6;
	double bucketHeight = ywidth / ((double)bucketeers) + 1e-6;
	if (bucketWidth < 1.0) bucketWidth = 1e-6;
	if (bucketHeight < 1.0) bucketHeight = 1e-6;
	//cout << "bucketWidth: " << std::fixed << bucketWidth << endl;
	//cout << "bucketHeight: " << std::fixed << bucketHeight << endl;

	std::vector<int>* buckets = new std::vector<int>[bucketeers * bucketeers];
	int* bucketLookupX = new int[N];
	int* bucketLookupY = new int[N];
	//int* bucketLookupZ = new int[N];
	for (int i = 0; i < N; i++) {
		int bucketX = (int)((vertices[i].x - minX) / bucketWidth);
		int bucketY = (int)((vertices[i].y - minY) / bucketHeight);
		int bucketCoord = bucketX + bucketeers * bucketY;
		//bucketLookupZ[i] = buckets[bucketCoord].size();
		buckets[bucketCoord].push_back(i);
		bucketLookupX[i] = bucketX;
		bucketLookupY[i] = bucketY;
	}

	std::vector<int>* bucketsLarge = new std::vector<int>[bucketeers * bucketeers];
	for (int i = 0; i < N; i++) {
		int bucketX = bucketLookupX[i];
		int bucketY = bucketLookupY[i];
		int bucketCoord = bucketX + bucketeers * bucketY;
		for (int a = -1; a <= 1; a++) {
			for (int b = -1; b <= 1; b++) {
				int bucketLargeIndex = bucketX + a + bucketeers * (bucketY + b);
				if (bucketX + a >= 0 && bucketX + a < bucketeers
					&&
					bucketY + b >= 0 && bucketY + b < bucketeers) {
					bucketsLarge[bucketLargeIndex].push_back(i);
				}
			}
		}
	}

	CircleSector_list* Cones = new CircleSector_list[N];

	// Call base algorithm
	ComputeGapSpannerBuckets<Heap>(T, edgeCount, myEdges, t, buckets, estimatedLongestEdgeLength, estimatedLongestTPath, bucketLookupX, bucketLookupY,
		bucketeers, Cones, bucketsLarge);

	int shortEdges = edgeCount;
	//cout << "Short Edges: " << edgeCount << endl;
	//cout << Timing::elapsed() << endl;

	AngleStorer* certificates = new AngleStorer[N];
	for (int i = 0; i < N; i++) {
		certificates[i].origin = vertices[i];

		CircleSector_list_node* v_item = Cones[i].first;
		CircleSector* cone;
		while (v_item != nullptr)
		{
			cone = v_item->data;
			certificates[i].Add(cone->startAngle, cone->endAngle);
			v_item = v_item->next;
		}
	}

	delete[] buckets;
	delete[] bucketLookupX;
	delete[] bucketLookupY;

	// Compute WSPD
	double s = 5.0;
	SplitTree<IsImplicit>* splitTree = SplitTree<IsImplicit>::createSplitTree(vertices);
	std::vector<WellSeparatedPair<IsImplicit>*> pairs;
	splitTree->ComputeWspdPairs(s, pairs, estimatedLongestEdgeLength);
	//cout << "Long pairs: " << pairs.size() << endl;

	std::vector<WellSeparatedPair<IsImplicit>*> pairs1;
	splitTree->ComputeWspdPairs(s, pairs1);

	// Discount pairs
	bool* discountedPairs = new bool[pairs.size()];
	bool** discounters = new bool* [pairs.size()];
	unsigned int discountedPairsCounter = 0;
	unsigned int discountersCounter = 0;
	unsigned int discountersTotalCounter = 0;
	for (unsigned int i = 0; i < pairs.size(); i++) {
		discountedPairs[i] = false;
		bool foundBadPoint = false;
		for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->first->begin(); it != pairs[i]->first->end(); ++it) {
			unsigned int a = *it;
			if (!certificates[a].Covers(pairs[i]->second->box)) {
				foundBadPoint = true;
				break;
			}
		}
		if (!foundBadPoint) {
			discountedPairs[i] = true;
			discountedPairsCounter++;
		}
		else {
			foundBadPoint = false;
			for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->second->begin(); it != pairs[i]->second->end(); ++it) {
				unsigned int b = *it;
				if (!certificates[b].Covers(pairs[i]->first->box)) {
					foundBadPoint = true;
					break;
				}
			}
			if (!foundBadPoint) {
				discountedPairs[i] = true;
				discountedPairsCounter++;
			}
			else {
				discounters[i] = new bool[pairs[i]->first->count()];
				discountersTotalCounter += pairs[i]->first->count();
				unsigned int index = 0;
				for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->first->begin(); it != pairs[i]->first->end(); ++it, index++) {
					unsigned int a = *it;
					bool discount = certificates[a].Covers(pairs[i]->second->box);
					discounters[i][index] = discount;
					if (discount) discountersCounter++;
				}
			}
		}
	}
	//cout << "Discounted " << discountedPairsCounter << "/" << pairs.size() << endl;
	//cout << "Discounted " << discountersCounter << "/" << discountersTotalCounter << endl;
	//char filename[80];
	//sprintf(filename, "XX1(%d-%f-%f).txt", N, theta, w);
	//std::ofstream myfile("XX1.txt", std::ofstream::out);

	std::vector<int>* hostPairs = new std::vector<int>[N];

	for (unsigned int i = 0; i < pairs.size(); i++) {
		if (discountedPairs[i] == true)
			continue;

		for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->first->begin(); it != pairs[i]->first->end(); ++it) {
			unsigned int a = *it;
			hostPairs[a].push_back(i);
		}
		for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->second->begin(); it != pairs[i]->second->end(); ++it) {
			unsigned int b = *it;
			hostPairs[b].push_back(i);
		}
	}

	unsigned int m = pairs.size();

	std::vector<std::pair<unsigned int, unsigned int>>* candidatePairs = new std::vector<std::pair<unsigned int, unsigned int>>[m];

	
	std::vector<WellSeparatedPairRepresentative*> QContent(m, nullptr);
	Heap Q(m, -std::numeric_limits<double>::infinity());
	bool* dirtyBits = new bool[m];
	for (unsigned int i = 0; i < m; i++)
		dirtyBits[i] = false;

	//cout << Timing::elapsed() << "\t\t";
	GapFillQueue4<Heap>(Q, QContent, pairs, vertices, myEdges, t, w, edgeCount, discountedPairs, dirtyBits, discounters, certificates, Cones, candidatePairs);
	//cout << Timing::elapsed() << endl;

	while (Q.getCount() > 0) {
		bool foundOne = false;
		std::pair<unsigned int, double> minQ;
		while (!foundOne) {
			minQ = Q.getMin();

			if (dirtyBits[minQ.first])
				GapUpdatePair4(minQ.first, Q, QContent, pairs[minQ.first], vertices, myEdges, t, w, dirtyBits, discounters, certificates, Cones, candidatePairs);
			else
				foundOne = true;
			if (Q.getCount() == 0)
				goto end;
		}
		Q.extractMin();
		WellSeparatedPairRepresentative* minRepresentative = QContent[minQ.first];

		std::pair<unsigned int, double> pair(minRepresentative->u, minRepresentative->dist);
		myEdges[minRepresentative->v].push_back(pair);
		pair.first = minRepresentative->v;
		myEdges[minRepresentative->u].push_back(pair);
		if (T->create_edge(edge(minRepresentative->u, minRepresentative->v))) return 0;
		if (T->create_edge(edge(minRepresentative->v, minRepresentative->u))) return 0;
		edgeCount += 2;





		double edgeDist = w * minQ.second;//distance(vertices[minRepresentative->u], vertices[minRepresentative->v]);
		edgeDist = edgeDist * edgeDist;
		double  angyx = atan3(vertices[minRepresentative->u], vertices[minRepresentative->v]);
		double  angxy = atan3(vertices[minRepresentative->v], vertices[minRepresentative->u]);

		for (unsigned int a = 0; a < N; a++) {
			if (NotSquaredDistance(vertices[a], vertices[minRepresentative->u]) <= edgeDist) {
				if (AddOrMergeCones(Cones, vertices, a, angxy, theta) != 5)
				{
					for (size_t s = 0; s < hostPairs[a].size(); s++)
					{
						dirtyBits[hostPairs[a][s]] = true;
					}
				}
			}
		}

		for (unsigned int b = 0; b < N; b++) {
			if (NotSquaredDistance(vertices[b], vertices[minRepresentative->v]) <= edgeDist) {
				if (AddOrMergeCones(Cones, vertices, b, angyx, theta) != 5)
				{
					for (size_t s = 0; s < hostPairs[b].size(); s++)
					{
						dirtyBits[hostPairs[b][s]] = true;
					}
				}
			}
		}



		GapReinsertPair4(minQ.first, Q, QContent, pairs[minQ.first], vertices, myEdges, t, w, dirtyBits, discounters, certificates, Cones, candidatePairs);
	}
end:
	delete[] myEdges;
	delete splitTree;
	for (unsigned int i = 0; i < m; i++) {
		if (!discountedPairs[i]) {
			delete[] discounters[i];
		}
	}
	delete[] discounters;
	delete[] discountedPairs;
	delete[] dirtyBits;
	delete[] bucketsLarge;
	delete[] Cones;
	delete[] hostPairs;



	return edgeCount;
}

// ***********************************************************
template<class Heap, bool IsImplicit>
void GapFillQueue4(Heap& Q,
	std::vector<WellSeparatedPairRepresentative*>& QContent,
	std::vector<WellSeparatedPair<IsImplicit>*>& pairs,
	const pointset& vertices,
	std::vector< std::pair<unsigned int, double> >* myEdges,
	double t,
	double w,
	int edgeCount,
	bool* discountedPairs,
	bool* dirtyBits,
	bool** discounters,
	AngleStorer* certificates,
	CircleSector_list* Cones,
	std::vector<std::pair<unsigned int, unsigned int>>* candidatePairs
) {
	for (unsigned int j = 0; j < pairs.size(); j++) {
		if (!discountedPairs[j]) {
			WellSeparatedPairRepresentative* p = GapClosestPair4<Heap, IsImplicit>(j, pairs[j], vertices, t, w, discounters, Cones, candidatePairs);
			if (p == nullptr)
				QContent[j] = nullptr;
			else {
				Q.insert(j, p->dist);
				QContent[j] = p;
			}
		}
	}
}

template<class Heap, bool IsImplicit>
WellSeparatedPairRepresentative* GapClosestPair4(unsigned int pairNo,
	WellSeparatedPair<IsImplicit>* pair,
	const pointset& vertices,
	double t,
	double w,
	bool** discounters,
	CircleSector_list* Cones,
	std::vector<std::pair<unsigned int, unsigned int>>* candidatePairs
) {
	int shortestPairInedx1, shortestPairInedx2;
	double bestDist = std::numeric_limits<double>::infinity();
	unsigned int j = 0;
	for (typename SplitTree<IsImplicit>::Iterator it = pair->first->begin(); it != pair->first->end(); ++it, j++) {
		unsigned int indexj = *it;
		if (discounters[pairNo][j] == false) {
			bool discounted = true;

			for (typename SplitTree<IsImplicit>::Iterator it1 = pair->second->begin(); it1 != pair->second->end(); ++it1) {
				unsigned int indexk = *it1;
				
				if (!IsInCircleSector(Cones, vertices, indexj, indexk) && !IsInCircleSector(Cones, vertices, indexk, indexj)) {
					double   dist = NotSquaredDistance(vertices[indexj], vertices[indexk]);//distance(vertices[indexj], vertices[indexk]);
					discounted = false;
					if (dist < bestDist)
					{
						bestDist = dist;
						shortestPairInedx1 = indexj;
						shortestPairInedx2 = indexk;
					}
					candidatePairs[pairNo].push_back(std::pair<unsigned int, unsigned int>(indexj, indexk));
				}
			}
			discounters[pairNo][j] = discounted;
		}
	}

	if (std::isinf(bestDist))
		return nullptr;
	else
		return new WellSeparatedPairRepresentative(shortestPairInedx1, shortestPairInedx2, sqrt(bestDist), pairNo);
}

template<class Heap, bool IsImplicit>
WellSeparatedPairRepresentative* GapClosestPair3(unsigned int pairNo,
	WellSeparatedPair<IsImplicit>* pair,
	const pointset& vertices,
	double t,
	double w,
	bool** discounters,
	CircleSector_list* Cones,
	std::vector<std::pair<unsigned int, unsigned int>>* candidatePairs
) {
	int shortestPairInedx1, shortestPairInedx2;
	double bestDist = std::numeric_limits<double>::infinity();
	std::vector<std::pair<unsigned int, unsigned int>> newCandidatePairs;

	for (unsigned int i = 0; i < candidatePairs[pairNo].size(); i++)
	{
		if (!IsInCircleSector(Cones, vertices, candidatePairs[pairNo][i].first, candidatePairs[pairNo][i].second)
			&& !IsInCircleSector(Cones, vertices, candidatePairs[pairNo][i].second, candidatePairs[pairNo][i].first))
		{
			double   dist = NotSquaredDistance(vertices[candidatePairs[pairNo][i].first], vertices[candidatePairs[pairNo][i].second]);
			if (dist < bestDist)
			{
				shortestPairInedx1 = candidatePairs[pairNo][i].first;
				shortestPairInedx2 = candidatePairs[pairNo][i].second;
				bestDist = dist;
			}
			
			newCandidatePairs.push_back(candidatePairs[pairNo][i]);
		}
	}
	candidatePairs[pairNo].clear();

	for (unsigned int i = 0; i < newCandidatePairs.size(); i++)
	{
		candidatePairs[pairNo].push_back(newCandidatePairs[i]);
	}

	if (std::isinf(bestDist))
		return nullptr;
	else
		return new WellSeparatedPairRepresentative(shortestPairInedx1, shortestPairInedx2, sqrt(bestDist), pairNo);
}

template<class Heap, bool IsImplicit>
void GapUpdatePair4(unsigned int j,
	Heap& Q,
	std::vector<WellSeparatedPairRepresentative*>& QContent,
	WellSeparatedPair<IsImplicit>* pair,
	const pointset& vertices,
	std::vector< std::pair<unsigned int, double> >* myEdges,
	double t,
	double w,
	bool* dirtyBits,
	bool** discounters,
	AngleStorer* certificates,
	CircleSector_list* Cones,
	std::vector<std::pair<unsigned int, unsigned int>>* candidatePairs
) {
	WellSeparatedPairRepresentative* p = GapClosestPair3<Heap, IsImplicit>(j, pair, vertices, t, w, discounters, Cones, candidatePairs);
	if (p == nullptr) {
		Q.remove(j);
		QContent[j] = nullptr;
	}
	else {
		Q.increaseKey(j, p->dist);
		QContent[j] = p;
	}
	dirtyBits[j] = false;
}

template<class Heap, bool IsImplicit>
void GapReinsertPair4(unsigned int j,
	Heap& Q,
	std::vector<WellSeparatedPairRepresentative*>& QContent,
	WellSeparatedPair<IsImplicit>* pair,
	const pointset& vertices,
	std::vector< std::pair<unsigned int, double> >* myEdges,
	double t,
	double w,
	bool* dirtyBits,
	bool** discounters,
	AngleStorer* certificates,
	CircleSector_list* Cones,
	std::vector<std::pair<unsigned int, unsigned int>>* candidatePairs
) {
	WellSeparatedPairRepresentative* p = GapClosestPair3<Heap, IsImplicit>(j, pair, vertices, t, w, discounters, Cones, candidatePairs);
	if (p != nullptr) {
		Q.insert(j, p->dist);
		QContent[j] = p;
	}
	//else
		//QContent[j] = nullptr;
	dirtyBits[j] = false;
}

///------------------------------------------------------------

template<class Heap, bool IsImplicit>
int GapSpannerBucketsFixed3(BaseTask* T) {
	//TODO:
	// Make ESA' local
	STARTLOG
		const pointset& vertices = T->p;
	unsigned int N = vertices.size();
	std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	unsigned int edgeCount = 0;

	// Bucketing stuff
	//double w_ratio = T->t;
	double theta = T->t;
	double w = T->var;
	double t = 1.0 / (cos(theta) - sin(theta) - 2 * w);
	double estimatedLongestEdgeLengthUnit = (log((double)N) / sqrt((double)N) / log(log((double)N)) / sqrt(sqrt(t - 1.0)));
	double minX = vertices[0].x, minY = vertices[0].y, maxX = vertices[0].x, maxY = vertices[0].y;
	for (int i = 1; i < N; i++) {
		if (vertices[i].x < minX) minX = vertices[i].x;
		if (vertices[i].x > maxX) maxX = vertices[i].x;
		if (vertices[i].y < minY) minY = vertices[i].y;
		if (vertices[i].y > maxY) maxY = vertices[i].y;
	}
	double xwidth = maxX - minX;
	double ywidth = maxY - minY;
	double maxwidth = xwidth > ywidth ? xwidth : ywidth;
	double estimatedLongestEdgeLength = estimatedLongestEdgeLengthUnit * maxwidth;
	double estimatedLongestTPath = estimatedLongestEdgeLength * (1 + t / 2.0);

	int bucketeers = (int)(1.0 / (estimatedLongestEdgeLengthUnit * (1 + t / 2.0)));

	if (bucketeers < 1) { bucketeers = 1; }

	double bucketWidth = xwidth / ((double)bucketeers) + 1e-6;
	double bucketHeight = ywidth / ((double)bucketeers) + 1e-6;
	if (bucketWidth < 1.0) bucketWidth = 1e-6;
	if (bucketHeight < 1.0) bucketHeight = 1e-6;

	std::vector<int>* buckets = new std::vector<int>[bucketeers * bucketeers];
	int* bucketLookupX = new int[N];
	int* bucketLookupY = new int[N];
	for (int i = 0; i < N; i++) {
		int bucketX = (int)((vertices[i].x - minX) / bucketWidth);
		int bucketY = (int)((vertices[i].y - minY) / bucketHeight);
		int bucketCoord = bucketX + bucketeers * bucketY;
		//bucketLookupZ[i] = buckets[bucketCoord].size();
		buckets[bucketCoord].push_back(i);
		bucketLookupX[i] = bucketX;
		bucketLookupY[i] = bucketY;
	}

	std::vector<int>* bucketsLarge = new std::vector<int>[bucketeers * bucketeers];
	for (int i = 0; i < N; i++) {
		int bucketX = bucketLookupX[i];
		int bucketY = bucketLookupY[i];
		int bucketCoord = bucketX + bucketeers * bucketY;
		for (int a = -1; a <= 1; a++) {
			for (int b = -1; b <= 1; b++) {
				int bucketLargeIndex = bucketX + a + bucketeers * (bucketY + b);
				if (bucketX + a >= 0 && bucketX + a < bucketeers
					&&
					bucketY + b >= 0 && bucketY + b < bucketeers) {
					bucketsLarge[bucketLargeIndex].push_back(i);
				}
			}
		}
	}

	CircleSector_list* Cones = new CircleSector_list[N];

	// Call base algorithm
	ComputeGapSpannerBuckets<Heap>(T, edgeCount, myEdges, t, buckets, estimatedLongestEdgeLength, estimatedLongestTPath, bucketLookupX, bucketLookupY,
		bucketeers, Cones, bucketsLarge);

	int shortEdges = edgeCount;
	//cout << "Short Edges: " << edgeCount << endl;
	//cout << Timing::elapsed() << endl;

	AngleStorer* certificates = new AngleStorer[N];
	for (int i = 0; i < N; i++) {
		certificates[i].origin = vertices[i];

		CircleSector_list_node* v_item = Cones[i].first;
		CircleSector* cone;
		while (v_item != nullptr)
		{
			cone = v_item->data;
			certificates[i].Add(cone->startAngle, cone->endAngle);
			v_item = v_item->next;
		}
	}

	delete[] buckets;
	delete[] bucketLookupX;
	delete[] bucketLookupY;

	// Compute WSPD
	double s = 5;//maximum((double)(4.0 / sin(theta / 2.0)), (double)(2.0 / w));//5.0;
	SplitTree<IsImplicit>* splitTree = SplitTree<IsImplicit>::createSplitTree(vertices);
	std::vector<WellSeparatedPair<IsImplicit>*> pairs;
	splitTree->ComputeWspdPairs(s, pairs, estimatedLongestEdgeLength);
	//cout << "Long pairs: " << pairs.size() << endl;

	std::vector<WellSeparatedPair<IsImplicit>*> pairs1;
	splitTree->ComputeWspdPairs(s, pairs1);

	//cout << Timing::elapsed() << "\t\t";
	std::sort(pairs.begin(), pairs.end(), WellSeparatedPairPointerComparer<IsImplicit>());

	HeapIndex heapIndex(pairs.size());
	// Discount pairs
	bool** discounters = new bool* [pairs.size()];
	
	for (size_t i = 0; i < pairs.size(); i++)
	{
		discounters[i] = new bool[pairs[i]->first->count()];
		for (size_t k = 0; k < pairs[i]->first->count(); k++)
		{
			discounters[i][k] = false;
		}
	}

	std::vector<int>* hostPairs = new std::vector<int>[N];

	for (unsigned int i = 0; i < pairs.size(); i++) {

		for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->first->begin(); it != pairs[i]->first->end(); ++it) {
			unsigned int a = *it;
			hostPairs[a].push_back(i);
		}
		for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->second->begin(); it != pairs[i]->second->end(); ++it) {
			unsigned int b = *it;
			hostPairs[b].push_back(i);
		}
	}
	//char filename[80];
	//sprintf(filename, "XX2.txt");
	//std::ofstream myfile("XX2.txt", std::ofstream::out);

	unsigned int	m = pairs.size(),
					index = 0;

	std::vector<WellSeparatedPairRepresentative*> QContent(m, nullptr);
	Heap Q(m, -std::numeric_limits<double>::infinity());
	bool* dirtyBits = new bool[m];
	for (size_t i = 0; i < m; i++)
	{
		dirtyBits[i] = false;
	}
	//cout << Timing::elapsed() << "\t\t";
	GapFillQueue5<Heap>(Q, QContent, pairs, vertices, index, theta, w, discounters, certificates, Cones, heapIndex);
	//cout << Timing::elapsed() << "\n";

	std::vector<HeapIndexNode*>	heapNodes;

	while (Q.getCount() > 0) {

		std::pair<unsigned int, double > minQ = Q.getMin();

		WellSeparatedPairRepresentative* minRepresentative = QContent[minQ.first];
		
		if (T->create_edge(edge(minRepresentative->u, minRepresentative->v))) return 0;
		if (T->create_edge(edge(minRepresentative->v, minRepresentative->u))) return 0;
		edgeCount += 2;



		double edgeDist = w * minQ.second;
		edgeDist = edgeDist * edgeDist;
		double  angyx = atan3(vertices[minRepresentative->u], vertices[minRepresentative->v]);
		double  angxy = atan3(vertices[minRepresentative->v], vertices[minRepresentative->u]);

		for (unsigned int a = 0; a < N; a++) {
			if (NotSquaredDistance(vertices[a], vertices[minRepresentative->u]) <= edgeDist) {
				AddOrMergeCones(Cones, vertices, a, angxy, theta);
				double	sa = angxy - theta,
						ea = angxy + theta;
				if (sa < 0)
					sa += 2.0 * M_PI;

				if (ea > 2.0 * M_PI)
					ea -= 2.0 * M_PI;
				
				certificates[a].Add(sa, ea);
				
				for (size_t s = 0; s < hostPairs[a].size(); s++)
				{
					if (!dirtyBits[hostPairs[a][s]])
					{
						HeapIndexNode* node = heapIndex.getNode(hostPairs[a][s]);
						if (node != nullptr)
						{
							dirtyBits[hostPairs[a][s]] = true;
							heapNodes.push_back(node);
						}
					}
				}
			}
		}

		for (unsigned int b = 0; b < N; b++) {
			if (NotSquaredDistance(vertices[b], vertices[minRepresentative->v]) <= edgeDist) {
				AddOrMergeCones(Cones, vertices, b, angyx, theta);
				double  sa = angyx - theta,
						ea = angyx + theta;
				
				if (sa < 0)
					sa += 2.0 * M_PI;

				if (ea > 2.0 * M_PI)
					ea -= 2.0 * M_PI;

				certificates[b].Add(sa, ea);
				
				for (size_t s = 0; s < hostPairs[b].size(); s++)
				{
					if (!dirtyBits[hostPairs[b][s]])
					{
						HeapIndexNode* node = heapIndex.getNode(hostPairs[b][s]);
						if (node != nullptr)
						{
							dirtyBits[hostPairs[b][s]] = true;
							heapNodes.push_back(node);
						}
					}
				}
			}
		}

		HeapIndexNode* node = heapIndex.first->next;
		
		for (int x = 0; x < heapNodes.size(); x++)
		{

		//while (node != nullptr)
		//{
			int j = heapNodes[x]->index; //node->index;
			//node = node->next;
			//if (!dirtyBits[j])
			//{
				//continue;
			//}
			dirtyBits[j] = false;
			
			if (!IsInCircleSector(Cones, vertices, QContent[j]->u, QContent[j]->v) && !IsInCircleSector(Cones, vertices, QContent[j]->v, QContent[j]->u))
			{
				continue;
			}

			bool foundBadPoint = false;
			unsigned int index = 0;
			for (typename SplitTree<IsImplicit>::Iterator it = pairs[j]->first->begin(); it != pairs[j]->first->end(); ++it, index++) {
				unsigned int a = *it;
				if (!discounters[j][index])
					if (!certificates[a].Covers(pairs[j]->second->box)) {
						foundBadPoint = true;
						break;
					}
					else
						discounters[j][index] = true;
			}
			if (!foundBadPoint) {
				Q.remove(j);
				QContent[j] = nullptr;
				heapIndex.Remove(j);
				continue;
			}
			else {
				foundBadPoint = false;
				for (typename SplitTree<IsImplicit>::Iterator it = pairs[j]->second->begin(); it != pairs[j]->second->end(); ++it) {
					unsigned int b = *it;
					if (!certificates[b].Covers(pairs[j]->first->box)) {
						foundBadPoint = true;
						break;
					}
				}
				if (!foundBadPoint) {
					Q.remove(j);
					QContent[j] = nullptr;
					heapIndex.Remove(j);
					continue;
				}
				else {
					index = 0;
					for (typename SplitTree<IsImplicit>::Iterator it = pairs[j]->first->begin(); it != pairs[j]->first->end(); ++it, index++) {
						if (!discounters[j][index]) {
							unsigned int a = *it;
							bool discount = certificates[a].Covers(pairs[j]->second->box);
							discounters[j][index] = discount;
						}
					}
				}
			}

			WellSeparatedPairRepresentative* p = GapClosestPair2<Heap, IsImplicit>(j, pairs[j], vertices, theta, w, discounters, Cones);
			if (p == nullptr) {
				Q.remove(j);
				QContent[j] = nullptr;
				heapIndex.Remove(j);
			}
			else {
				Q.increaseKey(j, p->dist);
				QContent[j] = p;
			}
		}
		GapFillQueue5<Heap>(Q, QContent, pairs, vertices, index, theta, w, discounters, certificates, Cones, heapIndex);
		heapNodes.clear();
		delete minRepresentative;
	}
end:
	delete[] myEdges;
	delete splitTree;
	for (unsigned int i = 0; i < m; i++) {
		//if (!discountedPairs[i]) {
			delete[] discounters[i];
		//}
	}
	delete[] discounters;
	delete[] dirtyBits;
	delete[] bucketsLarge;
	delete[] Cones;
	delete[] hostPairs;


	return edgeCount;
}

/// ************************************************************
template<class Heap, bool IsImplicit>
void GapFillQueue5(Heap& Q,
	std::vector<WellSeparatedPairRepresentative*>& QContent,
	std::vector<WellSeparatedPair<IsImplicit>*>& pairs,
	const pointset& vertices,
	unsigned int& i,
	double t,
	double w,
	bool** discounters,
	AngleStorer* certificates,
	CircleSector_list* Cones,
	HeapIndex& heapIndex
) {
	while (i < pairs.size() && (Q.getCount() == 0 || pairs[i]->minlength() <= Q.getMin().second)) {

		bool foundBadPoint = false;
		unsigned int index = 0;
		for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->first->begin(); it != pairs[i]->first->end(); ++it, index++) {
			unsigned int a = *it;
			if (!certificates[a].Covers(pairs[i]->second->box)) {
				foundBadPoint = true;
				break;
			}
			else
				discounters[i][index] = true;
		}
		if (!foundBadPoint) {
			i++;
			continue;
		}
		else {
			foundBadPoint = false;
			for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->second->begin(); it != pairs[i]->second->end(); ++it) {
				unsigned int b = *it;
				if (!certificates[b].Covers(pairs[i]->first->box)) {
					foundBadPoint = true;
					break;
				}
			}
			if (!foundBadPoint) {
				i++;
				continue;
			}
			else {
				index = 0;
				for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->first->begin(); it != pairs[i]->first->end(); ++it, index++) {
					if (!discounters[i][index]){
						unsigned int a = *it;
						bool discount = certificates[a].Covers(pairs[i]->second->box);
						discounters[i][index] = discount;
					}
				}
			}
		}

		WellSeparatedPairRepresentative* p = GapClosestPair2<Heap, IsImplicit>(i, pairs[i], vertices, t, w, discounters, Cones);
		//if (i % 10000000 == 0)
		//std::cout << i << "/" << pairs.size() << " " << Q.getCount() << " " << edgeCount << endl;
		if (p != nullptr) {
			Q.insert(i, p->dist);
			QContent[i] = p;
			heapIndex.Add(i);
		}
		i++;
	}
}

template<class Heap, bool IsImplicit>
WellSeparatedPairRepresentative* GapClosestPair5(unsigned int pairNo,
	WellSeparatedPair<IsImplicit>* pair,
	const pointset& vertices,
	double t,
	double w,
	bool** discounters,
	CircleSector_list* Cones,
	std::vector<std::pair<unsigned int, std::pair<unsigned int, double>>>* candidatePairs
) {
	unsigned int j = 0;

	for (typename SplitTree<IsImplicit>::Iterator it = pair->first->begin(); it != pair->first->end(); ++it, j++) {
		unsigned int indexj = *it;
		if (discounters[pairNo][j] == false) {
			for (typename SplitTree<IsImplicit>::Iterator it1 = pair->second->begin(); it1 != pair->second->end(); ++it1) {
				unsigned int indexk = *it1;

				if (!IsInCircleSector(Cones, vertices, indexj, indexk) && !IsInCircleSector(Cones, vertices, indexk, indexj)) {
					double   dist = NotSquaredDistance(vertices[indexj], vertices[indexk]);
					candidatePairs[pairNo].push_back(std::pair<unsigned int, std::pair<unsigned int, double>>
						(indexj, std::pair<unsigned int, double>(indexk, dist)));
				}
			}
		}
	}

	if (candidatePairs[pairNo].size() == 0)//(isinf(bestDist))
		return nullptr;
	else
	{
		std::sort(candidatePairs[pairNo].begin(), candidatePairs[pairNo].end(), PointPairComparer<IsImplicit>());
		//cout << "pairNo:" << pairNo << "\t" << pair->minlength() <<"\t"<< candidatePairs[pairNo].size() << std::endl;
		return new WellSeparatedPairRepresentative(candidatePairs[pairNo][0].first,
			candidatePairs[pairNo][0].second.first, sqrt(candidatePairs[pairNo][0].second.second), pairNo);
	}
}

template<class Heap, bool IsImplicit>
WellSeparatedPairRepresentative* GapClosestPair6(unsigned int pairNo,
	WellSeparatedPair<IsImplicit>* pair,
	const pointset& vertices,
	double t,
	double w,
	bool** discounters,
	CircleSector_list* Cones,
	std::vector<std::pair<unsigned int, std::pair<unsigned int, double>>>* candidatePairs
) {
	int shortestPairInedx1, shortestPairInedx2;
	double bestDist = std::numeric_limits<double>::infinity();
	//std::vector<std::pair<unsigned int, std::pair<unsigned int, double>>> newCandidatePairs;
	unsigned int i;

	for (i = 0; i < candidatePairs[pairNo].size(); i++)
	{
		if (!IsInCircleSector(Cones, vertices, candidatePairs[pairNo][i].first, candidatePairs[pairNo][i].second.first)
			&& !IsInCircleSector(Cones, vertices, candidatePairs[pairNo][i].second.first, candidatePairs[pairNo][i].first))
		{
			shortestPairInedx1 = candidatePairs[pairNo][i].first;
			shortestPairInedx2 = candidatePairs[pairNo][i].second.first;
			bestDist = candidatePairs[pairNo][i].second.second;;
			break;
		}
	}

	if (std::isinf(bestDist))
	{
		candidatePairs[pairNo].clear();
		return nullptr;
	}
	else
	{
		candidatePairs[pairNo].erase(candidatePairs[pairNo].begin(), candidatePairs[pairNo].begin() + i);
		//cout << "pairNo:" << pairNo << "\t" << candidatePairs[pairNo].size() << std::endl;
		return new WellSeparatedPairRepresentative(shortestPairInedx1, shortestPairInedx2, sqrt(bestDist), pairNo);
	}
}

template<class Heap, bool IsImplicit>
void GapUpdatePair5(unsigned int j,
	Heap& Q,
	std::vector<WellSeparatedPairRepresentative*>& QContent,
	WellSeparatedPair<IsImplicit>* pair,
	const pointset& vertices,
	std::vector< std::pair<unsigned int, double> >* myEdges,
	double t,
	double w,
	bool* dirtyBits,
	bool** discounters,
	AngleStorer* certificates,
	CircleSector_list* Cones,
	std::vector<std::pair<unsigned int, std::pair<unsigned int, double>>>* candidatePairs
) {
	//cout << "update("<<j<<")\t";
	WellSeparatedPairRepresentative* p = GapClosestPair6<Heap, IsImplicit>(j, pair, vertices, t, w, discounters, Cones, candidatePairs);
	if (p == nullptr) {
		Q.remove(j);
		QContent[j] = nullptr;
	}
	else {
		Q.increaseKey(j, p->dist);
		QContent[j] = p;
	}
	dirtyBits[j] = false;
}

template<class Heap, bool IsImplicit>
void GapReinsertPair5(unsigned int j,
	Heap& Q,
	std::vector<WellSeparatedPairRepresentative*>& QContent,
	WellSeparatedPair<IsImplicit>* pair,
	const pointset& vertices,
	std::vector< std::pair<unsigned int, double> >* myEdges,
	double t,
	double w,
	bool* dirtyBits,
	bool** discounters,
	AngleStorer* certificates,
	CircleSector_list* Cones,
	std::vector<std::pair<unsigned int, std::pair<unsigned int, double>>>* candidatePairs
) {
	//cout << "Reinst("<<j<<")\t";
	WellSeparatedPairRepresentative* p = GapClosestPair6<Heap, IsImplicit>(j, pair, vertices, t, w, discounters, Cones, candidatePairs);
	if (p != nullptr) {
		Q.insert(j, p->dist);
		QContent[j] = p;
	}

	dirtyBits[j] = false;
}


template<bool IsImplicit>
struct PointPairComparer {
	inline bool operator()(std::pair<unsigned int, std::pair<unsigned int, double>> const i, 
		std::pair<unsigned int, std::pair<unsigned int, double>> const j) const {
		
		return i.second.second < j.second.second;
	}
};

//***************************************************************
/*
الگوریتم زیر، استفاده از روش ایجاد مخروط عجولانه 
را در الگوریتم اصلی ساخت پوشش شکاف-حریصانه 
بکار خواهد گرفت
*/
//***************************************************************
template<class Heap, bool IsImplicit>
int GapGreedyOriginalWithCones(BaseTask* T) {
	STARTLOG
		const pointset& vertices = T->p;
	unsigned int N = vertices.size();
	unsigned int edgeCount = 0;
	double w = T->var, theta = T->t;

	std::vector<EdgeInfo> edgeList;
	edgeList.reserve(((N - 1) * N)/2);
	for (unsigned int x = 0; x < N; x++) {
		for (unsigned int y = x + 1; y < N; y++) {
			EdgeInfo e(x, y, distance(vertices[x], vertices[y]));
			edgeList.push_back(e);
		}
	}
	
	std::sort(edgeList.begin(), edgeList.end());

	std::vector<std::pair<unsigned int, unsigned int>> edges;
	CircleSector_list* Cones = new CircleSector_list[N];
	//CircleSector_list* ConesRev = new CircleSector_list[N];
	//AngleStorer* certificates = new AngleStorer[N];
	//AngleStorer* certificatesRev = new AngleStorer[N];

	

	for (unsigned int i = 0; i < edgeList.size(); i++)
	{
		if (!IsInCircleSector(Cones, vertices, edgeList[i].x, edgeList[i].y) &&
			!IsInCircleSector(Cones, vertices, edgeList[i].y, edgeList[i].x))
		{
			std::pair<unsigned int, unsigned int> e1(edgeList[i].x, edgeList[i].y);
			edges.push_back(e1);

			std::pair<unsigned int, unsigned int> e2(edgeList[i].y, edgeList[i].x);
			edges.push_back(e2);

			if (T->create_edge(edge(edgeList[i].x, edgeList[i].y))) return 0;
			if (T->create_edge(edge(edgeList[i].y, edgeList[i].x))) return 0;
			edgeCount += 2;

			double edgeDist = w * distance(vertices[edgeList[i].x], vertices[edgeList[i].y]);
			double  angxy = atan3(vertices[edgeList[i].y], vertices[edgeList[i].x]);
			double  angyx = atan3(vertices[edgeList[i].x], vertices[edgeList[i].y]);
			double edgeDistPow2 = edgeDist * edgeDist;

			for (unsigned int a = 0; a < N; a++) {
				if (NotSquaredDistance(vertices[a], vertices[edgeList[i].x]) <= edgeDistPow2) {
					AddOrMergeCones(Cones, vertices, a, angxy, theta);
				}
				if (NotSquaredDistance(vertices[a], vertices[edgeList[i].y]) <= edgeDistPow2) {
					AddOrMergeCones(Cones, vertices, a, angyx, theta);
				}
			}
		}
	}
	
	char filename[80];
	sprintf(filename, "IGG(%d-%f-%f).txt", N, theta, w);
	std::ofstream myfile(filename);
	PROCESS_MEMORY_COUNTERS_EX pmc;
	GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));

	//SIZE_T physMemUsedByMe = pmc.WorkingSetSize;
	myfile << pmc.WorkingSetSize;

	myfile.close();
	delete[] Cones;
	edgeList.clear();
	edges.clear();
//	delete[] ConesRev;

	return edgeCount;
}

//**********************************************************
int k;

template<class Heap, bool IsImplicit>
int BaseGapSpannerBucketsFixed(BaseTask* T) {
	//TODO:
	// Make ESA' local
	STARTLOG
		const pointset& vertices = T->p;
	unsigned int N = vertices.size();
	std::vector<std::pair<unsigned int, double> >* myEdges = new std::vector<std::pair<unsigned int, double> >[N];
	unsigned int edgeCount = 0;

	
	// Bucketing stuff
	//double w_ratio = T->t;
	double theta = T->t;
	double w = T->var;
	double t = 1.0 / (cos(theta) - sin(theta) - 2 * w);
	//double gapGreedyWeightRatio =  1 / w + 1 / 2 * theta;//(log((double)N)* log((double)N) / sqrt((double)N) / log(log((double)N)) / sqrt(sqrt(t - 1.0)));//(1 + 2 / w) * log((double)N) / (1 / pow((t - 1), 4));

	//double estimatedLongestEdgeLengthUnit = ((log((double)N)* gapGreedyWeightRatio) / sqrt((double)N) / log(log((double)N)) / sqrt(sqrt(t - 1.0)));
	double tavan = pow(2, k);

	double estimatedLongestEdgeLengthUnit = log((double)N) /sqrt((double)N) / pow(t - 1.0, tavan);

	if (estimatedLongestEdgeLengthUnit > 1)
		estimatedLongestEdgeLengthUnit = 1;

	double minX = vertices[0].x, minY = vertices[0].y, maxX = vertices[0].x, maxY = vertices[0].y;
	for (int i = 1; i < N; i++) {
		if (vertices[i].x < minX) minX = vertices[i].x;
		if (vertices[i].x > maxX) maxX = vertices[i].x;
		if (vertices[i].y < minY) minY = vertices[i].y;
		if (vertices[i].y > maxY) maxY = vertices[i].y;
	}
	double xwidth = maxX - minX;
	double ywidth = maxY - minY;
	double maxwidth = xwidth > ywidth ? xwidth : ywidth;
	double estimatedLongestEdgeLength = estimatedLongestEdgeLengthUnit * maxwidth;
	double estimatedLongestTPath = estimatedLongestEdgeLength * (1 + t / 2.0);

	/*cout << "log N: " << std::fixed << log((double)N)  << endl;
	cout << "Sqrt N: " << std::fixed << sqrt((double)N) << endl;
	cout << "log log N: " << std::fixed << log(log((double)N)) << endl;
	cout << "sqrt(sqrt(t - 1.0)): " << std::fixed << sqrt(sqrt(t - 1.0)) << endl;

	cout << "estimatedLongestEdgeLengthUnit: " << std::fixed << estimatedLongestEdgeLengthUnit << endl;
	cout << "estimatedLongestEdgeLength: " << std::fixed << estimatedLongestEdgeLength << endl;
	cout << "estimatedLongestTPath: " << std::fixed << estimatedLongestTPath << endl;*/

	int bucketeers = (int)(1.0 / (estimatedLongestEdgeLengthUnit * (1 + t / 2.0)));

	//cout << "xwidth: " << std::fixed << xwidth << endl;
	//cout << "ywidth: " << std::fixed << ywidth << endl;

	if (bucketeers < 1) { bucketeers = 1; }
	//cout << "bucketSide: " << bucketeers << endl;

	double bucketWidth = xwidth / ((double)bucketeers) + 1e-6;
	double bucketHeight = ywidth / ((double)bucketeers) + 1e-6;
	if (bucketWidth < 1.0) bucketWidth = 1e-6;
	if (bucketHeight < 1.0) bucketHeight = 1e-6;
	//cout << "bucketWidth: " << std::fixed << bucketWidth << endl;
	//cout << "bucketHeight: " << std::fixed << bucketHeight << endl;

	std::vector<int>* buckets = new std::vector<int>[bucketeers * bucketeers];
	int* bucketLookupX = new int[N];
	int* bucketLookupY = new int[N];
	//int* bucketLookupZ = new int[N];
	for (int i = 0; i < N; i++) {
		int bucketX = (int)((vertices[i].x - minX) / bucketWidth);
		int bucketY = (int)((vertices[i].y - minY) / bucketHeight);
		int bucketCoord = bucketX + bucketeers * bucketY;
		//bucketLookupZ[i] = buckets[bucketCoord].size();
		buckets[bucketCoord].push_back(i);
		bucketLookupX[i] = bucketX;
		bucketLookupY[i] = bucketY;
	}

	std::vector<int>* bucketsLarge = new std::vector<int>[bucketeers * bucketeers];
	for (int i = 0; i < N; i++) {
		int bucketX = bucketLookupX[i];
		int bucketY = bucketLookupY[i];
		int bucketCoord = bucketX + bucketeers * bucketY;
		for (int a = -1; a <= 1; a++) {
			for (int b = -1; b <= 1; b++) {
				int bucketLargeIndex = bucketX + a + bucketeers * (bucketY + b);
				if (bucketX + a >= 0 && bucketX + a < bucketeers
					&&
					bucketY + b >= 0 && bucketY + b < bucketeers) {
					bucketsLarge[bucketLargeIndex].push_back(i);
				}
			}
		}
	}

	CircleSector_list* Cones = new CircleSector_list[N];

	// Call base algorithm
	SIZE_T memUsed = ComputeGapSpannerBuckets<Heap>(T, edgeCount, myEdges, t, buckets, estimatedLongestEdgeLength, estimatedLongestTPath, bucketLookupX, bucketLookupY,
		bucketeers, Cones, bucketsLarge);

	int shortEdges = edgeCount;
	//cout << "Short Edges: " << edgeCount << endl;
	//cout << Timing::elapsed() << endl;

	AngleStorer* certificates = new AngleStorer[N];
	for (int i = 0; i < N; i++) {
		certificates[i].origin = vertices[i];

		CircleSector_list_node* v_item = Cones[i].first;
		CircleSector* cone;
		while (v_item != nullptr)
		{
			cone = v_item->data;
			certificates[i].Add(cone->startAngle, cone->endAngle);
			v_item = v_item->next;
		}
	}

	/*delete[] buckets;
	delete[] bucketLookupX;
	delete[] bucketLookupY;*/


		// Compute WSPD
		//double s = maximum((double)(4.0 / sin(theta / 2.0)), (double)(2.0 / w)); //4.0 * t / (t - 1.0);//5.0;
		double s = 5.0;
		SplitTree<IsImplicit>* splitTree = SplitTree<IsImplicit>::createSplitTree(vertices);
		std::vector<WellSeparatedPair<IsImplicit>*> pairs;
		splitTree->ComputeWspdPairs(s, pairs, estimatedLongestEdgeLength);
		//cout << "Long pairs: " << pairs.size() << endl;

		std::vector<WellSeparatedPair<IsImplicit>*> pairs1;
		splitTree->ComputeWspdPairs(s, pairs1);

		// Discount pairs
		bool* discountedPairs = new bool[pairs.size()];
		bool** discounters = new bool* [pairs.size()];
		unsigned int discountedPairsCounter = 0;
		unsigned int discountersCounter = 0;
		unsigned int discountersTotalCounter = 0;
		for (unsigned int i = 0; i < pairs.size(); i++) {
			discountedPairs[i] = false;
			bool foundBadPoint = false;
			for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->first->begin(); it != pairs[i]->first->end(); ++it) {
				unsigned int a = *it;
				if (!certificates[a].Covers(pairs[i]->second->box)) {
					foundBadPoint = true;
					break;
				}
			}
			if (!foundBadPoint) {
				discountedPairs[i] = true;
				discountedPairsCounter++;
			}
			else {
				foundBadPoint = false;
				for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->second->begin(); it != pairs[i]->second->end(); ++it) {
					unsigned int b = *it;
					if (!certificates[b].Covers(pairs[i]->first->box)) {
						foundBadPoint = true;
						break;
					}
				}
				if (!foundBadPoint) {
					discountedPairs[i] = true;
					discountedPairsCounter++;
				}
				else {
					discounters[i] = new bool[pairs[i]->first->count()];
					discountersTotalCounter += pairs[i]->first->count();
					unsigned int index = 0;
					for (typename SplitTree<IsImplicit>::Iterator it = pairs[i]->first->begin(); it != pairs[i]->first->end(); ++it, index++) {
						unsigned int a = *it;
						bool discount = certificates[a].Covers(pairs[i]->second->box);
						discounters[i][index] = discount;
						if (discount) discountersCounter++;
					}
				}
			}
		}
		//cout << "Discounted " << discountedPairsCounter << "/" << pairs.size() << endl;
		//cout << "Discounted " << discountersCounter << "/" << discountersTotalCounter << endl;
		//char filename[80];
		//sprintf(filename, "NMNM(%d-%f-%f).txt", N, theta, w);
		//std::ofstream myfile(filename, std::ofstream::out);

		unsigned int m = pairs.size();
		std::vector<WellSeparatedPairRepresentative*> QContent(m, nullptr);
		Heap Q(m, -std::numeric_limits<double>::infinity());
		bool* dirtyBits = new bool[m];
		for (unsigned int i = 0; i < m; i++)
			dirtyBits[i] = false;

		GapFillQueue3<Heap>(Q, QContent, pairs, vertices, myEdges, t, w, edgeCount, discountedPairs, dirtyBits, discounters, certificates, Cones);

		while (Q.getCount() > 0) {
			bool foundOne = false;
			std::pair<unsigned int, double> minQ;
			while (!foundOne) {
				minQ = Q.getMin();

				if (dirtyBits[minQ.first])
					GapUpdatePair2(minQ.first, Q, QContent, pairs[minQ.first], vertices, myEdges, t, w, dirtyBits, discounters, certificates, Cones);
				else
					foundOne = true;
				if (Q.getCount() == 0)
					goto end;
			}
			Q.extractMin();
			WellSeparatedPairRepresentative* minRepresentative = QContent[minQ.first];

			std::pair<unsigned int, double> pair(minRepresentative->u, minRepresentative->dist);
			myEdges[minRepresentative->v].push_back(pair);
			pair.first = minRepresentative->v;
			myEdges[minRepresentative->u].push_back(pair);
			if (T->create_edge(edge(minRepresentative->u, minRepresentative->v))) return 0;
			if (T->create_edge(edge(minRepresentative->v, minRepresentative->u))) return 0;
			edgeCount += 2;


			/*if (minRepresentative->v > minRepresentative->u)
			{
				myfile << "-" << minRepresentative->v << "&" << minRepresentative->u << "-" << "(" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ") ---> (" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ")" << endl;
				myfile << "-" << minRepresentative->u << "&" << minRepresentative->v << "-" << "(" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ") ---> (" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ")" << endl;
			}
			else
			{
				myfile << "-" << minRepresentative->u << "&" << minRepresentative->v << "-" << "(" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ") ---> (" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ")" << endl;
				myfile << "-" << minRepresentative->v << "&" << minRepresentative->u << "-" << "(" << vertices[minRepresentative->v].x << ", " << vertices[minRepresentative->v].y << ") ---> (" << vertices[minRepresentative->u].x << ", " << vertices[minRepresentative->u].y << ")" << endl;
			}*/


			double edgeDist = w * distance(vertices[minRepresentative->u], vertices[minRepresentative->v]);
			double  angyx = atan3(vertices[minRepresentative->u], vertices[minRepresentative->v]);
			double  angxy = atan3(vertices[minRepresentative->v], vertices[minRepresentative->u]);

			for (unsigned int a = 0; a < N; a++) {
				if (distance(vertices[a], vertices[minRepresentative->u]) <= edgeDist) {
					AddOrMergeCones(Cones, vertices, a, angxy, theta);
				}
			}

			for (unsigned int b = 0; b < N; b++) {
				if (distance(vertices[b], vertices[minRepresentative->v]) <= edgeDist) {
					AddOrMergeCones(Cones, vertices, b, angyx, theta);
				}
			}

			for (int j = 0; j < m; j++)
				if (QContent[j] != nullptr)
					dirtyBits[j] = true;

			GapReinsertPair2(minQ.first, Q, QContent, pairs[minQ.first], vertices, myEdges, t, w, dirtyBits, discounters, certificates, Cones);
		}
	end:

		char filename[80];
		sprintf(filename, "IAM%d(%d-%f-%f).txt", k + 3, N, theta, w);
		std::ofstream myfile(filename, std::ofstream::app);

		PROCESS_MEMORY_COUNTERS_EX pmc;
		GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
		SIZE_T physMemUsedByMe = pmc.WorkingSetSize;

		if (memUsed > physMemUsedByMe)
			physMemUsedByMe = memUsed;

		myfile << physMemUsedByMe;
		myfile << "\t" << pairs.size() << "\t" << discountedPairsCounter << "\t" << pairs1.size() << endl;
		myfile.close();


		delete[] buckets;
		delete[] bucketLookupX;
		delete[] bucketLookupY;

		delete[] myEdges;
		delete splitTree;
		for (unsigned int i = 0; i < m; i++) {
			if (!discountedPairs[i]) {
				delete[] discounters[i];
			}
		}
		delete[] discounters;
		delete[] discountedPairs;
		delete[] dirtyBits;

	
	delete[] bucketsLarge;
	delete[] Cones;

	/*
	cout <<endl<< "Long Pairs%\t\tDiscounted\t\t\tShort Edges" << endl;
	cout <<  (double)pairs.size() / (double)pairs1.size()<<"\t\t"
	<< discountedPairsCounter << "/" << pairs.size() << ", %=" << (double)discountedPairsCounter / (double)pairs.size() << "\t\t"
	<< (double)shortEdges / (double)edgeCount <<endl << endl;*/
	//cout << "Discounted " << discountersCounter << "/" << discountersTotalCounter << ", %=" << (double)discountersCounter / (double)discountersTotalCounter << endl;

	

	/*myfile << shortEdges << "\t" << edgeCount << "\t"
		<< (double)pairs.size() / (double)pairs1.size() << "\t"
		<< discountedPairsCounter << "/" << pairs.size() << ", %=" << (double)discountedPairsCounter / (double)pairs.size() << "\t"
		<< (double)shortEdges / (double)edgeCount << endl;*/

	myfile.close();

	return edgeCount;
}


template<class Heap, bool IsImplicit>
int DistributionSensitiveGapSpanner1(BaseTask* T) {
	k = -1;
	return BaseGapSpannerBucketsFixed<BinHeap<double>, false>(T);
	
}

template<class Heap, bool IsImplicit>
int DistributionSensitiveGapSpanner2(BaseTask* T) {
	k = 0;
	return BaseGapSpannerBucketsFixed<BinHeap<double>, false>(T);

}



#pragma endregion

#pragma endregion

#pragma endregion
//******************************************

Solver::type Solver::solvers[] = {
	// 0
    {"Greedy BinHeap Dijkstra", GreedySpanner<Dijkstra<BinHeap<double> > > },
	// 1
    {"Greedy BinHeap A-Star", GreedySpanner<AStar<BinHeap<double> > >},
	// 2
    {"Greedy ED", GreedySpannerOnCompleteGraph<EDGreedySpanner>},
	// 3
	{"Greedy approx", GreedySpannerOnCompleteGraph<ApproxGreedySpanner>},
	// 4
    {"Greedy Linspace", GreedyLinspace<BinHeap<double>, false>},
	// 5
    {"Greedy Linspace Implicit", GreedyLinspace<BinHeap<double>, true>},
	// 6 (Lazy Greedy Linspace - WSPD Greedy Lazy-)
    {"Greedy Linspace 2", GreedyLinspace2<BinHeap<double>, false>},
	// 7
    {"Greedy Linspace 2 Implicit", GreedyLinspace2<BinHeap<double>, true>},
	// 8 (LazyGreedy)
    {"Greedy Linspace 3", GreedyLinspace3<BinHeap<double> >},
	// 9
    {"Greedy FG BinHeap", GreedySpannerOnCompleteGraph<FGGreedySpanner<BinHeap<double> > >},
	// 10
    {"Greedy FGD BinHeap", GreedySpannerOnCompleteGraph<FGDGreedySpanner<BinHeap<double> > >},
	// 11
    {"Greedy preliminary", GreedySpannerPrelim<AStar<BinHeap<double> >, BinHeap<double> >},
	// 12
    {"Greedy new", GreedySpannerNew<BinHeap<double> >},
	// 13
    {"Greedy buckets", GreedySpannerBuckets<GreedySpannerBucketsWrapper<BinHeap<double> > >},
	// 14
    {"Greedy buckets 2", GreedySpannerBuckets<GreedySpannerBuckets2Wrapper<BinHeap<double> > >},
	// 15
    {"Greedy buckets fixed", GreedySpannerBucketsFixed<GreedySpannerBucketsWrapper<BinHeap<double> >, BinHeap<double>, false >},
	// 16
    {"Greedy buckets fixed Implicit", GreedySpannerBucketsFixed<GreedySpannerBucketsWrapper<BinHeap<double> >, BinHeap<double>, true >},
	// 17
    {"Greedy buckets 2 fixed", GreedySpannerBucketsFixed<GreedySpannerBuckets2Wrapper<BinHeap<double> >, BinHeap<double>, false >},
	// 18
    {"Greedy buckets 2 fixed Implicit", GreedySpannerBucketsFixed<GreedySpannerBuckets2Wrapper<BinHeap<double> >, BinHeap<double>, true >},
	// 19
    {"Hyperbola spanner", HyperbolaSpanner},
	// 20
    {"Hyperbola spanner fast", HyperbolaSpannerFast},
	// 21
    {"Symmetric hyperbola spanner", SymmetricHyperbolaSpanner},
	// 22
    {"Symmetric hyperbola spanner fast", SymmetricHyperbolaSpannerFast<BinHeap<double> >},
	// 23
    {"Theta Spanner", ThetaSpanner},
	// 24
    {"WSPD Spanner", WspdSpanner<false>},
	// 25
    {"WSPD Spanner Implicit", WspdSpanner<true>},
	// 26
    {"Hyperbole KD Spanner", HyperbolaKDSpanner},
	// 27
    {"Hybrid Theta Heur Spanner", HybridSpanner<ThetaEDHybrid>},
	// 28
    {"Hybrid Theta FG Spanner", HybridSpanner<ThetaFGHybrid<BinHeap<double> > >},
	// 29
	{"Hybrid Theta SymHyp Spanner", HybridSpanner<ThetaSymHypHybrid>},
	// 30
    {"Hybrid WSPD Heur Spanner", HybridSpanner<WspdEDHybrid<false> >},
	// 31
    {"Hybrid WSPD Heur Spanner Implicit", HybridSpanner<WspdEDHybrid<true> >},
	// 32
    {"Hybrid WSPD FG Spanner", HybridSpanner<WspdFGHybrid<BinHeap<double>, false > >},
	// 33
    {"Hybrid WSPD FG Spanner Implicit", HybridSpanner<WspdFGHybrid<BinHeap<double>, true > >},
	// 34
    {"Hybrid WSPD SymHyp Spanner", HybridSpanner<WspdSymHypHybrid<false > >},/*wspd!!*/
	// 35
    {"Hybrid WSPD SymHyp Spanner Implicit", HybridSpanner<WspdSymHypHybrid<true > >},
	// 36
	{"Theta SymHyp Sparsification Spanner", ThetaSymHypSparsificationSpanner},
	// 37
	{"Theta SymHyp Sparsification Spanner 2", ThetaSymHypSparsificationSpanner2},
	// 38
	{"Theta SymHyp Sparsification Spanner 3", ThetaSymHypSparsificationSpanner3},
	// 39
	{"WSPD Hyp Sparsification Spanner", WspdHypSparsificationSpanner<false>},
	// 40
	{"WSPD Hyp Sparsification Spanner Implicit", WspdHypSparsificationSpanner<true>},
	// 41
	{"WSPD SymHyp Sparsification Spanner", WspdSymHypSparsificationSpanner<false>},
	// 42
	{"WSPD SymHyp Sparsification Spanner Implicit", WspdSymHypSparsificationSpanner<true>},
	// 43
	{"WSPD SymHyp Sparsification Spanner 2", WspdSymHypSparsificationSpanner2<false>},
	// 44
	{"WSPD SymHyp Sparsification Spanner 2 Implicit", WspdSymHypSparsificationSpanner2<true>},
	// 45
    {"Yao Spanner", YaoSpanner},
	// 46
    {"Delaunay", DelaunaySpanner},
	// 47
    {"From file", FromFile},
	//48
	{ "IFGBN1", ComputeGreedySpannerIFGBN1<BinHeap<double>> },
	//49
	{ "IFGBN2", ComputeGreedySpannerIFGBN2<BinHeap<double>> },
	//50
	{ "IFGBN", ComputeGreedySpannerIFGBN<BinHeap<double>> },
	//51
	{ "IFGBC",ComputeGreedySpannerIFGBC<BinHeap<double>> },
	//52
	{ "Delta", DeltaGreedy<BinHeap<double>> },
		// 53 BakhsheshMethod1
	{ "BM1", BakhsheshMethod1<BinHeap<double>, false> },
	//54 BakhsheshMethod2
	{ "BM2", BakhsheshMethod2<BinHeap<double>, false> },
		// 55
	{ "ImprovedMethod1BakhsheshWithChangingCloseCriteria", ImprovedMethod1BakhsheshWithChangingCloseCriteria<BinHeap<double>, false> },
		// 56
	{ "MyGapGreedyLinspace", MyGapGreedyLinspace<BinHeap<double>, false> },
		// 57 GapGreedyOriginal
	{ "GapOrig", GapGreedyOriginal<BinHeap<double>, false> },
		// 58
	{ "AG_FC", NewAlphaGapGreedyFirstCandidate<BinHeap<double>, false> },
		// 59
	{ "AG_BC", NewAlphaGapGreedyBestCandidate<BinHeap<double>, false> },
		// 60
	{ "LAG_FC", NewLinearSpaceAlphaGapGreedyFirstCandidate<BinHeap<double>, false> },
		// 61
	{ "LAG_BC", NewLinearSpaceAlphaGapGreedyBestCandidate<BinHeap<double>, false> },
		// 62
	{ "MNM", MNM_Method<BinHeap<double>, false> },
		// 63 اطباق مستقیم روش بوچین برای ساخت پوشش شکاف حریصانه
	{ "AM", GapSpannerBucketsFixed<BinHeap<double>, false> },
		//64 (GapOrig+Cone)
	{ "IGG", GapGreedyOriginalWithCones<BinHeap<double>, false> },
		//65
	{ "MyMNM", GapSpannerBucketsFixed1<BinHeap<double>, false> },
		//66
	{ "OIAM1", GapSpannerBucketsFixed2<BinHeap<double>, false> },
		//67
	{ "OIAM2", GapSpannerBucketsFixed3<BinHeap<double>, false> },
		//68
		// نسخه انطباق مستقیم روش بوچین برای سخت پوشش شکاف حریصانه با این تفاوت که  
		// از وزن پوشش شکاف برای تعیین اندازه بزرگترین یال در پوشش استفاده کرده ایم
		// در واقع با توجه به اینکه وزن پوشش شکاف، 
		//log n
		//برابر وزن پوشش مسیر حریصانه است، لذا یک تخمین برای اندازه کبزرگترین یال ممکن در شکاف حریصانه نیز
		//از همین راه بدست می آِد
	{ "IAM1", DistributionSensitiveGapSpanner1<BinHeap<double>, false> },
	
		//69
		//همان الگوریتم بالایی فقط اندازه بزرگترین یال احتمالی، با تغییراتی در فرمول قبل، بدست می آِید
	{ "IAM2", DistributionSensitiveGapSpanner2<BinHeap<double>, false> },


    {nullptr,nullptr}, // Last element of list.
};

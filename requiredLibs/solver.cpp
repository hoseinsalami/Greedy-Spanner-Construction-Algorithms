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

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

using std::cout;
using std::endl;

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
void ModifiedUpdateDistances(pointset vertices, std::vector< std::pair<unsigned int, double> >* myEdges, unsigned int from, double* matrix, double t) {
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
			delete right;
		}
		if (isTop)
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
    cout << "estimatedLongestEdgeLength: " << std::fixed << estimatedLongestEdgeLength << endl;
    cout << "estimatedLongestTPath: " << std::fixed << estimatedLongestTPath << endl;
	int bucketeers = (int)(1.0 / (estimatedLongestEdgeLengthUnit * (1 + T->t / 2.0)));
	if(bucketeers<1){bucketeers=1;}
    cout << "bucketSide: " << bucketeers << endl;
	double bucketWidth = xwidth / ((double)bucketeers) + 1e-6;
	double bucketHeight = ywidth / ((double)bucketeers) + 1e-6;
	if (bucketWidth < 1.0) bucketWidth = 1e-6;
	if (bucketHeight < 1.0) bucketHeight = 1e-6;
    cout << "bucketWidth: " << std::fixed << bucketWidth << endl;
    cout << "bucketHeight: " << std::fixed << bucketHeight << endl;

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
    cout << "estimatedLongestEdgeLength: " << std::fixed << estimatedLongestEdgeLength << endl;
    cout << "estimatedLongestTPath: " << std::fixed << estimatedLongestTPath << endl;
	int bucketeers = (int)(1.0 / (estimatedLongestEdgeLengthUnit * (1 + t / 2.0)));
	if(bucketeers<1){bucketeers=1;}
    cout << "bucketSide: " << bucketeers << endl;
	double bucketWidth = xwidth / ((double)bucketeers) + 1e-6;
	double bucketHeight = ywidth / ((double)bucketeers) + 1e-6;
	if (bucketWidth < 1.0) bucketWidth = 1e-6;
	if (bucketHeight < 1.0) bucketHeight = 1e-6;
    cout << "bucketWidth: " << std::fixed << bucketWidth << endl;
    cout << "bucketHeight: " << std::fixed << bucketHeight << endl;

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
	cout << edgeCount << endl;
	cout << Timing::elapsed() << endl;

	// Find bridging points
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
	cout << "Long pairs: " << pairs.size() << endl;

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
	cout << "Discounted " << discountedPairsCounter << "/" << pairs.size() << endl;
	cout << "Discounted " << discountersCounter << "/" << discountersTotalCounter << endl;

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
	std::cout << "Cleaning up..." << std::endl;
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
double shortest_T_Path(std::vector<std::pair<unsigned int, double>>* myEdges, unsigned int N, double* matrix, int source, int target)
{
	double nearestNeighborToTargetDistance = INFINITY, d;

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

					}
				}
				else
				{
					std::pair<unsigned int, double> pair(edgeList[k].x, edgeList[k].distance);
					myEdges[edgeList[k].y].push_back(pair);
					pair.first = edgeList[k].y;
					myEdges[edgeList[k].x].push_back(pair);
					if (T->create_edge(edgeList[k])) { goto cleanup; }

					matrix[edgeList[k].y + N * edgeList[k].x] = matrix[edgeList[k].x + N * edgeList[k].y] = edgeList[k].distance;
				}
			}
		}
		else
			matrix[edgeList[k].y + N * edgeList[k].x] = matrix[edgeList[k].x + N * edgeList[k].y] = dist12<dist21 ? dist12 : dist21;
	}
cleanup:
	ENDLOG
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
const long double PI = 3.1415926535897932384626433832L;
class Cone
{
public:
	vertex apex;
	vertex bisector;
	vertex start;
	vertex end;
	double	angle;

	Cone(vertex& ap, vertex& bs, double ang)
	{
		apex = ap;
		bisector = bs;
		angle = ang;

		double  angleToRadian = (PI / 180.0) * angle;
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
void	AddCone(cone_list* list, pointset &vertices, int apex, int bisector, double	angle)
{
	Cone	*cone = new Cone(vertices[apex], vertices[bisector], angle);

	list[apex].add(cone);
}
//***************************************
bool IsInCone(cone_list*	Cones, pointset &vertices, int apex, int p)
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

int DeltaGreedy(BaseTask * T)
{
	STARTLOG
	pointset vertices = T->p;
	cone_list*	Cones = new cone_list[T->p.size()];
	double d, g, asn, tetha, delta = T->t;
	double * dist = new double[T->p.size()];
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

	for (int i = 0; i<(N*(N - 1)) / 2; i++)
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
				if (T->create_edge(edgeList[i])) {  goto cleanup; }
				d = 1;
			}
			g = d / (sqrt((double)2) * T->t);
			asn = asin((double)g) * (180.0 / PI);
			tetha = 45 - asn;

			AddCone(Cones, vertices, edgeList[i].x, edgeList[i].y, 2 * tetha);
			AddCone(Cones, vertices, edgeList[i].y, edgeList[i].x, 2 * tetha);
		}
	}
cleanup:
	ENDLOG

	delete[]dist;
	delete[]Cones;
	return c;
}
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
	// 6
    {"Greedy Linspace 2", GreedyLinspace2<BinHeap<double>, false>},
	// 7
    {"Greedy Linspace 2 Implicit", GreedyLinspace2<BinHeap<double>, true>},
	// 8
    {"Greedy Linspace 3", GreedyLinspace3<BinHeap<double> >},
	// 9
    {"Greedy FG BinHeap", GreedySpannerOnCompleteGraph<FGDGreedySpanner<BinHeap<double> > >},
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
    {"Hybrid WSPD SymHyp Spanner", HybridSpanner<WspdSymHypHybrid<false > >},
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

    {nullptr,nullptr}, // Last element of list.
};

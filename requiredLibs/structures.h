#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>
#define nullptr NULL //tmp fix omdat ik hier nog een oude gcc heb :P

struct vertex {
	double x,y;
	vertex() : x(0), y(0) {}
	vertex(double x, double y) : x(x), y(y) {}
	bool operator==(const vertex &v) const {return (v.x==x && v.y==y);}
};

struct edge {
	unsigned int x,y;
	edge() : x(0), y(0) {}
	edge(unsigned int x, unsigned int y) : x(x), y(y) {}
	bool operator==(const edge &e) const {return (e.x==x && e.y==y) || (e.x==y && e.y==x);}
    inline bool operator < (const edge &o) const {
        if (x != o.x) return x < o.x;
		else return y < o.y;
    }
};

typedef std::vector<vertex> pointset;
typedef std::vector<edge> edgelist;

#include <fstream>
#include <sstream>
#include <iomanip>
#include "timing.h"
class BaseTask {
public:
    BaseTask(pointset p, double t, double var, int landmarks, double alpha) : p(p), t(t), var(var), landmarks(landmarks),logGraph(false), alpha(alpha) {}
    virtual ~BaseTask(){logFile.close();}
    const pointset p;
    int seed;
    double t;
	double alpha;
    const double var;
    const int landmarks;
	std::string logName;
    std::ofstream logFile;
    std::ofstream graphLogFile;
    bool logGraph;
    void LogGraph(std::string txt){
    	if(!graphLogFile.is_open()){
    		std::stringstream ss;
    		ss << logName<< ".graph";
				graphLogFile.open(ss.str().c_str(),std::ios::app);
				graphLogFile<<"===== STARTING NEW GRAPH ("<<seed<<") ========\n";
			}
			graphLogFile<<txt;
		}
		void Log(std::string txt){
			if(!logFile.is_open()){
				logFile.open(logName.c_str(),std::ios::app);
			}
			logFile<<std::setprecision(20)<<Timing::elapsed()<<">>"<<txt<<"\n";
		}
		void Log(std::string name, long long value){
			if(!logFile.is_open()){
				logFile.open(logName.c_str(),std::ios::app);
			}
			logFile<<std::setprecision(20)<<Timing::elapsed()<<">>"<<name<<value<<"\n";
		}
		void LogDouble(std::string name, double value){
			if(!logFile.is_open()){
				logFile.open(logName.c_str(),std::ios::app);
			}
			logFile<<std::setprecision(20)<<Timing::elapsed()<<">>"<<name<<value<<"\n";
		}
		void Log(int edgeCount,int longestEdgeLength,int avgEdgeLength){
			std::stringstream ss;
			ss<<edgeCount<<":"<<longestEdgeLength<<":"<<avgEdgeLength;
			Log(ss.str().c_str());
		}
		virtual bool create_edge(edge e){
			if(logGraph){
				std::stringstream ss;
				ss<<e.x<<","<<e.y<<"\n";
				LogGraph(ss.str());
			}
            return false;
    }
};

#endif

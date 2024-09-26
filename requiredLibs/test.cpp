#include <cstdio>
#include <QHash>
#include <iostream>
#include <fstream>
#include <sstream>

#include "structures.h"
#include "stats.h"
#include "solver.h"
#include "generator.h"
#include "timing.h"
#include "verify.h"

unsigned int qHash(const edge &e) {
    return e.x*(e.x+31)+e.y*(e.y+31);
}

struct TestTask : public BaseTask {
    double origT;
    TestTask(pointset p, double t, double var, int landmarks) : BaseTask(p,t,var,landmarks), origT(t) {}
    edgelist E;
    map<edge,int> edges;
    virtual bool create_edge(edge e) {
        edges[e]++;
        if (edges[e]==1) {
            E.push_back(e);
        }
        return false;
    }
    void reset() {
        E.clear();
        t=origT;
    }
};

typedef pointset (*generator_fetcher)(int c, int seed, int &count, int max_coordinate, double shape, double &t);
pointset getNormalGenerator(int c, int seed, int& count, int max_coordinate, double shape, double& t) {
	return Generator::generators[c].function(seed,count,max_coordinate,shape,t,0);
}
pointset getChallengeGenerator(int c, int seed, int& count, int max_coordinate, double shape, double& t) {
	return Generator::generators[10].function(seed,count,max_coordinate,shape,t,c);
}

using namespace std;
const int repeats=1;
void DoTests(int N, int* algos, int algoCount, int* cases, int caseCount, double T, generator_fetcher fetcher, const char *_Filename, double alpha = 0.5) {
	std::stringstream outFileName;
		outFileName << _Filename << N << "-" << T << ".txt";
	ofstream myfile;
	myfile.open(outFileName.str());
    myfile << "# algorithm                                data generator     verts fix-t     time        +/- 95%%   t     +/-      edges     +/-    mdegr  +/-            weight      +/-  valid?      repeats" << endl;
	for (int ai = 0; ai < algoCount; ai++) {
		int a = algos[ai];
		for (int ci = 0; ci < caseCount; ci++) {
			int c = cases[ci];
			std::cerr << a << " " << c << endl;
            int n;
            double t=T;
            double shape=0;
            if (c==6 || c==7) shape=3;
            if (c==9) shape=.75;
            
            stats time, actt, edges, mdegr, weight;
            bool valid = true;
            int i;
            for (i=0; i<repeats; i++) {
                pointset p = fetcher(c, 100+i, N, 1000, shape, t);
                n = p.size();
                TestTask T(p,t,alpha,100);
                Timing::reset();
                Solver::solvers[a].function(&T);
                double elapsed = Timing::elapsed();
                time.add(elapsed);
                verification_info vi = verify(p, T.E);
                if (vi.t>t) valid = false;
                actt.add(vi.t);
                edges.add(T.E.count());
                mdegr.add(vi.maxdeg);
                weight.add(vi.weight);
                if (c==7 || c==8 || c==11) {i++; break;}
            }
            
			/*
            printf("%2d,\"%s\", %*s\"B%d.txt\", %*s%4d, %6.3lf, %12.7lf,          , %4.2lf,     , %8.1lf,       , %6.2lf,      , %15.0lf,        , %s, %2d\n",
                a,Solver::solvers[a].name, 36-(int)strlen(Solver::solvers[a].name), "",
                c, 14-6, "",
                n,
                t,
                time.mean(),
                actt.mean(),
                edges.mean(),
                mdegr.mean(),
                weight.mean(),
                n<=dilationVerificationBound?valid?"ok       ":"invalid! ":"unchecked",
                i
            );
			//*/
			myfile << a + 1 << ",";
			myfile << "\"" << Solver::solvers[a].name << "\"" << ",";
			myfile << "\"" << Generator::generators[c].name << "\"" << ",";
			myfile << n << ",";
			myfile << t << ",";
			myfile << time.mean() << ",";
			myfile << time.radius95p() << ",";
			myfile << actt.mean() << ",";
			myfile << actt.radius95p() << ",";
			myfile << edges.mean() << ",";
			myfile << edges.radius95p() << ",";
			myfile << mdegr.mean() << ",";
			myfile << mdegr.radius95p() << ",";
			myfile << weight.mean() << ",";
			myfile << weight.radius95p() << ",";
			myfile << (N <= dilationVerificationBound ? valid ? "ok" : "invalid!" : "unchecked") << ",";
			myfile << i << endl;
            /*
            if (i<=1)
            printf("%2d,\"%s\", %*s\"%s\", %*s%4d, %6.3lf, %12.7lf,          , %4.2lf,     , %8.1lf,       , %6.2lf,      , %15.0lf,        , %s, %2d\n",
                a,Solver::solvers[a].name, 36-(int)strlen(Solver::solvers[a].name), "",
                Generator::generators[c].name, 14-(int)strlen(Generator::generators[c].name), "",
                n,
                t,
                time.mean(),
                actt.mean(),
                edges.mean(),
                mdegr.mean(),
                weight.mean(),
                N<=dilationVerificationBound?valid?"ok       ":"invalid! ":"unchecked",
                i
            );
            else
            printf("%2d,\"%s\", %*s\"%s\", %*s%4d, %6.3lf, %12.7lf, %9.6lf, %4.2lf, %4.2lf, %8.1lf, %6.1lf, %6.2lf, %5.2lf, %15.0lf, %7.0lf, %s, %2d\n",
                a,Solver::solvers[a].name, 36-(int)strlen(Solver::solvers[a].name), "",
                Generator::generators[c].name, 14-(int)strlen(Generator::generators[c].name), "",
                n,
                t,
                time.mean(), time.radius95p(),
                actt.mean(), actt.radius95p(),
                edges.mean(), edges.radius95p(),
                mdegr.mean(), mdegr.radius95p(),
                weight.mean(), weight.radius95p(),
                N<=dilationVerificationBound?valid?"ok       ":"invalid! ":"unchecked",
                i
            );
            //*/
            fflush(stdout);
        }
    }
}

int main(int argc, char **argv) {
    //QVector<int> algos={0,1,2,3,4,5,6,7,9,10,11,12,14,15,16,17,18,19,20}; // n<=200
    //QVector<int> algos={2,3,4,6,9,10,11,12,14,15,16,17,18,20}; // n<=500
    //QVector<int> algos={2,4,9,10,11,12,14,16,17,20}; //n<=1000
    //QVector<int> algos={4,10,11,12,16,20}; //n<=2000
    //QVector<int> algos={11,12,16,20};
    //QVector<int> cases={1,2,4,5}; // for challenges.
	
	/*
	double T = 1.1;
    int algos1[] = {2,4,9,10,11,12,20,21};
    int cases1[] = {0,1,2,3,4,6,7,9,11};

	DoTests(100, algos1, sizeof(algos1) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, NULL);
	DoTests(200, algos1, sizeof(algos1) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, NULL);
	DoTests(500, algos1, sizeof(algos1) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, NULL);
	DoTests(1000, algos1, sizeof(algos1) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, NULL);
	
    int algos2[] = {11,12,20,21};

	DoTests(2000, algos2, sizeof(algos2) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, NULL);
	DoTests(5000, algos2, sizeof(algos2) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, NULL);
	DoTests(10000, algos2, sizeof(algos2) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, NULL);
	//*/

	///*
	double T = 2.0;
    int algos1[] = {0,1,2,3};
    int cases1[] = {0,1,2,3,4,6,7,9,11};

	DoTests(100, algos1, sizeof(algos1) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, "Greedy");
	DoTests(200, algos1, sizeof(algos1) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, "Greedy");
	
    int algos2[] = {2,3};

	//DoTests(500, algos2, sizeof(algos2) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, "Greedy");
	//DoTests(1000, algos2, sizeof(algos2) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, "Greedy");
	//DoTests(2000, algos2, sizeof(algos2) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, "Greedy");
	//*/
	
	/*
	double T = 1.1;
    int algos1[] = {14,15,16,17,18,19};
    int cases1[] = {0,1,2,3,4,6,7,9,11};

	DoTests(100, algos1, sizeof(algos1) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, "Hyb0.9", 0.9);
	DoTests(200, algos1, sizeof(algos1) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, "Hyb0.9", 0.9);
	DoTests(500, algos1, sizeof(algos1) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, "Hyb0.9", 0.9);
	DoTests(1000, algos1, sizeof(algos1) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, "Hyb0.9", 0.9);
	//*/
	
	/*
	double T = 2.0;
    int algos1[] = {2,4,9,10,11,12,20,21};
    int cases1[] = {1,2,3,4,5,6};

	DoTests(100, algos1, sizeof(algos1) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getChallengeGenerator, "challenge");
	//*/

	return 0;
}
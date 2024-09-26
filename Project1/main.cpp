#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <cstring>
#include <cmath>

#include "structures.h"
#include "stats.h"
#include "solver.h"
#include "generator.h"
#include "timing.h"
#include "verify.h"

#define	FirstCandidate	1
#define	BestCandidate	2

//candidatePolicy:
//1:FirstCandidate
//2:BestCandidate
inline double distance(const vertex &a, const double x, const double y)
{
	return hypot(a.x - x, a.y - y);
}

inline double distance(const vertex &a, const vertex &b)
{
	return distance(a, b.x, b.y);
}
struct TestTask : public BaseTask {
	double origT;
	double weight;
	double maxLenght;
	int *degs;
	

	TestTask(pointset p, double t, double var, int landmarks, double alpha) 
		: BaseTask(p, t, var, landmarks, alpha), origT(t)
	{
		weight = 0;
		maxLenght = 0;
		degs = new int[p.size()];
		for (size_t i = 0; i < p.size(); i++)
		{
			degs[i] = 0;
		}
	}
	edgelist E;
	std::map<edge, int> edges;

	virtual bool create_edge(edge e) {
		edges[e]++;
		if (edges[e] == 1) {
			E.push_back(e);
		}
		double d = distance(p[e.x], p[e.y]);
		if (maxLenght < d)
			maxLenght = d;
		
		weight += d;
		degs[e.x]++;
		degs[e.y]++;

		return false;
	}
	int MaxDeg()
	{
		int max = degs[0];
		for (size_t i = 1; i < p.size(); i++)
		{
			if (degs[i] > max)
				max = degs[i];
		}
		return max;
	}

	void reset() {
		E.clear();
		t = origT;
	}
};

typedef pointset(*generator_fetcher)(int c, int seed, int &count, int max_coordinate, double shape, double &t);
pointset getNormalGenerator(int c, int seed, int& count, int max_coordinate, double shape, double& t) {
	return Generator::generators[c].function(seed, count, max_coordinate, shape, t, 0);
}
pointset getChallengeGenerator(int c, int seed, int& count, int max_coordinate, double shape, double& t) {
	return Generator::generators[10].function(seed, count, max_coordinate, shape, t, c);
}

using namespace std;
const int repeats = 1;
void DoTests(pointset p, int N, int* algos, int algoCount, int distribution, double T, const char *_Filename, double alpha = 0.5) {
	std::stringstream outFileName;
	outFileName << _Filename << "(" << Generator::generators[distribution].name << "-" << N << "-" << T << ")" << ".txt";
	ofstream myfile;
	myfile.open(outFileName.str());
	
	myfile << "# algorithm            time        queryCount	edges     +/-    mdegr  +/-            weight      " << endl;
	for (int ai = 0; ai < algoCount; ai++) {
		int a = algos[ai];

		//std::cerr << a << " " << endl;
		double t = T;
		/*double shape = 0;
		if (c == 6 || c == 7) shape = 3;
		if (c == 9) shape = .75;*/

		stats time, actt, edges, mdegr, weight;
		bool valid = true;
		int i;
		for (i = 0; i < repeats; i++) {
			TestTask T(p, t, alpha, 100, 0);
			Timing::reset();
			qcount = 0;
			Solver::solvers[a].function(&T);
			double elapsed = Timing::elapsed();
			time.add(elapsed);
			verification_info vi = verify(p, T.E);
			if (vi.t > t) valid = false;
			actt.add(vi.t);
			edges.add(T.E.size());
			mdegr.add(vi.maxdeg);
			weight.add(vi.weight);
			char otherString[15]; 
			
			std::strncpy(otherString, Solver::solvers[a].name, 14);
			otherString[14] = '\0'; // place the null terminator

			cout << otherString << "\t\t\t" << T.E.size() << "\t\t\t\t" << time.mean() << "\n";
			//if (c == 7 || c == 8 || c == 11) { i++; break; }
		}

		myfile << a << ",";
		myfile << "\"" << Solver::solvers[a].name << "\"" << "\t\t";
		//myfile << n << ",";
		//myfile << t << ",";
		myfile << time.mean() << "\t\t";
		myfile << qcount << "\t";
		//myfile << time.radius95p() << ",";
		//myfile << actt.mean() << ",";
		//myfile << actt.radius95p() << ",";
		//myfile << edges.mean() << "\t";
		//myfile << edges.radius95p() << ",";
		//myfile << mdegr.mean() << "\t";
		//myfile << mdegr.radius95p() << ",";
		//myfile << weight.mean() << "\t";
		//myfile << weight.radius95p() << ",";
		//myfile << (N <= dilationVerificationBound ? valid ? "ok" : "invalid!" : "unchecked") << ",";
		//myfile << i;
		myfile << endl;

		std::fflush(stdout);
	}
}

string* split(char* str, int* count)
{
	int c = 1;
	for (size_t f = 0; f < std::strlen(str); f++)
	{
		if (str[f] == ',')
			c++;
	}
	
	string* arr = new string[c];
	int g = 0;
	for (size_t f = 0; f < c; f++)
	{
		while (str[g] != ',' && g < strlen(str))
		{
			string h = string(1, str[g]);
			arr[f].append(h);
			g++;
		}
		g++;
	}
	*count = c;
	return arr;
}

int*	toIntArray(string* arr, int count)
{
	int* result = new int[count];
	for (size_t k = 0; k < count; k++)
	{
		result[k] = std::atoi(arr[k].c_str());
	}
	return result;
}

double*	toDoubleArray(string* arr, int count)
{
	double* result = new double[count];
	for (size_t k = 0; k < count; k++)
	{
		result[k] = std::atof(arr[k].c_str());
	}
	return result;
}

void	greedySpannerTest(int argc, char **argv)
{
	//int algos1[] = { 4,/*12,*/13,19, 9, 10,34,35, 48,49,50, 51 };
	//int cases1[] = { 0,1,2,3,4,7,9/*,11*/ };
	/*argc = 10;
	argv = new char*[10];
	for (size_t i = 0; i < 10; i++)
	{
	argv[i] = new char[80];
	}
	strcpy_s(argv[1],80,  "n");
	strcpy_s(argv[2],80, "28");
	strcpy_s(argv[3],80, "d");
	strcpy_s(argv[4],80, "0");
	//strcpy_s(argv[5],80, "-f");
	strcpy_s(argv[5],80, "t");
	strcpy_s(argv[6],80, "2");
	strcpy_s(argv[7],80, "a");
	strcpy_s(argv[8],80, "12,66");*/
	int* algos = NULL, *N = NULL, *cases = NULL;
	double * T = NULL;
	int count, countT, countN, countA, countC;
	string* arr;

	if (argc != 9 && argc != 10 && argc != 6)
	{
		std::cout << "Error: Parameters Not Specified Correctly!";
		return;
	}

	bool useFile = false;

	for (size_t i = 0; i < argc; i++)
	{
		if (strcmp(argv[i], "N") == 0 || strcmp(argv[i], "n") == 0)
		{
			arr = split(argv[++i], &count);
			N = toIntArray(arr, count);
			countN = count;
		}
		else if (strcmp(argv[i], "T") == 0 || strcmp(argv[i], "t") == 0)
		{
			arr = split(argv[++i], &count);
			T = toDoubleArray(arr, count);
			countT = count;
		}
		else if (strcmp(argv[i], "A") == 0 || strcmp(argv[i], "a") == 0)
		{
			arr = split(argv[++i], &count);
			algos = toIntArray(arr, count);
			countA = count;
		}
		else if (strcmp(argv[i], "D") == 0 || strcmp(argv[i], "d") == 0)
		{
			arr = split(argv[++i], &count);
			cases = toIntArray(arr, count);
			countC = count;
		}
		else if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "-F") == 0)
			useFile = true;
	}

	//DoTests(500, algos1, sizeof(algos1) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, "Greedy");
	double dummyT;
	if (argc == 6)
	{
		using namespace std;
		char filename[80];
		generator_fetcher fetcher = getNormalGenerator;
		double x, y;

		for (size_t k = 0; k < countC; k++)
		{
			for (size_t i = 0; i < countN; i++)
			{
				pointset p = fetcher(cases[k], 100, N[i], 10000, 0, dummyT);// shape parameter is 0 for square and uniform distributions, dummyT is useless
				std::sprintf(filename, "%s-%d.txt", Generator::generators[k].name, N[i]);
				ofstream myfile(filename);

				for (size_t j = 0; j < p.size(); j++)
				{
					myfile << p[j].x << "\t";
					myfile << p[j].y << endl;
				}
				myfile.close();
			}
		}
	}
	else
	{
		pointset p;
		for (size_t k = 0; k < countC; k++)
		{
			for (size_t i = 0; i < countN; i++)
			{
				if (useFile)
				{
					char filename[80];
					std::sprintf(filename, "%s-%d.txt", Generator::generators[k].name, N[i]);
					ifstream myfile;
					double x, y;
					myfile.open(filename, ios::in);
					if (myfile.fail())
					{
						std::cout << "cannot open file" << "\'" << filename << "\'";
						return;
					}
					while (true)
					{
						myfile >> x;
						myfile >> y;
						if (myfile.fail())
							break;
						p.push_back(vertex(
							x,
							y
						));

						//cout << x << "," << y << endl;
					}
					myfile.close();
				}
				else
				{
					generator_fetcher fetcher = getNormalGenerator;
					p = fetcher(cases[k], 100, N[i], 1000, 0, dummyT);// shape parameter is 0 for square and uniform distributions, dummyT is useless
				}

				for (size_t j = 0; j < countT; j++)
				{
					DoTests(p, N[i], algos, countA, cases[k], T[j], "result");
				}
			}
		}
	}
	//DoTests(16000, algos1, sizeof(algos1) / sizeof(int), cases1, sizeof(cases1) / sizeof(int), T, getNormalGenerator, "Greedy");
}

void GapTests(pointset p, int N, int* algos, int algoCount, int distribution, double Theta, double w, const char *_Filename) {
	stats time, edges, mdegr, weight;
	std::stringstream k;
	k << "";
	std::cout << "N:" << N << endl;
	for (size_t i = 0; i < algoCount; i++)
	{
		TestTask T(p, Theta, w, 100, Theta);
		Timing::reset();
		qcount = 0;
		Solver::solvers[algos[i]].function(&T);
		double elapsed = Timing::elapsed();
		time.add(elapsed);
		edges.add(T.E.size());
		std::cout << "time:" << elapsed << ", edge:" << T.E.size() << ", query:" << qcount <<endl;
	}
	
	//cout << k.str() << endl;
	//getchar();
}

int gapSpannerTest(int argc, char **argv)
{
	int* algos = NULL, *N = NULL, *cases = NULL;
	double * tw = NULL;
	int count, countTW, countN, countA, countC;
	string* arr;

	/*argc = 9;
	argv = new char*[10];
	for (size_t i = 0; i < 10; i++)
	{
		argv[i] = new char[80];
	}
	strcpy_s(argv[1], 80, "n");
	strcpy_s(argv[2], 80, "1000,2000,4000,8000,16000");
	strcpy_s(argv[3], 80, "d");
	strcpy_s(argv[4], 80, "0");
	//strcpy_s(argv[5],80, "-f");
	strcpy_s(argv[5], 80, "tw");
	strcpy_s(argv[6], 80, "0.174533,0.1555");//, 0.139626,0.0088
	strcpy_s(argv[7], 80, "a");
	strcpy_s(argv[8], 80, "68,69");*/

	if (argc != 9)
	{
		std::cout << "Error: Parameters Not Specified Correctly!" << endl;
		return 0;
	}

	bool useFile = false;

	for (size_t i = 0; i < argc; i++)
	{
		if (strcmp(argv[i], "N") == 0 || strcmp(argv[i], "n") == 0)
		{
			arr = split(argv[++i], &count);
			N = toIntArray(arr, count);
			countN = count;
		}
		else if (strcmp(argv[i], "TW") == 0 || strcmp(argv[i], "tw") == 0)
		{
			arr = split(argv[++i], &count);
			if (count % 2 == 1)
			{
				std::cout << "Error: values of tw parameter must be even!" << endl;
				return 0;
			}
			tw = toDoubleArray(arr, count);
			countTW = count;
		}
		else if (strcmp(argv[i], "A") == 0 || strcmp(argv[i], "a") == 0)
		{
			arr = split(argv[++i], &count);
			algos = toIntArray(arr, count);
			countA = count;
		}
		else if (strcmp(argv[i], "D") == 0 || strcmp(argv[i], "d") == 0)
		{
			arr = split(argv[++i], &count);
			cases = toIntArray(arr, count);
			countC = count;
		}

	}
	std::time_t t = std::time(0);   // get time now
	//std::tm* now = std::localtime(&t);
	
	stats time, edges, mdegr, weight;
	generator_fetcher fetcher = getNormalGenerator;
	double dummyT, seed =/*100*/(double)t, theta, w, alpha;
	double	alphs[] = {1.1, 1.5, 2, 3};

	for (size_t i = 0; i < countC; i++)
	{
		for (size_t j = 0; j < countN; j++)
		{
			pointset ps = fetcher(cases[i], seed, N[j], 100000, 0, dummyT);

			for (size_t k = 0; k < countTW; k += 2)
			{
				theta = tw[k];
				w = tw[k + 1];
				
				//for (size_t h = 0; h < 4; h++)
				//{
					//alpha = pow(theta, alphs[h]);

					std::ostringstream strDouble;
					std::string twStr;
					strDouble << std::fixed;
					strDouble << std::setprecision(5);
					strDouble << theta;
					strDouble << "-" << w;
					//strDouble << "-" << alpha;
					twStr = strDouble.str();

					std::stringstream outFileName;

					//for (int a = 0; a < ps.size(); a++)
						///std::cout << "x:" << ps[a].x << ",y:" << ps[a].y <<endl;

					outFileName << Generator::generators[cases[i]].name << "-" << N[j] << "(" << twStr << ")" << ".txt";
					ofstream myfile;
					myfile.open(outFileName.str(), std::ios_base::app);

					//cout << "alpha:" << alpha << endl;
					for (size_t a = 0; a < countA; a++)
					{
						//if (h > 0 && (!(algos[a] >= 60 && algos[a] < 64)))
							//continue;

						TestTask T(ps, theta, w, 100, theta);
						Timing::reset();
						qcount = 0;
						Solver::solvers[algos[a]].function(&T);
						double elapsed = Timing::elapsed();
						time.add(elapsed);
						edges.add(T.E.size());
						std::cout << Solver::solvers[algos[a]].name
							<< "\t\ttime:" << elapsed
							//<< "\tQ:" << qcount 
							<< "\tE:" << T.E.size() << "\t\t"
							<< T.maxLenght << "\t\t"
							//<< T.maxLenght << "\t\t"
							//<< T.MaxDeg()
							<< endl;

						myfile << Solver::solvers[algos[a]].name << "\t\t";
						myfile << elapsed << "\t\t";
						//myfile << qcount<<"\t\t";
						myfile << T.maxLenght << "\t\t";
						//myfile << T.maxLenght << "\t\t";
						//myfile << T.MaxDeg();
						myfile << endl;

						std::fflush(stdout);

					}
					myfile.close();
				//}
			}
		}
	}

	/*
	int	al[] = { 57, 64};


	int n[] = { 125, 250, 500, 1000, 2000, 4000, 8000, 16000};//6
	double dummyT, seed = 5;//5
	generator_fetcher fetcher = getNormalGenerator;
	int distrib = 1;
	
	pointset p01;// = fetcher(distrib, seed, n[0], 1000000, 0, dummyT);// shape parameter is 0 for square and uniform distributions, dummyT is useless
	int algCount = 2;
	double theta = 0.174533, alpha = 0.139626, w = 0.1555;
	
	for (size_t i = 0; i < 5; i++)
	{
		p01 = fetcher(distrib, seed, n[i], 1000000, 0, dummyT);// shape parameter is 0 for square and uniform distributions, dummyT is useless
		GapTests(p01, n[i], al, algCount, distrib, theta, w, "aa");
	}
	pointset p01 = fetcher(distrib, seed, n[0], 1000000, 0, dummyT);// shape parameter is 0 for square and uniform distributions, dummyT is useless
	pointset p02 = fetcher(distrib, seed, n[1], 100000, 0, dummyT);// shape parameter is 0 for square and uniform distributions, dummyT is useless
	pointset p0 = fetcher(distrib, seed, n[2], 100000, 0, dummyT);// shape parameter is 0 for square and uniform distributions, dummyT is useless
	pointset p1 = fetcher(distrib, seed, n[3], 100000, 0, dummyT);// shape parameter is 0 for square and uniform distributions, dummyT is useless
	pointset p2 = fetcher(distrib, seed, n[4], 100000, 0, dummyT);// shape parameter is 0 for square and uniform distributions, dummyT is useless
	pointset p3 = fetcher(distrib, seed, n[5], 100000, 0, dummyT);// shape parameter is 0 for square and uniform distributions, dummyT is useless
	pointset p4 = fetcher(distrib, seed, n[6], 100000, 0, dummyT);// shape parameter is 0 for square and uniform distributions, dummyT is useless
	pointset p5 = fetcher(distrib, seed, n[7], 100000, 0, dummyT);// shape parameter is 0 for square and uniform distributions, dummyT is useless
	//int candidatePolicy = BestCandidate;
	
	//p01.clear();
	p01.push_back(vertex(5,5));
	p01.push_back(vertex(4.998, 3));
	p01.push_back(vertex(5, 7));
	p01.push_back(vertex(3, 10));
	p01.push_back(vertex(7, 11));

	pointset p01;
	p01.push_back(vertex(0,4));
	p01.push_back(vertex(16, 4));
	p01.push_back(vertex(20, 10));
	p01.push_back(vertex(30, 2));
	
	GapTests(p01, n[0], al, algCount, distrib, theta, w, "aa");
	GapTests(p02, n[1], al, algCount, distrib, theta, alpha, "aa");
	GapTests(p0, n[2], al, algCount, distrib, theta, alpha, "aa");
	GapTests(p1, n[3], al, algCount, distrib, theta, w, "aa");
	GapTests(p2, n[4], al, algCount, distrib, theta, w, "aa");
	GapTests(p3, n[5], al, algCount, distrib, theta, w,  "aa");
	GapTests(p4, n[6], al, algCount, distrib, theta, w,  "aa");
	GapTests(p5, n[7], al, algCount, distrib, theta, alpha, "aa");
	*/
	return 0;
	
}

void generatePoint(int argc, char **argv)
{
	int *N = NULL, *cases = NULL;

	int count, countN, countC;
	string* arr;

	bool useFile = false;

	for (size_t i = 0; i < argc; i++)
	{
		if (strcmp(argv[i], "N") == 0 || strcmp(argv[i], "n") == 0)
		{
			arr = split(argv[++i], &count);
			N = toIntArray(arr, count);
			countN = count;
		}
		else if (strcmp(argv[i], "D") == 0 || strcmp(argv[i], "d") == 0)
		{
			arr = split(argv[++i], &count);
			cases = toIntArray(arr, count);
			countC = count;
		}
	}

	

	for (size_t i = 0; i < countC; i++)
	{
		for (size_t j = 0; j < countN; j++)
		{
			std::time_t t = std::time(0);   // get time now
											//std::tm* now = std::localtime(&t);

			generator_fetcher fetcher = getNormalGenerator;
			double dummyT, seed = (double)t;
			pointset ps = fetcher(cases[i], seed, N[j], 100000, 0, dummyT);

			std::ostringstream strDouble;
			std::string twStr;
			strDouble << std::fixed;
			strDouble << std::setprecision(5);
			twStr = strDouble.str();

			std::stringstream outFileName;
			outFileName << Generator::generators[cases[i]].name << "-" << N[j] << ".txt";
			ofstream myfile;
			myfile.open(outFileName.str());

			for (int a = 0; a < ps.size(); a++)
				myfile << ps[a].x << "," << ps[a].y << endl;

			myfile.close();
		}
	}
}

int main(int argc, char **argv) 
{
	gapSpannerTest(argc, argv);
	//greedySpannerTest(argc, argv);
	//generatePoint(argc, argv);
	/*int* al = new int[3];
	al[0] = 54;
	al[1] = 62;
	al[2] = 57;

	pointset p01;
	p01.push_back(vertex(0, 4));
	p01.push_back(vertex(16, 4));
	p01.push_back(vertex(16, 10));
	p01.push_back(vertex(30, 2));

	GapTests(p01, 4, al, 3, 0, 0.261799, 0.2, "aa");*/

	std::cout << "Finish...";
	//getchar();
	
}
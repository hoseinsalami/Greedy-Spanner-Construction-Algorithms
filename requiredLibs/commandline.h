#ifndef COMMANDLINE_H
#define COMMANDLINE_H

//default number of test iterations
#define TEST_ITERATIONS 100
//if #iterations has been specified it can be atmost this:
#define MAX_TEST_ITERATIONS 100000
//default number of points for a test
#define TEST_NPOINTS 2000
#define TEST_GENID 0
#define TEST_MAX_COORD 10000
#define TEST_SHAPE 0.0
#define TEST_T 1.1
#define TEST_SOURCE 0
#define TEST_VAR 0.0
#define TEST_LANDMARKS 0

void PrintAlgos();
void PrintGenerators();
int ProcessCommand(int nParams, char ** command);
int Test(int genId, int algoId, int nPoints, int iterations, double t);
void Help();
int EnableLogging(int algoId, int genId, int nPoints, double t);
int CommandMain(int argc, char ** argv);

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

//inline uint qHash(const edge &e);
using namespace std;


#endif // COMMANDLINE_H


#include "commandline.h"

//uint qHash(const edge &e) {
//    return e.x*(e.x+31)+e.y*(e.y+31);
//}

struct CLTask : public BaseTask {
    double origT;
    CLTask(pointset p, double t, double var, int landmarks, int _seed) : BaseTask(p,t,var,landmarks), origT(t) {seed=_seed;}
    edgelist E;
    QHash<edge,int> edges;
    virtual bool create_edge(edge e) {
				BaseTask::create_edge(e); //for logging
            E.push_back(e);
        return false;
    }
    void reset() {
        E.clear();
        t=origT;
    }
};

int CommandMain(int argc, char ** argv){
	if(argc>=3){//argv: path, -c, command
		return ProcessCommand(argc-3,&argv[2]);
	}else{
		Help();
	}
	return 0;
}

int VerifyTestParams(int genId, int algoId, int iterations){
	int nGenerators=0,nAlgos=0;
	while(Generator::generators[nGenerators].function!=NULL){
		nGenerators++;
	}
	while(Solver::solvers[nAlgos].function!=NULL){
		nAlgos++;
	}

	if(!(genId>=0 && genId<nGenerators)){
		cerr<<"invalid generator id\n";
		return 1;
	}
	if(!(algoId>=0 && algoId<nAlgos)){
		cerr<<"invalid algorithm id\n";
		return 2;
	}
	if(!(iterations>0 && iterations<MAX_TEST_ITERATIONS)){
		cerr<<"invalid number of test iterations\n";
		return 3;
	}
	return 0; //valid params
}

int Test(int genId, int algoId, int nPoints, int iterations, double t){
	int ret=VerifyTestParams(genId,algoId,iterations);
	if(ret!=0){return ret;}
	stringstream ss;
	ss<<Solver::solvers[algoId].name<<"_"<<Generator::generators[genId].name<<"_";
	ss<<nPoints<<"_"<<t<<"_"<<iterations<<".log";
	string filename=ss.str();
	for(int i=0;i<filename.length();i++){
		if(filename[i]==' '){filename[i]='_';}
	}
	//truncate log file
	std::ofstream logFile;
	std::ofstream IPEFile;
	logFile.open(filename.c_str(),std::ios::trunc);
	logFile.close();

	stringstream ss2;
	ss2<<filename<<".graph";
	logFile.open(ss2.str().c_str(),std::ios::trunc);
	logFile.close();

	ss2<<".ipe";
	IPEFile.open(ss2.str().c_str(),std::ios::trunc);
	IPEFile.close();

	//Time to do the actual test
	cout<<"Starting tests:\n";
	cout<<"\tGenerator: "<<Generator::generators[genId].name<<" ("<<nPoints<<" points, t="<<t<<")\n";
	cout<<"\tAlgorithm: "<<Solver::solvers[algoId].name<<"\n";
	cout<<"\tIterations: "<<iterations<<"\n";

	for(int i=0;i<iterations;i++){
		//(int seed, int& count, int max_coordinate, double shape, double& t, int source)
		cout <<"Generating pointset "<<i<<"/"<<iterations<<"\n";
		int seed=rand();
		pointset p=Generator::generators[genId].function(seed, nPoints, TEST_MAX_COORD, TEST_SHAPE, t, TEST_SOURCE);
		CLTask* task=new CLTask(p,t,TEST_VAR,TEST_LANDMARKS,seed);
		task->logName=filename;
		cout <<"Running algorithm iteration "<<i<<"/"<<iterations<<"\n";
		Solver::solvers[algoId].function(task);
		delete task;
	}
	return 0;
};


int ProcessCommand(int nParams, char ** command){
	//listgen
	if(strcmp(command[0],"listgen")==0){PrintGenerators();return 0;}
	//listalgo
	if(strcmp(command[0],"listalgo")==0){PrintAlgos();return 0;}
	//test
	double t=TEST_T;
	int genId=TEST_GENID, algoId=0,nPoints=TEST_NPOINTS,iterations=TEST_ITERATIONS;
	if(nParams>=1 && strcmp(command[0],"test")==0){ //argv: path, -c, test, algoId, [t, genId, nPoints, iterations]
		stringstream ss;

		ss<<command[1];
		ss>>algoId;
		ss.clear();
		if(nParams>=2){
			ss<<command[2];
			ss>>t;
			ss.clear();
		}
		if(nParams>=3){
		ss<<command[3];
		ss>>genId;
		ss.clear();
		}
		if(nParams>=4){
			ss<<command[4];
			ss>>nPoints;
			ss.clear();
		}
		if(nParams==5){
			ss<<command[5];
			ss>>iterations;
			ss.clear();
		}
		return Test(genId,algoId,nPoints,iterations,t);
	}
	cerr<<"Invalid command or number of parameters.\n\n";
	Help();
	return 1;
}

void Help(){
    cout<<"Syntax: [command] [params]"<<endl;
    cout<<"Commands:"<<endl;
    cout<<"\ttest: performs tests, params: algoId, [t="<<TEST_T<<"], [generatorId="<<TEST_GENID<<"], [nPoints="<<TEST_NPOINTS<<"], [iterations="<<TEST_ITERATIONS<<"]"<<endl;;
    cout<<"\t\texample: test 2 1.2-> tests using:"<<endl;
    cout<<"\t\t\tAlgo: "<<Solver::solvers[2].name<<endl;
    cout<<"\t\t\tgenerator: "<<Generator::generators[TEST_GENID].name<<" ("<<TEST_NPOINTS<<" point with t="<<TEST_T<<")"<<endl<<endl;
    cout<<"\tlistgen: lists generators with their id."<<endl<<endl;
    cout<<"\tlistalgo: lists algorithms with their id."<<endl<<endl;
}

void PrintGenerators(){
	int i=0;
	Generator::type gType=Generator::generators[i];
	while(gType.function!=NULL){
		cout <<i<<": "<<gType.name<<endl;
		gType=Generator::generators[++i];
	}
}

void PrintAlgos(){
	int i=0;
	Solver::type sType=Solver::solvers[i];
	while(sType.function!=NULL){
		cout <<i<<": "<<sType.name<<endl;
		sType=Solver::solvers[++i];
	}
}

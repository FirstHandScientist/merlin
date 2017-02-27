/*
 * uai.cpp
 *
 *  Created on: Jun 6, 2014
 *      Author: radu
 */


#include "base.h"

#include "braobb.h"
#include "any_aaobf.h"

string GetEnv( const string & var ) {
     const char * val = ::getenv( var.c_str() );
     if ( val == 0 ) {
         return "";
     }
     else {
         return val;
     }
}

// UAI COMPETITION
int uai(int argc, char** argv) {

	// make sure we have the right number of arguments
	if (argc != 5) {
		std::cerr << "Incorrect number of arguments." << std::endl
				<< "  expected: <input> <evidence> <query> MMAP" << std::endl;
		return EXIT_FAILURE;
	}

	// get the arguments
	string strInputFileName( argv[1] );
	string strEvidenceFileName( argv[2] );
	string strQueryFileName( argv[3] );
	string strTaskName( argv[4] );

	size_t found = strInputFileName.find_last_of("/");
	string strProbName = (found != std::string::npos) ? strInputFileName.substr(found + 1) : strInputFileName;
	string strOutputFileName = "./" + strProbName + "." + strTaskName;
	std::cout << "Output filename: " << strOutputFileName << std::endl;

	// make sure we're running the MMAP query
	if (strTaskName.compare("MMAP") != 0) {
		std::cerr << "Unknown task name; expected MMAP!" << std::endl;
		return EXIT_FAILURE;
	}

	// init the options
	string algo, heur, prop;
	ProgramOptions* opt = new ProgramOptions;
	opt->executableName = argv[0];
	opt->timeLimit = NONE;
	opt->memoryLimit = 2.0; // 3GB of RAM for BRAOBB, 2GB for AAOBF
	opt->subprobOrder = 0;
	opt->ibound = 15; // 25, 35
	opt->heuristic = HEUR_WMB_MM; heur = "WMB-MM";
//	opt->algorithm = ALGO_BRAOBB; algo = "BRAOBB";
	opt->algorithm = ALGO_ANY_AAOBF; algo = "AOBF";
	opt->incremental = true;
	opt->propagation = CP_NONE; prop = "none"; // unit propagation
	opt->seed = 12345678;
	opt->verbose = true;
	opt->jglp = 100; // 10 iterations for WMB-JG
	opt->jglps = 10;
	opt->mplp = 100; // 100 iterations
	opt->mplps = 10; // 10 seconds timout
	opt->inputFile = strInputFileName;
	opt->outputFile = strOutputFileName;
	opt->evidenceFile = strEvidenceFileName;
	opt->mapFile = strQueryFileName;
	opt->format = FORMAT_UAI;
	opt->positive = true; // for numerical issues with WMB-MM
	opt->orderingIter = 1000;
	opt->orderingTime = -1;

	// set the time and memory limit
	string strTimeLimit = GetEnv("INF_TIME");
	string strMemoryLimit = GetEnv("INF_MEMORY");
	if (strTimeLimit.empty() == false) {
		opt->timeLimit = ::atof(strTimeLimit.c_str());
	}
	if (strMemoryLimit.empty() == false) {
		opt->memoryLimit = ::atof(strMemoryLimit.c_str());
	}

	// set the iterations and timeout for variable ordering
	if (opt->timeLimit == 60) {
		opt->orderingIter = 500;
		opt->orderingTime = 5; // 5 seconds
		opt->ibound = 15; // will be adjusted to fit the memory limit
		opt->mplp = 100; // 100 iterations MPLP for pure MAP assignment
		opt->mplps = 5; // 5 seconds timout
	} else if (opt->timeLimit == 1200) {
		opt->orderingIter = 10000;
		opt->orderingTime = 60; // 1 minute
		opt->ibound = 16;
		opt->ibound = 20; // will be adjusted to fit the memory limit
		opt->mplp = 1000; // 100 iterations MPLP for pure MAP assignment
		opt->mplps = 60; // 60 seconds timout
	} else if (opt->timeLimit >= 3600) {
		opt->orderingIter = 300000;
		opt->orderingTime = 180; // 3 minutes
		opt->ibound = 35; // will be adjusted to fit the memory limit
		opt->mplp = 10000; // 100 iterations MPLP for pure MAP assignment
		opt->mplps = 600; // 10 seconds timout
	}

	// set the seed
	if (opt->seed == NONE) {
		opt->seed = time(0);
	}
	rand::seed(opt->seed);

	std::ostringstream oss2;
	oss2 << "+ i-bound:\t\t" << opt->ibound << std::endl
		<< "+ c-bound:\t\t" << opt->cbound << std::endl
		<< "+ Suborder:\t\t" << opt->subprobOrder << " (" << subprob_order[opt->subprobOrder] << ")" << std::endl
		<< "+ Random seed:\t\t" << opt->seed << std::endl
		<< "+ Algorithm:\t\t" << algo << std::endl
		<< "+ Heuristic:\t\t" << heur << std::endl
		<< "+ Verbose:\t\t" << (opt->verbose ? "yes" : "no") << std::endl
		<< "+ Positive:\t\t" << (opt->positive ? "yes" : "no") << std::endl
		<< "+ Incremental:\t\t" << (opt->incremental ? "yes" : "no") << std::endl
		<< "+ Separator:\t\t" << (opt->separator ? "bounded" : "not bounded") << std::endl
		<< "+ UPBO only:\t\t" << (opt->upperBoundOnly ? "yes" : "no") << std::endl
		<< "+ Rotate limit:\t\t" << opt->rotateLimit << std::endl
		<< "+ JG iter:\t\t" << opt->jglp << std::endl
		<< "+ JG time:\t\t" << opt->jglps << " sec" << std::endl
		<< "+ BP iter:\t\t" << opt->mplp << std::endl
		<< "+ BP time:\t\t" << opt->mplps << std::endl
		<< "+ VO iter:\t\t" << opt->orderingIter << std::endl
		<< "+ VO time:\t\t" << opt->orderingTime << " sec" << std::endl
		<< "+ Propagation:\t\t" << prop << std::endl
		<< "+ Time limit:\t\t" << opt->timeLimit << " sec" << std::endl
		<< "+ Memory time:\t\t" << opt->memoryLimit << " GB" << std::endl
	;

	std::cout << oss2.str();

	// create the output file
	std::ofstream out;
	out.open(opt->outputFile.c_str(), std::ofstream::out);
	if (out.fail()) {
		std::cerr << "Cannot create the output file " << opt->outputFile << std::endl;
		exit(EXIT_FAILURE);
	}
	out << strTaskName << std::endl;
	out.close();

	// create the solver
	if (opt->algorithm == ALGO_BRAOBB) {
		BRAOBB s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_ANY_AAOBF) {
		AnyAAOBF s(opt);
		s.init();
		s.solve();
	}

	delete opt;
	return EXIT_SUCCESS;
}


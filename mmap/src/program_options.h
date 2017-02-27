// programoptions.h -- Command line options (using BOOST)

/*
 * Author: Radu Marinescu
 *
 * Copyright (c) IBM Research, 2012
 *
 * This program is under IBM license. Do not distribute.
 *
 */


#ifndef IBM_ANYTIME_PROGRAM_OPTIONS_H_
#define IBM_ANYTIME_PROGRAM_OPTIONS_H_

#include <boost/program_options.hpp>
#include <iostream>
#include <sstream>

#include "base.h"

namespace po = boost::program_options;

struct ProgramOptions {

	std::string executableName;		// program name
	double timeLimit;				// time limit (in seconds)
	double memoryLimit;				// memory limit (in Gigs)
	int ibound;						// mini-bucket i-bound
	int cbound;						// adaptive cache bound
	algorithm_t algorithm;			// algorithm
	std::string inputFile;			// input file
	std::string evidenceFile;		// evidence file
	std::string outputFile;			// output file
	heuristic_t heuristic;			// heuristic generator
	std::string orderingFile;		// constrained elimination ordering file
	std::string pseudotreeFile;		// pseudo tree file
	std::string mapFile;			// MAP variables
	int seed; 						// the seed for the random number generator
	int subprobOrder; 				// subproblem ordering, integers defined in base.h
	bool verbose;					// verbosity
	bool positive;					// force positive probabilities
	size_t rotateLimit;				// BRAOBB how many nodes to expand per subproblem before rotating: default 1000
	double percMapVars;				// Percentage of MAP variables (0.3 by default)
	int numInstances;				// Number of problem instances
	format_t format;				// Input file format
	bool incremental;				// Incremental heuristic updates
	bool separator; 				// Separator size is bounded by i-bound
	bool exactSubprob;				// Check if subproblem is solved exactly by MBE
	bool upperBoundOnly;			// Computes the initial upper bound
	int jglp;						// iterations for JGLP
	int jglps;						// time limit in seconds for JGLP
	int threads;					// number of threads for paralle processing
	int mplp;						// iterations for MPLP (preprocessor)
	int mplps;						// time limit in seconds for MPLP
	int propagation;				// constraint propagation (0 - none[default], 1 - unit, 2 - full)
	int orderingIter;
	int orderingTime;
	size_t cacheSize;				// cache size used by RBFAOO (in KB): default 1GB
	double overestimation;			// RBFAOO overestimation
	evaluation_t evaluation;		// MAP assignment evaluation (VE or AO)
	double weight;					// initial weight used by anytime algorithms
	weight_schedule_t weightSched;	// weight update schedule
	double weightParamK;			// weight update schedule param k
	double weightParamD;			// weight update schedule param d
	size_t cutoff;					// every 'cutoff' number of nodes call a dfs probe (hybrid dfs/bfs search)
	probing_t probing;				// probe type (for LAOBF, LRBFAOO)
	double probability;				// probability for selecting a tip node outside the current PST
	size_t afseStep;				// ASFE cardinality step size (default 1)

public:

	// default constructor
	ProgramOptions();
};

ProgramOptions* parseCommandLine(int argc, char** argv);

inline ProgramOptions::ProgramOptions() :
		timeLimit(NONE), memoryLimit(80), ibound(2), cbound(NONE), algorithm(ALGO_AOBB),
		heuristic(HEUR_MBE), seed(NONE), subprobOrder(NONE),
		verbose(false), positive(false), rotateLimit(1000),
		percMapVars(0.3), numInstances(1), format(FORMAT_UAI),
		incremental(false), separator(false), exactSubprob(false),
		      upperBoundOnly(false), jglp(-1), jglps(-1), threads(2), mplp(-1), mplps(-1), 
		propagation(CP_NONE),
		orderingIter(500), orderingTime(3), cacheSize(1024*1024),
		overestimation(1.0), evaluation(EVAL_VE),
		weight(1.0), weightSched(WEIGHT_SCHED_SQRT), weightParamK(1.0), weightParamD(1.0),
		cutoff(1000), probing(PROBE_BNB), probability(0.5), afseStep(1) {}

#endif /* IBM_MAP_PROGRAMOPTIONS_H_ */

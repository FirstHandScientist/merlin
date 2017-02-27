/*
 * main.cpp
 *
 *  Created on: Mar 29, 2013
 *      Author: radu
 */

#include "base.h"

#include "map_generator.h"
#include "map_evaluator.h"
#include "map_summation.h"

#include "ao.h"
#include "aobb.h"
#include "braobb.h"
#include "bb_park.h"
#include "bb_yuan.h"
#include "bbbt.h"

#include "astar.h"
#include "aobf.h"
#include "aobf2.h"
#include "rbfaoo.h"
#include "parallel_rbfaoo.h"

#include "any_bbbt.h"
#include "any_waobf.h"
#include "any_waobf2.h"
#include "any_wrbfaoo.h"
#include "any_wraobf.h"
#include "any_wraobf2.h"
#include "any_arastar.h"
#include "any_laobf.h"
#include "any_saobf.h"
#include "any_aaobf.h"
#include "any_wrbfs.h"
#include "any_wrbfs2.h"
#include "any_bbmb.h"
#include "any_ldfs.h"
#include "any_afse.h"

#include "any_ldfs_bounded.h"
#include "any_bbmb_bounded.h"

int main(int argc, char **argv) {

	// Prologue
	std::ostringstream oss;
	std::string version = "ibm-map ";
	version += VERSIONINFO;

	oss << "------------------------------------------------------------------" << std::endl
		<< version << std::endl
		<< "Using Boost "
		<< BOOST_VERSION / 100000     << "."  // major version
		<< BOOST_VERSION / 100 % 1000 << "."  // minior version
		<< BOOST_VERSION % 100                // patch level
		<< std::endl
		<< "Using minisat/zchaff"
		<< std::endl
		<< std::endl
		<< "  by Radu Marinescu, IBM Research - Ireland" << std::endl
		<< "  (C) Copyright IBM Corp. 2015, 2016" << std::endl
		<< "------------------------------------------------------------------" << std::endl;

	std::cout << oss.str();

	// Reprint command line
	for (int i = 0; i < argc; ++i) {
		std::cout << argv[i] << " ";
	}
	std::cout << endl;

#ifdef UAI_COMPETITION
	return uai(argc, argv);
#endif

	// parse command line
	ProgramOptions* opt = parseCommandLine(argc, argv);
	if (!opt) {
		return EXIT_FAILURE;
	}

	// set the seed
	if (opt->seed == NONE) {
		opt->seed = time(0);
	}
	rand::seed(opt->seed);

	// Get the algorithm name and heuristic
	std::string algo, heur, selection, propagation;
	switch (opt->heuristic) {
		case HEUR_MBE: heur = "mbe"; break;
		case HEUR_WMB_MM: heur = "wmb-mm"; break;
		case HEUR_WMB_JG: heur = "wmb-jg"; break;
		case HEUR_JOIN_TREE: heur = "join-tree"; break;
		case HEUR_MINI_CLUSTER_TREE: heur = "mini-cluster-tree"; break;
		default: heur = "UNKNOWN"; break;
	}
	switch (opt->algorithm) {
	case ALGO_BBBT: algo = "bbbt"; break;
	case ALGO_AOBB: algo = "aobb"; break;
	case ALGO_AOBF: algo = "aobf (AND contexts)"; break;
	case ALGO_AOBF2: algo = "aobf2 (OR contexts)"; break;
	case ALGO_BRAOBB: algo = "braobb"; break;
	case ALGO_PARK: algo = "park-bb"; break;
	case ALGO_YUAN: algo = "yuan-bb"; break;
	case ALGO_ASTAR: algo = "astar"; break;
	case ALGO_RBFAOO: algo = "rbfaoo"; break;
	case ALGO_SPRBFAOO: algo = "sprbfaoo"; break;
	case ALGO_ANY_BBBT: algo = "any-bbbt"; break;
	case ALGO_ANY_BBMB: algo = "any-bbmb"; break;
	case ALGO_ANY_WAOBF: algo = "any-waobf (restart)"; break;
	case ALGO_ANY_WAOBF2: algo = "any-waobf2 (restart)"; break;
	case ALGO_ANY_WRBFAOO: algo = "any-wrbfaoo (restart)"; break;
	case ALGO_ANY_WRAOBF: algo = "any-wraobf (repair)"; break;
	case ALGO_ANY_WRAOBF2: algo = "any-wraobf2 (repair)"; break;
	case ALGO_ANY_ARASTAR: algo = "any-arastar (repair)"; break;
	case ALGO_ANY_LAOBF: algo = "any-laobf"; break;
	case ALGO_ANY_SAOBF: algo = "any-saobf"; break;
	case ALGO_ANY_WRBFS: algo = "any-wrbfs"; break;
	case ALGO_ANY_WRBFS2: algo = "any-wrbfs2 (caching)"; break;
	case ALGO_ANY_AAOBF: algo = "any-aaobf"; break;
	case ALGO_ANY_LDFS: algo = "any-ldfs"; break;
	case ALGO_AFSE: algo = "afse"; break;
	case ALGO_ANY_LDFS_BOUNDED: algo = "any-ldfs-bounded"; break;
	case ALGO_ANY_BBMB_BOUNDED: algo = "any-bbmb-bounded"; break;
	case ALGO_GENERATOR: algo = "map-generator"; break;
	case ALGO_EVALUATOR: algo = "map-evaluator"; break;
	default: algo = "UNKNOWN"; break;
	}

	switch (opt->propagation) {
	case CP_UNIT: propagation = "unit"; break;
	case CP_FULL: propagation = "full"; break;
	default: propagation = "none"; break;
	}

	std::ostringstream oss2;
	oss2 << "+ i-bound:        " << opt->ibound << std::endl
		<< "+ c-bound:        " << opt->cbound << std::endl
		<< "+ Suborder:       " << opt->subprobOrder << " (" << subprob_order[opt->subprobOrder] << ")" << std::endl
		<< "+ Random seed:    " << opt->seed << std::endl
		<< "+ Algorithm:      " << algo << std::endl
		<< "+ Heuristic:      " << heur << std::endl
		<< "+ Verbose:        " << (opt->verbose ? "yes" : "no") << std::endl
		<< "+ Positive:       " << (opt->positive ? "yes" : "no") << std::endl
		<< "+ Incremental:    " << (opt->incremental ? "yes" : "no") << std::endl
		<< "+ Separator:      " << (opt->separator ? "bounded" : "not bounded") << std::endl
		<< "+ UPBO only:      " << (opt->upperBoundOnly ? "yes" : "no") << std::endl
		<< "+ Rotate limit:   " << opt->rotateLimit << std::endl
		<< "+ JG iter:        " << opt->jglp << std::endl
		<< "+ JG time:        " << opt->jglps << " sec" << std::endl
		<< "+ Propagation:    " << propagation << std::endl
		<< "+ Overestimation: " << opt->overestimation << std::endl
		<< "+ Cache size:     " << (opt->cacheSize/1024) << " MB" << std::endl
		<< "+ Schedule:       " << opt->weightSched << " (" << schedule_names[opt->weightSched] << ")" << std::endl
		<< "+   weight:       " << opt->weight << std::endl
		<< "+   param-k:      " << opt->weightParamK << std::endl
		<< "+   param-d:      " << opt->weightParamD << std::endl
		<< "+ Probing:        " << opt->probing << " (" << probing_names[opt->probing] << ")" << std::endl
		<< "+ Cutoff:         " << opt->cutoff << std::endl
		<< "+ Probability:    " << opt->probability << std::endl
		<< "+ Step (afse):    " << opt->afseStep << std::endl;
	;

	std::cout << oss2.str();

	// run the algorithm
	if (opt->algorithm == ALGO_AO) {
		AO s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_BBBT) {
		BBBT s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_AOBB) {
		AOBB s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_BRAOBB) {
		BRAOBB s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_PARK) {
		BBPark s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_YUAN) {
		BBYuan s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_ASTAR) {
		Astar s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_AOBF) {
		AOBF s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_AOBF2) {
		AOBF2 s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_RBFAOO) {
		RBFAOO s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_ANY_BBBT) {
		AnyBBBT s(opt);
		s.init();
		s.solve();
	}else if (opt->algorithm == ALGO_ANY_BBMB) {
		AnyBBMB s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_ANY_WAOBF) {
		AnyWAOBF s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_ANY_WAOBF2) {
		AnyWAOBF2 s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_ANY_WRBFAOO) {
		AnyWRBFAOO s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_ANY_WRAOBF) {
		AnyWRAOBF s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_ANY_WRAOBF2) {
		AnyWRAOBF2 s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_ANY_ARASTAR) {
		AnyARAstar s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_ANY_LAOBF) {
		AnyLAOBF s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_ANY_SAOBF) {
		AnySAOBF s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_ANY_WRBFS) {
		AnyWRBFS s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_ANY_WRBFS2) {
		AnyWRBFS2 s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_ANY_AAOBF) {
		AnyAAOBF s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_ANY_LDFS) {
		AnyLDFS s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_AFSE) {
		AnyAFSE s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_ANY_LDFS_BOUNDED) {
		AnyLDFS_Bounded s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_ANY_BBMB_BOUNDED) {
		AnyBBMB_Bounded s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_GENERATOR) {
		MapGenerator s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_EVALUATOR) {
		MapEvaluator s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_SPRBFAOO) {
		ParallelRBFAOO s(opt);
		s.init();
		s.solve();
	} else if (opt->algorithm == ALGO_SUMMATION) {
		MapSummation s(opt);
		s.init();
		s.solve();
	}

	// free memory
	delete opt;

	return EXIT_SUCCESS;
}




// programoptions.cpp -- Command line options (using BOOST)

/*
 * Author: Radu Marinescu
 *
 * Copyright (c) IBM Research, 2012
 *
 * This program is under IBM license. Do not distribute.
 *
 */

#include "program_options.h"

ProgramOptions* parseCommandLine(int argc, char** argv) {

	ProgramOptions* opt = new ProgramOptions;

	// executable name
	opt->executableName = argv[0];

	try {

		po::options_description desc("Valid options");
		desc.add_options()
			("input-file,f", po::value<std::string>(), "path to problem file (required)")
			("output-file,o", po::value<std::string>(), "path to output file (required)")
			("algorithm,a", po::value<std::string>(), "inference algorithm (required)")
			("heuristic,H", po::value<std::string>(), "lower bounding heuristic")
			("ibound,i", po::value<int>(), "mini-bucket i-bound")
			("cbound,c", po::value<int>()->default_value(NONE), "adaptive cache bound")
			("time-limit,L", po::value<int>(), "time limit in seconds")
			("ordering-file,E", po::value<std::string>(), "constrained elimination ordering file (input)")
			("pseudotree-file,P", po::value<std::string>(), "pseudo tree file (output)")
			("map-file,M", po::value<std::string>(), "map variables file (output)")
			("evidence-file,x", po::value<std::string>(), "evidence (input)")
			("suborder,r", po::value<int>()->default_value(0), "subproblem order (0:width-inc 1:width-dec 2:heur-inc 3:heur-dec)")
			("seed,s", po::value<int>(), "seed for random number generator, time() otherwise")
			("verbose,v", "enable verbose mode")
			("positive,p", "enable full positive distribution mode")
			("incremental,N", "incremental heuristic updates mode")
			("separator,S", "separator size bounded by i-bound")
			("exact-subprob,X", "subproblem solved exactly by MBE (subprob width <= i-bound)")
			("upper-bound,U", "computes only the initial upper bound")
			("rotate-limit,R", po::value<size_t>()->default_value(1000), "stack rotation limit")
			("threads,T", po::value<int>(), "number of worker threads (for parallel processing)")
			("map,Q", po::value<double>()->default_value(0.3), "percentage of random map variables")
			("instances,I", po::value<int>()->default_value(1), "number of random MAP problem instances")
			("jglp", po::value<int>(), "number of JGLP iterations")
			("jglps", po::value<int>(), "time limit in sec for JGLP")
			("propagation,G", po::value<std::string>(), "constraint propagation")
			("format,F", po::value<std::string>(), "input file format")
			("cache-size,C", po::value<std::string>()->default_value("1g"), "cache size (1m,10m,100m,1g,10g)")
			("overestimation,O", po::value<double>()->default_value(1.0), "overstimation")
			("evaluation,V", po::value<std::string>(), "evaluation method (ve or ao)")
			("weight,w", po::value<double>()->default_value(1), "initial weight")
			("weight-sched,S", po::value<std::string>(), "weight update schedule")
			("weight-paramk,K", po::value<double>(), "weight update parameter K")
			("weight-paramd,d", po::value<double>(), "weight update parameter D")
			("cutoff,u", po::value<size_t>(), "cutoff value (integer)")
			("probing", po::value<std::string>(), "DFS probing policy: dfs, bnb")
			("probability", po::value<double>()->default_value(0.5), "probability for selecting tip nodes outside PST")
			("step", po::value<size_t>(), "AFSE cardinality set step")
			("help,h", "produces this help message");

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		// parse help
		if (vm.count("help")) {
			std::cout << std::endl << desc << std::endl;
			delete opt;
			exit(0);
		}

		// parse verbose switch
		if (vm.count("verbose")) {
			opt->verbose = true;
		}

		// parse positive switch
		if (vm.count("positive")) {
			opt->positive = true;
		}

		// parse incremental switch
		if (vm.count("incremental")) {
			opt->incremental = true;
		}

		// parse separator switch
		if (vm.count("separator")) {
			opt->separator = true;
		}

		// parse exact subproblem switch
		if (vm.count("exact-subprob")) {
			opt->exactSubprob = true;
		}

		// parse upper bound only switch
		if (vm.count("upper-bound")) {
			opt->upperBoundOnly = true;
		}

		// parse input file
		if (!vm.count("input-file")) {
			std::cerr << "No or invalid arguments given, " << "call with '" << argv[0]
					<< " --help' for full list." << std::endl;
			delete opt;
			return NULL;
		}

		opt->inputFile = vm["input-file"].as<std::string>();

		// parse output file
		if (vm.count("output-file")) {
			opt->outputFile = vm["output-file"].as<std::string>();
		}

		// parse algorithm
		if (vm.count("algorithm")) {
			std::string alg = vm["algorithm"].as<std::string>();
			if (alg.compare("aobf") == 0) {
				opt->algorithm = ALGO_AOBF;
			} else if (alg.compare("aobf2") == 0) {
				opt->algorithm = ALGO_AOBF2;
			} else if (alg.compare("aobb") == 0) {
				opt->algorithm = ALGO_AOBB;
			} else if (alg.compare("braobb") == 0) {
				opt->algorithm = ALGO_BRAOBB;
			} else if (alg.compare("ao") == 0) {
				opt->algorithm = ALGO_AO;
			} else if (alg.compare("bb") == 0) {
				opt->algorithm = ALGO_BB;
			} else if (alg.compare("bbbt") == 0) {
				opt->algorithm = ALGO_BBBT;
			} else if (alg.compare("park") == 0) {
				opt->algorithm = ALGO_PARK;
			} else if (alg.compare("yuan") == 0) {
				opt->algorithm = ALGO_YUAN;
			} else if (alg.compare("astar") == 0) {
				opt->algorithm = ALGO_ASTAR;
			} else if (alg.compare("rbfaoo") == 0) {
				opt->algorithm = ALGO_RBFAOO;
			} else if (alg.compare("sprbfaoo") == 0) {
				opt->algorithm = ALGO_SPRBFAOO;
			} else if (alg.compare("any-bbbt") == 0) {
				opt->algorithm = ALGO_ANY_BBBT;
			} else if (alg.compare("any-bbmb") == 0) {
				opt->algorithm = ALGO_ANY_BBMB;
			} else if (alg.compare("any-waobf") == 0) {
				opt->algorithm = ALGO_ANY_WAOBF;
			} else if (alg.compare("any-waobf2") == 0) {
				opt->algorithm = ALGO_ANY_WAOBF2;
			} else if (alg.compare("any-wrbfaoo") == 0) {
				opt->algorithm = ALGO_ANY_WRBFAOO;
			} else if (alg.compare("any-wraobf") == 0) {
				opt->algorithm = ALGO_ANY_WRAOBF;
			} else if (alg.compare("any-wraobf2") == 0) {
				opt->algorithm = ALGO_ANY_WRAOBF2;
			} else if (alg.compare("any-arastar") == 0) {
				opt->algorithm = ALGO_ANY_ARASTAR;
			} else if (alg.compare("any-laobf") == 0) {
				opt->algorithm = ALGO_ANY_LAOBF;
			} else if (alg.compare("any-saobf") == 0) {
				opt->algorithm = ALGO_ANY_SAOBF;
			} else if (alg.compare("any-wrbfs") == 0) {
				opt->algorithm = ALGO_ANY_WRBFS;
			} else if (alg.compare("any-wrbfs2") == 0) {
				opt->algorithm = ALGO_ANY_WRBFS2;
			} else if (alg.compare("any-aaobf") == 0) {
				opt->algorithm = ALGO_ANY_AAOBF;
			} else if (alg.compare("any-ldfs") == 0) {
				opt->algorithm = ALGO_ANY_LDFS;
			} else if (alg.compare("afse") == 0) {
				opt->algorithm = ALGO_AFSE;
			} else if (alg.compare("any-ldfs-bounded") == 0) {
				opt->algorithm = ALGO_ANY_LDFS_BOUNDED;
			} else if (alg.compare("any-bbmb-bounded") == 0) {
				opt->algorithm = ALGO_ANY_BBMB_BOUNDED;
			} else if (alg.compare("generator") == 0) {
				opt->algorithm = ALGO_GENERATOR;
			} else if (alg.compare("evaluator") == 0) {
				opt->algorithm = ALGO_EVALUATOR;
			} else if (alg.compare("summation") == 0) {
				opt->algorithm = ALGO_SUMMATION;
			} else {
				std::cerr << "Unknown algorithm: " << alg << std::endl;
				delete opt;
				return NULL;
			}
		}

		// parse propagation
		if (vm.count("propagation")) {
			std::string cp = vm["propagation"].as<std::string>();
			if (cp.compare("none") == 0) {
				opt->propagation = CP_NONE;
			} else if (cp.compare("unit") == 0) {
				opt->propagation = CP_UNIT;
			} else if (cp.compare("full") == 0) {
				opt->propagation = CP_FULL;
			} else {
				std::cerr << "Unknown propagation: " << cp << std::endl;
				delete opt;
				return NULL;
			}
		}

		// parse mini-bucket i-bound
		if (vm.count("ibound")) {
			opt->ibound = vm["ibound"].as<int>();
		}

		// parse adaptive cache bound
		if (vm.count("cbound")) {
			opt->cbound = vm["cbound"].as<int>();
		}

		// parse the time limit
		if (vm.count("time-limit")) {
			opt->timeLimit = vm["time-limit"].as<int>();
		}

		// parse the rotate limit
		if (vm.count("rotate-limit")) {
			opt->rotateLimit = vm["rotate-limit"].as<size_t>();
		}

		// parse heuristic generator
		if (vm.count("heuristic")) {
			std::string heur = vm["heuristic"].as<std::string>();
			if (heur.compare("mbe") == 0) {
				opt->heuristic = HEUR_MBE;
			} else if (heur.compare("wmb-mm") == 0) {
				opt->heuristic = HEUR_WMB_MM;
			} else if (heur.compare("wmb-jg") == 0) {
				opt->heuristic = HEUR_WMB_JG;
			} else if (heur.compare("jt") == 0) {
				opt->heuristic = HEUR_JOIN_TREE;
			} else if (heur.compare("mcte") == 0) {
				opt->heuristic = HEUR_MINI_CLUSTER_TREE;
			} else {
				std::cerr << "Unknown heuristic generator." << std::endl;
				delete opt;
				return NULL;
			}
		}

		// parse file format
		if (vm.count("format")) {
			std::string fmt = vm["format"].as<std::string>();
			if (fmt.compare("uai") == 0) {
				opt->format = FORMAT_UAI;
			} else if (fmt.compare("ergo") == 0) {
				opt->format = FORMAT_ERGO;
			} else {
				std::cerr << "Unknown input file format." << std::endl;
				delete opt;
				return NULL;
			}
		}

		// parse the elimination ordering file (constrained; MAP)
		if (vm.count("ordering-file")) {
			opt->orderingFile = vm["ordering-file"].as<std::string>();
		}

		// parse the pseudo tree file
		if (vm.count("pseudotree-file")) {
			opt->pseudotreeFile = vm["pseudotree-file"].as<std::string>();
		}

		// parse the MAP variables file
		if (vm.count("map-file")) {
			opt->mapFile = vm["map-file"].as<std::string>();
		}

		// parse the evidence file
		if (vm.count("evidence-file")) {
			opt->evidenceFile = vm["evidence-file"].as<std::string>();
		}

		// parse the subproblem order
		if (vm.count("suborder")) {
			opt->subprobOrder = vm["suborder"].as<int>();
			if (opt->subprobOrder < 0 || opt->subprobOrder > 3) {
				std::cout << std::endl << desc << std::endl;
				exit(0);
			}
		}

		// parse the random generator seed
		if (vm.count("seed")) {
			opt->seed = vm["seed"].as<int>();
		}

		// parse the number of random MAP variables
		if (vm.count("map")) {
			opt->percMapVars = vm["map"].as<double>();
		}

		// parse the number of random MAP problem instances
		if (vm.count("instances")) {
			opt->numInstances = vm["instances"].as<int>();
		}

		// parse the number of JGLP iterations
		if (vm.count("jglp")) {
			opt->jglp = vm["jglp"].as<int>();
		}

		// parse threads
		if (vm.count("threads")) {
			opt->threads = vm["threads"].as<int>();
		}

		// parse the number of random MAP problem instances
		if (vm.count("jglps")) {
			opt->jglps = vm["jglps"].as<int>();
		}

		// parse the cache size (RBFAOO)
		if (vm.count("cache-size")) {
			std::string str = vm["cache-size"].as<std::string>();
			size_t pos_mega = str.find('m');
			size_t pos_giga = str.find('g');
			if (pos_mega != std::string::npos) { // cache size in MB
				assert( pos_giga == std::string::npos );
				int num = std::atoi(str.substr(0, pos_mega).c_str());
				opt->cacheSize = num * 1024;
			} else if (pos_giga != std::string::npos) { // cache size in GB
				assert( pos_mega == std::string::npos );
				int num = std::atoi(str.substr(0, pos_giga).c_str());
				opt->cacheSize = num * 1024 * 1024;
			} else { // unknown cache size
				std::cerr << "Unknown cache size (use MB or GB)." << std::endl;
				delete opt;
				return NULL;
			}
		}

		// parse the overestimation (RBFAOO)
		if (vm.count("overestimation")) {
			opt->overestimation = vm["overestimation"].as<double>();
		}

		// parse the evaluation
		if (vm.count("evaluation")) {
			std::string eval = vm["evaluation"].as<std::string>();
			if (eval.compare("ve") == 0) {
				opt->evaluation = EVAL_VE;
			} else if (eval.compare("ao") == 0) {
				opt->evaluation = EVAL_AO;
			} else {
				std::cerr << "Unknown MAP assignment evaluation method (use ve or ao)." << std::endl;
				delete opt;
				return NULL;
			}
		}

		// parse the weight
		if (vm.count("weight")) {
			opt->weight = vm["weight"].as<double>();
		}

		// parse the weight delta
		if (vm.count("weight-paramk")) {
			opt->weightParamK = vm["weight-paramk"].as<double>();
		}

		// parse the weight denominator
		if (vm.count("weight-paramd")) {
			opt->weightParamD = vm["weight-paramd"].as<double>();
		}

		// parse weight update schedule
		if (vm.count("weight-sched")) {
			std::string sched = vm["weight-sched"].as<std::string>();
			if (sched.compare("subtract") == 0) {
				opt->weightSched = WEIGHT_SCHED_SUBTRACT;
			} else if (sched.compare("divide") == 0) {
				opt->weightSched = WEIGHT_SCHED_DIVIDE;
			} else if (sched.compare("sqrt") == 0) {
				opt->weightSched = WEIGHT_SCHED_SQRT;
			} else if (sched.compare("inverse") == 0) {
				opt->weightSched = WEIGHT_SCHED_INVERSE;
			} else if (sched.compare("piecewise") == 0) {
				opt->weightSched = WEIGHT_SCHED_PIECEWISE;
			} else {
				std::cerr << "Unknown weight update schedule: " << sched << std::endl;
				delete opt;
				return NULL;
			}
		}

		// parse the cutoff parameter
		if (vm.count("cutoff")) {
			opt->cutoff = vm["cutoff"].as<size_t>();
		}

		// parse probing policy
		if (vm.count("probing")) {

			std::string la = vm["probing"].as<std::string>();
			if (la.compare("dfs") == 0) {
				opt->probing = PROBE_DFS;
			} else if (la.compare("bnb") == 0) {
				opt->probing = PROBE_BNB;
			}
		}

		// parse the probability parameter
		if (vm.count("probability")) {
			opt->probability = vm["probability"].as<double>();
		}

		// parse the AFSE cardinality step size
		if (vm.count("step")) {
			opt->afseStep = vm["step"].as<size_t>();
		}

	} catch (std::exception& e) {
		std::cerr << "error: " << e.what() << std::endl;
		delete opt;
		return NULL;
	} catch (...) {
		std::cerr << "Exception of unknown type!" << std::endl;
		delete opt;
		exit(1);
	}

	return opt;
}


/*
 ============================================================================
 Name        : demo.cpp
 Author      : Radu Marinescu
 Version     :
 Copyright   : Copyright (c) IBM Corp. 2015
 Description : Uses shared library to print greeting
 To run the resulting executable the LD_LIBRARY_PATH must be
 set to ${project_loc}/merlin/.libs
 Alternatively, libtool creates a wrapper shell script in the
 build directory of this program which can be used to run it.
 Here the script will be called exampleProgram.
 ============================================================================
 */
#include <string.h>
#include <string>
#include <iostream>
#include "merlin.h"
#include "util.h"

// Debugging only
void demo_debug() {

	// Init parameters
	unsigned int ibound = 10;
	unsigned int iterations = 100;
	unsigned int samples = 1000;
	const char* inputFile = "/home/radu/workspace/ibm-map/data/50-12-5.uai";
	const char* evidenceFile = "/home/radu/git/merlin/example/simple5.evid";
	const char* queryFile = "/home/radu/workspace/ibm-map/data/50-12-5.uai.map.N";
	const char* outputFile = "/home/radu/git/merlin/example/simple5.out";

//	const char* inputFile = "/home/radu/Downloads/chain.xmlbif.UAI";
//	const char* evidenceFile = "/home/radu/Downloads/chain.xmlbif.UAI.EVID";
//	const char* queryFile = "/home/radu/Downloads/chain.xmlbif.UAI.QUERY";

	// MMAP task
//	run(inputFile, evidenceFile, queryFile, outputFileMMAP, "MAR", ibound, iterations);

	// Initialize the Merlin engine
	Merlin eng;
	eng.set_param_ibound(ibound);
	eng.set_param_iterations(iterations);
	eng.set_param_samples(samples);
	eng.read_model(inputFile);
	eng.read_evidence(evidenceFile);
//	eng.read_query(queryFile);
	eng.set_task(MERLIN_TASK_PR);
	eng.set_algorithm(MERLIN_ALGO_WMB);
	eng.run();
}

// Demo the black-box run
void demo_run() {

	// Init parameters
	unsigned int ibound = 4;
	unsigned int iterations = 300;
	const char* inputFile = "pedigree1.uai";
	const char* evidenceFile = "pedigree1.evid";
	const char* queryFile = "pedigree1.map";
	const char* outputFileMAR = "pedigree1.MAR.out";
	const char* outputFileMAP = "pedigree1.MAP.out";
	const char* outputFileMMAP = "pedigree1.MMAP.out";

	// MAR task
	run(inputFile, evidenceFile, "", outputFileMAR, "MAR", ibound, iterations);

	// MAP task
	run(inputFile, evidenceFile, "", outputFileMAP, "MAP", ibound, iterations);

	// MMAP task
	run(inputFile, evidenceFile, queryFile, outputFileMMAP, "MMAP", ibound, iterations);

}

// Demo the Merlin API
void demo_api() {

	// Init parameters
	unsigned int ibound = 4;
	unsigned int iterations = 300;
	const char* model_file = "pedigree1.uai";
	const char* evid_file = "pedigree1.evid";
	const char* query_file = "pedigree1.map";


	// Initialize the Merlin engine
	Merlin eng;
	eng.set_param_ibound(8);
	eng.set_param_iterations(300);
	eng.read_model(model_file);
	eng.read_evidence(evid_file);

	// Solve a MAR task
	eng.set_task(MERLIN_TASK_MAR);
	eng.set_algorithm(MERLIN_ALGO_WMB);
	eng.run();

	// Solve a MAP task
	eng.set_task(MERLIN_TASK_MAP);
	eng.set_algorithm(MERLIN_ALGO_WMB);
	eng.run();

	// Solve a MMAP task
	eng.read_query(query_file);
	eng.set_task(MERLIN_TASK_MMAP);
	eng.run();
}

// Demo: convert a weigted cnf file into a factor graph (UAI format)
void demo_wcnf2uai(const char* file_name) {
	Merlin eng;

	std::string f = std::string(file_name) + ".uai";

	eng.read_model(file_name, MERLIN_INPUT_WCNF);
	std::cout << "Read wcnf file: " << file_name << std::endl;
	eng.write_model(f.c_str());
	std::cout << "Wrote uai file: " << f << std::endl;
}

// Run the MAR solver for the UAI competition
// ./solver <input-model-file> <input-evidence-file> <input-query-file> <PR|MPE|MAR|MMAP>
void uai(int argc, char** argv) {
	if (argc != 5) {
		std::cout << "UAI-2016 Inference evaluation format is:" << std::endl;
		std::cout << "./solver <input-model-file> <input-evidence-file> <input-query-file> <PR|MPE|MAR|MMAP>" << std::endl;
		exit(1);
	}

	const char* input_model_file = argv[1];		// input file
	const char* input_evid_file = argv[2]; 		// input evidence file
	const char* input_query_file = argv[3];		// input query file (for MMAP)
	const char* input_task = argv[4];			// input task name

	// Select the inference task
	size_t task = MERLIN_TASK_PR;
	size_t START_IBOUND = 2; 					// default ibound is 2
	size_t END_IBOUND = 30;						// MAX ibound to run
	size_t iterations = 20;						// default iterations is 20
	if (strcmp(input_task, "PR") == 0) {
		task = MERLIN_TASK_PR;
	} else if (strcmp(input_task, "MAR") == 0) {
		task = MERLIN_TASK_MAR;
	} else if (strcmp(input_task, "MPE") == 0) {
		task = MERLIN_TASK_MAP;
	} else if (strcmp(input_task, "MMAP") == 0) {
		task = MERLIN_TASK_MMAP;
	}

	// Setup Merlin engine
	Merlin eng;
	eng.set_task(task);
	eng.set_algorithm(MERLIN_ALGO_WMB);
	eng.read_model(input_model_file);
	eng.read_evidence(input_evid_file);
	if (task == MERLIN_TASK_MMAP) {
		eng.read_query(input_query_file);
	}

	for (size_t ib = START_IBOUND; ib <= END_IBOUND; ++ib) {
		eng.set_param_ibound(ib);
		eng.set_param_iterations(iterations);
		eng.run(); // should "append" instead of overwrite the output file.
	}
}

// Main
int main(int argc, char** argv) {

	// Call the "uai" solver
//	uai(argc, argv);

	// Call the 'run' function
	//demo_run();

	// Call Merlin API
	//demo_api();

	// Call the 'debug' function
	demo_debug();

//	demo_convert(argv[1]);

//	demo_fileformat();

//	demo_merlin(argv[1]);

//	demo_wcnf2uai(argv[1]);

//	demo_wcnf2net(argv[1]);

	return 0;
}



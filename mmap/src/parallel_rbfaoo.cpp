/*
 * paralllel_rbfaoo.cpp
 *
 */

#include "parallel_rbfaoo.h"
#include "parallel_rbfaoo_worker.h"
#include "mini_bucket_heur.h"
#include "mini_bucket_mm_heur.h"
#include "mini_bucket_jglp_heur.h"
#include <boost/thread.hpp>

class ParallelRBFAOOWorkerThread {

private:
	ParallelRBFAOOWorker *m_worker;

public:
	ParallelRBFAOOWorkerThread(ParallelRBFAOOWorker *worker) : m_worker(worker) {}
    
	void operator()() {
#if 111
		cpu_set_t mask;
		CPU_ZERO(&mask);
		int id = m_worker->getThreadId();
		CPU_SET(id, &mask);
		if (sched_setaffinity(0, sizeof(cpu_set_t), &mask) == -1) {
		  cerr << "Failed to set CPU affinity for thread " << id << endl;;
		  exit(1);
		}
#endif
		m_worker->start();
	}
};


//Kishi: define memory pool using boost to reduce overhead of new/delete
//Kishi: change implementation of TT to support constant-size memory allocation
ParallelRBFAOO::ParallelRBFAOO(ProgramOptions* opt) : Search(opt) {
	m_num_threads = opt->threads;
	m_cache = NULL;
}

ParallelRBFAOO::~ParallelRBFAOO() {
	if (m_cache)
	  delete m_cache;
}

// Kishi: Copied from Radu's code
// initialize
void ParallelRBFAOO::init() {

	// init problem instance
	Search::init();

	Timer tm;
	tm.start();

	// init heuristic generator
	if (m_options->heuristic == HEUR_MBE) {
		std::cout << "Computing MBE heuristic (i=" << m_options->ibound << ") ..." << std::endl;
		m_heuristic.reset(new MiniBucketHeur(m_problem.get(), m_pseudotree.get(),
				m_options, m_options->ibound));
	} else 	if (m_options->heuristic == HEUR_WMB_MM) {
		std::cout << "Computing WMB-MM heuristic (i=" << m_options->ibound << ") ..." << std::endl;
		m_heuristic.reset(new MiniBucketMMHeur(m_problem.get(), m_pseudotree.get(),
				m_options, m_options->ibound));
	} else 	if (m_options->heuristic == HEUR_WMB_JG) {
		std::cout << "Computing WMB-JG heuristic (i=" << m_options->ibound << ") ..." << std::endl;
		m_heuristic.reset(new MiniBucketJglpHeur(m_problem.get(), m_pseudotree.get(),
				m_options, m_options->ibound));
	}

	size_t sz = m_heuristic->build(&m_assignment, true);
	tm.stop();
	m_tmHeuristic = tm.elapsed();

	std::cout << "\tMini bucket finished in " << m_tmHeuristic
			<< " seconds" << std::endl;
	std::cout << "\tUsing " << (sz / (1024*1024.0)) * sizeof(double)
			<< " MBytes of RAM" << std::endl;

	// heuristic might have changed problem functions, pseudotree needs remapping
	m_pseudotree->resetFunctionInfo(m_problem->getFunctions());

	// check if heuristic is accurate
	if (m_heuristic->isAccurate()) {
		std::cout << "Heuristic is exact!" << std::endl;
		m_solutionCost = (m_heuristic->getGlobalUB());
	    m_solved = true;
	}

	// Record time after initialization
	cout << "Initialization complete: " << (m_tmLoad + m_tmHeuristic)
			<< " seconds" << std::endl;

	// check for upper bound computation only
	if (m_options->upperBoundOnly) {
		std::cout << "Done." << std::endl;
		exit(EXIT_SUCCESS);
	}
}

//Kishi: DONE
int ParallelRBFAOO::solve() {

	// check if solved during initialization
	if (m_solved) {
		std::cout << "--------- Solved during initialization ---------" << std::endl;
		std::cout << "Problem name:\t" << m_problem->getName() << std::endl;
		std::cout << "Status:\t\t" << solver_status[0] << std::endl;
		std::cout << "OR nodes:\t" << 0 << std::endl;
		std::cout << "AND nodes:\t" << 0 << std::endl;
		std::cout << "Time elapsed:\t" << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
		std::cout << "Preprocessing:\t" << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
		//std::cout << "Solution:\t" << std::setiosflags(std::ios::fixed) << std::setprecision(4) << m_solutionCost << std::endl;
		std::cout << "Solution:            " << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
		std::cout << "-------------------------------" << std::endl;

		return SEARCH_SUCCESS;
	}

	Zobrist::initOnce(*(m_problem.get()));

	// init search space
	Timer tm;
	tm.start();
	assert( m_options->cacheSize > 0 ); // in kilo bytes
	size_t cache_kilobytes = m_options->cacheSize;
	if (!m_cache) {
		m_cache = new ParallelRbfaooCacheTable;
		m_cache->init(m_num_threads, cache_kilobytes);
	}
	tm.stop();
	std::cout << "Cache allocation complete: " << tm.elapsed() << " seconds" << std::endl;

	// prologue
	std::cout << "--- Starting search ---" << std::endl;
	m_tmSolve = 0;
	m_timer.reset(m_options->timeLimit);
	m_timer.start();

	vector<ParallelRBFAOOWorker *> workers;
	{
		boost::thread_group thr_gp;
	  for (size_t i = 0; i < m_num_threads; i ++) {
	    ParallelRBFAOOWorker *worker = new ParallelRBFAOOWorker(i, &(*m_problem), &(*m_pseudotree), m_options, m_cache, &(*m_heuristic));
		workers.push_back(worker);
		thr_gp.add_thread(new boost::thread(ParallelRBFAOOWorkerThread(worker)));
	  }
	  thr_gp.join_all();
	}
	int res = workers[0]->getSearchResult();
	m_solutionCost = workers[0]->getSolutionCost();
	size_t total_or_node_expansion = 0, total_and_node_expansion = 0;
	for (size_t i = 0; i < m_num_threads; i ++) {
		size_t or_node_expansion = workers[i]->getNumExpandedOR();
		size_t and_node_expansion = workers[i]->getNumExpandedAND()
;		std::cout << "Thread " << i <<": OR node expansion = " << or_node_expansion << std::endl;
		std::cout << "Thread " << i <<": AND node expansion = " << and_node_expansion << std::endl;
		delete workers[i];
		total_or_node_expansion += or_node_expansion;
		total_and_node_expansion += and_node_expansion;
	}
	// stop timer
	m_timer.stop();
	m_tmSolve = m_timer.elapsed();

	// output solution (if found)
	std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::endl;
	std::cout << "--- Search done ---" << std::endl;
	std::cout << "Problem name:\t" << m_problem->getName() << std::endl;
	std::cout << "Status:\t\t" << solver_status[res] << std::endl;
	std::cout << "Total OR nodes:\t" << total_or_node_expansion << std::endl;
	std::cout << "Total AND nodes:\t" << total_and_node_expansion << std::endl;
	std::cout << "Time elapsed:\t" << (m_tmLoad + m_tmHeuristic + m_tmSolve / m_num_threads) << " seconds" << std::endl;
	std::cout << "Preprocessing:\t" << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
	//std::cout << "Solution:\t" << std::setiosflags(std::ios::fixed) << std::setprecision(4) << m_solutionCost << std::endl;
	std::cout << "Solution:            " << 1/ELEM_DECODE(m_solutionCost) << " (" << -m_solutionCost  << ")" << std::endl;

	std::cout << "-------------------------------" << std::endl;

	// clean up
	Zobrist::finishOnce();

	return res;
}

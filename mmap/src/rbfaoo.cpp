/*
 * rbfaoo.cpp
 *
 *  Created on: Apr 2, 2014
 *      Author: radu
 */



#include "rbfaoo.h"
#include "mini_bucket_heur.h"
#include "mini_bucket_mm_heur.h"
#include "mini_bucket_jglp_heur.h"

const float epsilon = 1e-10;

// Radu: DONE
// initialize
void RBFAOO::init() {

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

// Radu: DONE
// Kishi: define memory pool using boost to reduce overhead of new/delete
// Kishi: change implementation of TT to support constant-size memory allocation
RBFAOO::RBFAOO(ProgramOptions* opt) : Search(opt) {

	// Preallocate space for expansion vector.
	m_expand.reserve(65536);
	m_num_expanded_or = 0;
	m_num_expanded_and = 0;
	m_num_expanded = 0;
	m_num_expanded_or_map = 0;
	m_num_expanded_and_map = 0;
	m_prevTimeCheckpoint = 0;
	m_prevNodeCheckpoint = 0;
	m_overestimation = opt->overestimation;
	m_epsilon = opt->weight;
	m_numSumEvals = 0;
}

// Radu: DONE
int RBFAOO::solve() {

	// check if solved during initialization
	if (m_solved) {
		std::cout << "--------- Solved during initialization ---------" << std::endl;
		std::cout << "Problem name:        " << m_problem->getName() << std::endl;
		std::cout << "Status:              " << solver_status[0] << std::endl;
		std::cout << "OR nodes:            " << 0 << std::endl;
		std::cout << "AND nodes:           " << 0 << std::endl;
		std::cout << "OR nodes (MAP):      " << 0 << std::endl;
		std::cout << "AND nodes (MAP):     " << 0 << std::endl;
		std::cout << "SUM evaluations:     " << 0 << std::endl;
		std::cout << "Time elapsed:        " << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
		std::cout << "Preprocessing:       " << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
		std::cout << "Solution:            " << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
		std::cout << "-------------------------------" << std::endl;

		return SEARCH_SUCCESS;
	}

	Zobrist::initOnce(*(m_problem.get()));

	// init search space
	Timer tm;
	tm.start();
	m_space.reset(new RbfaooSearchSpace(m_pseudotree.get(), m_options));
	assert( m_options->cacheSize > 0 ); // in kilo bytes
	size_t cache_kilobytes = (m_options->cacheSize);
	if (!m_space->dfpncache) {
		m_space->dfpncache = new RbfaooCacheTable;
		m_space->dfpncache->init(cache_kilobytes);
	}
	tm.stop();
	std::cout << "Cache allocation complete: " << tm.elapsed() << " seconds" << std::endl;

	// prologue
	std::cout << "--- Starting search ---" << std::endl;
	int res = SEARCH_SUCCESS;
	char* EmergencyMemory = new char[16384]; // in order to handle bad alloc

	m_tmSolve = 0;
	m_timer.reset(m_options->timeLimit);
	m_timer.start();
	m_solutionCost = ELEM_NAN;
	m_assignment.resize(m_problem->getN(), UNKNOWN);

	try {
		// init the AO search space
//		ao_space.reset(new AobbSearchSpace());
//		ao_propagator.reset(new BoundPropagator(m_problem.get(), ao_space.get(), m_pseudotree.get(), m_options));
//		if (m_options->propagation == CP_UNIT) {
//			ao_propagator->setSatSolver(&m_zchaff);
//		} else {
//			ao_propagator->setSatSolver(NULL);
//		}
//		ao_space->cache = new AobbCacheTable(m_problem->getN());

		// do the search
		res = rbfs();

	} catch (std::bad_alloc& ba) {
		delete[] EmergencyMemory;
		EmergencyMemory = NULL;
		std::cerr << "out of memory";
		res = SEARCH_FAILURE;
	} catch (...) {
		std::cout << "unexpected error: ";
		res = SEARCH_FAILURE;
	}

	// stop timer
	m_timer.stop();
	m_tmSolve = m_timer.elapsed();

	// output solution (if found)
	std::cout << std::endl << std::endl;
	std::cout << "--- Search done ---" << std::endl;
	std::cout << "Problem name:        " << m_problem->getName() << std::endl;
	std::cout << "Status:              " << solver_status[res] << std::endl;
	std::cout << "OR nodes:            " << m_num_expanded_or << std::endl;
	std::cout << "AND nodes:           " << m_num_expanded_and << std::endl;
	std::cout << "OR nodes (MAP):      " << m_num_expanded_or_map << std::endl;
	std::cout << "AND nodes (MAP):     " << m_num_expanded_and_map << std::endl;
	std::cout << "SUM evaluations:     " << m_numSumEvals << std::endl;
	std::cout << "Time elapsed:        " << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
	std::cout << "Preprocessing:       " << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
//	std::cout << "Solution:            " << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
	std::cout << "Solution:            " << 1/ELEM_DECODE(m_solutionCost) << " (" << -m_solutionCost  << ")" << std::endl;
	std::cout << "-------------------------------" << std::endl;

	// clean up
	Zobrist::finishOnce();
	if (EmergencyMemory) delete[] EmergencyMemory;
	return res;
}

// Radu: DONE
// Recursive best first search; assumes that the cache was already initialized
int RBFAOO::rbfs() {

	// Root of the search space
	RbfaooSearchNode* root = NULL;
	m_solutionCost = ELEM_NAN;
	m_assignment.resize(m_problem->getN(), NONE);
	int res = SEARCH_SUCCESS;
	char* EmergencyMemory = new char[16384]; // in order to handle bad alloc

	try {

		// Add initial set of dummy nodes.
		PseudotreeNode* ptroot = m_pseudotree->getRoot();

		// create dummy AND node (domain size 1) with global constant as label
		int var = ptroot->getVar();
		double label = -ELEM_ENCODE( m_problem->getGlobalConstant() );
		root = new RbfaooSearchNode(NODE_AND, var, 0, -1, label);

		Zobrist z;
		const set<int> &full_context = ptroot->getFullContext();
		assert(root->getCacheContext().size() == 0);
		z.encodeOR(full_context, m_assignment);
		root->getZobrist().encodeAND(z, 0);

//		std::cout << "Root is " << *root << std::endl;

		multipleIterativeDeepeningAND(*root, 0, ELEM_INF, DNINF);

		// Search done; retrieve the optimal solution from the root node
		RbfaooCacheElementPtr p = m_space->dfpncache->read(root->getZobrist(), var,
			0, root->getCacheContext());

		assert(p);
		m_solutionCost = p->result;

	} catch (std::bad_alloc& ba) {
		delete[] EmergencyMemory;
		EmergencyMemory = NULL;
		std::cerr << "out of memory";
		res = SEARCH_FAILURE;
	} catch (int) {
		res = SEARCH_TIMEOUT;
	}

	// clean-up
	if (root) delete root;
	if (EmergencyMemory) delete[] EmergencyMemory;

	return res;
} // rbfs()

// Radu: DONE
bool RBFAOO::doProcess(RbfaooSearchNode* node) {
	assert(node);
	if (node->getType() == NODE_AND) {
		int var = node->getVar();
		int val = node->getVal();
		m_assignment[var] = val; // record assignment
	} else { // NODE_OR
		// nothing
	}

	return false; // default
} // doProcess

// Kishi: DONE
void RBFAOO::multipleIterativeDeepeningOR(RbfaooSearchNode &n,
		size_t table_index, double th_value, int dn_threshold) {
	bool children_generated_flag = false;

	int num_children;
	int var = n.getVar();
	context_t ctxt = n.getCacheContext();
	unsigned int increased_subtree_size = m_num_expanded;
	double value;
	int dn, best_index, solved_flag;
	for (;;) {
		double second_best_value, cutoff_value;

		if (!children_generated_flag) {
//			std::cout << "Expanding " << n << " with threshold " << th_value << std::endl;
			num_children = generateChildrenOR(n, second_best_value, best_index,
					solved_flag, cutoff_value);
			m_num_expanded++;
			m_num_expanded_or++;
			if (m_problem->isMap(var)) m_num_expanded_or_map++;
		} else {
//			std::cout << "Evaluating " << n << std::endl;
			calculateNodeValueandDNOR(n, table_index, num_children,
					second_best_value, best_index, solved_flag, cutoff_value);
		}

		//cerr << "OR: children_generated_flag = " << children_generated_flag << " d=" << n.getDepth() << " num_exp = " << m_num_expanded << " mapvarf = " << m_problem->isMap(n.getVar()) << " z = " << n.getZobrist().getC0() << " value = " << n.getNodeValue() << endl;
		children_generated_flag = true;
		value = n.getNodeValue();
		dn = n.getDN();

		if (value > th_value + epsilon || dn == DNINF || solved_flag) {
//			cout << "****OR threshold violation****" << endl;
//			cout << "  - value:       " << (value) << endl
//				 << "  - threshold:   " << th_value << endl
//				 << "  - solved flag: " << solved_flag << endl;
			break;
		}

		assert(best_index >= 0 && best_index < num_children);
		RbfaooSearchNode *c = m_expand[table_index + best_index];
		// Kishi: !! overestmation technique --- later use avarage like my AAAI 2005 paper !!
		//double new_th_value = max(th_value, second_best_value OP_DIVIDE exp(-n.getDepth()));
		double new_th_value;
		if (m_problem->isMap(c->getVar())) {
			new_th_value = min(th_value, second_best_value + m_overestimation);
			new_th_value = min(new_th_value, cutoff_value);
			new_th_value -= c->getLabel();
		} else {
			new_th_value = ELEM_INF;
		}
		int new_dn_threshold = min(static_cast<int>(DNINF),
						dn_threshold - dn + c->getDN());
		size_t new_table_index = table_index + num_children;
//		std::cout << "New (OR) threshold is " << new_th_value << std::endl;

		// check for time limit violation
		if (m_timer.timeout()) {
			throw SEARCH_TIMEOUT;
		}

		// output intermediate results if verbose mode is enabled
		if (m_options->verbose) {
			// logging
#ifdef SEARCH_CHECKPOINT_NODE
			size_t expansions = m_space->getAndNodes();
			if (expansions % SEARCH_CHECKPOINT_NODE == 0 &&
				expansions > m_prevNodeCheckpoint) {

				double timestamp = m_timer.elapsed();
				m_prevNodeCheckpoint = expansions;
				std::cout << "[*"
					<< std::setw(8) << timestamp << "] w "
					<< std::setw(8) << m_weight << " "
					<< std::setw(12) << m_space->getOrNodes() << " "
					<< std::setw(12) << m_space->getAndNodes() << " "
					<< std::setw(12) << m_upperBound
					<< std::endl;
			}
#else
			double elapsed = m_timer.elapsed();
			if (elapsed - m_prevTimeCheckpoint >= SEARCH_CHECKPOINT_TIME) {

				double timestamp = elapsed;
				m_prevTimeCheckpoint = elapsed;
				std::cout << "[*"
					<< std::setw(8) << timestamp << "] w "
					<< std::setw(8) << m_epsilon << " "
					<< std::setw(12) << m_num_expanded_or_map << " "
					<< std::setw(12) << m_num_expanded_and_map << " "
					<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ")"
					<< std::endl;
			}
#endif
		}

//		if (m_num_expanded % 1000000 == 0)
//			m_space->dfpncache->printWriteStats();

		doProcess(c);
		multipleIterativeDeepeningAND(*c, new_table_index, new_th_value,
				new_dn_threshold);
	}
	increased_subtree_size = m_num_expanded - increased_subtree_size;
	for (size_t i = table_index; i < table_index + num_children; i++)
		delete m_expand[i];

	m_expand.resize(table_index);
	m_space->dfpncache->write(n.getZobrist(), var, ctxt, 1, best_index, value,
			dn, solved_flag, increased_subtree_size);
	if (m_pseudotree->isSumRoot(n.getVar()))
		m_numSumEvals++;
}

// Kishi: DONE
void RBFAOO::multipleIterativeDeepeningAND(RbfaooSearchNode &n,
		size_t table_index, double th_value, int dn_threshold) {
	bool children_generated_flag = false;

	int num_children, best_index;
	int var = n.getVar();
	context_t ctxt = n.getCacheContext();
	unsigned int increased_subtree_size = m_num_expanded;
	double value;
	int dn, solved_flag;
	for (;;) {
		int second_best_dn;

		if (!children_generated_flag) {
//			std::cout << "Expanding " << n << " with threshold " << th_value << endl;
			num_children = generateChildrenAND(n, second_best_dn, best_index,
					solved_flag);
			m_num_expanded++;
			m_num_expanded_and++;
			if (m_problem->isMap(var)) m_num_expanded_and_map++;
		} else {
//			std::cout << "Evaluating " << n << endl;
			calculateNodeValueandDNAND(n, table_index, num_children,
					second_best_dn, best_index, solved_flag);
		}

		//cerr << "AND: children_generated_flag = " << children_generated_flag << " d=" << n.getDepth() << " num_exp = " << m_num_expanded << " z = " << n.getZobrist().getC0() << " value = " << n.getNodeValue() << endl;
		children_generated_flag = true;

		value = n.getNodeValue();
		dn = n.getDN();

		if (value > th_value + epsilon || dn == DNINF || solved_flag) {
//			cout << "****AND threshold violation****" << endl;
//			cout << "  - value:       " << (value) << endl
//				 << "  - threshold:   " << th_value << endl
//				 << "  - solved flag: " << solved_flag << endl;
			break;
		}

		assert(best_index >= 0 && best_index < num_children);
		RbfaooSearchNode *c = m_expand[table_index + best_index];
		double c_value = c->getNodeValue();
		double new_th_value;
		if (m_problem->isMap(c->getVar())) {
			new_th_value = (th_value - value) + c_value;
		} else {
			new_th_value = ELEM_INF;
		}
		int new_dn_threshold = min(sumDN(second_best_dn, 1), dn_threshold);
		//; assert(new_dn_threshold > 0);
		size_t new_table_index = table_index + num_children;

		// check for time limit violation
		if (m_timer.timeout()) {
			throw SEARCH_TIMEOUT;
		}

		// output intermediate results if verbose mode is enabled
		if (m_options->verbose) {
			// logging
#ifdef SEARCH_CHECKPOINT_NODE
			size_t expansions = m_space->getAndNodes();
			if (expansions % SEARCH_CHECKPOINT_NODE == 0 &&
				expansions > m_prevNodeCheckpoint) {

				double timestamp = m_timer.elapsed();
				m_prevNodeCheckpoint = expansions;
				std::cout << "[*"
					<< std::setw(8) << timestamp << "] w "
					<< std::setw(8) << m_weight << " "
					<< std::setw(12) << m_space->getOrNodes() << " "
					<< std::setw(12) << m_space->getAndNodes() << " "
					<< std::setw(12) << m_upperBound
					<< std::endl;
			}
#else
			double elapsed = m_timer.elapsed();
			if (elapsed - m_prevTimeCheckpoint >= SEARCH_CHECKPOINT_TIME) {

				double timestamp = elapsed;
				m_prevTimeCheckpoint = elapsed;
				std::cout << "[*"
					<< std::setw(8) << timestamp << "] w "
					<< std::setw(8) << m_epsilon << " "
					<< std::setw(12) << m_num_expanded_or_map << " "
					<< std::setw(12) << m_num_expanded_and_map << " "
					<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ")"
					<< std::endl;
			}
#endif
		}

//		if (m_num_expanded % 1000000 == 0)
//			m_space->dfpncache->printWriteStats();
		multipleIterativeDeepeningOR(*c, new_table_index, new_th_value,
				new_dn_threshold);
	}
	increased_subtree_size = m_num_expanded - increased_subtree_size;
	for (size_t i = table_index; i < table_index + num_children; i++)
		delete m_expand[i];
	m_expand.resize(table_index);
	m_space->dfpncache->write(n.getZobrist(), var, ctxt, 0, best_index, value,
			dn, solved_flag, increased_subtree_size);
}

// Kishi: DONE
int RBFAOO::generateChildrenAND(RbfaooSearchNode &n, int &second_best_dn,
		int &best_index, int &solved_flag) {

	assert(n.getType() == NODE_AND);

	//m_space->nodesAND += 1;
	m_space->incNodesExpanded(NODE_AND);

	int var = n.getVar();
	PseudotreeNode* ptnode = m_pseudotree->getNode(var);
	int depth = ptnode->getDepth();

	double value = ELEM_ZERO; //ELEM_ONE;
	int dn = DNINF;

	second_best_dn = DNINF;
	best_index = UNDETERMINED;
	solved_flag = 0;

	int num_children = 0;
	// create new OR children (going in reverse due to reversal on stack),
	for (vector<PseudotreeNode*>::const_iterator it =
			ptnode->getChildren().begin(); it != ptnode->getChildren().end();
			++it, num_children++) {
		int vChild = (*it)->getVar();
		PseudotreeNode* ptnode = m_pseudotree->getNode(vChild);
		const set<int> &full_context = ptnode->getFullContext();
		RbfaooSearchNode* c = new RbfaooSearchNode(NODE_OR, vChild, depth + 1);
		for (set<int>::const_iterator it = full_context.begin();
				it != full_context.end(); ++it) {
			assert(*it != vChild);
			c->addCacheContext(m_assignment[*it]);
		}
		Zobrist &z = c->getZobrist();
		z.encodeOR(full_context, m_assignment);
		RbfaooCacheElementPtr p = m_space->dfpncache->read(z, vChild, 1,
				c->getCacheContext());
		double child_value;
		int child_dn;
		if (p == NULL) {
			// Compute and set heuristic estimate, includes child labels
			// Kishi: perhaps we need to use a different heuristic
			child_value = heuristicOR(*c);
			if (m_problem->isSum(vChild)) child_value = ELEM_ZERO;
			child_dn = 1;
		} else {
			//!!KISHI: ADHOC FIX! NEED TO THINK ABOUT A WAY TO REMOVE THIS !!
			heuristicOR(*c);
			// edge cost is zero from AND node to OR child
			child_value = p->result;
			child_dn = p->dn;
		}
		c->setNodeValue(child_value);
		c->setDN(child_dn);
		value += child_value;
		if (child_dn < dn) {
			second_best_dn = dn;
			dn = child_dn;
			best_index = num_children;
		} else if (child_dn < second_best_dn)
			second_best_dn = child_dn;
		m_expand.push_back(c);
//		cout << "  - child " << *c << endl;
	} // for loop over new OR children

	n.setNodeValue(value);
	n.setDN(dn);

	if (dn == DNINF)
		solved_flag = 1;

#if 000
	if (num_children == 0)
	cout << "terminal AND: m_num_expanded = " << m_num_expanded << "\n";
#endif

	return num_children;
} // RBFAOO::generateChildrenAND

// Kishi: DONE
void RBFAOO::calculateNodeValueandDNAND(RbfaooSearchNode &n, size_t table_index,
		int num_children, int &second_best_dn, int &best_index,
		int &solved_flag) {

	assert(n.getType() == NODE_AND);

	double value = ELEM_ZERO; //ELEM_ONE;
	int dn = DNINF;
	second_best_dn = DNINF;
	best_index = UNDETERMINED;

	solved_flag = 0;
	// create new OR children (going in reverse due to reversal on stack),
	for (size_t i = 0; i < static_cast<size_t>(num_children); i++) {
		RbfaooSearchNode *c = m_expand[table_index + i];
		Zobrist &z = c->getZobrist();
		int vChild = c->getVar();
		RbfaooCacheElementPtr p = m_space->dfpncache->read(z, vChild, 1,
				c->getCacheContext());
		double child_value;
		int child_dn;
		if (p == NULL) {
			// Compute and set heuristic estimate, includes child labels
			// Kishi: perhaps we need to use a different heuristic
			child_value = c->getHeur();
			child_dn = 1;
		} else {
			// edge cost is zero from AND node to OR child
			child_value = p->result;
			child_dn = p->dn;
		}
		c->setNodeValue(child_value);
		c->setDN(child_dn);
		value += child_value;
		if (child_dn < dn) {
			second_best_dn = dn;
			dn = child_dn;
			best_index = i;
		} else if (child_dn < second_best_dn)
			second_best_dn = child_dn;
	} // for loop over new OR children

	n.setNodeValue(value);
	n.setDN(dn);

	if (dn == DNINF)
		solved_flag = 1;
} // RBFAOO::calculateNodeValueandDNAND

// Kishi: DONE
int RBFAOO::generateChildrenOR(RbfaooSearchNode &n, double &second_best_value,
		int &best_index, int &solved_flag, double &cutoff_value) {

	assert(n.getType() == NODE_OR);

	//m_space->nodesOR +=1;
	m_space->incNodesExpanded(NODE_OR);

	int var = n.getVar();
	bool mapVar = m_problem->isMap(var);
	PseudotreeNode* ptnode = m_pseudotree->getNode(var);
	int depth = ptnode->getDepth();
	solved_flag = 0;
	if (!mapVar) solved_flag = 1;

	// retrieve precomputed labels and heuristic values
	double* heur = n.getHeurCache();
	;
	assert(heur);
	double value = (mapVar) ? ELEM_INF : ELEM_ZERO;
	cutoff_value = second_best_value = ELEM_INF;
	int dn = 0, dn_with_ignoring_proven_child = 0;
	int num_children = 0;

	for (val_t i = 0; i < m_problem->getDomainSize(var); i++, num_children++) {
		RbfaooSearchNode *c = new RbfaooSearchNode(NODE_AND, var, i, depth,
				heur[2 * i + 1]); // uses cached label
		Zobrist &z = c->getZobrist();
		z.encodeAND(n.getZobrist(), i);
		c->setCacheContext(n.getCacheContext());
		context_t &context = c->getCacheContext();
		RbfaooCacheElementPtr p = m_space->dfpncache->read(z, var, 0, context);
		double child_value;
		int child_dn, child_solved_flag;
		c->setHeur(heur[2 * i]);
		if (p == NULL) {
			child_value = heur[2 * i];
			child_dn = 1;
			child_solved_flag = 0;
		} else {
			child_value = p->result + c->getLabel();
			child_dn = p->dn;
			child_solved_flag = p->solved_flag;
		}
		if (heur[2 * i] == ELEM_INF) {
			; assert(child_value == ELEM_INF);
			child_dn = 0;
		}
		if (mapVar) {
			if (child_value < value) {
				best_index = i;
				second_best_value = value;
				value = child_value;
				solved_flag = child_solved_flag;
			} else if (child_value == value && child_solved_flag) {
				best_index = i;
				second_best_value = value;
				solved_flag = child_solved_flag;
			} else if (child_value < second_best_value)
				second_best_value = child_value;
			if (child_dn != DNINF) {
				dn_with_ignoring_proven_child = sumDN(dn_with_ignoring_proven_child,
						child_dn);
			} else {
				;
				assert(child_solved_flag);
				cutoff_value = min(cutoff_value, child_value); // Radu: was 'max'
			}
			dn = sumDN(dn, child_dn);
		} else {
			value += 1/ELEM_DECODE(child_value);
			best_index = i;
			second_best_value = value;
			solved_flag &= child_solved_flag;
			if (child_dn != DNINF) {
				dn_with_ignoring_proven_child = sumDN(dn_with_ignoring_proven_child,
						child_dn);
			} else {
				;
				assert(child_solved_flag);
				cutoff_value = ELEM_INF;
			}
			dn = sumDN(dn, child_dn);
		}

		c->setNodeValue(child_value);
		c->setDN(child_dn);
		m_expand.push_back(c);
//		cout << "  - child " << *c << endl;
	}

	if (!solved_flag && dn == DNINF) {
		//need to recalculate dn
		;
		assert(dn_with_ignoring_proven_child != 0);
		dn = max(1, dn_with_ignoring_proven_child);
	}
	if (!mapVar) value = -ELEM_ENCODE(value);

	n.setNodeValue(value);
	n.setDN(dn);

#if 000
	if (num_children == 0)
	cout << "terminal OR: m_num_expanded = " << m_num_expanded << "\n";
#endif

	return num_children;
} // RBFAOO::generateChildrenOR

// Kishi: DONE
void RBFAOO::calculateNodeValueandDNOR(RbfaooSearchNode &n, size_t table_index,
		int num_children, double &second_best_value, int &best_index,
		int &solved_flag, double &cutoff_value) {

	assert(n.getType() == NODE_OR);

	solved_flag = 0;

	bool mapVar = m_problem->isMap(n.getVar());
	if (!mapVar) solved_flag = 1;

	double value = (mapVar) ? ELEM_INF : ELEM_ZERO;
	double* heur = n.getHeurCache();
	cutoff_value = second_best_value = ELEM_INF;
	int dn = 0, dn_with_ignoring_proven_child = 0;

	for (int i = 0; i < num_children; i++) {
		RbfaooSearchNode *c = m_expand[table_index + i];
		Zobrist &z = c->getZobrist();
		context_t &context = c->getCacheContext();
		int vChild = c->getVar();
		RbfaooCacheElementPtr p = m_space->dfpncache->read(z, vChild, 0, context);
		double child_value;
		int child_dn, child_solved_flag;
		if (p == NULL) {
			if (heur[2 * i] == ELEM_INF) {
				child_value = ELEM_INF;
				child_dn = 0;
			} else {
				child_value = c->getHeur();
				child_dn = 1;
			}
			child_solved_flag = 0;
		} else {
			child_value = p->result + c->getLabel();
			child_dn = p->dn;
			child_solved_flag = p->solved_flag;
		}

		if (mapVar) {
			if (child_value < value) {
				best_index = i;
				second_best_value = value;
				value = child_value;
				solved_flag = child_solved_flag;
			} else if (child_value == value) {
				if (child_solved_flag) {
					best_index = i;
					second_best_value = value;
					solved_flag = child_solved_flag;
				}
			} else if (child_value < second_best_value)
				second_best_value = child_value;
			if (child_dn != DNINF) {
				dn_with_ignoring_proven_child = sumDN(dn_with_ignoring_proven_child,
						child_dn);
			} else {
				;
				assert(child_solved_flag);
				cutoff_value = min(cutoff_value, child_value); // Radu: was 'max'
			}
			dn = sumDN(dn, child_dn);
		} else {
			value += 1/ELEM_DECODE(child_value);
			second_best_value = value;
			solved_flag &= child_solved_flag;
			if (!child_solved_flag) best_index = i;

			if (child_dn != DNINF) {
				dn_with_ignoring_proven_child = sumDN(dn_with_ignoring_proven_child,
						child_dn);
			} else {
				;
				assert(child_solved_flag);
				cutoff_value = ELEM_INF;
			}
			dn = sumDN(dn, child_dn);
		}

		c->setNodeValue(child_value);
		c->setDN(child_dn);
	}

	if (!solved_flag && dn == DNINF) {
		//need to recalculate dn
		;
		assert(dn_with_ignoring_proven_child != 0);
		dn = max(1, dn_with_ignoring_proven_child);
	}

	if (!mapVar) value = -ELEM_ENCODE(value);
	n.setNodeValue(value);
	n.setDN(dn);
} // RBFAOO::generateChildrenOR

// Kishi: DONE (we have a Marginal MAP problem : MAX-PROD)
double RBFAOO::heuristicOR(RbfaooSearchNode &n) {

	int var = n.getVar();
	bool mapVar = m_problem->isMap(var);
	int oldValue = m_assignment[var];
	double* dv = new double[m_problem->getDomainSize(var) * 2];
	double h = (mapVar) ? INFINITY : ELEM_ZERO; // the new OR nodes h value
	const list<Function*>& funs = m_pseudotree->getFunctions(var);
	for (val_t i = 0; i < m_problem->getDomainSize(var); ++i) {
		m_assignment[var] = i;

		// compute (weighted) heuristic value
		dv[2 * i] = m_epsilon * -ELEM_ENCODE(m_heuristic->getHeur(var, m_assignment));

		// precompute label value
		double d = ELEM_ONE;
		for (std::list<Function*>::const_iterator it = funs.begin();
				it != funs.end(); ++it) {
			d *= (*it)->getValue(m_assignment);
		}

		// store label and heuristic into cache table
		dv[2 * i + 1] = -ELEM_ENCODE(d); // label
		dv[2 * i] += -ELEM_ENCODE(d); // heuristic

		if (mapVar) {
			if (dv[2 * i] < h) {
				h = dv[2 * i]; // keep min. for OR node heuristic (MAP var)
			}
		} else {
			//h += dv[2 * i]; // keep sum. for OR node heuristic (SUM var)
			dv[2 * i] = ELEM_ZERO;
		}
	}

	n.setHeur(h);
	n.setHeurCache(dv);
	m_assignment[var] = oldValue;

	return h;

} // RBFAOO::heuristicOR


//////////////////////////////////////////////////////////////////////////////////
//// solve the problem
//double RBFAOO::_ao(int var, const std::vector<int>& assignment) {
//
//	// a temporary buffer
//	double result = ELEM_NAN;
//
//	// reinit the search space
//	ao_space->root = NULL;
//	m_assignment = assignment;
//
//	// init root of the search space
//	_initAoSearch(var);
//
//	AobbSearchNode* n = _nextLeaf();
//	while (n) { // throws timeout
//		ao_propagator->propagate(n, m_timer, m_currentDomains, false); // true = report solutions
//		n = _nextLeaf();
//	}
//
//	// Conditioned SUM problem solved
//	result = ao_space->root->getValue();
//	delete ao_space->root;
//	ao_space->root = NULL;
//
//	return result;
//}
//
//// DONE: radu
//void RBFAOO::_initAoSearch(int var) {
//
//	// create root OR node
//	AobbSearchNode* node = new AobbSearchNodeOR(NULL, var, 0);
//	ao_space->root = node;
//	_heuristicOR(node);
//
//	assert(ao_stack.empty());
//	ao_stack.push(node);
//}
//
//// DONE: radu
//void RBFAOO::_heuristicOR(AobbSearchNode* n) {
//
//	// Safety checks
//	assert(n->getType() == NODE_OR);
//
//	int v = n->getVar();
//	double d;
//	double* dv = new double[m_problem->getDomainSize(v) * 2];
//	const list<Function*>& funs = m_pseudotree->getFunctions(v);
//
//	for (val_t i = 0; i < m_problem->getDomainSize(v); ++i) {
//		m_assignment[v] = i;
//
//		// compute heuristic value
//		dv[2 * i] = ELEM_NAN;
//
//		// precompute label value
//		d = ELEM_ONE;
//		for (std::list<Function*>::const_iterator it = funs.begin();
//				it != funs.end(); ++it) {
//			d *= (*it)->getValue(m_assignment);
//		}
//
//		// store label and heuristic into cache table
//		dv[2 * i + 1] = d; // label
//		dv[2 * i] *= d; // heuristic (includes label)
//
//	}
//
//	n->setHeurCache(dv);
//
//} // AOBB::heuristicOR
//
//// DONE: radu
//bool RBFAOO::_doProcess(AobbSearchNode* node) {
//	assert(node);
//	if (node->getType() == NODE_AND) {
//		int var = node->getVar();
//		int val = node->getVal();
//		m_assignment[var] = val; // record assignment
//
//	} else { // NODE_OR
//		// do nothing
//	}
//
//	return false; // default
//}
//
//// DONE: radu
//void RBFAOO::_addCacheContext(AobbSearchNode* node, const set<int>& ctxt) const {
//
//	context_t sig;
//	sig.reserve(ctxt.size());
//	for (set<int>::const_iterator itC = ctxt.begin(); itC != ctxt.end(); ++itC) {
//		sig.push_back(m_assignment[*itC]);
//	}
//
//	node->setCacheContext(sig);
//}
//
//// DONE: radu
//bool RBFAOO::_doCaching(AobbSearchNode* node) {
//
//	assert(node);
//	int var = node->getVar();
//	PseudotreeNode* ptnode = m_pseudotree->getNode(var);
//
//	if (node->getType() == NODE_AND) { // AND node -> reset associated adaptive cache tables
//
//		// do nothing at AND nodes
//	    const list<int>& resetList = ptnode->getCacheReset();
//	    for (list<int>::const_iterator it=resetList.begin(); it!=resetList.end(); ++it)
//	    	ao_space->cache->reset(*it);
//
//	} else { // OR node, try actual caching
//
//		if (!ptnode->getParent())
//			return false;
//
//		if (ptnode->getFullContext().size() <= ptnode->getParent()->getFullContext().size()) {
//			// add cache context information
//			//addCacheContext(node, ptnode->getFullContext());
//			_addCacheContext(node, ptnode->getCacheContext());
//
//			// try to get value from cache
//			try {
//				// will throw int(UNKNOWN) if not found
//				double entry = ao_space->cache->read(var, node->getCacheContext());
//				node->setValue(entry); // set value
//				node->setLeaf(); // mark as leaf
//
//				return true;
//			} catch (...) { // cache lookup failed
//				node->setCachable(); // mark for caching later
//			}
//		}
//	} // if on node type
//
//	return false; // default, no caching applied
//
//} // BBBT::doCaching
//
//// DONE: radu
//bool RBFAOO::_doPruning(AobbSearchNode* node) {
//
//	assert(node);
//
//	if (_canBePruned(node)) {
//		node->setLeaf();
//		node->setPruned();
//		node->setValue(ELEM_ZERO); // dead end
//
//		return true;
//	}
//
//	return false; // default false
//} // AOBB::doPruning
//
//bool RBFAOO::_canBePruned(AobbSearchNode* n) {
//
//	// propagate current assignment (AND nodes only)
//	if (n->getType() == NODE_AND) {
//
////		if (m_options->propagation == CP_FULL) {
////
////			// propagate the assignment (as unit propagation)
////			std::vector<int> vars = m_assumptions;
////			getPathAssignment(n, vars);
////
////			return propagate(vars); // 'true' if conflict
////
////		} else if (m_options->propagation == CP_UNIT) {
////
////			// Propagate the assignment
////			m_zchaff.dlevel()++;
////			if (!lookahead(n->getVar(), n->getVal(), n->changes(),
////					m_pseudotree->getNode(n->getVar())->getSubprobVars())) {
////				// found inconsistency
////				return true; // 'true' if conflict
////			}
////		}
//	}
//
//	// default, no pruning possible
//	return false;
//
//} // BBBT::canBePruned
//
//// DONE: radu
//AobbSearchNode* RBFAOO::_nextLeaf() {
//
//	AobbSearchNode* node = _nextNode();
//	while (node) {
//
//		// check for time limit violation
//		if (m_timer.timeout()) {
//			throw SEARCH_TIMEOUT;
//		}
//
//		if (_doProcess(node)) { // initial processing
//			return node;
//		}
//		if (_doCaching(node)) { // caching?
//			return node;
//		}
//		if (_doPruning(node)) { // pruning?
//			return node;
//		}
//		if (_doExpand(node)) { // node expansion
//			return node;
//		}
//		node = _nextNode();
//	}
//	return NULL;
//}
//
//// DONE: radu
//bool RBFAOO::_generateChildrenAND(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {
//
//	assert(n && n->getType() == NODE_AND);
//	//assert(n->getChildren() == NULL);
//
//	if (n->getChildren()) {  // node previously expanded
//		if (n->getChildCountAct() == 0) {
//			n->clearChildren(); // previously expanded, but all children deleted
//		} else {
//			for (size_t i = 0; i < n->getChildCountFull(); ++i) {
//				if (n->getChildren()[i])
//					chi.push_back(n->getChildren()[i]);
//			}
//			return false;
//		}
//	}
//
//	int var = n->getVar();
//	PseudotreeNode* ptnode = m_pseudotree->getNode(var);
//	int depth = ptnode->getDepth();
//
//	// increase AND node expansions
//	ao_space->nodesAND += 1;
//
//	// create new OR children (going in reverse due to reversal on stack)
//	for (vector<PseudotreeNode*>::const_iterator it =
//			ptnode->getChildren().begin(); it != ptnode->getChildren().end();
//			++it) {
//		int vChild = (*it)->getVar();
//		AobbSearchNodeOR* c = new AobbSearchNodeOR(n, vChild, depth + 1);
//		// Compute and set heuristic estimate, includes child labels
//		_heuristicOR(c);
//		chi.push_back(c);
//
//	} // for loop over new OR children
//
//	if (chi.empty()) {
//		n->setLeaf(); // terminal node
//		n->setValue(ELEM_ONE);
//		return true; // no children
//	}
//
//	// (use reverse iterator due to stack reversal)
//	n->addChildren(chi);
//
//	return false; // default
//
//} // AO::generateChildrenAND
//
//// DONE: radu
//bool RBFAOO::_generateChildrenOR(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {
//
//	assert(n && n->getType() == NODE_OR);
//	//assert(n->getChildren() == NULL);
//
//	if (n->getChildren()) {  // node previously expanded
//		if (n->getChildCountAct() == 0) {
//			n->clearChildren(); // previously expanded, but all children deleted
//		} else {
//			for (size_t i = 0; i < n->getChildCountFull(); ++i) {
//				if (n->getChildren()[i])
//					chi.push_back(n->getChildren()[i]);
//			}
//			return false;
//		}
//	}
//
//
//	int var = n->getVar();
//
//	// increase OR node expansions
//	ao_space->nodesOR += 1;
//
//	// retrieve precomputed labels
//	double* heur = n->getHeurCache();
//	for (val_t i = 0 ; i < m_problem->getDomainSize(var); ++i) {
//
//		// consistent value
//		if (m_currentDomains.empty() == false &&
//			m_currentDomains[var][i] == false) {
//			continue;
//		}
//
//		double label = heur[2 * i + 1];
//		AobbSearchNodeAND* c = new AobbSearchNodeAND(n, i, label); // uses cached label
//		// set cached heur. value (includes the weight)
//		c->setHeur(heur[2 * i]);
//		chi.push_back(c);
//	}
//
//	if (chi.empty()) { // deadend
//		n->setLeaf();
//		n->setValue(ELEM_ZERO);
//		return true; // no children
//	}
//
//	// (use reverse iterator due to stack reversal)
//	n->addChildren(chi);
//
//	return false; // default
//
//} // AO::generateChildrenOR
//
//// DONE: radu
//AobbSearchNode* RBFAOO::_nextNode() {
//	AobbSearchNode* n = NULL;
//	if (!n && ao_stack.size()) {
//		n = ao_stack.top();
//		ao_stack.pop();
//	}
//	return n;
//}
//
//// DONE: radu
//bool RBFAOO::_doExpand(AobbSearchNode* n) {
//	assert(n);
//	ao_expand.clear();
//
//	if (n->getType() == NODE_AND) {  // AND node
//
//		if (_generateChildrenAND(n, ao_expand))
//			return true; // no children
//		for (vector<AobbSearchNode*>::reverse_iterator it = ao_expand.rbegin();
//				it != ao_expand.rend(); ++it)
//			ao_stack.push(*it);
//
//	} else {  // OR node
//
//		if (_generateChildrenOR(n, ao_expand))
//			return true; // no children
//		for (vector<AobbSearchNode*>::reverse_iterator it = ao_expand.rbegin();
//				it != ao_expand.rend(); ++it) {
//			ao_stack.push(*it);
//		} // for loop
//
//	} // if over node type
//
//	return false; // default false
//}


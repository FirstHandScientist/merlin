/*
 * any_wrbfs2.cpp
 *
 *  Created on: 18 Nov 2015
 *      Author: radu
 */


#include "any_wrbfs2.h"
#include "mini_bucket_heur.h"
#include "mini_bucket_mm_heur.h"
#include "mini_bucket_jglp_heur.h"

const float epsilon = 1e-10;

// Radu: DONE
// initialize
void AnyWRBFS2::init() {

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

	// Now, we need to rearrange the pseudo tree such that the MAP variables
	// form a chain at the top of the pseudo tree (and they're searched in the
	// order suggested by the bucket tree heuristic ...
	std::list<int> order;
	m_pseudotree->dfs(order, true);
	m_pseudotree->forceMapChain(order, m_options->cbound);
	m_lastMapVar = order.back();
	assert( m_problem->isMap(m_lastMapVar) );

	// heuristic might have changed problem functions, pseudotree needs remapping
	m_pseudotree->resetFunctionInfo(m_problem->getFunctions());

	// check if heuristic is accurate
	if (m_heuristic->isAccurate()) {
		std::cout << "Heuristic is exact!" << std::endl;
		m_solutionCost = (m_heuristic->getGlobalUB());
	    m_solved = true;
	}

	// Record time after initialization
	cout << "Last MAP variable: " << m_lastMapVar << std::endl;
	cout << "Initialization complete: " << (m_tmLoad + m_tmHeuristic)
			<< " seconds" << std::endl;

	// output the pseudo tree
	m_pseudotree->outputToDot("pseudotree.dot");

	// check for upper bound computation only
	if (m_options->upperBoundOnly) {
		std::cout << "Done." << std::endl;
		exit(EXIT_SUCCESS);
	}
}

// Radu: DONE
// Kishi: define memory pool using boost to reduce overhead of new/delete
// Kishi: change implementation of TT to support constant-size memory allocation
AnyWRBFS2::AnyWRBFS2(ProgramOptions* opt) : Search(opt) {

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
	m_lastMapVar = UNKNOWN;
}

// Radu: DONE
int AnyWRBFS2::solve() {

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
int AnyWRBFS2::rbfs() {

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

		multipleIterativeDeepeningAND(*root, 0, ELEM_INF, DNINF, label);

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
bool AnyWRBFS2::doProcess(RbfaooSearchNode* node) {
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
void AnyWRBFS2::multipleIterativeDeepeningOR(RbfaooSearchNode &n,
		size_t table_index, double th_value, int dn_threshold, double g_value) {
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
		double new_th_value, new_g_value;
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
		new_g_value = g_value + c->getLabel();
//		std::cout << "New (OR) threshold is " << new_th_value << std::endl;

		// check for time limit violation
		if (m_timer.timeout()) {
			throw SEARCH_TIMEOUT;
		}

		// output intermediate results if verbose mode is enabled
		if (m_options->verbose) {
			// logging
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
		}

//		if (m_num_expanded % 1000000 == 0)
//			m_space->dfpncache->printWriteStats();

		doProcess(c);
		multipleIterativeDeepeningAND(*c, new_table_index, new_th_value,
				new_dn_threshold, new_g_value);
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
void AnyWRBFS2::multipleIterativeDeepeningAND(RbfaooSearchNode &n,
		size_t table_index, double th_value, int dn_threshold, double g_value) {
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
		}

		multipleIterativeDeepeningOR(*c, new_table_index, new_th_value,
				new_dn_threshold, g_value);
	}
	increased_subtree_size = m_num_expanded - increased_subtree_size;
	for (size_t i = table_index; i < table_index + num_children; i++)
		delete m_expand[i];
	m_expand.resize(table_index);
	m_space->dfpncache->write(n.getZobrist(), var, ctxt, 0, best_index, value,
			dn, solved_flag, increased_subtree_size);
}

// Kishi: DONE
int AnyWRBFS2::generateChildrenAND(RbfaooSearchNode &n, int &second_best_dn,
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

	if (dn == DNINF) solved_flag = 1;

	return num_children;
} // RBFAOO::generateChildrenAND

// Kishi: DONE
void AnyWRBFS2::calculateNodeValueandDNAND(RbfaooSearchNode &n, size_t table_index,
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
int AnyWRBFS2::generateChildrenOR(RbfaooSearchNode &n, double &second_best_value,
		int &best_index, int &solved_flag, double &cutoff_value) {

	assert(n.getType() == NODE_OR);

	//m_space->nodesOR +=1;
	m_space->incNodesExpanded(NODE_OR);

	int var = n.getVar();
	bool mapVar = m_problem->isMap(var);
	PseudotreeNode* ptnode = m_pseudotree->getNode(var);
	int depth = ptnode->getDepth();
//	solved_flag = 0; // MAP var (exactly one child must be solved)
//	if (!mapVar) solved_flag = 1; // SUM var (all children must be solved)

	solved_flag = 1; // MAP or SUM var (all children must be solved)
					 // need to ensure that search continues

	// retrieve precomputed labels and heuristic values
	// (the heuristic value is already inflated by the weight)
	double* heur = n.getHeurCache();
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
			value += child_value;
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
	}

	if (!solved_flag && dn == DNINF) {
		//need to recalculate dn
		;
		assert(dn_with_ignoring_proven_child != 0);
		dn = max(1, dn_with_ignoring_proven_child);
	}

	n.setNodeValue(value);
	n.setDN(dn);

	return num_children;
} // RBFAOO::generateChildrenOR

// Kishi: DONE
void AnyWRBFS2::calculateNodeValueandDNOR(RbfaooSearchNode &n, size_t table_index,
		int num_children, double &second_best_value, int &best_index,
		int &solved_flag, double &cutoff_value) {

	// Safety checks
	assert(n.getType() == NODE_OR);

	bool mapVar = m_problem->isMap(n.getVar());
//	solved_flag = 0; // MAP var (only one AND child (best) must be solved)
//	if (!mapVar) solved_flag = 1; // SUM var (all AND children must be solved)

	solved_flag = 1; // all AND children must be solved (or pruned)
					 // need this to be able to continue search

	double value = (mapVar) ? ELEM_INF : ELEM_ZERO;
	double* heur = n.getHeurCache();
	cutoff_value = second_best_value = ELEM_INF;
	int dn = 0, dn_with_ignoring_proven_child = 0;

	std::vector<RbfaooSearchNode*> succ_nodes; // only for MAP variables
	std::vector<bool> succ_solved;
	std::vector<int> succ_dn;

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

		if (mapVar) { // this is a MAP variable
			succ_nodes.push_back(c);
			succ_solved.push_back(child_solved_flag);
			succ_dn.push_back(child_dn);

//			if (child_value < value) {
//				best_index = i;
//				second_best_value = value;
//				value = child_value;
//				solved_flag = child_solved_flag;
//			} else if (child_value == value) {
//				if (child_solved_flag) {
//					best_index = i;
//					second_best_value = value;
//					solved_flag = child_solved_flag;
//				}
//			} else if (child_value < second_best_value)
//				second_best_value = child_value;
//			if (child_dn != DNINF) {
//				dn_with_ignoring_proven_child = sumDN(dn_with_ignoring_proven_child,
//						child_dn);
//			} else {
//				;
//				assert(child_solved_flag);
//				cutoff_value = min(cutoff_value, child_value); // Radu: was 'max'
//			}
//			dn = sumDN(dn, child_dn);
		} else { // this is a SUM variable
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
			c->setNodeValue(child_value);
			c->setDN(child_dn);
		}
	}

	if (mapVar) { // MAP vars: find best value/index and second best
		for (int i = 0; i < num_children; ++i) {
			RbfaooSearchNode* c = succ_nodes[i];
			double child_value = c->getNodeValue();
			bool child_solved_flag = succ_solved[i];
			int child_dn = succ_dn[i];
			if (child_solved_flag) {
				cutoff_value = std::min(cutoff_value, child_value);
			} else {
				if (child_value < value) {
					value = child_value;
					best_index = i;
				}
			}

		}

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
double AnyWRBFS2::heuristicOR(RbfaooSearchNode &n) {

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




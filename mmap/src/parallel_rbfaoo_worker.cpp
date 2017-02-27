/*
 * parallel_rbfaoo_worker.cpp
 * Programmed by Akihiro Kishimoto
 */

#include "parallel_rbfaoo_worker.h"

//Kishi: define memory pool using boost to reduce overhead of new/delete
//Kishi: change implementation of TT to support constant-size memory allocation
ParallelRBFAOOWorker::ParallelRBFAOOWorker(int thread_id, Problem *problem, Pseudotree *pseudotree, ProgramOptions* opt, 
					   ParallelRbfaooCacheTable *cache, Heuristic *heuristic) : 
  m_thread_id(thread_id), m_depth(0), m_problem(problem), m_pseudotree(pseudotree), m_options(opt), m_cache(cache), 
	m_heuristic(heuristic) {

	// Preallocate space for expansion vector.
	m_expand.reserve(65536);
	m_num_expanded_or = 0;
	m_num_expanded_and = 0;
	m_num_expanded = 0;
	m_num_expanded_or_map = 0;
	m_num_expanded_and_map = 0;
	m_prevTimeCheckpoint = 0;
	m_prevNodeCheckpoint = 0;
	m_overestimation = 1.0;
	m_weight = 1.0; // should be 1.0 by default
	m_thread_id = thread_id;
	m_num_threads = cache->getNumThreads();
}

//Kishi: DONE
void ParallelRBFAOOWorker::start() {

	// init search space
	// prologue
	std::cout << "--- Starting search ---" << std::endl;
	m_timer.reset(m_options->timeLimit * m_options->threads);
	m_timer.start();

	// do the search
	rbfs();

	// stop timer
	m_timer.stop();

	// output solution (if found)
	std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(6) << std::endl;
	std::cout << "Thread " << m_thread_id << ": --- Search done ---" << std::endl;
	std::cout << "Thread " << m_thread_id << ": Problem name:\t" << m_problem->getName() << std::endl;
	std::cout << "Thread " << m_thread_id << ": OR nodes:\t" << m_num_expanded_or <<  std::endl;
	std::cout << "Thread " << m_thread_id << ": AND nodes:\t" << m_num_expanded_and <<  std::endl;
	std::cout << "Thread " << m_thread_id << ": -------------------------------" << std::endl;
}

// Recursive best first search; assumes that the cache was already initialized
int ParallelRBFAOOWorker::rbfs() {

	// Root of the search space
	ParallelRbfaooSearchNode* root = NULL;
	m_solutionCost = ELEM_NAN;
	m_assignment.resize(m_problem->getN(), NONE);
	int res = SEARCH_SUCCESS;
	char* EmergencyMemory = new char[16384]; // in order to handle bad alloc

	try {

		// Add initial set of dummy nodes.
		PseudotreeNode* ptroot = m_pseudotree->getRoot();

		// create dummy AND node (domain size 1) with global constant as label
		int var = ptroot->getVar();
		root = new ParallelRbfaooSearchNode(NODE_AND, var, 0, -1, -ELEM_ENCODE(m_problem->getGlobalConstant()));

		Zobrist z;
		const set<int> &full_context = ptroot->getFullContext();
		assert(root->getCacheContext().size() == 0);
		z.encodeOR(full_context, m_assignment);
		root->getZobrist().encodeAND(z, 0);
		//Value to use does not matter 
		context_t ctxt = root->getCacheContext();
		m_cache->enter(m_thread_id, root->getZobrist(), var,
			       m_problem->isMap(root->getVar()), ctxt, 0, ELEM_ZERO);
		multipleIterativeDeepeningAND(*root, 0, ELEM_INF, ELEM_INF);
		m_cache->setFirstThread(m_thread_id);

		if (m_thread_id == 0) {
		  double value, upper_bound, pseudo_value;
			int dn;
			unsigned int solved_flag;
			bool flag  = m_cache->read(root->getZobrist(), var,
						   0, root->getCacheContext(), value, upper_bound, pseudo_value, dn, solved_flag);
			assert(flag);
			m_solutionCost = upper_bound;
		}

	} catch (std::bad_alloc& ba) {
		delete[] EmergencyMemory;
		EmergencyMemory = NULL;
		std::cerr << "out of memory";
		res = SEARCH_FAILURE;
	} catch (int) {
		res = SEARCH_TIMEOUT;
	}
	m_result = res;

	if (root) delete root;
	if (EmergencyMemory) delete[] EmergencyMemory;

	return res;
}


bool ParallelRBFAOOWorker::doProcess(ParallelRbfaooSearchNode* node) {
	assert(node);
	if (node->getType() == NODE_AND) {
		int var = node->getVar();
		int val = node->getVal();
		m_assignment[var] = val; // record assignment

	} else { // NODE_OR
		// nothing
	}

	return false; // default

}

// Kishi: DONE
void ParallelRBFAOOWorker::multipleIterativeDeepeningOR(ParallelRbfaooSearchNode &n,
							size_t table_index, double th_value, double upper_bound) {

	bool children_generated_flag = false;

	int num_children;
	int var = n.getVar();
	context_t ctxt = n.getCacheContext();
	unsigned int increased_subtree_size = m_num_expanded;
	double value; 
	int dn, pseudo_best_index;
	unsigned int solved_flag;
#if 000
	for (size_t i = 0; i < m_depth; i ++) 
	  cerr << "|";
	cerr << "//<[OR]\n";
	cerr << "O:e:" << m_num_expanded << "d:=" << m_depth << ":(" << n.getZobrist().getC0() << "," << n.getZobrist().getC1() << "):th=" << th_value << ":ub=" << upper_bound << endl;
#endif
	
	for (;;) {
		double pseudo_second_best_value, new_upper_bound;

		if (!children_generated_flag)
			num_children = generateChildrenOR(n, pseudo_second_best_value, pseudo_best_index,
							  solved_flag);
		else
			calculateNodeValueandDNOR(n, table_index, num_children,
						  pseudo_second_best_value, pseudo_best_index, solved_flag);

		//cerr << "OR: children_generated_flag = " << children_generated_flag << " d=" << n.getDepth() << " num_exp = " << m_num_expanded << " mapvarf = " << m_problem->isMap(n.getVar()) << " z =" << n.getZobrist().getC0() << " value = " << n.getNodeValue() << endl;
		children_generated_flag = true;
		value = n.getNodeValue();
		new_upper_bound = n.getUpperbound();
		upper_bound = min(upper_bound, new_upper_bound);
		double pseudo_value = n.getPseudoValue();
		dn = n.getDN();
		if (dn == DNINF) {
			; assert(solved_flag);
		}
		else {
			; assert(!solved_flag);
		}
		if (m_cache->getFirstThread() >= 0)
			break;
		if (solved_flag || (value + RBFAOO_EPS) >= upper_bound || pseudo_value > th_value + RBFAOO_EPS) {
			break;
		}

		assert(pseudo_best_index >= 0 && pseudo_best_index < num_children);
		ParallelRbfaooSearchNode *c = m_expand[table_index + pseudo_best_index];
		// Kishi: !!!!!!!!!!!!!!!!!!!!! overestmation technique --- later use avarage like my AAAI 2005 paper !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//double new_th_value = max(th_value, second_best_value OP_DIVIDE exp(-n.getDepth()));
		double new_th_value;
		if (m_problem->isMap(c->getVar())) {
			new_th_value = min(th_value, pseudo_second_best_value + m_overestimation);
			new_th_value -= c->getLabel();
			new_upper_bound = upper_bound - c->getLabel();
		}
		else {
			new_th_value = new_upper_bound = ELEM_INF;
		}
		size_t new_table_index = table_index + num_children;
		m_num_expanded++;
		m_num_expanded_or++;
		if (m_problem->isMap(var)) m_num_expanded_or_map++;

		// check for time limit violation
		if (m_timer.timeout()) {
			throw SEARCH_TIMEOUT;
		}

		// output intermediate results if verbose mode is enabled
		if (m_options->verbose) {
			// logging
#ifdef SEARCH_CHECKPOINT_NODE
			if (m_num_expanded_and % SEARCH_CHECKPOINT_NODE == 0 &&
				m_num_expanded_and > m_prevNodeCheckpoint) {

				double timestamp = m_timer.elapsed();
				m_prevNodeCheckpoint = m_num_expanded_and;
				std::cout << "[*" << std::setiosflags(std::ios::fixed) << std::setprecision(2)
					<< std::setw(8) << timestamp << "] w " << std::setiosflags(std::ios::fixed) << std::setprecision(4)
					<< std::setw(8) << m_weight << " "
					<< std::setw(12) << m_space->getOrNodes() << " "
					<< std::setw(12) << m_space->getAndNodes() << " " << std::setiosflags(std::ios::fixed) << std::setprecision(4)
					<< std::setw(12) << m_upperBound
					<< std::endl;
			}
#else
			double elapsed = m_timer.elapsed();
			if (elapsed - m_prevTimeCheckpoint >= SEARCH_CHECKPOINT_TIME) {

				double timestamp = elapsed;
				m_prevTimeCheckpoint = elapsed;
				std::cout << "[*" << std::setiosflags(std::ios::fixed) << std::setprecision(2)
					<< std::setw(8) << timestamp << "] w " << std::setiosflags(std::ios::fixed) << std::setprecision(4)
					<< std::setw(8) << m_weight << " "
					<< std::setw(12) << m_num_expanded_or << " "
					<< std::setw(12) << m_num_expanded_and << " " << std::setiosflags(std::ios::fixed) << std::setprecision(4)
					<< std::endl;
			}
#endif
		}

//		if (m_num_expanded % 1000000 == 0)
//			m_space->dfpncache->printWriteStats();

		doProcess(c);
		context_t child_ctxt = c->getCacheContext();
		if (m_problem->isMap(n.getVar()))
			m_cache->enter(m_thread_id, c->getZobrist(), c->getVar(), 
				       m_problem->isMap(c->getVar()), child_ctxt, 0, value - c->getLabel());
		else
			m_cache->enter(m_thread_id, c->getZobrist(), c->getVar(), 
				       m_problem->isMap(c->getVar()), child_ctxt, 0, ELEM_ZERO);
		m_depth ++;
		multipleIterativeDeepeningAND(*c, new_table_index, new_th_value, new_upper_bound);
		m_depth --;
	}
	increased_subtree_size = m_num_expanded - increased_subtree_size;
	for (size_t i = table_index; i < table_index + num_children; i++)
		delete m_expand[i];
	m_expand.resize(table_index);
	if (solved_flag) 
		n.setUpperbound(n.getNodeValue());
	m_cache->writeandexit(m_thread_id, n.getZobrist(), var, m_problem->isMap(n.getVar()), 
			      ctxt, 1, pseudo_best_index, value, n.getUpperbound(), 
			      dn, solved_flag, (num_children == 0), increased_subtree_size);
#if 000
	for (size_t i = 0; i < m_depth; i ++) 
	  cerr << "*";
	cerr << "O:d=" << m_depth << ":(" << n.getZobrist().getC0() << "," << n.getZobrist().getC1() << "):pv=" << n.getPseudoValue() << ":v=" << value 
	     << ":ub=" << n.getUpperbound() << (solved_flag ? ":t" : ":f") << endl;
	cerr << "//>\n";
#endif
}

// Kishi: DONE
void ParallelRBFAOOWorker::multipleIterativeDeepeningAND(ParallelRbfaooSearchNode &n,
							 size_t table_index, double th_value, double upper_bound) {

	bool children_generated_flag = false;

	int num_children, pseudo_best_index;
	int var = n.getVar();
	context_t ctxt = n.getCacheContext();
	unsigned int increased_subtree_size = m_num_expanded;
	double value;
	int dn;
	unsigned int solved_flag;
#if 000
	for (size_t i = 0; i < m_depth; i ++) 
	  cerr << "|";
	cerr << "//<[AND]\n";
	cerr << "A:e=" << m_num_expanded << ":d=" << m_depth << ":(" << n.getZobrist().getC0() << "," << n.getZobrist().getC1() << "):th=" << th_value << ":ub=" << upper_bound << endl;
#endif
	for (;;) {
		int second_best_dn;

		if (!children_generated_flag)
			num_children = generateChildrenAND(n, second_best_dn, pseudo_best_index,
							   solved_flag);
		else
			calculateNodeValueandDNAND(n, table_index, num_children,
					second_best_dn, pseudo_best_index, solved_flag);

		//cerr << "AND: children_generated_flag = " << children_generated_flag << " d=" << n.getDepth() << " num_exp = "  << m_num_expanded << " mapvarf = " << m_problem->isMap(n.getVar()) << " z =" << n.getZobrist().getC0() << " value = " << n.getNodeValue() << endl;
		children_generated_flag = true;

		value = n.getNodeValue();
		double pseudo_value = n.getPseudoValue();
		double new_upper_bound = n.getUpperbound();
		upper_bound = min(new_upper_bound, upper_bound);
		dn = n.getDN();
		if (dn == DNINF) {
			; assert(solved_flag);
		}
		else {
			; assert(!solved_flag);
		}
		if (m_cache->getFirstThread() >= 0)
			break;
		if (solved_flag || (value + RBFAOO_EPS) >= upper_bound || pseudo_value > th_value + RBFAOO_EPS)
		  break;

		ParallelRbfaooSearchNode *c = m_expand[table_index + pseudo_best_index];
		double c_pseudo_value = c->getPseudoValue();
		double c_value = c->getNodeValue();

		double new_th_value;
		if (m_problem->isMap(c->getVar())) {
			new_th_value = (th_value - pseudo_value) + c_pseudo_value;
			new_upper_bound  = upper_bound - value + c_value;
		}
		else {
			new_th_value = new_upper_bound = ELEM_INF;
		}
		; assert(pseudo_best_index >= 0 && pseudo_best_index < num_children);
		size_t new_table_index = table_index + num_children;
		m_num_expanded++;
		m_num_expanded_and++;
		if (m_problem->isMap(var)) m_num_expanded_and_map++;

		// check for time limit violation
		if (m_timer.timeout()) {
			throw SEARCH_TIMEOUT;
		}

		// output intermediate results if verbose mode is enabled
		if (m_options->verbose) {
			// logging
#ifdef SEARCH_CHECKPOINT_NODE
			if (m_num_expanded_and % SEARCH_CHECKPOINT_NODE == 0 &&
				m_num_expanded_and > m_prevNodeCheckpoint) {

				double timestamp = m_timer.elapsed();
				m_prevNodeCheckpoint = m_num_expanded_and;
				std::cout << "[*" << std::setiosflags(std::ios::fixed) << std::setprecision(2)
					<< std::setw(8) << timestamp << "] w " << std::setiosflags(std::ios::fixed) << std::setprecision(4)
					<< std::setw(8) << m_weight << " "
					<< std::setw(12) << m_num_expanded_or << " "
					<< std::setw(12) << m_num_expanded_and << " " << std::setiosflags(std::ios::fixed) << std::setprecision(4)
					<< std::endl;
			}
#else
			double elapsed = m_timer.elapsed();
			if (elapsed - m_prevTimeCheckpoint >= SEARCH_CHECKPOINT_TIME) {

				double timestamp = elapsed;
				m_prevTimeCheckpoint = elapsed;
				std::cout << "[*" << std::setiosflags(std::ios::fixed) << std::setprecision(2)
					<< std::setw(8) << timestamp << "] w " << std::setiosflags(std::ios::fixed) << std::setprecision(4)
					<< std::setw(8) << m_weight << " "
					<< std::setw(12) << m_num_expanded_or << " "
					<< std::setw(12) << m_num_expanded_and << " " << std::setiosflags(std::ios::fixed) << std::setprecision(4)
					<< std::endl;
			}
#endif
		}


//		if (m_num_expanded % 1000000 == 0)
//			m_space->dfpncache->printWriteStats();
		context_t child_ctxt = c->getCacheContext();
		m_cache->enter(m_thread_id, c->getZobrist(), c->getVar(), 
			       m_problem->isMap(c->getVar()), child_ctxt, 1, c_value);
		m_depth ++;
		multipleIterativeDeepeningOR(*c, new_table_index, new_th_value, new_upper_bound);
		m_depth --;
	}
	increased_subtree_size = m_num_expanded - increased_subtree_size;
	for (size_t i = table_index; i < table_index + num_children; i++)
		delete m_expand[i];
	m_expand.resize(table_index);
	if (solved_flag) 
		n.setUpperbound(n.getNodeValue());
	m_cache->writeandexit(m_thread_id, n.getZobrist(), var,  m_problem->isMap(n.getVar()), 
			      ctxt, 0, pseudo_best_index, value, n.getUpperbound(), 
			      dn, solved_flag, (num_children == 0), increased_subtree_size);
#if 000
	for (size_t i = 0; i < m_depth; i ++) 
	  cerr << "*";
	cerr << "A:d=" << m_depth << ":(" << n.getZobrist().getC0() << "," << n.getZobrist().getC1() << "):pv=" << n.getPseudoValue() << ":v=" << value 
	     << ":ub=" << n.getUpperbound() << (solved_flag ? ":t" : ":f") << endl;
	cerr << "//>\n";
#endif

}

// Kishi: DONE
int ParallelRBFAOOWorker::generateChildrenAND(ParallelRbfaooSearchNode &n, int &second_best_dn,
					      int &best_index, unsigned int &solved_flag) {

	assert(n.getType() == NODE_AND);

	//m_space->nodesAND += 1;
	//KISHI: COMPUTE THIS LOCALLY TO AVOID TRUE SHARING!
	//m_space->incNodesExpanded(NODE_AND);

	int var = n.getVar();
	PseudotreeNode* ptnode = m_pseudotree->getNode(var);
	int depth = ptnode->getDepth();

	double value = ELEM_ZERO, pseudo_value = ELEM_ZERO, upper_bound = ELEM_ZERO; //ELEM_ONE;
	int dn = DNINF;

	second_best_dn = DNINF;
	best_index = UNDETERMINED;
	solved_flag = 0;

	int num_children = 0;
	bool child_map_var_flag = m_problem->isMap(n.getVar());
	// create new OR children (going in reverse due to reversal on stack),
	for (vector<PseudotreeNode*>::const_iterator it =
			ptnode->getChildren().begin(); it != ptnode->getChildren().end();
			++it, num_children++) {
		int vChild = (*it)->getVar();
		PseudotreeNode* ptnode = m_pseudotree->getNode(vChild);
		const set<int> &full_context = ptnode->getFullContext();
		ParallelRbfaooSearchNode* c = new ParallelRbfaooSearchNode(NODE_OR, vChild, depth + 1);
		for (set<int>::const_iterator it = full_context.begin();
				it != full_context.end(); ++it) {
			assert(*it != vChild);
			c->addCacheContext(m_assignment[*it]);
		}
		Zobrist &z = c->getZobrist();
		z.encodeOR(full_context, m_assignment);
		//RbfaooCacheElementPtr p = m_space->dfpncache->read(z, vChild, 1, c->getCacheContext());
		double child_value, child_ubound, child_pseudo_value;
		int child_dn;
		unsigned int child_solved_flag;
		bool entry_flag = m_cache->read(z, vChild, 1, c->getCacheContext(), child_value, child_ubound, child_pseudo_value, child_dn, child_solved_flag);
		child_map_var_flag = m_problem->isMap(vChild);
		if (!entry_flag || (!m_problem->isMap(vChild) && !child_solved_flag)) {
			// Compute and set heuristic estimate, includes child labels
			// Kishi: perhaps we need to use a different heuristic
			child_value = heuristicOR(*c);
			if (m_problem->isSum(vChild)) child_value = ELEM_ZERO;
			child_pseudo_value = child_value;
			child_ubound = ELEM_INF;
			child_dn = 1;
			child_solved_flag = 0;
		} else 
			//!!!!!!!!!!!!!!!!!!!!!!!!!KISHI: ADHOC FIX! NEED TO THINK ABOUT A WAY TO REMOVE THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			heuristicOR(*c);
			// edge cost is zero from AND node to OR child
		c->setNodeValue(child_value);
		c->setDN(child_dn);
		c->setUpperbound(child_ubound);
		c->setPseudoValue(child_pseudo_value);
		value += child_value;
		upper_bound += child_ubound;
		pseudo_value += child_pseudo_value;
		if (child_dn < dn) {
			second_best_dn = dn;
			dn = child_dn;
			best_index = num_children;
		} else if (child_dn < second_best_dn)
			second_best_dn = child_dn;
		m_expand.push_back(c);
	} // for loop over new OR children

	n.setNodeValue(value);
	n.setDN(dn);
	n.setUpperbound(upper_bound);
	n.setPseudoValue(pseudo_value);

	if (dn == DNINF)
		solved_flag = 1;
	if (!child_map_var_flag && n.getNodeValue() == ELEM_INF) {
		dn = 0;
		solved_flag = 1;
	}
#if 000
	if (num_children == 0)
	cout << "terminal AND: m_num_expanded = " << m_num_expanded << "\n";
#endif

	return num_children;
} // ParallelRBFAOOWorker::generateChildrenAND

// Kishi: DONE
void ParallelRBFAOOWorker::calculateNodeValueandDNAND(ParallelRbfaooSearchNode &n, size_t table_index,
		int num_children, int &second_best_dn, int &best_index,
		unsigned int &solved_flag) {


	assert(n.getType() == NODE_AND);

	double value = ELEM_ZERO, pseudo_value = ELEM_ZERO, upper_bound = ELEM_ZERO; //ELEM_ONE;
	int dn = DNINF;
	second_best_dn = DNINF;
	best_index = UNDETERMINED;

	solved_flag = 0;
	bool child_map_var_flag = m_problem->isMap(n.getVar());
	// create new OR children (going in reverse due to reversal on stack),
	for (size_t i = 0; i < static_cast<size_t>(num_children); i++) {
		ParallelRbfaooSearchNode *c = m_expand[table_index + i];
		Zobrist &z = c->getZobrist();
		int vChild = c->getVar();
		//RbfaooCacheElementPtr p = m_space->dfpncache->read(z, vChild, 1, c->getCacheContext());
		double child_value, child_ubound, child_pseudo_value;
		int child_dn;
		unsigned int child_solved_flag;
		bool entry_flag = m_cache->read(z, vChild, 1, c->getCacheContext(), child_value, child_ubound, child_pseudo_value, child_dn, child_solved_flag);
		child_map_var_flag = m_problem->isMap(vChild);
		if (!entry_flag || (!m_problem->isMap(vChild) && !child_solved_flag)) {
			// Compute and set heuristic estimate, includes child labels
			// Kishi: perhaps we need to use a different heuristic
			child_value = c->getHeur();
			child_dn = 1;
			child_pseudo_value = child_value;
			child_ubound = ELEM_INF;
			child_dn = 1;
			child_solved_flag = 0;
		}
		c->setNodeValue(child_value);
		c->setDN(child_dn);
		c->setUpperbound(child_ubound);
		c->setPseudoValue(child_pseudo_value);
		value += child_value;
		upper_bound += child_ubound;
		pseudo_value += child_pseudo_value;
		if (child_dn < dn) {
			second_best_dn = dn;
			dn = child_dn;
			best_index = i;
		} else if (child_dn < second_best_dn)
			second_best_dn = child_dn;
	} // for loop over new OR children

	n.setNodeValue(value);
	n.setDN(dn);
	n.setUpperbound(upper_bound);
	n.setPseudoValue(pseudo_value);

	if (dn == DNINF)
		solved_flag = 1;
	if (!child_map_var_flag && n.getNodeValue() == ELEM_INF) {
		dn = 0;
		solved_flag = 1;
	}
} // ParallelRBFAOOWorker::calculateNodeValueandDNAND

// Kishi: DONE
int ParallelRBFAOOWorker::generateChildrenOR(ParallelRbfaooSearchNode &n, double &pseudo_second_best_value, 
					     int &pseudo_best_index, unsigned int &solved_flag) {

	assert(n.getType() == NODE_OR);

	//m_space->nodesOR +=1;
	//KISHI: COMPUTE THIS LOCALLY TO AVOID TRUE SHARING!
	//m_space->incNodesExpanded(NODE_OR);

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
	double pseudo_value, value, upper_bound = ELEM_INF; //ELEM_ZERO;
	if (mapVar) {
		pseudo_value = value = ELEM_INF; //ELEM_ZERO;
	}
	else {
		pseudo_value = value = ELEM_ZERO;
	}
	int dn = 0, dn_with_ignoring_proven_child = 0;
	int num_children = 0, best_index = UNDETERMINED;
	unsigned int min_threads = m_num_threads + 1;

	pseudo_best_index = UNDETERMINED;
	pseudo_second_best_value = ELEM_INF;

	for (val_t i = 0; i < m_problem->getDomainSize(var); i++, num_children++) {
		ParallelRbfaooSearchNode *c = new ParallelRbfaooSearchNode(NODE_AND, var, i, depth,
				heur[2 * i + 1]); // uses cached label
		Zobrist &z = c->getZobrist();
		z.encodeAND(n.getZobrist(), i);
		c->setCacheContext(n.getCacheContext());
		context_t &context = c->getCacheContext();
		//RbfaooCacheElementPtr p = m_space->dfpncache->read(z, var, 0, context);
		double child_value, child_ubound, child_pseudo_value;
		int child_dn;
		unsigned int child_solved_flag, num_threads;
		bool entry_flag = m_cache->read(z, var, 0, context, child_value, child_ubound, child_pseudo_value, child_dn, child_solved_flag, num_threads);
		c->setHeur(heur[2 * i]);
		if (!entry_flag) {
			child_value = heur[2 * i]; 
			child_pseudo_value = child_value;
			child_dn = 1;
			child_ubound = ELEM_INF;
			child_solved_flag = 0;
			num_threads = 0;
		} else {
			child_value += c->getLabel();
			child_pseudo_value += c->getLabel();
			child_ubound +=  c->getLabel();
		}
		if (heur[2 * i] == ELEM_INF) {
			; assert(child_value == ELEM_INF);
			child_dn = 0;
		}
		if (mapVar) {
			//best_index in terms of the pseudo best
			if (child_pseudo_value < pseudo_value) {
				pseudo_best_index = i;
				pseudo_second_best_value = pseudo_value;
				pseudo_value = child_pseudo_value;
			} else if (child_pseudo_value <= pseudo_second_best_value)
			  pseudo_second_best_value = child_pseudo_value;
			//best value itself
			if (child_value < value) {
			  best_index = i;
			  value = child_value;
			  solved_flag = child_solved_flag;
			} else if (child_value == value && child_solved_flag) {
			  best_index = i;
			  solved_flag = child_solved_flag;
			}
			if (child_dn != DNINF) 
			  dn_with_ignoring_proven_child = sumDN(dn_with_ignoring_proven_child, child_dn);
			upper_bound = min(child_ubound, upper_bound);
		}
		else {
#if 000
			//CONFIRM THIS WITH RADU TOMORROW!!!
			pseudo_value += child_pseudo_value;
			value += child_value;
#else
			if (!child_solved_flag) 
				child_value = child_pseudo_value = heur[2 * i];
			pseudo_value += 1/ELEM_DECODE(child_pseudo_value);
			value += 1/ELEM_DECODE(child_value);
#endif
			solved_flag &= child_solved_flag;
			if (!child_solved_flag) {
			  if (!entry_flag || num_threads < min_threads) {
				pseudo_best_index = i;
				min_threads = num_threads;
			  }
			}
			if (child_dn != DNINF) 
				dn_with_ignoring_proven_child = sumDN(dn_with_ignoring_proven_child, child_dn);
		}
		dn = sumDN(dn, child_dn);
		c->setNodeValue(child_value);
		c->setDN(child_dn);
		c->setUpperbound(child_ubound);
		c->setPseudoValue(child_pseudo_value);
		m_expand.push_back(c);
	}
	//ASK RADU IF ADDING THE FOLLOWING OP IS CORRECT
	if (!mapVar) {
	  value = -ELEM_ENCODE(value);
	  pseudo_value = -ELEM_ENCODE(pseudo_value);
	}
	if (solved_flag) 
		pseudo_best_index = best_index;
	if (!solved_flag && dn == DNINF) {
	  //need to recalculate dn
		; assert(dn_with_ignoring_proven_child != 0);
		dn = max(1, dn_with_ignoring_proven_child);
	}
	n.setNodeValue(value);
	n.setDN(dn);
	n.setUpperbound(upper_bound);
	n.setPseudoValue(pseudo_value);
#if 000
	if (num_children == 0)
	cout << "terminal OR: m_num_expanded = " << m_num_expanded << "\n";
#endif

	return num_children;
} // ParallelRBFAOOWorker::generateChildrenOR

// Kishi: DONE
void ParallelRBFAOOWorker::calculateNodeValueandDNOR(ParallelRbfaooSearchNode &n, size_t table_index, 
						     int num_children, double &pseudo_second_best_value, int &pseudo_best_index,
						     unsigned int &solved_flag) {

	bool mapVar = m_problem->isMap(n.getVar());
	assert(n.getType() == NODE_OR);

	solved_flag = 0;
	if (!mapVar) solved_flag = 1;

	double* heur = n.getHeurCache();
	pseudo_second_best_value = ELEM_INF; //ELEM_ZERO;
	double pseudo_value, value, upper_bound = ELEM_INF; //ELEM_ZERO;
	if (mapVar) {
		pseudo_value = value = ELEM_INF; //ELEM_ZERO;
	}
	else {
		pseudo_value = value = ELEM_ZERO;
	}
	int dn = 0, dn_with_ignoring_proven_child = 0, best_index = UNDETERMINED;
	unsigned int min_threads = m_num_threads + 1;
	pseudo_best_index = UNDETERMINED;

	for (int i = 0; i < num_children; i++) {
		ParallelRbfaooSearchNode *c = m_expand[table_index + i];
		Zobrist &z = c->getZobrist();
		context_t &context = c->getCacheContext();
		int vChild = c->getVar();
		//RbfaooCacheElementPtr p = m_space->dfpncache->read(z, vChild, 0, context);

		double child_value, child_ubound, child_pseudo_value;
		int child_dn;
		unsigned int child_solved_flag, num_threads;
		bool entry_flag = m_cache->read(z, vChild, 0, context, child_value, child_ubound, child_pseudo_value, child_dn, child_solved_flag, num_threads);
		if (!entry_flag) {
			child_solved_flag = 0;
			child_ubound = ELEM_INF;
			child_value = c->getHeur();
			if (heur[2 * i] == ELEM_INF) {
			  ; assert(child_value == ELEM_INF);
				child_dn = 0;
			} else 
				child_dn = 1;
			child_pseudo_value = child_value;
			num_threads = 0;
		} else {
			child_value += c->getLabel();
			child_pseudo_value += c->getLabel();
			child_ubound += c->getLabel();
		}

		if (mapVar) {
			if (!child_solved_flag && child_pseudo_value < pseudo_value) {
				pseudo_best_index = i;
				pseudo_second_best_value = pseudo_value;
				pseudo_value = child_pseudo_value;
			} else if (!child_solved_flag && child_pseudo_value <= pseudo_second_best_value) 
				pseudo_second_best_value = child_pseudo_value;
			//best value itself
			if (child_value < value) {
				best_index = i;
				value = child_value;
				solved_flag = child_solved_flag;
			} else if (child_value == value && child_solved_flag) {
				best_index = i;
				solved_flag = child_solved_flag;
			}
			if (child_dn != DNINF) 
				dn_with_ignoring_proven_child = sumDN(dn_with_ignoring_proven_child,
								      child_dn);
			upper_bound = min(child_ubound, upper_bound);
		}
		else {
			if (!child_solved_flag) 
				child_value = child_pseudo_value = c->getHeur();
			pseudo_value += 1/ELEM_DECODE(child_pseudo_value);
			value += 1/ELEM_DECODE(child_value);
			solved_flag &= child_solved_flag;
			if (!child_solved_flag) {
			  if (!entry_flag || num_threads < min_threads) {
				pseudo_best_index = i;
				min_threads = num_threads;
			  }
			}
			if (child_dn != DNINF) 
				dn_with_ignoring_proven_child = sumDN(dn_with_ignoring_proven_child, child_dn);
		}
		dn = sumDN(dn, child_dn);
		c->setNodeValue(child_value);
		c->setDN(child_dn);
		c->setUpperbound(child_ubound);
		c->setPseudoValue(child_pseudo_value);
	}
	if (solved_flag) 
		pseudo_best_index = best_index;
	if (!solved_flag && dn == DNINF) {
		//need to recalculate dn
		; assert(dn_with_ignoring_proven_child != 0);
		dn = max(1, dn_with_ignoring_proven_child);
	}

	if (!mapVar) {
		value = -ELEM_ENCODE(value);
		pseudo_value = -ELEM_ENCODE(pseudo_value);
	}
	n.setNodeValue(value);
	n.setDN(dn);
	n.setPseudoValue(pseudo_value);
	n.setUpperbound(upper_bound);
} // ParallelRBFAOOWorker::generateChildrenOR

// Kishi: DONE (we have a MIN-SUM problem)
double ParallelRBFAOOWorker::heuristicOR(ParallelRbfaooSearchNode &n) {

	int var = n.getVar();
	bool mapVar = m_problem->isMap(var);
	int oldValue = m_assignment[var];
	double* dv = new double[m_problem->getDomainSize(var) * 2];
	double h = (mapVar) ? INFINITY : ELEM_ZERO; // the new OR nodes h value
	const list<Function*>& funs = m_pseudotree->getFunctions(var);
	for (val_t i = 0; i < m_problem->getDomainSize(var); ++i) {
		m_assignment[var] = i;
		double m_epsilon = 1.0;
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

} // ParallelRBFAOOWorker::heuristicOR

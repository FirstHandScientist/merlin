/*
 * pseudo_tree.cpp
 *
 *  Created on: Apr 2, 2013
 *      Author: radu
 */

#include "pseudo_tree.h"

#undef DEBUG

void Pseudotree::resetFunctionInfo(const std::vector<Function*>& fns) {
	// reset existing function mapping (if any)
	for (std::vector<PseudotreeNode*>::iterator it = m_nodes.begin();
			it != m_nodes.end(); ++it)
		(*it)->resetFunctions();
	// generate new function mapping
	for (std::vector<Function*>::const_iterator itF = fns.begin();
			itF != fns.end(); ++itF) {
		const std::set<int>& scope = (*itF)->getScope();
		if (scope.size() == 0) {
			m_nodes[m_elimOrder.back()]->addFunction(*itF);
			continue;
		}
		std::vector<int>::const_iterator it = m_elimOrder.begin();
		for (;; ++it) {
			if (scope.find(*it) != scope.end()) {
				m_nodes[*it]->addFunction(*itF);
				break;
			}
		}
	}
}

/* builds the pseudo tree according to order 'elim' */
void Pseudotree::build(Graph G, const std::vector<int>& elim, const int cachelimit) {

	if (m_height != UNKNOWN) {
		this->reset();
	}

	const int n = G.getNumNodes();
	assert(n == (int) m_nodes.size());
	assert(n == (int) elim.size());

	std::list<PseudotreeNode*> roots;

	Graph GG(G);
	GG.triangulate(elim);

	// build actual pseudo tree(s)
	for (std::vector<int>::const_iterator it = elim.begin(); it != elim.end(); ++it) {
		const set<int>& N = G.getNeighbors(*it);
		m_width = max(m_width, (int) N.size());
		insertNewNode(*it, N, roots);
		G.addClique(N);
		G.removeNode(*it);
	}

	std::cout << "Number of disconnected sub pseudo trees: " << roots.size() << std::endl;

	// compute contexts for adaptive caching
	for (vector<PseudotreeNode*>::iterator it = m_nodes.begin();
			it != m_nodes.end(); ++it) {
		const set<int>& ctxt = (*it)->getFullContext();
		if (cachelimit == NONE || cachelimit >= (int) ctxt.size()) {
			(*it)->setCacheContext(ctxt);
			continue; // proceed to next node
		}
		int j = cachelimit;
		PseudotreeNode* p = (*it)->getParent();
		while (j--) { // note: cachelimit < ctxt.size()
			while (ctxt.find(p->getVar()) == ctxt.end())
				p = p->getParent();
			(*it)->addCacheContext(p->getVar());
			p = p->getParent();
		}
		// find reset variable for current node's cache table
		while (ctxt.find(p->getVar()) == ctxt.end())
			p = p->getParent();
		p->addCacheReset((*it)->getVar());
		//    cout << "AC for var. " << (*it)->getVar() << " context " << (*it)->getCacheContext()
		//         << " out of " << (*it)->getFullContext() << ", reset at " << p->getVar() << endl;
	}

	// add artificial root node to connect disconnected components
	int bogusIdx = elim.size();
	PseudotreeNode* p = new PseudotreeNode(this, bogusIdx, std::set<int>());
	for (std::list<PseudotreeNode*>::iterator it = roots.begin();
			it != roots.end(); ++it) {
		p->addChild(*it);
		(*it)->setParent(p);
		++m_components; // increase component count
	}
	m_nodes.push_back(p);
	m_root = p;

	// remember the elim. order
	m_elimOrder = elim;
	m_elimOrder.push_back(bogusIdx); // add dummy variable as last node in ordering

	// initiate depth/height computation for tree and its nodes (bogus variable
	// gets depth -1), then need to subtract 1 from height for bogus as well
	m_height = m_root->updateDepthHeight(-1) - 1;

	// initiate subproblem width computation (recursive)
	m_root->updateSubWidth();

//  const vector<Pseudotree*> rootC = m_root->getChildren();
//  for (vector<Pseudotree*>::iterator it = rootC.begin(); it != rootC.end(); ++it) {
//    cout << (*it)->getSubWidth();
//  }

	// reorder each list of child nodes by complexity (width/height)
	for (std::vector<PseudotreeNode*>::iterator it = m_nodes.begin();
			it != m_nodes.end(); ++it)
		(*it)->orderChildren(m_subOrder);

	// initiate subproblem variables computation (recursive)
	m_root->updateSubprobVars(m_nodes.size());  // includes dummy
	m_size = m_root->getSubprobSize() - 1;  // -1 for dummy

	m_sumRoots.resize(m_size+1);
	for (std::vector<PseudotreeNode*>::iterator it = m_nodes.begin();
			it != m_nodes.end(); ++it) {

		PseudotreeNode *n = (*it);
		PseudotreeNode *p = n->getParent();
		if (p != NULL && (p->getVar() == bogusIdx || m_problem->isMap(p->getVar()))) {
			if (m_problem->isSum(n->getVar()))
				m_sumRoots[n->getVar()] = true;
		}
	}

	// create the AND contexts
	for (std::vector<PseudotreeNode*>::iterator it = m_nodes.begin();
			it != m_nodes.end(); ++it) {

		PseudotreeNode *n = (*it);
		if (n->getVar() == bogusIdx) {
			continue; // skip dummy variable
		}

		std::set<int> ancestors; // this is the and context
		const std::set<int>& descendants = n->getSubprobVars();
		ancestors.insert(n->getVar());

		// find the ancestors of the node that are connected to descendats of the node
		PseudotreeNode *p = n->getParent();
		while (p != NULL) {

			bool found = false;
			for (std::set<int>::const_iterator si = descendants.begin();
					si != descendants.end(); ++si) {
				int desc = (*si);
				if (desc == n->getVar()) {
					continue; // skip current node variable
				}

				if (GG.hasEdge(p->getVar(), desc)) {
					found = true;
					break;
				}
			}

			// found an ancestor connected to a descendant in the induced graph
			if (found) {
				ancestors.insert(p->getVar());
			}

			p = p->getParent();
		}

		// set the AND context
		n->setAndContext(ancestors);
	}

#ifdef DEBUG
	// dump the contexts (OR)
	std::cout << "OR contexts:" << std::endl;
	for (std::vector<PseudotreeNode*>::iterator it = m_nodes.begin();
			it != m_nodes.end(); ++it) {
		(*it)->dumpFullContext(std::cout);
	}

	// dump the contexts (AND)
	std::cout << "AND contexts:" << std::endl;
	for (std::vector<PseudotreeNode*>::iterator it = m_nodes.begin();
			it != m_nodes.end(); ++it) {
		(*it)->dumpAndContext(std::cout);
	}
#endif

	return;
}

/* builds the pseudo tree according to order 'elim' */
void Pseudotree::forceMapChain(const std::list<int>& mapSearchOrder,
		const int cachelimit) {

	// safety checks
	assert( m_problem->hasDummy() );

	std::cout << "Re-arranging the MAP pseudo tree (as chain)" << std::endl;
	std::cout << "  MAP search order:";
	std::copy(mapSearchOrder.begin(), mapSearchOrder.end(),
			std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;

	// create the primal graph
	Graph G(m_problem->getN());
	for (int var = 0; var < m_problem->getN(); ++var) {
		G.addNode(var);
	}
	const std::vector<Function*>& fns = m_problem->getFunctions();
	for (std::vector<Function*>::const_iterator it = fns.begin();
			it != fns.end(); ++it) {
		Function* f = (*it);
		G.addClique(f->getScope());
	}

	// triangulate the graph according to the constrained elimination order
	std::vector<int> elim = m_elimOrder;
	G.triangulate(elim);

	// collect the SUM nodes whose parents are MAP variables
	std::list<PseudotreeNode*> roots;
	for (std::vector<PseudotreeNode*>::iterator it = m_nodes.begin();
			it != m_nodes.end(); ++it) {
		PseudotreeNode* n = (*it);
		int var = n->getVar();
		if ( !m_problem->isMap(var) ) { // SUM variable
			PseudotreeNode* p = n->getParent();
			assert(p != NULL); // SUM variables cannot be root
			if ( m_problem->isMap(p->getVar()) ) { // has MAP parent
				roots.push_back(n);
			}
		}
	}

	std::cout << "  Number of SUM sub pseudo trees: " << roots.size() << std::endl;
	m_sumRoots.resize(m_size+1);
	for (std::list<PseudotreeNode*>::iterator li = roots.begin();
			li != roots.end(); ++li) {
		PseudotreeNode* r = *li;
		m_sumRoots[ r->getVar() ] = true;
	}

	// MAP variables form the start pseudo tree, arrange them as chain
	std::list<int>::const_iterator li = mapSearchOrder.begin();
	m_root = m_nodes[*li];
	m_root->setParent(NULL);
	PseudotreeNode* prev = m_root;
	li++;
	while ( li != mapSearchOrder.end() ){
		int var = (*li);
		m_nodes[var]->setParent(prev);
		prev->setChild(m_nodes[var]);
		prev = m_nodes[var];
		li++;
	}

	bool firstInList = true;
	PseudotreeNode* last = prev;
	for (std::list<PseudotreeNode*>::iterator li = roots.begin();
			li != roots.end(); ++li) {
		PseudotreeNode* n = (*li);
		n->setParent(last);
		if (firstInList) {
			last->setChild(n);
			firstInList = false;
		} else {
			last->addChild(n);
		}
	}

	// initiate depth/height computation for tree and its nodes
	m_height = m_root->updateDepthHeight(0);

	// initiate subproblem width computation (recursive)
	m_root->updateSubWidth();

	// reorder each list of child nodes by complexity (width/height)
	for (std::vector<PseudotreeNode*>::iterator it = m_nodes.begin();
			it != m_nodes.end(); ++it)
		(*it)->orderChildren(m_subOrder);

	// initiate subproblem variables computation (recursive)
	m_root->updateSubprobVars(m_nodes.size());  // includes dummy
	m_size = m_root->getSubprobSize();  // -1 for dummy

	// recompute the OR contexts
	for (std::vector<PseudotreeNode*>::iterator it = m_nodes.begin();
			it != m_nodes.end(); ++it) {

		PseudotreeNode *n = (*it);
		if (n->getParent() == NULL) {
			continue; // skip dummy variable
		}

		std::set<int> ancestors; // this is the context
		const std::set<int>& descendants = n->getSubprobVars();

		// find the ancestors of the node that are connected to it or to its descendants
		PseudotreeNode *p = n->getParent();
		while (p != NULL) {

			bool found = false;
			for (std::set<int>::const_iterator si = descendants.begin();
					si != descendants.end(); ++si) {
				int desc = (*si); // includes current node as well

				if (G.hasEdge(p->getVar(), desc)) {
					found = true;
					break;
				}
			}

			// found an ancestor connected to a descendant in the induced graph
			if (found) {
				ancestors.insert(p->getVar());
			}

			p = p->getParent();
		}

		// set the context
		n->setFullContext(ancestors);
	}

	// compute contexts for adaptive caching
	for (vector<PseudotreeNode*>::iterator it = m_nodes.begin();
			it != m_nodes.end(); ++it) {
		const set<int>& ctxt = (*it)->getFullContext();
		if (cachelimit == NONE || cachelimit >= (int) ctxt.size()) {
			(*it)->setCacheContext(ctxt);
			continue; // proceed to next node
		}
		int j = cachelimit;
		PseudotreeNode* p = (*it)->getParent();
		while (j--) { // note: cachelimit < ctxt.size()
			while (ctxt.find(p->getVar()) == ctxt.end())
				p = p->getParent();
			(*it)->addCacheContext(p->getVar());
			p = p->getParent();
		}
		// find reset variable for current node's cache table
		while (ctxt.find(p->getVar()) == ctxt.end())
			p = p->getParent();
		p->addCacheReset((*it)->getVar());
		//    cout << "AC for var. " << (*it)->getVar() << " context " << (*it)->getCacheContext()
		//         << " out of " << (*it)->getFullContext() << ", reset at " << p->getVar() << endl;
	}

	// create the AND contexts
	for (std::vector<PseudotreeNode*>::iterator it = m_nodes.begin();
			it != m_nodes.end(); ++it) {

		PseudotreeNode *n = (*it);

		std::set<int> ancestors; // this is the and context
		const std::set<int>& descendants = n->getSubprobVars();
		ancestors.insert(n->getVar());

		// find the ancestors of the node that are connected to descendats of the node
		PseudotreeNode *p = n->getParent();
		while (p != NULL) {

			bool found = false;
			for (std::set<int>::const_iterator si = descendants.begin();
					si != descendants.end(); ++si) {
				int desc = (*si);
				if (desc == n->getVar()) {
					continue; // skip current node variable
				}

				if (G.hasEdge(p->getVar(), desc)) {
					found = true;
					break;
				}
			}

			// found an ancestor connected to a descendant in the induced graph
			if (found) {
				ancestors.insert(p->getVar());
			}

			p = p->getParent();
		}

		// set the AND context
		n->setAndContext(ancestors);
	}

	// remap the elimination order (dfs)
	std::list<int> order;
	std::stack<PseudotreeNode*> dfs;
	dfs.push(m_root);
	while ( !dfs.empty() ) {
		PseudotreeNode* n = dfs.top();
		dfs.pop();
		order.push_front(n->getVar());
		for (std::vector<PseudotreeNode*>::const_iterator it = n->getChildren().begin();
				it != n->getChildren().end(); ++it) {
			dfs.push( (*it) );
		}
	}

	m_elimOrder.clear();
	std::copy(order.begin(), order.end(), std::back_inserter(m_elimOrder));

#ifdef DEBUG
	// dump the contexts (OR)
	std::cout << "OR contexts:" << std::endl;
	for (std::vector<PseudotreeNode*>::iterator it = m_nodes.begin();
			it != m_nodes.end(); ++it) {
		(*it)->dumpFullContext(std::cout);
	}

	// dump the contexts (AND)
	std::cout << "AND contexts:" << std::endl;
	for (std::vector<PseudotreeNode*>::iterator it = m_nodes.begin();
			it != m_nodes.end(); ++it) {
		(*it)->dumpAndContext(std::cout);
	}
#endif

	return;
}

Pseudotree::Pseudotree(const Pseudotree& pt) {

	m_problem = pt.m_problem;
	m_size = pt.m_size;
	m_subOrder = pt.m_subOrder;

	m_nodes.reserve(m_size + 1);
	m_nodes.resize(m_size + 1, NULL);

	m_elimOrder = pt.m_elimOrder;

	// clone PseudotreeNode structure
	std::stack<PseudotreeNode*> stack;

	PseudotreeNode *ptnNew = NULL, *ptnPar = NULL;
	// start with bogus variable, always has highest index
	m_nodes.at(m_size) = new PseudotreeNode(this, m_size, set<int>());
	m_root = m_nodes.at(m_size);
	stack.push(m_root);

	while (!stack.empty()) {
		ptnPar = stack.top();
		stack.pop();
		int var = ptnPar->getVar();
		const std::vector<PseudotreeNode*>& childOrg =
				pt.m_nodes.at(var)->getChildren();
		for (std::vector<PseudotreeNode*>::const_iterator it = childOrg.begin();
				it != childOrg.end(); ++it) {
			ptnNew = new PseudotreeNode(this, (*it)->getVar(),
					(*it)->getFullContext());
			m_nodes.at((*it)->getVar()) = ptnNew;
			ptnPar->addChild(ptnNew);
			ptnNew->setParent(ptnPar);
			ptnNew->setFunctions((*it)->getFunctions());
			ptnNew->setCacheContext( (*it)->getCacheContext() );
			ptnNew->setCacheReset( (*it)->getCacheReset() );

			stack.push(ptnNew);
		}
	}

	m_height = m_root->updateDepthHeight(-1) - 1;
	m_root->updateSubWidth();
	m_root->updateSubprobVars(m_nodes.size());
	m_size = m_root->getSubprobSize() - 1; // -1 for bogus

	m_width = pt.m_width;
	m_components = pt.m_components;

	// reorder each list of child nodes by complexity (width/height)
	for (std::vector<PseudotreeNode*>::iterator it = m_nodes.begin();
			it != m_nodes.end(); ++it) {
		(*it)->orderChildren(m_subOrder);
	}

}

/*
 * orders the sub pseudo trees of a node ("less" -> ascending by w*)
 */
void PseudotreeNode::orderChildren(int subOrder) {
	if (subOrder == SUBPROB_WIDTH_INC) {
		sort(m_children.begin(), m_children.end(), PseudotreeNode::compLess);
	} else if (subOrder == SUBPROB_WIDTH_DEC) {
		sort(m_children.begin(), m_children.end(), PseudotreeNode::compGreater);
	}
}

/* compares two pseudotree nodes, returns true if subtree below a is more complex
 * than below b (looks at width, then height) */
bool PseudotreeNode::compGreater(PseudotreeNode* a, PseudotreeNode* b) {
	assert(a && b);
	if (a->getSubWidth() > b->getSubWidth())
		return true;
	if (a->getSubWidth() < b->getSubWidth())
		return false;
	if (a->getSubHeight() > b->getSubHeight())
		return true;
	return false;
}

bool PseudotreeNode::compLess(PseudotreeNode* a, PseudotreeNode* b) {
	assert(a && b);
	if (a->getSubWidth() < b->getSubWidth())
		return true;
	if (a->getSubWidth() > b->getSubWidth())
		return false;
	if (a->getSubHeight() < b->getSubHeight())
		return true;
	return false;
}

/* updates a single node's depth and height, recursively updating the child nodes.
 * return value is height of node's subtree */
int PseudotreeNode::updateDepthHeight(int d) {

	m_depth = d;

	if (m_children.empty()) {
		m_subHeight = 0;
	} else {
		int m = 0;
		for (std::vector<PseudotreeNode*>::iterator it = m_children.begin();
				it != m_children.end(); ++it)
			m = max(m, (*it)->updateDepthHeight(m_depth + 1));
		m_subHeight = m + 1;
	}

	return m_subHeight;

}

/* recursively finds and returns the max. width in this node's subproblem */
int PseudotreeNode::updateSubWidth() {

	m_subWidth = m_context.size();

	for (std::vector<PseudotreeNode*>::iterator it = m_children.begin();
			it != m_children.end(); ++it)
		m_subWidth = max(m_subWidth, (*it)->updateSubWidth());

	return m_subWidth;

}

/* recursively updates the set of variables in the current subproblem */
const set<int>& PseudotreeNode::updateSubprobVars(int numVars) {

	// clear current set
	m_subproblemVars.clear();
	// add self
	m_subproblemVars.insert(m_var);

	// iterate over children and collect their subproblem variables
	for (std::vector<PseudotreeNode*>::iterator it = m_children.begin();
			it != m_children.end(); ++it) {
		const std::set<int>& childVars = (*it)->updateSubprobVars(numVars);
		for (std::set<int>::const_iterator itC = childVars.begin();
				itC != childVars.end(); ++itC)
			m_subproblemVars.insert(*itC);
	}

	m_subproblemVarMap.clear();
	m_subproblemVarMap.resize(numVars, NONE);
	size_t i = 0;
	for (std::set<int>::const_iterator it = m_subproblemVars.begin();
			it != m_subproblemVars.end(); ++it, ++i) {
		m_subproblemVarMap[*it] = i;
	}

#ifdef DEBUG
//  cout << "Subproblem at var. " << m_var << ": " << m_subproblemVars << endl;
#endif

	// return a const reference
	return m_subproblemVars;
}

void Pseudotree::insertNewNode(const int i, const std::set<int>& N,
		std::list<PseudotreeNode*>& roots) {
	// create new node in pseudo tree
	PseudotreeNode* p = new PseudotreeNode(this, i, N);
	m_nodes[i] = p;

	// incorporate new pseudo tree node
	std::list<PseudotreeNode*>::iterator it = roots.begin();
	while (it != roots.end()) {
		const set<int>& context = (*it)->getFullContext();
		if (context.find(i) != context.end()) {
			p->addChild(*it); // add child to current node
			(*it)->setParent(p); // set parent of previous node
			it = roots.erase(it); // remove previous node from roots list
		} else {
			++it;
		}
	}
	roots.push_back(p); // add current node to roots list

}

void Pseudotree::outputToFileNode(const PseudotreeNode* node,
		std::ostringstream& oss) const {
	oss << "(" << node->getVar();
	for (std::vector<PseudotreeNode*>::const_iterator it =
			node->getChildren().begin(); it != node->getChildren().end(); ++it)
		outputToFileNode(*it, oss);
	oss << ")";
}

void Pseudotree::outputToFile(std::string of_name) const {
	assert(m_root);
	std::ostringstream oss;
	outputToFileNode(m_root, oss);  // recursive, call on root node

	std::ofstream of;
	of.open(of_name.c_str(), std::ios_base::out);
	of << oss.str() << endl;
	of.close();
}

// select top M variables (start pseudo tree)
// -- need to make sure that the elimination order is still valid
void Pseudotree::bfs(std::vector<int>& start) {

	std::cout << " - Selecting start pseudo tree:" << std::endl;
	std::cout << " - (dummy) root depth: " << m_root->getDepth() << std::endl;

	start.clear();

	std::deque<PseudotreeNode*> bfs;
	for (std::vector<PseudotreeNode*>::const_iterator it = m_root->getChildren().begin();
			it != m_root->getChildren().end(); ++it) {
		bfs.push_back( *it );
	}

	while (!bfs.empty()) {
		PseudotreeNode* n = bfs.front();
		bfs.pop_front();
		start.push_back(n->getVar());
		for (std::vector<PseudotreeNode*>::const_iterator it =
				n->getChildren().begin(); it != n->getChildren().end(); ++it) {
			bfs.push_back(*it);
		}
	}

}

// Return a dfs order of the nodes in the pseudo tree. If mapOnly is true, then
// return a dfs order of the MAP variables only.
void Pseudotree::dfs(std::list<int>& order, const bool mapOnly) {

	order.clear();

	std::stack<PseudotreeNode*> dfs;
	dfs.push(m_root);
	while (!dfs.empty()) {
		PseudotreeNode* n = dfs.top();
		dfs.pop();
		int var = n->getVar();
		if (mapOnly) {
			if (m_problem->isMap(var)) {
				order.push_back(n->getVar());
			}
		} else {
			order.push_back(var);
		}

		for (std::vector<PseudotreeNode*>::const_iterator it =
				n->getChildren().begin(); it != n->getChildren().end(); ++it) {
			dfs.push(*it);
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
void Pseudotree::outputToDot(const char *dotfile ) {

	std::ofstream outfile( dotfile );
	if (outfile.fail()) {
		std::cerr << "Cannot open Graphviz file.";
		exit(EXIT_FAILURE);
	}

	outfile << "digraph g {\n";
	outfile << "node [shape = record];\n";
	outfile << "size = \"10, 7.5\";\n";
	outfile << "rotate = \"90\";\n";
	outfile << "ratio = \"fill\";\n";
	drawNodesForDot( outfile );
	drawLinesForDot( outfile );
	outfile << "}\n";
}

void Pseudotree::drawNodesForDot( std::ofstream& outfile ) {

	std::queue<PseudotreeNode*> q;

	q.push(m_root);
	while (!q.empty()) {
		PseudotreeNode* node = q.front();
		q.pop();

		if (m_problem->isSum(node->getVar())) {
			outfile << "node" << node->getVar()
				<< "[ shape=ellipse, label = \"" << node->getVar() << "\"];\n";
		} else {
			outfile << "node" << node->getVar()
				<< "[ shape=box, color=gold, label = \"" << node->getVar() << "\"];\n";
		}

		const std::vector<PseudotreeNode*>& ch = node->getChildren();
		for (size_t i = 0; i < ch.size(); ++i) {
			PseudotreeNode* child = ch[i];
			q.push(child);
		}
	}
}

void Pseudotree::drawLinesForDot( std::ofstream& outfile ) {
	std::queue<PseudotreeNode*> q;

	q.push(m_root);
	while (!q.empty()) {
		PseudotreeNode* node = q.front();
		q.pop();

		const std::vector<PseudotreeNode*>& ch = node->getChildren();
		for (size_t i = 0; i < ch.size(); ++i) {

			PseudotreeNode* child = ch[i];

			outfile << "node" << node->getVar()
				<< " -> node" << child->getVar() << ";\n";

			q.push(child);
		}
	}
}


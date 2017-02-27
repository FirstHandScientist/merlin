/*
 * join_tree.cpp
 *
 *  Created on: Sep 12, 2013
 *      Author: radu
 */

#include "join_tree.h"
#include "pseudo_tree.h"


#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

using namespace boost;
typedef adjacency_list<vecS, vecS, undirectedS, no_property,
		property<edge_weight_t, int> > BoostGraph;
typedef graph_traits<BoostGraph>::edge_descriptor BoostEdge;
typedef graph_traits<BoostGraph>::vertex_descriptor BoostVertex;
typedef std::pair<int, int> E;

//#define DEBUG

void PathEdge::updateMessage(const std::vector<int>& assignment) {
	assert(edge != NULL);
	assert(activeMessage == 1 || activeMessage == 2);
	if (activeMessage == 1) {
		edge->updateMessage1to2(assignment);
	} else {
		edge->updateMessage2to1(assignment);
	}
}

void PathEdge::updateMessage(int ibound, const std::vector<int>& assignment) {
	assert(edge != NULL);
	assert(activeMessage == 1 || activeMessage == 2);
	if (activeMessage == 1) {
		edge->updateMessage1to2(ibound, assignment);
	} else {
		edge->updateMessage2to1(ibound, assignment);
	}
}

// check if A dominates B (set B is included in set A)
bool JTree::dominates(const std::set<int>& A, const std::set<int>& B) {
	std::set<int> temp;
	std::set_intersection(A.begin(), A.end(),
			B.begin(), B.end(),
			std::inserter(temp, temp.begin()));

	return (temp == B);
}

// Kruskal stuff
void JTree::buildJunctionTree(bool mc) {

	// safety checks
	assert( m_root != NULL );

	int n = (int)m_clusters.size();
	std::vector<std::vector<int> > adjacencies;
	adjacencies.resize(n);
	for (std::vector<std::vector<int> >::iterator it = adjacencies.begin();
			it != adjacencies.end(); ++it) {
		it->resize(n);
	}

	// fill in the weights
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i != j) {
				JTreeNode* I = m_id2cluster[i];
				JTreeNode* J = m_id2cluster[j];

				std::set<int> sep;
				std::set_intersection(I->vars().begin(), I->vars().end(),
						J->vars().begin(), J->vars().end(),
						std::inserter(sep, sep.begin()));

				if (sep.empty()) {
					adjacencies[i][j] = (100000);
				} else {
					adjacencies[i][j] = -(sep.size());
				}
			}
		}
	}

	// create a boost weighted graph and run Kruskal on it
	const int num_nodes = n;
	std::vector<E> temp1;
	std::vector<int> temp2;
	for (int i = 0; i < num_nodes - 1; ++i) {
		for (int j = i + 1; j < num_nodes; ++j) {
			if (i != j) {
				temp1.push_back( E(i, j) );
				temp2.push_back( adjacencies[i][j] );
			}
		}
	}

	E* edge_array = new E[temp1.size()];
	int* weights = new int[temp1.size()];
	for (size_t i = 0; i < temp1.size(); ++i) {
		edge_array[i] = temp1[i];
		weights[i] = temp2[i];
	}

	std::size_t num_edges = temp1.size();

	BoostGraph g(edge_array, edge_array + num_edges, weights, num_nodes);
	std::vector<BoostEdge> spanning_tree;

	kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

//#ifdef PRINT_JUNCTION_TREE
//	property_map<BoostGraph, edge_weight_t>::type weight = get(edge_weight, g);
//	std::cout << "Undirected junction tree (MST):" << std::endl;
//	for (std::vector<BoostEdge>::iterator ei = spanning_tree.begin();
//			ei != spanning_tree.end(); ++ei) {
//		std::cout << source(*ei, g) << " <--> " << target(*ei, g)
//				<< " with weight of " << weight[*ei] << "\n";
//	}
//
//#endif

	// select the root and redirect edges outwards
	int root = m_root->id();
	std::set<int> visited;
	std::stack<int> dfs;
	dfs.push(root);
	while ( !dfs.empty() ) {

		int n = dfs.top();
		dfs.pop();

		// get all unvisited children and direct them towards n
		for (std::vector<BoostEdge>::iterator ei = spanning_tree.begin();
				ei != spanning_tree.end(); ++ei) {
			int src = source(*ei, g);
			int trg = target(*ei, g);
			JTreeNode* from = NULL, *to = NULL;

			if (src == n &&
				std::find(visited.begin(), visited.end(), trg) == visited.end()) {

				from = m_id2cluster[trg]; // from
				to = m_id2cluster[src]; // to
				dfs.push(trg);
			}

			if (trg == n &&
				std::find(visited.begin(), visited.end(), src) == visited.end()) {

				from = m_id2cluster[src]; // from
				to = m_id2cluster[trg]; // to
				dfs.push(src);
			}

			// create the directed edge
			if (from != NULL && to != NULL) {
				std::set<int> sep;
				std::set_intersection(from->vars().begin(), from->vars().end(),
						to->vars().begin(), to->vars().end(),
						std::inserter(sep, sep.begin()));
				JTreeEdge* e = new JTreeEdge(from, to, sep);
				e->setProblem(m_problem, m_options, mc);
				m_edges.push_back(e);
				from->edges().push_back(e);
				to->edges().push_back(e);
				to->addChild(from);
				from->setParent(to);
			}

		}

		visited.insert(n);
	}

	// find the message order (schedule) by a bfs of the junction tree
	std::queue<JTreeNode*> bfs;
	bfs.push(m_root);
	while (!bfs.empty()) {
		JTreeNode* n = bfs.front();
		bfs.pop();

		JTreeNode* p = n->parent();
		if (p != NULL) {
			std::vector<JTreeEdge*>& elist = p->edges();
			for (std::vector<JTreeEdge*>::iterator it = elist.begin();
					it != elist.end(); ++it) {
				JTreeEdge* e = (*it);
				if (e->getNode1() == n && e->getNode2() == p) {
					m_messageOrder.push_back(e);
					break;
				}
			}
		}

		for (std::vector<JTreeNode*>::iterator it = n->children().begin();
				it != n->children().end(); ++it) {
			bfs.push(*it);
		}
	}

	assert(m_messageOrder.size() == m_edges.size());

	// reverse the order to start from the leaves
	std::reverse(m_messageOrder.begin(), m_messageOrder.end());

	// clean up
	delete[] weights;
	delete[] edge_array;

#ifdef PRINT_JUNCTION_TREE
	std::ofstream fout("/home/radu/tmp/kruskal.dot");
	fout << "graph JT {\n" << " rankdir=LR\n" << " size=\"3,3\"\n"
			<< " ratio=\"filled\"\n" << " edge[style=\"bold\"]\n"
			<< " node[shape=\"circle\"]\n";
	// nodes
	for (size_t i = 0; i < m_clusters.size(); ++i) {
		fout << "node" << m_clusters[i]->id()
			<< "[ label = \"" << i << "\"];\n";
	}

	// edges
	for (size_t i = 0; i < m_edges.size(); ++i) {
		JTreeEdge* e = m_edges[i];
		fout << "node" << e->getNode1()->id() << " -> " << "node" << e->getNode2()->id() << "\n";
	}

	fout << "}\n";
#endif
}


// build the join tree (Judea Pearl)
// - assume the dummy variable has been added to the model
// - assume the dummy functions were added to the model as well
// - input graph G should be connected (and contain dummy functions)
// - input elim order has the dummy variable at the end
void JTree::build(const Graph& G, const std::vector<int>& elim, const bool mc) {

	// safety checks
	assert( m_problem->hasDummy() );
	assert( (int)elim.size() == m_problem->getN() );

	// check if GG is disconnected
	std::set<int> vars;
	Graph GGG(G);
	for (int i = 0; i < m_problem->getN(); ++i) vars.insert(i);
	std::map<int, std::set<int> > comps1 = GGG.connectedComponents(vars);
	std::cout << "JT: Input graph connected components: " << comps1.size();
	std::cout << std::endl;

	// triangulate the graph
	Graph GG(G);
	GG.triangulate(elim);

	std::map<int, std::set<int> > comps2 = GG.connectedComponents(vars);
	std::cout << "JT: Triangulated graph connected components: " << comps2.size();
	std::cout << std::endl;

	// get the elimination cliques
	int lastVar = elim.back();
	std::vector<std::set<int> > elimCliques;
	for (size_t i = 0; i < elim.size(); ++i) {
		int var = elim[i];
		const set<int>& N = GG.getNeighbors(var);
		m_treewidth = std::max(m_treewidth, (int) N.size());

		// form the clique
		std::set<int> cand;
		cand.insert(var);
		cand.insert(N.begin(), N.end());
		elimCliques.push_back(cand);

		GG.removeNode(var);
	}

	std::cout << "  Found " << elimCliques.size() << " elimination cliques" << std::endl;

	// remove dominated cliques
	// : A dominates B ( A > B ) iff scope(A) includes scope(B)
	std::vector<std::set<int> > clusters;
	for (std::vector<set<int> >::iterator it1 = elimCliques.begin();
			it1 != elimCliques.end(); ++it1) {
		std::set<int>& A = *it1;

		bool found = false;
		for (std::vector<std::set<int> >::iterator it2 = clusters.begin();
				it2 != clusters.end(); ++it2) {

			std::set<int>& B = *it2;

			if ( dominates(B, A) ) {
				found = true;
				break;
			}
		}

		if (!found) {
			for (std::vector<std::set<int> >::iterator it3 = clusters.begin();
					it3 != clusters.end();) {

				std::set<int>& B = (*it3);
				if ( dominates(A, B) ) {
					clusters.erase(it3);
				} else {
					++it3;
				}
			}

			clusters.push_back(A);
		}
	}


	std::cout << "  Found " << clusters.size() << " clusters (maximal cliques))" << std::endl;
#ifdef PRINT_JUNCTION_TREE
	for (size_t i = 0; i < clusters.size(); ++i) {
		std::cout << "( ";
		std::copy(clusters[i].begin(), clusters[i].end(),
				std::ostream_iterator<int>(std::cout, " "));
		std::cout << ")" << std::endl;
	}
#endif

	// create the nodes of the join tree
	int id = 0;
	m_var2cluster.resize(elim.size(), NULL);
	for (std::vector<std::set<int> >::iterator it = clusters.begin();
			it != clusters.end(); ++it) {

		JTreeNode* cluster = new JTreeNode(id++, *it);
		cluster->setProblem(m_problem, m_options);
		m_id2cluster[ cluster->id() ] = cluster;
		for (std::set<int>::iterator si = it->begin();
				si != it->end(); ++si) {
			int var = *si;
			if (m_var2cluster[var] == NULL) {
				m_var2cluster[var] = cluster;
			}

			if (m_problem->isMap(var)) {
				cluster->addMaxVar(var);
			} else {
				cluster->addSumVar(var);
			}
		}

		m_clusters.push_back(cluster);

		// check if the cluster contains the last variable
		if (cluster->contains(lastVar)) {
			m_root = cluster;
		}
	}

	// root shouldn't be null
	assert(m_root != NULL);

	// put the functions in the appropriate node (cluster)
	double globalConstant = 1.0;
	const std::vector<Function*>& functions = m_problem->getFunctions();
	int counter = (int)functions.size();
	std::vector<Function*>::const_iterator fi = functions.begin();
	for (; fi != functions.end(); ++fi) {

		Function* f = (*fi);
		if ( f->getScope().empty() ) {
			globalConstant *= f->getTable()[0];
			--counter;
			continue; // add constants to the last node
		}

		// find a cluster node that contains the function scope
		const std::set<int>& scope = f->getScope();
		for (size_t j = 0; j < m_clusters.size(); ++j) {
			const std::set<int>& clique = m_clusters[j]->vars();
			std::vector<int> temp;
			if (dominates(clique, scope)) {
				if (mc) {
					m_clusters[j]->add(f->asFactor());
				} else {
					m_clusters[j]->factor() *= f->asFactor();
				}

				--counter;
				break;
			}
		}
	}

	// safety checks: all functions allocated
	assert (counter == 0);

	// create the junction tree (using Kruscal on weighted graph)
	buildJunctionTree(mc);

	// collect the MAP variables (including dummy which is assumed MAX)
	m_varsMax.clear();
	for (int var = 0; var < m_problem->getN(); ++var) {
		if ( m_problem->getVartypes()[var] == VAR_MAX ) {
			m_varsMax.insert( var );
		}
	}

	// max cluster size
	for (size_t i = 0; i < m_clusters.size(); ++i) {
		m_clusterSize = std::max(m_clusterSize, m_clusters[i]->vars().size());
	}

	// max separator size
	for (size_t i = 0; i < m_edges.size(); ++i) {
		m_separatorSize = std::max(m_separatorSize, m_edges[i]->separator().size());
	}

	// report some stats
	std::cout << "  Join tree created (" << m_treewidth << "/" << m_clusters.size() << ")"<< std::endl;
	std::cout << "    Mini-Clustering: " << (mc ? "yes" : "no") << std::endl;
	std::cout << "    Global constant: " << globalConstant << std::endl;
	std::cout << "    Max cluster size:   " << m_clusterSize << std::endl;
	std::cout << "    Max separator size: " << m_separatorSize << std::endl;
	std::cout << "  Unconstrained elimination order: ";
	std::copy(elim.begin(), elim.end(), std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;

	// Find the static MAP search order
	// From now on, we're pre-processing the bucket tree for incremental
	// upper bounds computation by:
	// - fixing the MAP search order: post-order of the bucket tree
	// - precomputing the paths between succesive clusters containing MAP vars
	// - mark edges that need to be cached (for efficient backtracking)

	// create the MAP search order (static): post order traversal of the tree
	findPostOrder();

	// get the paths between succesive clusters
	for (int i = 0; i < (int)m_searchClusters.size()-1; ++i) {

		JTreeNode* from = m_searchClusters[i];
		JTreeNode* to = m_searchClusters[i+1];
		std::vector<PathEdge*> path;
		findPath(from, to, path);
		m_paths[ from->id() ] = path;
	}

	// Ouput the join tree
#ifdef PRINT_JUNCTION_TREE

	std::cout << "DEBUG Join tree created (" << m_treewidth << "/" << m_clusters.size() << ")"<< std::endl;
	std::cout << "DEBUG JTNodes:" << std::endl;
	for (size_t i = 0; i < m_clusters.size(); ++i) {
		std::cout << m_clusters[i]->toString() << std::endl;
		std::cout << "   : " << m_clusters[i]->factor() << std::endl;
	}

	std::cout << "DEBUG JTEdges:" << std::endl;
	for (size_t i = 0; i < m_edges.size(); ++i) {
		std::cout << "  Edge " << i << std::endl;
		std::cout << "   .from(1): " << m_edges[i]->getNode1()->id() << std::endl;
		std::cout << "   .to(2)  : " << m_edges[i]->getNode2()->id() << std::endl;
		std::cout << "   .sep    : " << "[ ";
		std::copy(m_edges[i]->separator().begin(),
				m_edges[i]->separator().end(),
				std::ostream_iterator<int>(std::cout, " "));
		std::cout << "]" << std::endl;
	}

	std::cout << "DEBUG  Propagation schedule:" << std::endl;
	for (size_t i = 0; i < m_messageOrder.size(); ++i) {
		JTreeEdge* edge = m_messageOrder[i];
		std::cout << " " << edge->getNode1()->id() << " --> "
				<< edge->getNode2()->id() << std::endl;
	}

#endif

}

// find the post-order for searching the MAP variables
// - ignore the dummy variable
void JTree::findPostOrder() {
	m_searchOrder.clear();
	m_searchClusters.clear();
	if (m_varsMax.empty())
		return; // nothing to do

	std::set<int> visited;
	std::set<int> mapVars;
	std::stack<JTreeNode*> dfs;
	dfs.push(m_root);
	while ( !dfs.empty() ) {
		JTreeNode* cluster = dfs.top();
		std::vector<JTreeNode*>& children = cluster->children();
		bool flag = false;
		for (std::vector<JTreeNode*>::iterator it = children.begin();
				it != children.end(); ++it) {
			JTreeNode* child = *it;
			int id = child->id();
			if (visited.find(id) == visited.end()) {
				dfs.push(child);
				flag = true;
			}
		}

		if (flag) continue;

		// update the search order and search clusters
		bool marked = false;
		for (std::set<int>::const_iterator it = cluster->vars().begin();
				it != cluster->vars().end(); ++it) {
			int var = *it;
			if (!m_problem->isMap(var)) continue; // ignore SUM variables

			// check if current cluster still has new MAP variables
			if (mapVars.find(var) == mapVars.end()) {
				m_searchOrder.push_back(var);
				mapVars.insert(var);
				marked = true;
			}
		}

		if (marked) {
			m_searchClusters.push_back(cluster);
		}

		visited.insert(cluster->id());
		dfs.pop();
	}

	// remove dominated clusters (wrt MAP variables)
	std::vector<JTreeNode*> temp;
	for (std::vector<JTreeNode*>::iterator it1 = m_searchClusters.begin();
			it1 != m_searchClusters.end(); ++it1) {
		JTreeNode* u = (*it1);
		const std::set<int>& A = u->varsMax();

		bool found = false;
		for (std::vector<JTreeNode*>::iterator it2 = temp.begin();
				it2 != temp.end(); ++it2) {

			JTreeNode* v = (*it2);
			const std::set<int>& B = v->varsMax();

			if (dominates(B, A)) {
				found = true;
				break;
			}
		}

		if (!found) {
			for (std::vector<JTreeNode*>::iterator it3 = temp.begin();
					it3 != temp.end();) {

				JTreeNode* v = (*it3);
				const std::set<int>& B = v->varsMax();

				if (dominates(A, B)) {
					temp.erase(it3);
				} else {
					++it3;
				}
			}

			temp.push_back(u);
		}
	}

	m_searchClusters = temp;

	// set the search cluster corresponding to each of the MAP variables in order
	for (std::vector<JTreeNode*>::iterator it = m_searchClusters.begin();
			it != m_searchClusters.end(); ++it) {
		JTreeNode* cluster = (*it);
		for (std::set<int>::const_iterator si = cluster->varsMax().begin();
				si != cluster->varsMax().end(); ++si) {
			int var = (*si);
			if (m_map2cluster.find(var) == m_map2cluster.end()) {
				m_map2cluster[var] = cluster;
			}
		}
	}

	// print the MAP search order
	std::cout << "Static MAP search order (" << m_searchOrder.size() << "): ";
	std::copy(m_searchOrder.begin(), m_searchOrder.end(),
			std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;
	std::cout << "Static MAP cluster order (" << m_searchClusters.size() << "): ";
	for (std::vector<JTreeNode*>::iterator it = m_searchClusters.begin();
			it != m_searchClusters.end(); ++it) {
		std::cout << " C" << (*it)->id();
		std::cout << ":( ";
		std::copy((*it)->varsMax().begin(), (*it)->varsMax().end(),
				std::ostream_iterator<int>(std::cout, " "));
		std::cout << ") ";
	}
	std::cout << std::endl;
#ifdef PRINT_JUNCTION_TREE
	std::cout << "Assignment MAP var to JT cluster:" << std::endl;
	for (std::list<int>::iterator li = m_searchOrder.begin();
			li != m_searchOrder.end(); ++li) {
		std::cout << "  " << *li << "\t-- C" << m_map2cluster[*li]->id()
				<< " " << m_map2cluster[*li]->toString() << std::endl;
	}
#endif

}

// find a path between two nodes in the join tree
void JTree::findPath(JTreeNode* from, JTreeNode* to, std::vector<PathEdge*>& path) {

	// safety checks
	assert(from != NULL);
	assert(to != NULL);
	path.clear();

	// check if same nodes
	if (from == to)
		return; // nothing to do

	// find the LCA of the two nodes
	std::list<int> path1, path2;
	path1.push_front(from->id());
	JTreeNode* n = from->parent();
	while (n != NULL) {
		path1.push_front(n->id());
		n = n->parent();
	}

	path2.push_front(to->id());
	n = to->parent();
	while (n != NULL) {
		path2.push_front(n->id());
		n = n->parent();
	}

	// find the last element that's common to the two lists
	int lca = NONE;
	std::list<int>::iterator it1 = path1.begin();
	std::list<int>::iterator it2 = path2.begin();
	while (it1 != path1.end() && it2 != path2.end()) {
		if (*it1 == *it2) {
			lca = *it1;
			++it1;
			++it2;
		} else {
			break;
		}
	}

	// assemble the path (the actual message schedule)
	assert(lca != NONE);
	JTreeNode* c = from;
	while ( c->id() != lca ) {
		JTreeNode* p = c->parent();
		for (std::vector<JTreeEdge*>::iterator it = c->edges().begin();
				it != c->edges().end(); ++it) {

			JTreeEdge* edge = (*it);
			if (edge->getNode1()->id() == c->id() &&
				edge->getNode2()->id() == p->id()) {

				PathEdge* e = new PathEdge(c->id(), p->id());
				e->edge = edge;
				e->activeMessage = 1; // send 1-2 message;
				path.push_back(e);
				break;
			}
		}

		c = p;
		p = c->parent();
	}

	c = to;
	std::list<PathEdge*> temp;
	while (c->id() != lca) {
		JTreeNode* p = c->parent();
		for (std::vector<JTreeEdge*>::iterator it = c->edges().begin();
				it != c->edges().end(); ++it) {

			JTreeEdge* edge = (*it);
			if (edge->getNode1()->id() == c->id() &&
				edge->getNode2()->id() == p->id()) {

				PathEdge* e = new PathEdge(p->id(), c->id());
				e->edge = edge;
				e->activeMessage = 2; // send 2-1 message
				temp.push_front(e);
				break;
			}
		}

		c = p;
		p = c->parent();
	}

	// update the path
	for (std::list<PathEdge*>::iterator it = temp.begin();
			it != temp.end(); ++it) {
		path.push_back(*it);
	}
}

void JTree::clean() {
	for (std::vector<JTreeEdge*>::iterator it = m_edges.begin();
			it != m_edges.end(); ++it) {
		JTreeEdge* e = (*it);
		e->clean();
	}
}

// propagate messages
void JTree::propagate(const std::vector<int>& assignment) {

	// clean previous messages
	clean();

	// initialize the nodes
	for (size_t i = 0; i < m_clusters.size(); ++i) {
		m_clusters[i]->init();
	}

	// collect phase
#ifdef PRINT_JUNCTION_TREE
	std::cout << "---- Collect messages (leaves-root) ----" << std::endl;
#endif
	for (std::vector<JTreeEdge*>::iterator it = m_messageOrder.begin();
			it != m_messageOrder.end(); ++it) {
		JTreeEdge* edge = *it;
#ifdef PRINT_JUNCTION_TREE
		std::cout << " msg from C" << edge->getNode1()->id()
				<< " to C" << edge->getNode2()->id() << std::endl;
#endif
		edge->sendMessage1to2(assignment);
	}

	// distribute phase
#ifdef PRINT_JUNCTION_TREE
	std::cout << "---- Distribute messages (root-leaves) ----" << std::endl;
#endif
	for (std::vector<JTreeEdge*>::reverse_iterator it = m_messageOrder.rbegin();
			it != m_messageOrder.rend(); ++it) {
		JTreeEdge* edge = *it;
#ifdef PRINT_JUNCTION_TREE
		std::cout << " msg from C" << edge->getNode2()->id()
				<< " to C" << edge->getNode1()->id() << std::endl;
#endif
		edge->sendMessage2to1(assignment);
	}

	// compute the marginals of the unassigned MAP variables (if any)
	m_maxmarginals.clear();
	for (std::set<int>::iterator si = m_varsMax.begin();
			si != m_varsMax.end(); ++si) {
		int var = *si;
		if (assignment[var] == NONE) {
			JTreeNode* cl = m_var2cluster[var];
			m_maxmarginals[var] = cl->marginal(var, assignment);
		}
	}

	// get the (partial) MAP upper bound
	if (!m_maxmarginals.empty()) {
		mex::Factor& F = (*m_maxmarginals.begin()).second;
		m_upperBound = F.max();
	} else { // full MAP assignment
		m_upperBound = m_root->upperBound(assignment);
	}
}

// propagate messages (mini-clustering)
void JTree::propagate(int ibound, const std::vector<int>& assignment) {

	// clean previous messages
	clean();

	// initialize the nodes
	for (size_t i = 0; i < m_clusters.size(); ++i) {
		m_clusters[i]->init();
	}

	// collect phase
#ifdef PRINT_JUNCTION_TREE
	std::cout << "---- Collect messages (leaves-root) ----" << std::endl;
#endif
	for (std::vector<JTreeEdge*>::iterator it = m_messageOrder.begin();
			it != m_messageOrder.end(); ++it) {
		JTreeEdge* edge = *it;
		edge->sendMessage1to2(ibound, assignment);

#ifdef PRINT_JUNCTION_TREE
		std::cout << " msg from C" << edge->getNode1()->id()
				<< " to C" << edge->getNode2()->id() << std::endl;
		std::vector<mex::Factor>& vmsg = edge->getVectorMessage1();
		for (std::vector<mex::Factor>::iterator mi = vmsg.begin(); mi != vmsg.end(); ++mi) {
			mex::Factor& ff = (*mi);
			std::cout << ff << "; ";
		}
		std::cout << std::endl;
#endif
	}

	// distribute phase
#ifdef PRINT_JUNCTION_TREE
	std::cout << "---- Distribute messages (root-leaves) ----" << std::endl;
#endif
	for (std::vector<JTreeEdge*>::reverse_iterator it = m_messageOrder.rbegin();
			it != m_messageOrder.rend(); ++it) {
		JTreeEdge* edge = *it;
		edge->sendMessage2to1(ibound, assignment);

#ifdef PRINT_JUNCTION_TREE
		std::cout << " msg from C" << edge->getNode2()->id()
				<< " to C" << edge->getNode1()->id() << std::endl;
		std::vector<mex::Factor>& vmsg = edge->getVectorMessage1();
		for (std::vector<mex::Factor>::iterator mi = vmsg.begin(); mi != vmsg.end(); ++mi) {
			mex::Factor& ff = (*mi);
			std::cout << ff << "; ";
		}
		std::cout << std::endl;
#endif

	}

	// compute the marginals of the unassigned MAP variables (if any)
	m_maxmarginals.clear();
	for (std::set<int>::iterator si = m_varsMax.begin();
			si != m_varsMax.end(); ++si) {
		int var = *si;
		if (assignment[var] == NONE) {
			JTreeNode* cl = m_var2cluster[var];
			m_maxmarginals[var] = cl->marginal(var, ibound, assignment);
		}
	}

	// get the (partial) MAP upper bound
	if (!m_maxmarginals.empty()) {
		mex::Factor& F = (*m_maxmarginals.begin()).second;
		m_upperBound = F.max();
	} else { // full MAP assignment
		m_upperBound = m_root->upperBound(ibound, assignment);
	}
}


// get an overall upper bound on the full MAP assignment
double JTreeNode::upperBound(const std::vector<int>& assignment) {

	// collect all original functions of this node
	mex::Factor F = this->factor();

	// project the current assignment
	mex::VarSet vs = F.vars();
	for (mex::VarSet::const_iterator it = vs.begin(); it != vs.end(); ++it) {
		if (assignment[it->label()] != NONE) {
			mex::VarSet vc(*it);
			F = F.condition(vc, assignment[it->label()]);
		}
	}

	// collect all incoming messages to this node
	for (std::vector<JTreeEdge*>::iterator it = m_edges.begin();
			it != m_edges.end(); ++it) {
		JTreeEdge* currEdge = *it;
		if (currEdge->getNode2()->id() == this->id()) { // incoming edge
			F *= currEdge->getMessage1();
		} else { // outgoing edge
			F *= currEdge->getMessage2();
		}
	}

	// current scope of F should only contain SUM vars
	return F.sum();
}

typedef std::pair<size_t, size_t> fscope;
typedef std::vector<fscope> mini_cluster;

// check if a mini-cluster allows a new function
bool allows_factor(mini_cluster& mb,
		fscope& f,
		int ibound,
		const std::vector<mex::Factor>& factors) {

	mex::VarSet jointScope;
	for (mini_cluster::iterator it = mb.begin(); it != mb.end(); ++it) {
		fscope& tmp = (*it);
		jointScope |= factors[tmp.first].vars();
	}

	jointScope |= factors[f.first].vars();
	if ((int)jointScope.size() <= ibound) return true;
	else return false;
}

// assumes that all input factors have been projected on the current assignment
size_t make_partition(int ibound, const std::vector<mex::Factor>& factors,
		std::vector<mex::Factor>& partition) {

	// initialize the scopes of the factors
	std::vector<fscope> scopes;
	for (size_t i = 0; i < factors.size(); ++i) {
		scopes.push_back( std::make_pair(i, factors[i].vars().size()) );
	}

	// sort the scopes in decreasing order of the scopes
	std::sort(scopes.begin(), scopes.end(), scopeIsLargerSet);

	// create the partition
	std::vector<mini_cluster> part;
	for (std::vector<fscope>::iterator itS = scopes.begin();
			itS != scopes.end(); ++itS) {

		fscope& s = (*itS);
		bool placed = false;
		for (std::vector<mini_cluster>::iterator itB = part.begin();
				!placed && itB != part.end(); ++itB) {

			mini_cluster& mb = (*itB);
			if (allows_factor(mb, s, ibound, factors)) { // checks if function fits into bucket
				mb.push_back(s);
				placed = true;
			}
		}

		if (!placed) { // no fit, need to create new bucket
			mini_cluster mb;
			mb.push_back(s);
			part.push_back(mb);
		}
	}

	// create the join factors corresponding to the partitioning
	for (std::vector<mini_cluster>::iterator itB = part.begin();
			itB != part.end(); ++itB) {

		mini_cluster& mb = (*itB);
		mex::Factor F;
		for (mini_cluster::iterator it = mb.begin(); it != mb.end(); ++it) {
			fscope& f = (*it);
			F *= factors[f.first];
		}

		partition.push_back(F);
	}

	return partition.size();
}

// get an overall upper bound on the full MAP assignment
double JTreeNode::upperBound(int ibound, const std::vector<int>& assignment) {

	// collect all original functions of this node
	std::vector<mex::Factor> factors;
	std::vector<mex::Factor>& FF = this->factors();

	// project the current assignment
	for (std::vector<mex::Factor>::iterator fi = FF.begin();
			fi != FF.end(); ++fi) {

		mex::Factor F = *fi;
		mex::VarSet vs = F.vars();
		for (mex::VarSet::const_iterator it = vs.begin(); it != vs.end(); ++it) {
			if (assignment[it->label()] != NONE) {
				mex::VarSet vc(*it);
				F = F.condition(vc, assignment[it->label()]);
			}
		}

		factors.push_back(F);
	}

	// collect all incoming messages to this node
	for (std::vector<JTreeEdge*>::iterator it = m_edges.begin();
			it != m_edges.end(); ++it) {
		JTreeEdge* currEdge = *it;
		if (currEdge->getNode2()->id() == this->id()) { // incoming edge
			std::copy(currEdge->getVectorMessage1().begin(),
					currEdge->getVectorMessage1().end(),
					std::back_inserter(factors));
		} else { // outgoing edge
			std::copy(currEdge->getVectorMessage2().begin(),
					currEdge->getVectorMessage2().end(),
					std::back_inserter(factors));
		}
	}

	// partition the factors into mini-clusters
	std::vector<mex::Factor> partition;
	size_t num = make_partition(ibound, factors, partition);
	double upbo = 1;
	for (size_t i = 0; i < num; ++i) {
		if (i == 0) upbo *= partition[i].sum();
		else upbo *= partition[i].max();
	}

	return upbo;
}

// get an overall upper bound on the full MAP assignment
double JTreeNode::upperBoundIncr(const std::vector<int>& assignment) {

	// safety check
	for (size_t i = 0; i < assignment.size(); ++i) {
		if (m_problem->isMap(i)) {
			assert(assignment[i] != NONE);
		}
	}

	// collect all original functions of this node
	mex::Factor F = this->factor();

	// project the current assignment
	mex::VarSet vs = F.vars();
	for (mex::VarSet::const_iterator it = vs.begin(); it != vs.end(); ++it) {
		if (assignment[it->label()] != NONE) {
			mex::VarSet vc(*it);
			F = F.condition(vc, assignment[it->label()]);
		}
	}

	// collect all incoming messages to this node
	for (std::vector<JTreeEdge*>::iterator it = m_edges.begin();
			it != m_edges.end(); ++it) {
		JTreeEdge* currEdge = *it;
		if (currEdge->getNode2()->id() == this->id()) { // incoming edge

			mex::Factor tmp = currEdge->getMessage1();

			// project the current assignment
			mex::VarSet vars = tmp.vars();
			for (mex::VarSet::const_iterator it = vars.begin(); it != vars.end(); ++it) {
				if (assignment[it->label()] != NONE) {
					mex::VarSet vc(*it);
					tmp = tmp.condition(vc, assignment[it->label()]);
				}
			}

			F *= tmp;

		} else { // outgoing edge

			mex::Factor tmp = currEdge->getMessage2();

			// project the current assignment
			mex::VarSet vars = tmp.vars();
			for (mex::VarSet::const_iterator it = vars.begin(); it != vars.end(); ++it) {
				if (assignment[it->label()] != NONE) {
					mex::VarSet vc(*it);
					tmp = tmp.condition(vc, assignment[it->label()]);
				}
			}

			F *= tmp;
		}
	}

	// current scope of F should only contain SUM vars
	return F.sum();
}

// get an overall upper bound on the full MAP assignment
double JTreeNode::upperBoundIncr(int ibound, const std::vector<int>& assignment) {

	// safety check
	for (size_t i = 0; i < assignment.size(); ++i) {
		if (m_problem->isMap(i)) {
			assert(assignment[i] != NONE);
		}
	}

	// collect all original functions of this node
	std::vector<mex::Factor> factors;
	std::vector<mex::Factor>& FF = this->factors();

	// project the current assignment
	for (std::vector<mex::Factor>::iterator fi = FF.begin();
			fi != FF.end(); ++fi) {

		mex::Factor F = *fi;
		mex::VarSet vs = F.vars();
		for (mex::VarSet::const_iterator it = vs.begin(); it != vs.end(); ++it) {
			if (assignment[it->label()] != NONE) {
				mex::VarSet vc(*it);
				F = F.condition(vc, assignment[it->label()]);
			}
		}

		factors.push_back(F);
	}

	// collect all incoming messages to this node
	for (std::vector<JTreeEdge*>::iterator it = m_edges.begin();
			it != m_edges.end(); ++it) {
		JTreeEdge* currEdge = *it;
		if (currEdge->getNode2()->id() == this->id()) { // incoming edge

			std::vector<mex::Factor>& vtmp = currEdge->getVectorMessage1();
			for (size_t i = 0; i < vtmp.size(); ++i) {
				mex::Factor tmp = vtmp[i];

				// project the current assignment
				mex::VarSet vars = tmp.vars();
				for (mex::VarSet::const_iterator it = vars.begin(); it != vars.end(); ++it) {
					if (assignment[it->label()] != NONE) {
						mex::VarSet vc(*it);
						tmp = tmp.condition(vc, assignment[it->label()]);
					}
				}

				factors.push_back(tmp);
			}
		} else { // outgoing edge

			std::vector<mex::Factor>& vtmp = currEdge->getVectorMessage2();
			for (size_t i = 0; i < vtmp.size(); ++i) {
				mex::Factor tmp = vtmp[i];

				// project the current assignment
				mex::VarSet vars = tmp.vars();
				for (mex::VarSet::const_iterator it = vars.begin(); it != vars.end(); ++it) {
					if (assignment[it->label()] != NONE) {
						mex::VarSet vc(*it);
						tmp = tmp.condition(vc, assignment[it->label()]);
					}
				}

				factors.push_back(tmp);
			}
		}
	}

	// partition the factors into mini-clusters
	std::vector<mex::Factor> partition;
	size_t num = make_partition(ibound, factors, partition);
	double upbo = 1;
	for (size_t i = 0; i < num; ++i) {
		if (i == 0) upbo *= partition[i].sum();
		else upbo *= partition[i].max();
	}

	return upbo;
}

// compute the marginal of a variable
mex::Factor JTreeNode::marginal(int var, const std::vector<int>& assignment) {

	// safety checks
	assert( m_vars.find(var) != m_vars.end() );
	assert( assignment[var] == NONE );

	// collect all original functions of this node
	mex::Factor F = this->factor();

	// project the current assignment
	mex::VarSet vs = F.vars();
	for (mex::VarSet::const_iterator it = vs.begin(); it != vs.end(); ++it) {
		if (assignment[it->label()] != NONE) {
			mex::VarSet vc(*it);
			F = F.condition(vc, assignment[it->label()]);
		}
	}

	// collect all incoming messages to this node (should already be projected)
	for (std::vector<JTreeEdge*>::iterator it = m_edges.begin();
			it != m_edges.end(); ++it) {
		JTreeEdge* currEdge = *it;
		if (currEdge->getNode2()->id() == this->id()) { // incoming edge

			F *= currEdge->getMessage1();

			// safety checks
			mex::VarSet vars = currEdge->getMessage1().vars();
			for (mex::VarSet::const_iterator it = vars.begin();
					it != vars.end(); ++it) {
				const mex::Var& v = (*it);
				assert( assignment[v.label()] == NONE );
			}
		} else { // outgoing edge
			assert( currEdge->getNode1()->id() == this->id() );

			F *= currEdge->getMessage2();

			// safety checks
			mex::VarSet vars = currEdge->getMessage2().vars();
			for (mex::VarSet::const_iterator it = vars.begin();
					it != vars.end(); ++it) {
				const mex::Var& v = (*it);
				assert( assignment[v.label()] == NONE );
			}
		}
	}

	// safety check
	mex::VarSet vars = F.vars();
	for (mex::VarSet::const_iterator it = vars.begin();
			it != vars.end(); ++it) {
		const mex::Var& v = (*it);
		assert( assignment[v.label()] == NONE );
	}


	mex::VarSet sumOut;
	for (std::set<int>::const_iterator it = m_varsSum.begin();
			it != m_varsSum.end(); ++it) {
		int x = (*it);
		if (assignment[x] != NONE) continue;
		if ( var != x ) {
			sumOut.insert(mex::Var(x, m_problem->getDomainSize(x)));
		}
	}

	if (sumOut.size() > 0) F = F.sum(sumOut);

	mex::VarSet maxOut;
	for (std::set<int>::const_iterator it = m_varsMax.begin();
			it != m_varsMax.end(); ++it) {
		int x = (*it);
		if (assignment[x] != NONE) continue;
		if ( var != x ) {
			maxOut.insert(mex::Var(x, m_problem->getDomainSize(x)));
		}
	}

	if (maxOut.size() > 0) F = F.max(maxOut);

	return F;
}

// compute the marginal of a variable
mex::Factor JTreeNode::marginal(int var, int ibound, const std::vector<int>& assignment) {

	// safety checks
	assert( m_vars.find(var) != m_vars.end() );
	assert( assignment[var] == NONE );

	// collect all original functions of this node
	std::vector<mex::Factor> factors;
	std::vector<mex::Factor>& FF = this->factors();

	// project the current assignment
	for (std::vector<mex::Factor>::iterator fi = FF.begin();
			fi != FF.end(); ++fi) {

		mex::Factor F = *fi;
		mex::VarSet vs = F.vars();
		for (mex::VarSet::const_iterator it = vs.begin(); it != vs.end(); ++it) {
			if (assignment[it->label()] != NONE) {
				mex::VarSet vc(*it);
				F = F.condition(vc, assignment[it->label()]);
			}
		}

		factors.push_back(F);
	}

	// collect all incoming messages to this node (should already be projected)
	for (std::vector<JTreeEdge*>::iterator it = m_edges.begin();
			it != m_edges.end(); ++it) {
		JTreeEdge* currEdge = *it;
		if (currEdge->getNode2()->id() == this->id()) { // incoming edge

			std::vector<mex::Factor>& vtmp = currEdge->getVectorMessage1();
			for (size_t i = 0; i < vtmp.size(); ++i) {
				factors.push_back(vtmp[i]);
			}
		} else { // outgoing edge
			assert( currEdge->getNode1()->id() == this->id() );

			std::vector<mex::Factor>& vtmp = currEdge->getVectorMessage2();
			for (size_t i = 0; i < vtmp.size(); ++i) {
				factors.push_back(vtmp[i]);
			}
		}
	}

	// partition the factors into mini-clusters
	std::vector<mex::Factor> partition;
	size_t num = make_partition(ibound, factors, partition);
	mex::Factor marg;
	for (size_t i = 0; i < num; ++i) {
		mex::Factor& mb = partition[i];
		mex::Var x(var, m_problem->getDomainSize(var));
		if (i == 0) {
			mex::VarSet sumOut, maxOut;
			for (mex::VarSet::const_iterator it = mb.vars().begin();
					it != mb.vars().end(); ++it) {
				if (m_problem->isMap( it->label() )) maxOut += (*it);
				else sumOut += (*it);
			}
			maxOut -= x;

			if (sumOut.size() > 0) mb = mb.sum(sumOut);
			if (maxOut.size() > 0) mb = mb.max(maxOut);

			marg *= mb;
		} else {
			mex::VarSet maxOut = mb.vars();
			maxOut -= x;

			marg *= mb.max(maxOut);
		}
	}

	return marg;
}

// compute the marginal of a variable (incremental update)
mex::Factor JTreeNode::marginalIncr(int var, const std::vector<int>& assignment) {

	// safety checks
	assert( m_vars.find(var) != m_vars.end() );
	assert( assignment[var] == NONE );

	// collect all original functions of this node
	mex::Factor F = this->factor();

	// project the current assignment
	mex::VarSet vs = F.vars();
	for (mex::VarSet::const_iterator it = vs.begin(); it != vs.end(); ++it) {
		if (assignment[it->label()] != NONE) {
			mex::VarSet vc(*it);
			F = F.condition(vc, assignment[it->label()]);
		}
	}

	// collect all incoming messages to this node (should already be projected)
	for (std::vector<JTreeEdge*>::iterator it = m_edges.begin();
			it != m_edges.end(); ++it) {
		JTreeEdge* currEdge = *it;
		if (currEdge->getNode2()->id() == this->id()) { // incoming edge

			mex::Factor tmp = currEdge->getMessage1();

			// project the current assignment
			mex::VarSet vars = tmp.vars();
			for (mex::VarSet::const_iterator it = vars.begin(); it != vars.end(); ++it) {
				if (assignment[it->label()] != NONE) {
					mex::VarSet vc(*it);
					tmp = tmp.condition(vc, assignment[it->label()]);
				}
			}

			F *= tmp;

		} else { // outgoing edge
			assert( currEdge->getNode1()->id() == this->id() );

			mex::Factor tmp = currEdge->getMessage2();
			// project the current assignment
			mex::VarSet vars = tmp.vars();
			for (mex::VarSet::const_iterator it = vars.begin(); it != vars.end(); ++it) {
				if (assignment[it->label()] != NONE) {
					mex::VarSet vc(*it);
					tmp = tmp.condition(vc, assignment[it->label()]);
				}
			}

			F *= tmp;
		}
	}

	// safety check
	mex::VarSet vars = F.vars();
	for (mex::VarSet::const_iterator it = vars.begin();
			it != vars.end(); ++it) {
		const mex::Var& v = (*it);
		assert( assignment[v.label()] == NONE );
	}


	mex::VarSet sumOut;
	for (std::set<int>::const_iterator it = m_varsSum.begin();
			it != m_varsSum.end(); ++it) {
		int x = (*it);
		if (assignment[x] != NONE) continue;
		if ( var != x ) {
			sumOut.insert(mex::Var(x, m_problem->getDomainSize(x)));
		}
	}

	if (sumOut.size() > 0) F = F.sum(sumOut);

	mex::VarSet maxOut;
	for (std::set<int>::const_iterator it = m_varsMax.begin();
			it != m_varsMax.end(); ++it) {
		int x = (*it);
		if (assignment[x] != NONE) continue;
		if ( var != x ) {
			maxOut.insert(mex::Var(x, m_problem->getDomainSize(x)));
		}
	}

	if (maxOut.size() > 0) F = F.max(maxOut);

	return F;
}

// compute the marginal of a variable (incremental update)
mex::Factor JTreeNode::marginalIncr(int var, int ibound, const std::vector<int>& assignment) {

	// safety checks
	assert( m_vars.find(var) != m_vars.end() );
	assert( assignment[var] == NONE );

	// collect all original functions of this node
	std::vector<mex::Factor> factors;
	std::vector<mex::Factor>& FF = this->factors();

	// project the current assignment
	for (std::vector<mex::Factor>::iterator fi = FF.begin();
			fi != FF.end(); ++fi) {

		mex::Factor F = *fi;
		mex::VarSet vs = F.vars();
		for (mex::VarSet::const_iterator it = vs.begin(); it != vs.end(); ++it) {
			if (assignment[it->label()] != NONE) {
				mex::VarSet vc(*it);
				F = F.condition(vc, assignment[it->label()]);
			}
		}

		factors.push_back(F);
	}

	// collect all incoming messages to this node (should already be projected)
	for (std::vector<JTreeEdge*>::iterator it = m_edges.begin();
			it != m_edges.end(); ++it) {
		JTreeEdge* currEdge = *it;
		if (currEdge->getNode2()->id() == this->id()) { // incoming edge

			std::vector<mex::Factor>& vtmp = currEdge->getVectorMessage1();
			for (size_t i = 0; i < vtmp.size(); ++i) {
				mex::Factor tmp = vtmp[i];

				// project the current assignment
				mex::VarSet vars = tmp.vars();
				for (mex::VarSet::const_iterator it = vars.begin(); it != vars.end(); ++it) {
					if (assignment[it->label()] != NONE) {
						mex::VarSet vc(*it);
						tmp = tmp.condition(vc, assignment[it->label()]);
					}
				}

				factors.push_back(tmp);
			}

		} else { // outgoing edge
			assert( currEdge->getNode1()->id() == this->id() );

			std::vector<mex::Factor>& vtmp = currEdge->getVectorMessage2();
			for (size_t i = 0; i < vtmp.size(); ++i) {
				mex::Factor tmp = vtmp[i];

				// project the current assignment
				mex::VarSet vars = tmp.vars();
				for (mex::VarSet::const_iterator it = vars.begin(); it != vars.end(); ++it) {
					if (assignment[it->label()] != NONE) {
						mex::VarSet vc(*it);
						tmp = tmp.condition(vc, assignment[it->label()]);
					}
				}

				factors.push_back(tmp);
			}
		}
	}


	// partition the factors into mini-clusters
	std::vector<mex::Factor> partition;
	size_t num = make_partition(ibound, factors, partition);
	mex::Factor marg;
	for (size_t i = 0; i < num; ++i) {
		mex::Factor& mb = partition[i];
		mex::Var x(var, m_problem->getDomainSize(var));
		if (i == 0) {
			mex::VarSet sumOut, maxOut;
			for (mex::VarSet::const_iterator it = mb.vars().begin();
					it != mb.vars().end(); ++it) {
				if (m_problem->isMap( it->label() )) maxOut += (*it);
				else sumOut += (*it);
			}
			maxOut -= x;

			if (sumOut.size() > 0) mb = mb.sum(sumOut);
			if (maxOut.size() > 0) mb = mb.max(maxOut);

			marg *= mb;
		} else {
			mex::VarSet maxOut = mb.vars();
			maxOut -= x;

			marg *= mb.max(maxOut);
		}
	}

	return marg;
}

// send message 1 to 2 (used for full message passing)
void JTreeEdge::sendMessage1to2(const std::vector<int>& assignment) {

	// collect the original functions
	mex::Factor F = m_node1->factor();

	// project the current assignment
	mex::VarSet vs = F.vars();
	for (mex::VarSet::const_iterator it = vs.begin(); it != vs.end(); ++it) {
		if (assignment[it->label()] != NONE) {
			mex::VarSet vc(*it);
			F = F.condition(vc, assignment[it->label()]);
		}
	}

	// collect the incoming messages (should already be projected)
	for (size_t i = 0; i < m_node1->edges().size(); ++i) {
		JTreeEdge* currEdge = m_node1->edges()[i];
		if (currEdge == this) {
			continue; // skip current edge
		}
		if (currEdge->getNode1()->id() == m_node1->id()) {

			F *= currEdge->getMessage2();

			// safety checks
			mex::VarSet vars = currEdge->getMessage2().vars();
			for (mex::VarSet::const_iterator it = vars.begin();
					it != vars.end(); ++it) {
				const mex::Var& v = (*it);
				assert( assignment[v.label()] == NONE );
			}
		} else if (currEdge->getNode2()->id() == m_node1->id()) {

			F *= currEdge->getMessage1();

			// safety checks
			mex::VarSet vars = currEdge->getMessage1().vars();
			for (mex::VarSet::const_iterator it = vars.begin();
					it != vars.end(); ++it) {
				const mex::Var& v = (*it);
				assert( assignment[v.label()] == NONE );
			}
		}
	}

	// safety checks
	mex::VarSet vars = F.vars();
	for (mex::VarSet::const_iterator it = vars.begin();
			it != vars.end(); ++it) {
		const mex::Var& v = (*it);
		assert( assignment[v.label()] == NONE );
	}


	// marginalize the eliminator (SUM first, then MAX)
	mex::VarSet sumOut, maxOut;
	for (std::set<int>::const_iterator it = m_node1->varsSum().begin();
			it != m_node1->varsSum().end(); ++it) {
		int x = (*it);
		if (assignment[x] != NONE) continue;
		if (m_separator.find( x ) == m_separator.end()) {
			sumOut.insert(mex::Var(x, m_problem->getDomainSize(x)));
		}
	}

	if (sumOut.size() > 0) F = F.sum(sumOut);

	for (std::set<int>::const_iterator it = m_node1->varsMax().begin();
			it != m_node1->varsMax().end(); ++it) {
		int x = (*it);
		if (assignment[x] != NONE) continue;
		if (m_separator.find( x ) == m_separator.end()) {
			maxOut.insert(mex::Var(x, m_problem->getDomainSize(x)));
		}
	}

	if (maxOut.size() > 0) F = F.max(maxOut);
	m_message1to2 = F;
}

// send message 1 to 2 (used for full message passing)
void JTreeEdge::sendMessage1to2(int ibound, const std::vector<int>& assignment) {

	// collect all original functions
	std::vector<mex::Factor> factors;
	std::vector<mex::Factor>& FF = m_node1->factors();

	// project the current assignment
	for (std::vector<mex::Factor>::iterator fi = FF.begin();
			fi != FF.end(); ++fi) {

		mex::Factor F = *fi;
		mex::VarSet vs = F.vars();
		for (mex::VarSet::const_iterator it = vs.begin(); it != vs.end(); ++it) {
			if (assignment[it->label()] != NONE) {
				mex::VarSet vc(*it);
				F = F.condition(vc, assignment[it->label()]);
			}
		}

		factors.push_back(F);
	}

	// collect the incoming messages (should already be projected)
	for (size_t i = 0; i < m_node1->edges().size(); ++i) {
		JTreeEdge* currEdge = m_node1->edges()[i];
		if (currEdge == this) {
			continue; // skip current edge
		}
		if (currEdge->getNode1()->id() == m_node1->id()) {

			std::vector<mex::Factor>& vtmp = currEdge->getVectorMessage2();
			for (size_t i = 0; i < vtmp.size(); ++i) {
				factors.push_back( vtmp[i] );
			}

		} else if (currEdge->getNode2()->id() == m_node1->id()) {

			std::vector<mex::Factor>& vtmp = currEdge->getVectorMessage1();
			for (size_t i = 0; i < vtmp.size(); ++i) {
				factors.push_back( vtmp[i] );
			}

		}
	}


	// partition the factors into mini-clusters
	std::vector<mex::Factor> partition;
	size_t num = make_partition(ibound, factors, partition);
	std::vector<mex::Factor> vmsg;
	for (size_t i = 0; i < num; ++i) {
		mex::Factor& mb = partition[i];

		if (i == 0) {
			mex::VarSet sumOut, maxOut;
			for (mex::VarSet::const_iterator it = mb.vars().begin();
					it != mb.vars().end(); ++it) {

				int x = it->label();
				if (m_separator.find(x) == m_separator.end()) {
					if (m_problem->isMap(x))
						maxOut += mex::Var(x, m_problem->getDomainSize(x));
					else
						sumOut += mex::Var(x, m_problem->getDomainSize(x));
				}
			}

			if (sumOut.size() > 0) mb = mb.sum(sumOut);
			if (maxOut.size() > 0) mb = mb.max(maxOut);

			vmsg.push_back(mb);

		} else {
			mex::VarSet maxOut;
			for (mex::VarSet::const_iterator it = mb.vars().begin();
					it != mb.vars().end(); ++it) {

				int x = it->label();
				if (m_separator.find(x) == m_separator.end()) {
					maxOut += mex::Var(x, m_problem->getDomainSize(x));
				}
			}

			if (maxOut.size() > 0) mb = mb.max(maxOut);

			vmsg.push_back(mb);
		}
	}

	m_mcMessage1to2 = vmsg;
}

// send message 2 to 1 (used for full message passing)
void JTreeEdge::sendMessage2to1(const std::vector<int>& assignment) {

	// collect the original functions
	mex::Factor F = m_node2->factor();

	// project the current assignment
	mex::VarSet vs = F.vars();
	for (mex::VarSet::const_iterator it = vs.begin(); it != vs.end(); ++it) {
		if (assignment[it->label()] != NONE) {
			mex::VarSet cond(*it);
			F = F.condition(cond, assignment[it->label()]);
		}
	}

	// collect the incoming messages
	for (size_t i = 0; i < m_node2->edges().size(); ++i) {
		JTreeEdge* currEdge = m_node2->edges()[i];
		if (currEdge == this) {
			continue; // skip current edge
		}
		if (currEdge->getNode1()->id() == m_node2->id()) {

			F *= currEdge->getMessage2();

			// safety checks
			mex::VarSet vars = currEdge->getMessage2().vars();
			for (mex::VarSet::const_iterator it = vars.begin();
					it != vars.end(); ++it) {
				const mex::Var& v = (*it);
				assert( assignment[v.label()] == NONE );
			}
		} else if (currEdge->getNode2()->id() == m_node2->id()) {

			F *= currEdge->getMessage1();

			// safety checks
			mex::VarSet vars = currEdge->getMessage1().vars();
			for (mex::VarSet::const_iterator it = vars.begin();
					it != vars.end(); ++it) {
				const mex::Var& v = (*it);
				assert( assignment[v.label()] == NONE );
			}
		}
	}

	// safety checks
	mex::VarSet vars = F.vars();
	for (mex::VarSet::const_iterator it = vars.begin();
			it != vars.end(); ++it) {
		const mex::Var& v = (*it);
		assert( assignment[v.label()] == NONE );
	}

	// marginalize the eliminator (SUM first, then MAX)
	mex::VarSet sumOut, maxOut;
	for (std::set<int>::const_iterator it = m_node2->varsSum().begin();
			it != m_node2->varsSum().end(); ++it) {
		int x = (*it);
		if (assignment[x] != NONE) continue;
		if (m_separator.find( x ) == m_separator.end()) {
			sumOut.insert(mex::Var(x, m_problem->getDomainSize(x)));
		}
	}

	if (sumOut.size() > 0) F = F.sum(sumOut);

	for (std::set<int>::const_iterator it = m_node2->varsMax().begin();
			it != m_node2->varsMax().end(); ++it) {
		int x = (*it);
		if (assignment[x] != NONE) continue;
		if (m_separator.find( x ) == m_separator.end()) {
			maxOut.insert(mex::Var(x, m_problem->getDomainSize(x)));
		}
	}

	if (maxOut.size() > 0) F = F.max(maxOut);
	m_message2to1 = F;
}

// send message 2 to 1 (used for full message passing)
void JTreeEdge::sendMessage2to1(int ibound, const std::vector<int>& assignment) {

	// collect all original functions
	std::vector<mex::Factor> factors;
	std::vector<mex::Factor>& FF = m_node2->factors();

	// project the current assignment
	for (std::vector<mex::Factor>::iterator fi = FF.begin();
			fi != FF.end(); ++fi) {

		mex::Factor F = *fi;
		mex::VarSet vs = F.vars();
		for (mex::VarSet::const_iterator it = vs.begin(); it != vs.end(); ++it) {
			if (assignment[it->label()] != NONE) {
				mex::VarSet vc(*it);
				F = F.condition(vc, assignment[it->label()]);
			}
		}

		factors.push_back(F);
	}

	// collect the incoming messages (should already be projected on the assignment)
	for (size_t i = 0; i < m_node2->edges().size(); ++i) {
		JTreeEdge* currEdge = m_node2->edges()[i];
		if (currEdge == this) {
			continue; // skip current edge
		}
		if (currEdge->getNode1()->id() == m_node2->id()) {

			std::vector<mex::Factor>& vtmp = currEdge->getVectorMessage2();
			for (size_t i = 0; i < vtmp.size(); ++i) {
				factors.push_back( vtmp[i] );
			}

		} else if (currEdge->getNode2()->id() == m_node2->id()) {

			std::vector<mex::Factor>& vtmp = currEdge->getVectorMessage1();
			for (size_t i = 0; i < vtmp.size(); ++i) {
				factors.push_back( vtmp[i] );
			}
		}
	}

	// partition the factors into mini-clusters
	std::vector<mex::Factor> partition;
	size_t num = make_partition(ibound, factors, partition);
	std::vector<mex::Factor> vmsg;
	for (size_t i = 0; i < num; ++i) {
		mex::Factor& mb = partition[i];

		if (i == 0) {
			mex::VarSet sumOut, maxOut;
			for (mex::VarSet::const_iterator it = mb.vars().begin();
					it != mb.vars().end(); ++it) {

				int x = it->label();
				if (m_separator.find(x) == m_separator.end()) {
					if (m_problem->isMap(x))
						maxOut += mex::Var(x, m_problem->getDomainSize(x));
					else
						sumOut += mex::Var(x, m_problem->getDomainSize(x));
				}
			}

			if (sumOut.size() > 0) mb = mb.sum(sumOut);
			if (maxOut.size() > 0) mb = mb.max(maxOut);

			vmsg.push_back(mb);

		} else {
			mex::VarSet maxOut;
			for (mex::VarSet::const_iterator it = mb.vars().begin();
					it != mb.vars().end(); ++it) {

				int x = it->label();
				if (m_separator.find(x) == m_separator.end()) {
					maxOut += mex::Var(x, m_problem->getDomainSize(x));
				}
			}

			if (maxOut.size() > 0) mb = mb.max(maxOut);

			vmsg.push_back(mb);
		}
	}

	m_mcMessage2to1 = vmsg;
}

// send message 1 to 2 (used for incremental message passing)
void JTreeEdge::updateMessage1to2(const std::vector<int>& assignment) {

	// collect the original functions
	mex::Factor F = m_node1->factor();

	// project the current assignment
	mex::VarSet vs = F.vars();
	for (mex::VarSet::const_iterator it = vs.begin(); it != vs.end(); ++it) {
		if (assignment[it->label()] != NONE) {
			mex::VarSet vc(*it);
			F = F.condition(vc, assignment[it->label()]);
		}
	}

	// collect the incoming messages and project the assignment
	for (size_t i = 0; i < m_node1->edges().size(); ++i) {
		JTreeEdge* currEdge = m_node1->edges()[i];
		if (currEdge == this) {
			continue; // skip current edge
		}
		if (currEdge->getNode1()->id() == m_node1->id()) {

			mex::Factor tmp = currEdge->getMessage2();

			// project the current assignment
			mex::VarSet vars = tmp.vars();
			for (mex::VarSet::const_iterator it = vars.begin(); it != vars.end(); ++it) {
				if (assignment[it->label()] != NONE) {
					mex::VarSet vc(*it);
					tmp = tmp.condition(vc, assignment[it->label()]);
				}
			}

			F *= tmp;

		} else if (currEdge->getNode2()->id() == m_node1->id()) {

			mex::Factor tmp = currEdge->getMessage1();

			// project the current assignment
			mex::VarSet vars = tmp.vars();
			for (mex::VarSet::const_iterator it = vars.begin(); it != vars.end(); ++it) {
				if (assignment[it->label()] != NONE) {
					mex::VarSet vc(*it);
					tmp = tmp.condition(vc, assignment[it->label()]);
				}
			}

			F *= tmp;

		}
	}

	// safety checks
	mex::VarSet vars = F.vars();
	for (mex::VarSet::const_iterator it = vars.begin();
			it != vars.end(); ++it) {
		const mex::Var& v = (*it);
		assert( assignment[v.label()] == NONE );
	}


	// marginalize the unassigned eliminator (SUM first, then MAX)
	mex::VarSet sumOut, maxOut;
	for (std::set<int>::const_iterator it = m_node1->varsSum().begin();
			it != m_node1->varsSum().end(); ++it) {
		int x = (*it);
		if (assignment[x] != NONE) continue;
		if (m_separator.find( x ) == m_separator.end()) {
			sumOut.insert(mex::Var(x, m_problem->getDomainSize(x)));
		}
	}

	if (sumOut.size() > 0) F = F.sum(sumOut);

	for (std::set<int>::const_iterator it = m_node1->varsMax().begin();
			it != m_node1->varsMax().end(); ++it) {
		int x = (*it);
		if (assignment[x] != NONE) continue;
		if (m_separator.find( x ) == m_separator.end()) {
			maxOut.insert(mex::Var(x, m_problem->getDomainSize(x)));
		}
	}

	if (maxOut.size() > 0) F = F.max(maxOut);
	m_message1to2 = F;
}

// send message 1 to 2 (used for incremental message passing)
void JTreeEdge::updateMessage1to2(int ibound, const std::vector<int>& assignment) {

	// collect all original functions
	std::vector<mex::Factor> factors;
	std::vector<mex::Factor>& FF = m_node1->factors();

	// project the current assignment
	for (std::vector<mex::Factor>::iterator fi = FF.begin();
			fi != FF.end(); ++fi) {

		mex::Factor F = *fi;
		mex::VarSet vs = F.vars();
		for (mex::VarSet::const_iterator it = vs.begin(); it != vs.end(); ++it) {
			if (assignment[it->label()] != NONE) {
				mex::VarSet vc(*it);
				F = F.condition(vc, assignment[it->label()]);
			}
		}

		factors.push_back(F);
	}

	// collect the incoming messages and project the assignment
	for (size_t i = 0; i < m_node1->edges().size(); ++i) {
		JTreeEdge* currEdge = m_node1->edges()[i];
		if (currEdge == this) {
			continue; // skip current edge
		}
		if (currEdge->getNode1()->id() == m_node1->id()) {

			std::vector<mex::Factor>& vtmp = currEdge->getVectorMessage2();
			for (size_t i = 0; i < vtmp.size(); ++i) {
				mex::Factor tmp = vtmp[i];

				// project the current assignment
				mex::VarSet vars = tmp.vars();
				for (mex::VarSet::const_iterator it = vars.begin(); it != vars.end(); ++it) {
					if (assignment[it->label()] != NONE) {
						mex::VarSet vc(*it);
						tmp = tmp.condition(vc, assignment[it->label()]);
					}
				}

				factors.push_back(tmp);
			}

		} else if (currEdge->getNode2()->id() == m_node1->id()) {

			std::vector<mex::Factor>& vtmp = currEdge->getVectorMessage1();
			for (size_t i = 0; i < vtmp.size(); ++i) {
				mex::Factor tmp = vtmp[i];

				// project the current assignment
				mex::VarSet vars = tmp.vars();
				for (mex::VarSet::const_iterator it = vars.begin(); it != vars.end(); ++it) {
					if (assignment[it->label()] != NONE) {
						mex::VarSet vc(*it);
						tmp = tmp.condition(vc, assignment[it->label()]);
					}
				}

				factors.push_back(tmp);
			}
		}
	}

	// partition the factors into mini-clusters
	std::vector<mex::Factor> partition;
	size_t num = make_partition(ibound, factors, partition);
	std::vector<mex::Factor> vmsg;
	for (size_t i = 0; i < num; ++i) {
		mex::Factor& mb = partition[i];

		if (i == 0) {
			mex::VarSet sumOut, maxOut;
			for (mex::VarSet::const_iterator it = mb.vars().begin();
					it != mb.vars().end(); ++it) {

				int x = it->label();
				if (m_separator.find(x) == m_separator.end()) {
					if (m_problem->isMap(x))
						maxOut += mex::Var(x, m_problem->getDomainSize(x));
					else
						sumOut += mex::Var(x, m_problem->getDomainSize(x));
				}
			}

			if (sumOut.size() > 0) mb = mb.sum(sumOut);
			if (maxOut.size() > 0) mb = mb.max(maxOut);

			vmsg.push_back(mb);

		} else {
			mex::VarSet maxOut;
			for (mex::VarSet::const_iterator it = mb.vars().begin();
					it != mb.vars().end(); ++it) {

				int x = it->label();
				if (m_separator.find(x) == m_separator.end()) {
					maxOut += mex::Var(x, m_problem->getDomainSize(x));
				}
			}

			if (maxOut.size() > 0) mb = mb.max(maxOut);

			vmsg.push_back(mb);
		}
	}

	m_mcMessage1to2 = vmsg;
}


// send message 2 to 1 (used for incremental message passing)
void JTreeEdge::updateMessage2to1(const std::vector<int>& assignment) {

	// collect the original functions
	mex::Factor F = m_node2->factor();

	// project the current assignment
	mex::VarSet vs = F.vars();
	for (mex::VarSet::const_iterator it = vs.begin(); it != vs.end(); ++it) {
		if (assignment[it->label()] != NONE) {
			mex::VarSet cond(*it);
			F = F.condition(cond, assignment[it->label()]);
		}
	}

	// collect the incoming messages
	for (size_t i = 0; i < m_node2->edges().size(); ++i) {
		JTreeEdge* currEdge = m_node2->edges()[i];
		if (currEdge == this) {
			continue; // skip current edge
		}
		if (currEdge->getNode1()->id() == m_node2->id()) {

			mex::Factor tmp = currEdge->getMessage2();

			// project the current assignment
			mex::VarSet vars = tmp.vars();
			for (mex::VarSet::const_iterator it = vars.begin(); it != vars.end(); ++it) {
				if (assignment[it->label()] != NONE) {
					mex::VarSet vc(*it);
					tmp = tmp.condition(vc, assignment[it->label()]);
				}
			}

			F *= tmp;

		} else if (currEdge->getNode2()->id() == m_node2->id()) {

			mex::Factor tmp = currEdge->getMessage1();

			// project the current assignment
			mex::VarSet vars = tmp.vars();
			for (mex::VarSet::const_iterator it = vars.begin(); it != vars.end(); ++it) {
				if (assignment[it->label()] != NONE) {
					mex::VarSet vc(*it);
					tmp = tmp.condition(vc, assignment[it->label()]);
				}
			}

			F *= tmp;

		}
	}

	// safety checks
	mex::VarSet vars = F.vars();
	for (mex::VarSet::const_iterator it = vars.begin();
			it != vars.end(); ++it) {
		const mex::Var& v = (*it);
		assert( assignment[v.label()] == NONE );
	}

	// marginalize the eliminator (SUM first, then MAX)
	mex::VarSet sumOut, maxOut;
	for (std::set<int>::const_iterator it = m_node2->varsSum().begin();
			it != m_node2->varsSum().end(); ++it) {
		int x = (*it);
		if (assignment[x] != NONE) continue;
		if (m_separator.find( x ) == m_separator.end()) {
			sumOut.insert(mex::Var(x, m_problem->getDomainSize(x)));
		}
	}

	if (sumOut.size() > 0) F = F.sum(sumOut);

	for (std::set<int>::const_iterator it = m_node2->varsMax().begin();
			it != m_node2->varsMax().end(); ++it) {
		int x = (*it);
		if (assignment[x] != NONE) continue;
		if (m_separator.find( x ) == m_separator.end()) {
			maxOut.insert(mex::Var(x, m_problem->getDomainSize(x)));
		}
	}

	if (maxOut.size() > 0) F = F.max(maxOut);
	m_message2to1 = F;
}

// send message 2 to 1 (used for incremental message passing)
void JTreeEdge::updateMessage2to1(int ibound, const std::vector<int>& assignment) {

	// collect all original functions
	std::vector<mex::Factor> factors;
	std::vector<mex::Factor>& FF = m_node2->factors();

	// project the current assignment
	for (std::vector<mex::Factor>::iterator fi = FF.begin();
			fi != FF.end(); ++fi) {

		mex::Factor F = *fi;
		mex::VarSet vs = F.vars();
		for (mex::VarSet::const_iterator it = vs.begin(); it != vs.end(); ++it) {
			if (assignment[it->label()] != NONE) {
				mex::VarSet vc(*it);
				F = F.condition(vc, assignment[it->label()]);
			}
		}

		factors.push_back(F);
	}

	// collect the incoming messages
	for (size_t i = 0; i < m_node2->edges().size(); ++i) {
		JTreeEdge* currEdge = m_node2->edges()[i];
		if (currEdge == this) {
			continue; // skip current edge
		}
		if (currEdge->getNode1()->id() == m_node2->id()) {

			std::vector<mex::Factor>& vtmp = currEdge->getVectorMessage2();
			for (size_t i = 0; i < vtmp.size(); ++i) {
				mex::Factor tmp = vtmp[i];

				// project the current assignment
				mex::VarSet vars = tmp.vars();
				for (mex::VarSet::const_iterator it = vars.begin(); it != vars.end(); ++it) {
					if (assignment[it->label()] != NONE) {
						mex::VarSet vc(*it);
						tmp = tmp.condition(vc, assignment[it->label()]);
					}
				}

				factors.push_back(tmp);
			}

		} else if (currEdge->getNode2()->id() == m_node2->id()) {

			std::vector<mex::Factor>& vtmp = currEdge->getVectorMessage1();
			for (size_t i = 0; i < vtmp.size(); ++i) {
				mex::Factor tmp = vtmp[i];

				// project the current assignment
				mex::VarSet vars = tmp.vars();
				for (mex::VarSet::const_iterator it = vars.begin(); it != vars.end(); ++it) {
					if (assignment[it->label()] != NONE) {
						mex::VarSet vc(*it);
						tmp = tmp.condition(vc, assignment[it->label()]);
					}
				}

				factors.push_back(tmp);
			}

		}
	}

	// partition the factors into mini-clusters
	std::vector<mex::Factor> partition;
	size_t num = make_partition(ibound, factors, partition);
	std::vector<mex::Factor> vmsg;
	for (size_t i = 0; i < num; ++i) {
		mex::Factor& mb = partition[i];

		if (i == 0) {
			mex::VarSet sumOut, maxOut;
			for (mex::VarSet::const_iterator it = mb.vars().begin();
					it != mb.vars().end(); ++it) {

				int x = it->label();
				if (m_separator.find(x) == m_separator.end()) {
					if (m_problem->isMap(x))
						maxOut += mex::Var(x, m_problem->getDomainSize(x));
					else
						sumOut += mex::Var(x, m_problem->getDomainSize(x));
				}
			}

			if (sumOut.size() > 0) mb = mb.sum(sumOut);
			if (maxOut.size() > 0) mb = mb.max(maxOut);

			vmsg.push_back(mb);

		} else {
			mex::VarSet maxOut;
			for (mex::VarSet::const_iterator it = mb.vars().begin();
					it != mb.vars().end(); ++it) {

				int x = it->label();
				if (m_separator.find(x) == m_separator.end()) {
					maxOut += mex::Var(x, m_problem->getDomainSize(x));
				}
			}

			if (maxOut.size() > 0) mb = mb.max(maxOut);

			vmsg.push_back(mb);
		}
	}

	m_mcMessage2to1 = vmsg;
}

// output join tree to the Graphviz file format
void JTree::toDot(const std::string& fileName) {

	std::ofstream out( fileName.c_str() );
	if (out.fail()) {
		std::cerr << "Cannot open Graphviz file.";
		exit(EXIT_FAILURE);
	}

	out << "digraph g {\n";
	out << "node [shape = record];\n";
	out << "size = \"10, 7.5\";\n";
	out << "rotate = \"90\";\n";
	out << "ratio = \"fill\";\n";

	// nodes
	for (std::vector<JTreeNode*>::iterator it = m_clusters.begin();
			it != m_clusters.end(); ++it) {
		JTreeNode* n = (*it);

		out << "node" << n->id()
			<< "[ label = \"" << n->id() << " (" << n->toString() << ")" << "\"];\n";
	}

	// edges
	for (std::vector<JTreeEdge*>::iterator it = m_edges.begin();
			it != m_edges.end(); ++it) {
		JTreeEdge* e = (*it);

		out << "node" << e->getNode1()->id()
			<< " -> node" << e->getNode2()->id() << ";\n";
	}

	out << "}\n";
}


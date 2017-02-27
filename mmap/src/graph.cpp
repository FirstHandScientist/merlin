/*
 * graph.cpp
 *
 *  Created on: Mar 29, 2013
 *      Author: radu
 */

#include "graph.h"

nCost Graph::scoreMinfill(const int& i) {
	std::map<int, std::set<int> >::iterator iti = m_neighbors.find(i);
	assert(iti != m_neighbors.end());

	std::set<int>& S = iti->second;
	nCost c = 0;
	std::set<int>::iterator it1, it2;
	for (it1 = S.begin(); it1 != S.end(); ++it1) {
		it2 = it1;
		while (++it2 != S.end()) {
			if (!hasEdge(*it1, *it2))
				++c;
		}
	}

	return c;
}

// Extract the connected components
std::map<int, std::set<int> > Graph::connectedComponents(const std::set<int>& s) {

	std::map<int, std::set<int> > comps; // will hold the components

	// return zero components for empty argument set
	if (s.size() == 0) {
		return comps;
	}

	// keeps track of which node has been visited
	std::set<int> nodes(s.begin(), s.end());

	// stack for DFS search
	std::stack<int> dfs;

	std::set<int> comp;
	int c = 0; // current component and index
	for (; !nodes.empty(); ++c) {

		comp.clear();

		dfs.push(*nodes.begin());
		nodes.erase(dfs.top());

		// depth first search
		while (!dfs.empty()) {
			std::set<int>& neighbors = m_neighbors[dfs.top()];
			comp.insert(dfs.top());
			dfs.pop();
			for (std::set<int>::const_iterator it = neighbors.begin();
					it != neighbors.end(); ++it) {
				if (nodes.find(*it) != nodes.end()) {
					dfs.push(*it);
					nodes.erase(*it);
				}
			}
		}

		// record component
		comps[c] = comp;
	}

	return comps;
}

// Triangulate the graph; elimOrder[0] is the first variable to eliminate
int Graph::triangulate(const std::vector<int>& elimOrder) {

	int w = 0;

	// declare variables
	int n = (int)elimOrder.size();
	std::vector<int> position(n);

	// set the position map
	for (int i = 0; i < n; ++i) {
		position[elimOrder[i]] = i;
	}

	// eliminate variables and create induced edges
	for (int i = 0; i < n; ++i) {
		int var = elimOrder[i];
		int pos = position[var];

		const std::set<int>& ns = getNeighbors(var);
		std::vector<int> tmp;
		for (std::set<int>::const_iterator si = ns.begin(); si != ns.end(); ++si) {
			if (position[*si] > pos) {
				tmp.push_back(*si);
			}
		}

		// connect the earlier neighbors
		int m = (int)tmp.size();
		w = std::max(w, m);
		for (int j = 0; j < m-1; ++j) {
			for (int k = j+1; k < m; ++k) {
				if (!hasEdge(tmp[j], tmp[k])) {
					addEdge(tmp[j], tmp[k]);
				}
			}
		}
	}

	return w;
}

/* computes an (unconstrained) elimination order into 'elim' and returns its tree width */
int Graph::eliminate(std::vector<int>& elim) {

	int width = UNKNOWN;

	Graph G(*this);
	int n = G.getNumNodes();

	if (elim.size())
		elim.clear();
	elim.reserve(n);

	const std::set<int>& nodes = G.getNodes();

	// keeps track of node scores
	std::vector<nCost> scores;

	// have initial scores already been computed? then reuse!
	if (m_initialScores.size()) {
		scores = m_initialScores;
	} else {
		scores.resize(n);
		for (std::set<int>::const_iterator it = nodes.begin();
				it != nodes.end(); ++it) {
			scores[*it] = G.scoreMinfill(*it);
		}
		m_initialScores = scores;
	}

	int nextNode = NONE;
	nCost minScore = UNKNOWN;

	// eliminate nodes until all gone
	while (G.getNumNodes() != 0) {

		// keeps track of minimal score nodes
		std::vector<int> minCand; // minimal score of 1 or higher
		std::vector<int> simplicial; // simplicial nodes (score 0)

		minScore = std::numeric_limits<nCost>::max();

		// find node to eliminate
		for (int i = 0; i < n; ++i) {
			if (scores[i] == 0) { // score 0
				simplicial.push_back(i);
			} else if (scores[i] < minScore) { // new, lower score (but greater 0)
				minScore = scores[i];
				minCand.clear();
				minCand.push_back(i);
			} else if (scores[i] == minScore) { // current min. found again
				minCand.push_back(i);
			}
		}

		// eliminate all nodes with score=0 -> no edges will have to be added
		for (std::vector<int>::iterator it = simplicial.begin();
				it != simplicial.end(); ++it) {
			elim.push_back(*it);
			width = std::max(width, (int) G.getNeighbors(*it).size());
			G.removeNode(*it);
			scores[*it] = std::numeric_limits<nCost>::max();
		}

		// anything left to eliminate? If not, we are done!
		if (minScore == std::numeric_limits<nCost>::max())
			return width;

		// Pick one of the minimal score nodes (with score >= 1),
		// breaking ties randomly
		nextNode = minCand[rand::next(minCand.size())];
		elim.push_back(nextNode);

		// remember it's neighbors, to be used later
		const std::set<int>& neighbors = G.getNeighbors(nextNode);

		// update width of implied tree decomposition
		width = std::max(width, (int) neighbors.size());

		// connect neighbors in primal graph
		G.addClique(neighbors);

		// compute candidates for score update (node's neighbors and their neighbors)
		std::set<int> updateCand(neighbors);
		for (std::set<int>::const_iterator it = neighbors.begin();
				it != neighbors.end(); ++it) {
			const std::set<int>& X = G.getNeighbors(*it);
			updateCand.insert(X.begin(), X.end());
		}
		updateCand.erase(nextNode);

		// remove node from primal graph
		G.removeNode(nextNode);
		scores[nextNode] = std::numeric_limits<nCost>::max(); // tag score

		// update scores in primal graph (candidate nodes computed earlier)
		for (std::set<int>::const_iterator it = updateCand.begin();
				it != updateCand.end(); ++it) {
			scores[*it] = G.scoreMinfill(*it);
		}
	}

	return width;
}


// Minfill (constrained) elimination order into 'elim' and returns its tree width
int Graph::eliminate(std::set<int> mapVars, std::vector<int>& elim) {

	// treewidth
	int width = UNKNOWN;

	Graph G(*this);
	int n = G.getNumNodes();

	if (elim.size())
		elim.clear();
	elim.reserve(n);

	const std::set<int>& nodes = G.getNodes();
	std::set<int> sumVars(nodes);
	for (std::set<int>::const_iterator si = mapVars.begin();
			si != mapVars.end(); ++si) {
		sumVars.erase( *si );
	}

	// keeps track of node scores
	std::vector<nCost> scores;

	// have initial scores already been computed? then reuse!
	if (m_initialScores.size()) {
		scores = m_initialScores;
	} else {
		scores.resize(n);
		for (std::set<int>::const_iterator it = nodes.begin();
				it != nodes.end(); ++it) {
			scores[*it] = G.scoreMinfill(*it);
		}
		m_initialScores = scores;
	}

	int nextNode = NONE;
	nCost minScore = UNKNOWN;

	// eliminate nodes until all gone
	while (G.getNumNodes() != 0) {

		// keeps track of minimal score nodes
		std::vector<int> minCand; // minimal score of 1 or higher
		std::vector<int> simplicial; // simplicial nodes (score 0)

		minScore = std::numeric_limits<nCost>::max();

		// find node to eliminate
		std::set<int>::iterator start, stop;
		start = (sumVars.empty() ? mapVars.begin() : sumVars.begin());
		stop = (sumVars.empty() ? mapVars.end(): sumVars.end());
		for (std::set<int>::iterator si = start; si != stop; ++si) {
			int var = *si;
			if (scores[var] == 0) {
				simplicial.push_back(var);
			} else if (scores[var] < minScore) { // new, lower score (but greater 0)
				minScore = scores[var];
				minCand.clear();
				minCand.push_back(var);
			} else if (scores[var] == minScore) { // current min. found again
				minCand.push_back(var);
			}
		}

		// eliminate all nodes with score=0 -> no edges will have to be added
		for (std::vector<int>::iterator it = simplicial.begin();
				it != simplicial.end(); ++it) {
			int var = *it;
			elim.push_back(var);
			width = std::max(width, (int) G.getNeighbors(var).size());
			G.removeNode(var);
			scores[var] = std::numeric_limits<nCost>::max();
			// remove the node from either the SUM of MAP variables
			if (sumVars.find(var) != sumVars.end()) {
				sumVars.erase(var);
			} else if (mapVars.find(var) != mapVars.end()) {
				mapVars.erase(var);
			}
		}

		// anything left to eliminate? If not, we are done!
		if (minScore == std::numeric_limits<nCost>::max()) {
			continue;
		}

		// Pick one of the minimal score nodes (with score >= 1),
		// breaking ties randomly
		nextNode = minCand[rand::next(minCand.size())];
		elim.push_back(nextNode);

		// remember it's neighbors, to be used later
		const std::set<int>& neighbors = G.getNeighbors(nextNode);

		// update width of implied tree decomposition
		width = std::max(width, (int) neighbors.size());

		// connect neighbors in primal graph
		G.addClique(neighbors);

		// compute candidates for score update (node's neighbors and their neighbors)
		std::set<int> updateCand(neighbors);
		for (std::set<int>::const_iterator it = neighbors.begin();
				it != neighbors.end(); ++it) {
			const std::set<int>& X = G.getNeighbors(*it);
			updateCand.insert(X.begin(), X.end());
		}
		updateCand.erase(nextNode);

		// remove node from primal graph
		G.removeNode(nextNode);
		scores[nextNode] = std::numeric_limits<nCost>::max(); // tag score

		// update scores in primal graph (candidate nodes computed earlier)
		for (std::set<int>::const_iterator it = updateCand.begin();
				it != updateCand.end(); ++it) {
			scores[*it] = G.scoreMinfill(*it);
		}

		// remove the node from either the SUM of MAP variables
		if (sumVars.find(nextNode) != sumVars.end()) {
			sumVars.erase(nextNode);
		} else if (mapVars.find(nextNode) != mapVars.end()) {
			mapVars.erase(nextNode);
		}
	}

	return width;
}

/*
// Minfill (constrained) elimination order into 'elim' and returns its tree width
int Graph::eliminateS(std::set<int> mapVars, std::vector<int>& elim) {

	// treewidth
	int width = UNKNOWN;
	int widthS = UNKNOWN;

	Graph G(*this);
	int n = G.getNumNodes();

	if (elim.size())
		elim.clear();
	elim.reserve(n);

	const std::set<int>& nodes = G.getNodes();
	std::set<int> sumVars(nodes);
	for (std::set<int>::const_iterator si = mapVars.begin();
			si != mapVars.end(); ++si) {
		sumVars.erase( *si );
	}

	// keeps track of node scores
	std::vector<nCost> scores;

	// have initial scores already been computed? then reuse!
	if (m_initialScores.size()) {
		scores = m_initialScores;
	} else {
		scores.resize(n);
		for (std::set<int>::const_iterator it = nodes.begin();
				it != nodes.end(); ++it) {
			scores[*it] = G.scoreMinfill(*it);
		}
		m_initialScores = scores;
	}

	int nextNode = NONE;
	nCost minScore = UNKNOWN;

	// eliminate nodes until all gone
	while (G.getNumNodes() != 0) {

		// keeps track of minimal score nodes
		std::vector<int> minCand; // minimal score of 1 or higher
		std::vector<int> simplicial; // simplicial nodes (score 0)

		minScore = std::numeric_limits<nCost>::max();

		// find node to eliminate
		std::set<int>::iterator start, stop;
		start = (sumVars.empty() ? mapVars.begin() : sumVars.begin());
		stop = (sumVars.empty() ? mapVars.end(): sumVars.end());
		for (std::set<int>::iterator si = start; si != stop; ++si) {
			int var = *si;
			if (scores[var] == 0) {
				simplicial.push_back(var);
			} else if (scores[var] < minScore) { // new, lower score (but greater 0)
				minScore = scores[var];
				minCand.clear();
				minCand.push_back(var);
			} else if (scores[var] == minScore) { // current min. found again
				minCand.push_back(var);
			}
		}

		// eliminate all nodes with score=0 -> no edges will have to be added
		for (std::vector<int>::iterator it = simplicial.begin();
				it != simplicial.end(); ++it) {
			int var = *it;
			elim.push_back(var);
			width = std::max(width, (int) G.getNeighbors(var).size());
			G.removeNode(var);
			scores[var] = std::numeric_limits<nCost>::max();
			// remove the node from either the SUM of MAP variables
			if (sumVars.find(var) != sumVars.end()) {
				sumVars.erase(var);
			} else if (mapVars.find(var) != mapVars.end()) {
				mapVars.erase(var);
			}
		}

		// anything left to eliminate? If not, we are done!
		if (minScore == std::numeric_limits<nCost>::max()) {
			continue;
		}

		// Pick one of the minimal score nodes (with score >= 1),
		// breaking ties randomly
		nextNode = minCand[rand::next(minCand.size())];
		elim.push_back(nextNode);

		// remember it's neighbors, to be used later
		const std::set<int>& neighbors = G.getNeighbors(nextNode);

		// update width of implied tree decomposition
		width = std::max(width, (int) neighbors.size());

		// conunt only the SUM neighbors
		int svars = 0;
		for (std::set<int>::const_iterator ci = neighbors.begin();
				ci != neighbors.end(); ++ci) {
			if (mapVars.find(*ci) == mapVars.end())
				svars++;
		}
		widthS = std::max(widthS, svars);

		// connect neighbors in primal graph
		G.addClique(neighbors);

		// compute candidates for score update (node's neighbors and their neighbors)
		std::set<int> updateCand(neighbors);
		for (std::set<int>::const_iterator it = neighbors.begin();
				it != neighbors.end(); ++it) {
			const std::set<int>& X = G.getNeighbors(*it);
			updateCand.insert(X.begin(), X.end());
		}
		updateCand.erase(nextNode);

		// remove node from primal graph
		G.removeNode(nextNode);
		scores[nextNode] = std::numeric_limits<nCost>::max(); // tag score

		// update scores in primal graph (candidate nodes computed earlier)
		for (std::set<int>::const_iterator it = updateCand.begin();
				it != updateCand.end(); ++it) {
			scores[*it] = G.scoreMinfill(*it);
		}

		// remove the node from either the SUM of MAP variables
		if (sumVars.find(nextNode) != sumVars.end()) {
			sumVars.erase(nextNode);
		} else if (mapVars.find(nextNode) != mapVars.end()) {
			mapVars.erase(nextNode);
		}
	}

	return widthS;
}
*/

// write the graph to the output file (graphviz format)
void Graph::outputToFile(const char* output) {
	std::ostringstream oss;
	oss << "graph {\n";
	//oss << "node [shape = record];\n";
	oss << "size = \"10, 7.5\";\n";
	oss << "rotate = \"90\";\n";
	oss << "ratio = \"fill\";\n";

	std::set<int> nodes = getNodes();
	std::vector<int> temp;
	std::copy(nodes.begin(), nodes.end(), std::back_inserter(temp));
	int n = (int)nodes.size();
	for (int i = 0; i < n - 1; ++i) {
		for (int j = i + 1; j < n; ++j) {
			int u = temp[i];
			int v = temp[j];
			if (hasEdge(u, v)) {
				oss << "x" << u << " -- x" << v << ";" << std::endl;
			}
		}
	}
	oss << "}\n";

	std::ofstream of;
	of.open(output, std::ios_base::out);
	of << oss.str() << std::endl;
	of.close();
}

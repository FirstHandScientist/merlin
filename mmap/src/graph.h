/*
 * graph.h
 *
 *  Created on: Mar 28, 2013
 *      Author: radu
 */

#ifndef IBM_ANYTIME_GRAPH_H_
#define IBM_ANYTIME_GRAPH_H_

#include "base.h"

typedef int nCost;

/*
 * Implements an undirected graph. Internal representation
 * is an adjacency list, which stores the neighbors of each
 * node. It could could thus represent a directed graph, but
 * in the implementation edges and their reverse edges are
 * kept in sync.
 */
class Graph {

protected:

	// Adjacency lists
	std::map<int, std::set<int> > m_neighbors;

    // No. of vertices
	size_t m_numNodes;

    // No. of edges
	size_t m_numEdges;

public:

	// Constructor and destructor
	Graph(const int& n);
	~Graph() {};

public:

	// Get the neighbors of a node
	const std::set<int>& getNeighbors(const int& var);

	// Get the nodes in the graph
	std::set<int> getNodes();

	// Get the number of nodes
	inline size_t getNumNodes() {
		return m_numNodes;
	}

	// Get the number of edges
	inline size_t getNumEdges() {
		return m_numEdges;
	}

	// Get the graph density
	double getDensity();

public:

	// Add a node or an edge to the graph
	void addNode(const int& i);
	void addEdge(const int& i, const int& j);

	// Remove a node or an edge from the graph
	void removeNode(const int& i);
	void removeEdge(const int& i, const int& j);

	// Check if a node or an edge is in the graph
	bool hasNode(const int& i);
	bool hasEdge(const int& i, const int& j);

public:

	nCost scoreMinfill(const int& i);

protected:

	// Add a new adjacency in the graph
	void addAdjacency(const int& i, const int& j);

	// Remove an adjency from the graph
	void removeAdjacency(const int& i, const int& j);

public:

	// Add a new clique to the graph
	void addClique(const std::set<int>& s);

	// Finds the connected components of the graph
	std::map<int, std::set<int> > connectedComponents(const std::set<int>& s);

	// Triangulate the graph
	int triangulate(const std::vector<int>& elimOrder);

	// Minfill elimination order (unconstrained); returns the treewidth
	int eliminate(std::vector<int>& elim);

	// Minfill elimination order (constrained); returns the treewidth
	int eliminate(std::set<int> mapVars, std::vector<int>& elim);

	// Output to file (graphviz format)
	void outputToFile(const char* output);

private:

	// initial minfill scores
	std::vector<nCost> m_initialScores;
};

/****************************
 *  Inline implementations  *
 ****************************/

/* Constructor */
inline Graph::Graph(const int& n) :
		m_numNodes(0), m_numEdges(0) {
}

/* Returns a node's neighbors */
inline const std::set<int>& Graph::getNeighbors(const int& var) {
	std::map<int, std::set<int> >::const_iterator it = m_neighbors.find(var);
	assert(it != m_neighbors.end());
	return it->second;
}

/* returns the set of graph nodes */
inline std::set<int> Graph::getNodes() {
	std::set<int> s;
	for (std::map<int, std::set<int> >::iterator it = m_neighbors.begin();
			it != m_neighbors.end(); ++it) {
		s.insert(it->first);
	}
	return s;
}

/* Computes the graph density*/
inline double Graph::getDensity() {
	if (m_numNodes) {
		return 2.0 * m_numEdges / (m_numNodes) / (m_numNodes - 1);
	} else {
		return 0.0;
	}
}

/* Adds a node to the graph */
inline void Graph::addNode(const int& i) {
	if (m_neighbors.find(i) == m_neighbors.end()) {
		m_neighbors.insert(std::make_pair(i, std::set<int>()));
		++m_numNodes;
	}
}

/* Adds the edge (i,j) to the graph, also adds the reversed edge. */
inline void Graph::addEdge(const int& i, const int& j) {
	addAdjacency(i, j);
	addAdjacency(j, i);
	++m_numEdges;
}

/* removes the node and all related edges from the graph */
inline void Graph::removeNode(const int& i) {
	std::map<int, std::set<int> >::iterator iti = m_neighbors.find(i);
	if (iti != m_neighbors.end()) {
		std::set<int>& S = iti->second;
		for (std::set<int>::iterator it = S.begin(); it != S.end(); ++it) {
			removeAdjacency(*it, i);
			--m_numEdges;
		}
		m_neighbors.erase(iti);
		--m_numNodes;
	}
}

/* removes a single (undirected) edge from the graph */
inline void Graph::removeEdge(const int& i, const int& j) {
	removeAdjacency(i, j);
	removeAdjacency(j, i);
	--m_numEdges;
}

/* returns TRUE iff node i is in the graph */
inline bool Graph::hasNode(const int& i) {
	std::map<int, std::set<int> >::iterator iti = m_neighbors.find(i);
	return (iti != m_neighbors.end());
}

/* returns TRUE iff edge between nodes i and j exists */
inline bool Graph::hasEdge(const int& i, const int& j) {
	std::map<int, std::set<int> >::iterator iti = m_neighbors.find(i);
	if (iti != m_neighbors.end()) {
		std::set<int>::iterator itj = iti->second.find(j);
		return ( itj != iti->second.end() );
	} else {
		return false;
	}
}

/* Adds the edge (i,j) to the adjacency list. Does not add the reversed edge. */
inline void Graph::addAdjacency(const int& i, const int& j) {
	// Make sure node exists
	addNode(i);
	// Add edge (i,j) to the graph.
	std::map<int, std::set<int> >::iterator it = m_neighbors.find(i); // guaranteed to find node
	std::set<int>& s = it->second;
	s.insert(j);
}

/* Removes the adjacency entry (i,j) but NOT the reverse (j,i) */
inline void Graph::removeAdjacency(const int& i, const int& j) {
	std::map<int, std::set<int> >::iterator iti = m_neighbors.find(i);
	if (iti != m_neighbors.end()) {
		iti->second.erase(j);
	}
}

/* Adds the nodes in s to the graph and fully connects them */
inline void Graph::addClique(const std::set<int>& s) {
	for (std::set<int>::const_iterator it = s.begin(); it != s.end(); ++it) {
		addNode(*it); // insert the node
		std::set<int>::const_iterator it2 = it;
		while (++it2 != s.end()) {
			addEdge(*it, *it2);
		}
	}
}

#endif /* GRAPH_H_ */

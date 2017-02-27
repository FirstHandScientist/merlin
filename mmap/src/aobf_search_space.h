/*
 * cache_table.h
 *
 *  Created on: Apr 19, 2013
 *      Author: radu
 */

#ifndef IBM_MAP_AOBF_SEARCH_SPACE_H_
#define IBM_MAP_AOBF_SEARCH_SPACE_H_

#include "base.h"
#include "aobf_search_node.h"
#include "hash_murmur.h"

// Using the boost hash table implementation
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>

// A search state (corresponding to a search node)
struct AobfSearchState {
	int type;
	std::string context;
	bool operator==(const AobfSearchState& a) const {
		return (type == a.type && context == a.context);
	}
	AobfSearchState(int t, const std::string& s) : type(t), context(s) {}
};

// Hasher for search states
struct HashSearchState {
	size_t operator()(const AobfSearchState& a) const {
#ifdef HASH_MURMUR
		size_t hash = 0;
		size_t seed = 0x9e3779b9;
		int len = (int)a.context.size();
		register unsigned char *key = (unsigned char*)a.context.c_str();
		MurmurHash3_x64_64(key, len, seed, &hash);
		return hash;
#elif defined(HASH_DEFAULT)
		return boost::hash<std::string>()(a.context);
#endif
	}
};


/* The context minimal AND/OR search graph */
typedef boost::unordered_map<AobfSearchState, AobfSearchNode*, HashSearchState> AOGRAPH;
typedef std::vector<std::list<AobfSearchNode*> > LAYERS;

class AobfSearchSpace {

protected:

	// Hash table index by state
	AOGRAPH m_nodes;

	// Layers used for repairing (anytime weighted best-first search)
	LAYERS m_layers;

	// Number of AND nodes expanded
	size_t m_andNodes;

	// Number of OR nodes expanded
	size_t m_orNodes;

	// Root of the search space
	AobfSearchNode* m_root;

public:

	// add a new node to the graph
	void add(const AobfSearchState& state, AobfSearchNode* node);

	// find a node in the graph
	bool find(const AobfSearchState& state) const;

	// returns a node from the graph
	AobfSearchNode* get(const AobfSearchState& state);

	// erase a node from the graph
	void erase(const AobfSearchState& state);

	// clear the search space
	void clear();

	// returns the root
	inline AobfSearchNode* getRoot() {
		return m_root;
	}

	// set the root node
	inline void setRoot(AobfSearchNode* node) {
		m_root = node;
	}

	// returns the number of AND nodes expanded
	inline size_t getAndNodes() const {
		return m_andNodes;
	}

	// returns the number of OR nodes expanded
	inline size_t getOrNodes() const {
		return m_orNodes;
	}

	// increment the number of nodes expanded
	inline void incNodesExpanded(int nodeType) {
		assert(nodeType == NODE_AND || nodeType == NODE_OR);
		if (nodeType == NODE_AND) {
			m_andNodes += 1;
		} else {
			m_orNodes += 1;
		}
	}

	// get the nodes
	const AOGRAPH& getNodes() const {
		return m_nodes;
	}

	// get the nodes
	AOGRAPH& getNodes() {
		return m_nodes;
	}

	// get the layers (const reference)
	const LAYERS& getLayers() const {
		return m_layers;
	}

	// get the layers
	LAYERS& getLayers() {
		return m_layers;
	}

	// reserve the number of layers
	void initLayers(int n) {
		assert(n > 0);
		m_layers.resize(n);
	}

public:

	// constructor
	AobfSearchSpace();

	// destructor
	~AobfSearchSpace();

private:
	// prevents copy construction
	AobfSearchSpace(const AobfSearchSpace&);
};

/* inline declarations */

// constructor
inline AobfSearchSpace::AobfSearchSpace() : m_andNodes(0), m_orNodes(0), m_root(NULL) {

}

// destructor
inline AobfSearchSpace::~AobfSearchSpace() {
	clear();
}

// clears all nodes in the graph
inline void AobfSearchSpace::clear() {
	AOGRAPH::iterator mi = m_nodes.begin();
	for (; mi != m_nodes.end(); ++mi) {
		if (mi->second != NULL) {
			delete mi->second;
		}
	}
	m_nodes.clear();
}

// adds a node to the graph
inline void AobfSearchSpace::add(const AobfSearchState& state, AobfSearchNode* node) {
	assert(node != NULL);
	m_nodes.insert(std::make_pair(state, node));
	if ( !m_layers.empty() ) {
		int d = node->getDepth();
		assert( d >= 0 && d < (int)m_layers.size() );
		m_layers[d].push_front(node);
	}
}

// finds a node in the graph
inline bool AobfSearchSpace::find(const AobfSearchState& state) const {
	AOGRAPH::const_iterator mi = m_nodes.find(state);
	return (mi != m_nodes.end());
}

// returns a node from the graph
inline AobfSearchNode* AobfSearchSpace::get(const AobfSearchState& state) {
	AOGRAPH::iterator mi = m_nodes.find(state);
	return (mi != m_nodes.end() ? mi->second : NULL);
}

// erase a node from the graph
inline void AobfSearchSpace::erase(const AobfSearchState& state) {
	AOGRAPH::iterator mi = m_nodes.find(state);
	if (mi != m_nodes.end()) {
		AobfSearchNode* node = mi->second;
		m_nodes.erase(mi);
		delete node;
	}
}

////////////////////////////////////////////////////////////////////////////////

class AobfSearchSpace2 {

protected:

	// Hash table index by state
	std::vector<AOGRAPH> m_nodes;

	// Layers used for repairing (anytime weighted best-first search)
	LAYERS m_layers;

	// Number of AND nodes expanded
	size_t m_andNodes;

	// Number of OR nodes expanded
	size_t m_orNodes;

	// Root of the search space
	AobfSearchNode* m_root;

public:

	// add a new node to the graph
	void add(const AobfSearchState& state, AobfSearchNode* node);

	// find a node in the graph
	bool find(int var, const AobfSearchState& state) const;

	// returns a node from the graph
	AobfSearchNode* get(int var, const AobfSearchState& state);

	// erase a node from the graph
	void erase(int var, const AobfSearchState& state);

	// clear the search space
	void clear();

	// returns the root
	inline AobfSearchNode* getRoot() {
		return m_root;
	}

	// set the root node
	inline void setRoot(AobfSearchNode* node) {
		m_root = node;
	}

	// returns the number of AND nodes expanded
	inline size_t getAndNodes() const {
		return m_andNodes;
	}

	// returns the number of OR nodes expanded
	inline size_t getOrNodes() const {
		return m_orNodes;
	}

	// increment the number of nodes expanded
	inline void incNodesExpanded(int nodeType) {
		assert(nodeType == NODE_AND || nodeType == NODE_OR);
		if (nodeType == NODE_AND) {
			m_andNodes += 1;
		} else {
			m_orNodes += 1;
		}
	}

	// get the nodes
	const std::vector<AOGRAPH>& getNodes() const {
		return m_nodes;
	}

	// get the nodes
	std::vector<AOGRAPH>& getNodes() {
		return m_nodes;
	}

	// get the layers (const reference)
	const LAYERS& getLayers() const {
		return m_layers;
	}

	// get the layers
	LAYERS& getLayers() {
		return m_layers;
	}

	// reserve the number of layers
	void initLayers(int n) {
		assert(n > 0);
		m_layers.resize(n);
	}

public:

	// constructor
	AobfSearchSpace2(int size);

	// destructor
	~AobfSearchSpace2();

private:
	// prevents copy construction
	AobfSearchSpace2(const AobfSearchSpace2&);
};

/* inline declarations */

// constructor
inline AobfSearchSpace2::AobfSearchSpace2(int size) : m_andNodes(0), m_orNodes(0), m_root(NULL) {
	m_nodes.resize(size);
}

// destructor
inline AobfSearchSpace2::~AobfSearchSpace2() {
	clear();
}

// clears all nodes in the graph
inline void AobfSearchSpace2::clear() {
	for (size_t i = 0; i < m_nodes.size(); ++i) {
		AOGRAPH::iterator mi = m_nodes[i].begin();
		for (; mi != m_nodes[i].end(); ++mi) {
			if (mi->second != NULL) {
				delete mi->second;
			}
		}
	}

	m_nodes.clear();
}

// adds a node to the graph
inline void AobfSearchSpace2::add(const AobfSearchState& state, AobfSearchNode* node) {
	assert(node != NULL);
	int var = node->getVar();
	m_nodes[var].insert(std::make_pair(state, node));
	if ( !m_layers.empty() ) {
		int d = node->getDepth();
		assert( d >= 0 && d < (int)m_layers.size() );
		m_layers[d].push_front(node);
	}
}

// finds a node in the graph
inline bool AobfSearchSpace2::find(int var, const AobfSearchState& state) const {
	AOGRAPH::const_iterator mi = m_nodes[var].find(state);
	return (mi != m_nodes[var].end());
}

// returns a node from the graph
inline AobfSearchNode* AobfSearchSpace2::get(int var, const AobfSearchState& state) {
	AOGRAPH::iterator mi = m_nodes[var].find(state);
	return (mi != m_nodes[var].end() ? mi->second : NULL);
}

// erase a node from the graph
inline void AobfSearchSpace2::erase(int var, const AobfSearchState& state) {
	AOGRAPH::iterator mi = m_nodes[var].find(state);
	if (mi != m_nodes[var].end()) {
		AobfSearchNode* node = mi->second;
		m_nodes[var].erase(mi);
		delete node;
	}
}



#endif /* CACHE_TABLE_H_ */

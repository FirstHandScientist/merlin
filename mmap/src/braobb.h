/*
 * braobb.h
 *
 *  Created on: Jun 7, 2013
 *      Author: radu
 */

#ifndef IBM_ANYTIME_BRAOBB_H_
#define IBM_ANYTIME_BRAOBB_H_

#include "aobb.h"


// Breadth Rotating AND/OR Branch and Bound (BRAOBB) - [Otten and Dechter, 2012]

class BRAOBB : public AOBB {

protected:

	size_t m_stackCount;        // counter for current stack
	size_t m_stackLimit;        // expansion limit for stack rotation
	MyStack* m_rootStack;       // the root stack
	queue<MyStack*> m_stacks;   // the queue of active stacks

protected:
	void reset(AobbSearchNode* p);
	bool isDone() const;
	bool doExpand(AobbSearchNode* n);
	AobbSearchNode* nextNode();

public:
	void setStackLimit(size_t s) {
		m_stackLimit = s;
	}

	// Main solver routine
	int solve();

public:
	BRAOBB(ProgramOptions* opt);
	virtual ~BRAOBB();
};

/* inline definitions */

inline BRAOBB::BRAOBB(ProgramOptions* opt) : AOBB(opt) {
	m_stackCount = 0;
	m_stackLimit = opt->rotateLimit; // default: 1000 node expansions per stack
	m_rootStack = NULL;
}

inline BRAOBB::~BRAOBB() {

}

inline bool BRAOBB::isDone() const {
	assert(false); // TODO
	return false;
}

inline void BRAOBB::reset(AobbSearchNode* p) {
	assert(p);
	while (m_stacks.size()) {
		delete m_stacks.front();
		m_stacks.pop();
	}
	m_rootStack = new MyStack(NULL);
	m_rootStack->push(p);
	m_stacks.push(m_rootStack);
}

#endif /* BRAOBB_H_ */

/*
 * message.h
 *
 *  Created on: Oct 15, 2013
 *      Author: radu
 */

#ifndef IBM_MAP_MESSAGE_H_
#define IBM_MAP_MESSAGE_H_


#include "mini_bucket.h"


/**
 * Represents a "message" between two clusters in a join/bucket-tree and
 * can be used for full join/bucket-tree message passing as well as
 * mini-bucket/mini-cluset tree elimination (in which case it is represented
 * as a collection of smaller arity functions)
 *
 */

class Message {

protected:

	// message as a set of functions
	std::vector<Function*> m_message;

	// mini-bucket i-bound
	int m_ibound;

public:

	// add a function to the list
	inline void addFunction(Function* f) {
		assert(f != NULL);
		m_message.push_back(f);
	}

	// returns the list of mini-buckets
	inline std::vector<Function*>& getMessage() {
		return m_message;
	}

	// checks if exact message
	inline bool isExact() {
		return (m_ibound == NONE);
	}

	inline void setIBound(int ib) {
		m_ibound = ib;
	}

	inline int getIBound() {
		return m_ibound;
	}

	// clone a message
	Message* clone();

public:

	// constructor
	Message(int ibound = NONE);
	Message(Function*);

	// destructor
	~Message();

};

#endif /* MESSAGE_H_ */

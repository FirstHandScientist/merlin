/*
 * message.cpp
 *
 *  Created on: Oct 15, 2013
 *      Author: radu
 */

#include "message.h"

// constructor
Message::Message(int ibound) {
	m_ibound = ibound;
}

Message::Message(Function* f) {
	assert(f != NULL);
	m_message.push_back(f);
	m_ibound = NONE;
}

// destructor
Message::~Message() {

	for (std::vector<Function*>::iterator fi = m_message.begin();
			fi != m_message.end(); ++fi) {
		Function* f = (*fi);
		if ( !f->isOriginal() ) {
			delete f;
		}
	}

}

// clone a message
Message* Message::clone() {

	Message* cl = new Message();
	cl->setIBound(m_ibound);
	for (std::vector<Function*>::iterator fi = m_message.begin();
			fi != m_message.end(); ++fi) {
		cl->addFunction( (*fi)->clone() );
	}

	return cl;
}

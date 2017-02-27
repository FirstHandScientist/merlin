/*
 * DfpnpluCacheTable.h
 *
 */

#ifndef DFPNPLUSCACHETABLE_H_
#define DFPNPLUSCACHETABLE_H_

#include <malloc.h>

#include "base.h"
#include "utils.h"
#include "hash_zobrist.h"

struct RbfaooCacheElement {
	double result;
	RbfaooCacheElement *next;
	context_t context;
	Zobrist zobkey;
	int var;
	int best_index;
	int dn;
	unsigned int subtree_size :30;
	unsigned int or_node_flag :1;
	unsigned int solved_flag :1;
};

typedef RbfaooCacheElement * RbfaooCacheElementPtr;

class RbfaooCacheTable {

private:

	enum {
		GC_STEP_SIZE = 1
	};
	RbfaooCacheElement **m_table;
	RbfaooCacheElement *m_entries;
	size_t m_table_size;

	size_t m_gc_subtree_size;

	size_t m_num_saved;

	size_t m_num_newly_saved;

	RbfaooCacheElement *m_free_list;

	void smallTreeGC();

	size_t tableIndex(const Zobrist &zob) const {
		return zob.getC0() % m_table_size;
	}

public:
	enum {
		TERMINAL_INDEX = 0xffffffff
	};

	RbfaooCacheTable() :
			m_table(NULL), m_entries(NULL), m_table_size(0),
			m_gc_subtree_size(0), m_num_saved(0), m_num_newly_saved(0),
			m_free_list(NULL) {
	}

	~RbfaooCacheTable();

	void printWriteStats();
	void init(size_t size);
	void clear();

	RbfaooCacheElementPtr read(const Zobrist &zob, int var, int or_flag,
			const context_t &ctxt);

	RbfaooCacheElementPtr write(const Zobrist &zob, int var, context_t &ctxt,
			int or_flag, unsigned int best_index, double result, int dn,
			int solved_flag, unsigned int increased_subtree_size);

};

#endif /* CACHETABLE_H_ */

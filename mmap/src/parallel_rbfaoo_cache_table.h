/*
 * parallel_rbfaoo_cache_table.h
 * Programmed by Akihiro Kishimoto
 *
 */

#ifndef PARALLEL_RBFAOO_CACHE_TABLE_H
#define PARALLEL_RBFAOO_CACHE_TABLE_H

#include <malloc.h>
#include <vector>


#include "base.h"
#include "utils.h"
#include "hash_zobrist.h"
#include "light_mutex.h"
//#include <boost/thread.hpp>
#include "rbfaoo_cache_table.h"

//typedef boost::mutex LightMutex;

#define RBFAOO_EPS (1e-10)

struct ParallelRbfaooCacheElement : RbfaooCacheElement {
  double upperbound;
  double pseudo_result;
  short  num_threads;
  ParallelRbfaooCacheElement *next;
};

typedef ParallelRbfaooCacheElement * ParallelRbfaooCacheElementPtr;

struct RbfaooCacheLocalStruct {

  ParallelRbfaooCacheElement *m_entries;
  size_t m_gc_subtree_size;
  size_t m_num_saved;
  size_t m_num_newly_saved;
  
  ParallelRbfaooCacheElement *m_free_list;
};

class ParallelRbfaooCacheTable {

private:

	enum {
		GC_STEP_SIZE = 1, TERMINAL_INDEX = 0xffffffff
	};
	ParallelRbfaooCacheElement **m_table;
	size_t m_table_size;
	std::vector<RbfaooCacheLocalStruct *> m_thread_local;

	size_t m_num_threads;

        LightMutex *m_mutex_table;

 	LightMutex m_first_thread_lock; 
 
 	volatile int m_first_thread;

	void smallTreeGC(size_t thread_id);

	size_t tableIndex(const Zobrist &zob) const {
		return zob.getC0() % m_table_size;
	}

public:

	ParallelRbfaooCacheTable() : m_table(NULL), m_table_size(0), m_num_threads(0) { }

	~ParallelRbfaooCacheTable();

 	bool setFirstThread(int thread_id); 
 
 	int  getFirstThread() const { 
 		return m_first_thread;
 	}

	unsigned int getNumThreads() const {
		return m_num_threads;
	}

	void printWriteStats();
	void init(size_t num_threads, size_t size);

	bool read(const Zobrist &zob, int var,
		  int or_flag, const context_t &ctxt, double &value, double &upperbound, 
		  double &pseudo_value, int &dn, unsigned int &solved_flag);

	bool read(const Zobrist &zob, int var,
		  int or_flag, const context_t &ctxt, double &value, double &upperbound, 
		  double &pseudo_value, int &dn, unsigned int &solved_flag, 
		  unsigned int &num_threads);

	void writeandexit(size_t thread_id, const Zobrist &zob, int var, bool mapvar_flag, 
			  context_t &ctxt, int or_flag, unsigned int best_index, double result, 
			  double upperbound, int dn, int solved_flag, bool terminal_node_flag, 
			  unsigned int increased_subtree_size);

	ParallelRbfaooCacheElement *enter(size_t thread_id, const Zobrist &zob, int var, bool mapvar_flag, 
					  context_t &ctxt, int or_flag, double result);

};

#endif /* PARALLEL_RBFAOO_CACHE_TABLE_H */

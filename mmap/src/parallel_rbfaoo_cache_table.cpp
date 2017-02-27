/*
 * ParallelRbfaooCacheTable.cpp
 *
 */

#include <assert.h>

#include <algorithm>
#include <iostream>
#include "rbfaoo_search_node.h"
#include "parallel_rbfaoo_cache_table.h"

//const double virtual_value = 0.001;
const double virtual_value = 0.01;

using namespace std;

typedef LightMutex::scoped_lock scoped_lock_t;
//typedef boost::mutex::scoped_lock scoped_lock_t;

ParallelRbfaooCacheTable::~ParallelRbfaooCacheTable() {
	printWriteStats();
	if (m_table) {
		delete[] m_table;
		delete[] m_mutex_table;
		for  (size_t i = 0; i < m_num_threads; i ++) {
			delete [] m_thread_local[i]->m_entries;
			delete m_thread_local[i];
		}
	}
}

void ParallelRbfaooCacheTable::printWriteStats() {
	size_t total_saved = 0, total_newly_saved = 0;
	for (size_t i = 0; i < m_num_threads; i ++) {
		cout << "Thread " << i << ": # of times written in TT = " << m_thread_local[i]->m_num_saved << endl;
		cout << "Thread " << i << ": # of nodes newly saved in TT = " << m_thread_local[i]->m_num_newly_saved << endl;
		total_saved += m_thread_local[i]->m_num_saved;
		total_newly_saved += m_thread_local[i]->m_num_newly_saved;
	}
	cout << "Total # of times written in TT = " << total_saved << endl;
	cout << "Total # of nodes newly saved in TT = " << total_newly_saved << endl;
}

void ParallelRbfaooCacheTable::init(size_t num_threads, size_t size) {
	if (m_table) {
		delete[] m_table;
		for  (size_t i = 0; i < m_num_threads; i ++) {
			delete [] m_thread_local[i]->m_entries;
			delete m_thread_local[i];
		}
		m_thread_local.resize(0);
	}

	// kishi: need to think about how much size is used for context
	m_table_size = size * 1024
			/ (sizeof(ParallelRbfaooCacheElement) + sizeof(ParallelRbfaooCacheElementPtr));
	m_table = new ParallelRbfaooCacheElementPtr[m_table_size];

	m_mutex_table = new LightMutex[m_table_size];

	cout << "Num of cache entries = " << m_table_size << "\n";
	
	m_num_threads = num_threads;
	size_t entry_size = m_table_size / m_num_threads;
	for (size_t i = 0; i < m_num_threads; i ++) {
		m_thread_local.push_back(new RbfaooCacheLocalStruct());
		m_thread_local[i]->m_entries = new ParallelRbfaooCacheElement[entry_size];
		m_thread_local[i]->m_free_list = m_thread_local[i]->m_entries;
		for (size_t j = 0; j < entry_size - 1; j++)
		  m_thread_local[i]->m_entries[j].next = &m_thread_local[i]->m_entries[j + 1];
		m_thread_local[i]->m_entries[entry_size - 1].next = NULL;
		m_thread_local[i]->m_gc_subtree_size = 0; // Kishi: I made up a number
	}
	fill(&m_table[0], &m_table[m_table_size],
			reinterpret_cast<ParallelRbfaooCacheElementPtr>(NULL));

	m_first_thread = -1;
}

void ParallelRbfaooCacheTable::smallTreeGC(size_t thread_id) {

	const double TABLE_OCCUPACY = 0.70;
	size_t num_collected = 0;
	double occupacy = 1.0;
	context_t empty_context;

	//cout << "Thread " << thread_id << ": Garbage collecting...\n";
	size_t range = m_table_size / m_num_threads;
	size_t start = range * thread_id;
	do {
		unsigned int smallest_subtree = 0x7fffffff;
		size_t i = start;
		for (size_t num = 0; num < m_table_size; num ++, i = (i + 1) % m_table_size) {
			ParallelRbfaooCacheElementPtr prev, element;
			size_t mutex_index = i % m_table_size;
			scoped_lock_t lock(m_mutex_table[mutex_index]);
			for (prev = NULL, element = m_table[i]; element != NULL;) {
				if (element->num_threads == 0 && element->subtree_size <= m_thread_local[thread_id]->m_gc_subtree_size) {
					ParallelRbfaooCacheElementPtr tmp;

					element->context.swap(empty_context);
					if (prev == NULL) {
						assert(element == m_table[i]);
						m_table[i] = m_table[i]->next;
					} else
						prev->next = element->next;
					num_collected ++;
					tmp = element;
					element = element->next;
					tmp->next = m_thread_local[thread_id]->m_free_list;
					m_thread_local[thread_id]->m_free_list = tmp;
				} else {
					if (element->num_threads == 0) 					
						smallest_subtree = min(smallest_subtree, element->subtree_size);
					prev = element;
					element = element->next;
				}
			}
			occupacy = 1.0 - static_cast<double>(num_collected) * m_num_threads / m_table_size;
			if (occupacy < TABLE_OCCUPACY)
				break;
			if (m_first_thread >= 0 && m_thread_local[thread_id]->m_free_list != NULL)
				return;
		}
		//occupacy = 1.0 - static_cast<double>(num_collected) * m_num_threads / m_table_size;
		//occupacy = 1.0 - static_cast<double>(num_collected) / num_total;
#ifdef DEBUG
		cout << "! Thread" << thread_id << ": Garbage collection: Subtree Size = " << std::setw(4) << m_thread_local[thread_id]->m_gc_subtree_size
		     << " " << num_collected << " collected " << " : Occupacy = " << std::setiosflags(std::ios::fixed)
		     << std::setprecision(4) << occupacy << " % !\n";
#endif
		m_thread_local[thread_id]->m_gc_subtree_size = smallest_subtree;
	} while (occupacy >= TABLE_OCCUPACY);
}

bool ParallelRbfaooCacheTable::read(const Zobrist &zob, int var,
				    int or_flag, const context_t &ctxt, double &value, double &upperbound, double &pseudo_value, int &dn, 
				    unsigned int &solved_flag, unsigned int &num_threads) {
	size_t index = tableIndex(zob);
	//KISHI: SAFE APPROACH BUT MIGHT NEED TO THINK ABOUT HOW TO MAKE THIS LOCK-FREE
	size_t mutex_index = index % m_table_size;
	scoped_lock_t lock(m_mutex_table[mutex_index]);
	for (ParallelRbfaooCacheElement *element = m_table[index]; element != NULL;
			element = element->next) {
		if (zob == element->zobkey && var == element->var
			&& element->or_node_flag == or_flag && element->context == ctxt) {
			value = element->result;
			upperbound = element->upperbound;
			dn = element->dn;
			solved_flag = element->solved_flag;
			pseudo_value = element->pseudo_result;
			num_threads = element->num_threads;
			if (upperbound == value && value != ELEM_INF) {
			  ; assert(solved_flag);
			}
			//Kishi: adding diversity at AND node
			if (dn != 0 && dn != DNINF && or_flag) {
				dn += element->num_threads;
			}
			return true;
		}
	}
	return false;
}

bool ParallelRbfaooCacheTable::read(const Zobrist &zob, int var,
				    int or_flag, const context_t &ctxt, double &value, double &upperbound, double &pseudo_value, int &dn, unsigned int &solved_flag) {
	unsigned int dummy;
	return read(zob, var, or_flag, ctxt, value, upperbound, pseudo_value, dn, solved_flag, dummy);
}

void ParallelRbfaooCacheTable::writeandexit(size_t thread_id, const Zobrist &zob, int var, bool mapvar_flag, 
				    context_t &ctxt, int or_flag, unsigned int best_index, double result, 
				    double upperbound, int dn, int solved_flag, bool terminal_node_flag, 
				    unsigned int increased_subtree_size) {
	size_t index = tableIndex(zob);

	m_thread_local[thread_id]->m_num_saved++;

	//if (zob.getC0() == 3377473968147770893 && zob.getC1() == 1277974532) {
	//cerr << "BREAK: r=" << result << "ub=" << upperbound << "sf=" << solved_flag << endl;
	//}

	//KISHI: SAFE APPROACH BUT MIGHT NEED TO THINK ABOUT HOW TO MAKE THIS LOCK-FREE
	{
		size_t mutex_index = index % m_table_size;
		scoped_lock_t lock(m_mutex_table[mutex_index]);
		for (ParallelRbfaooCacheElement *element = m_table[index]; element != NULL;
		     element = element->next) {
			if (zob == element->zobkey && var == element->var
					&& element->or_node_flag == or_flag
					&& element->context == ctxt) {
				if (mapvar_flag) {
					element->result = max(element->result, result);
					element->pseudo_result = max(result, element->pseudo_result);
					if (element->dn != 0 && element->dn != DNINF) {
						element->dn = dn;
						element->solved_flag = solved_flag;
						element->best_index = best_index;
					}
				}
				else {
				  if (solved_flag) {
					element->result = result;
					element->pseudo_result = result;
					element->dn = dn;
					element->solved_flag = solved_flag;
					element->best_index = best_index;
				  }
				}
				element->subtree_size += increased_subtree_size;
				element->upperbound = min(element->upperbound, upperbound);
				//WHY DO WE NEED TO CHECK TERMINAL_NODE_FLAG
#if 000
				if (!terminal_node_flag) {
				  //cerr << "Thread " << thread_id << ": write (" << zob.getC0() 
				  //<< "," << zob.getC1() << "): nth=" << element->num_threads << "\n"; 
					; assert(element->num_threads > 0);
					element->num_threads --;
				}
				else 
					element->best_index = TERMINAL_INDEX;
#else
				; assert(element->num_threads > 0);
				element->num_threads --;
				if (terminal_node_flag) 
					element->best_index = TERMINAL_INDEX;

#endif
				if (mapvar_flag) {
					if (!solved_flag && element->result != ELEM_INF && element->result + RBFAOO_EPS >= element->upperbound) {
						element->solved_flag = 1;
						element->dn = DNINF;
					}
					if (element->solved_flag || element->num_threads == 0)
						element->pseudo_result = result;
					if (element->result == element->upperbound && element->result != ELEM_INF) {
						; assert(element->solved_flag);
						; assert(element->dn == DNINF);
					}
				}
				return;
			}
		}
	}
	if (m_thread_local[thread_id]->m_free_list == NULL)
		smallTreeGC(thread_id);

	; assert(0);
	ParallelRbfaooCacheElementPtr element = m_thread_local[thread_id]->m_free_list;
	m_thread_local[thread_id]->m_num_newly_saved++;
	m_thread_local[thread_id]->m_free_list = m_thread_local[thread_id]->m_free_list->next;
	element->zobkey = zob;
	element->subtree_size = increased_subtree_size;
	element->var = var;
	element->result = result;
	element->pseudo_result = result;
	element->dn = dn;
	element->solved_flag = solved_flag;
	element->context.swap(ctxt);
	if (!terminal_node_flag) 
		element->best_index = best_index;
	else 
		element->best_index = TERMINAL_INDEX;
	element->or_node_flag = or_flag;
	element->upperbound = upperbound;
	element->num_threads = 0; 
	if (element->result == element->upperbound && element->result != ELEM_INF) {
	  ; assert(element->solved_flag);
	}
	bool written_by_another_flag = false;
	{
		size_t mutex_index = index % m_table_size;
		scoped_lock_t lock(m_mutex_table[mutex_index]);
		// Need to check again if the same entry is added when the lock was released. 
		for (ParallelRbfaooCacheElement *e = m_table[index]; e != NULL; e = e->next) {
			if (zob == e->zobkey && var == e->var
					&& e->or_node_flag == or_flag
					&& e->context == element->context) {
				if (mapvar_flag) {
					e->result = max(element->result, result);
					e->pseudo_result = max(result, e->pseudo_result);
					if (e->dn != 0 && e->dn != DNINF) {
						e->dn = dn;
						e->solved_flag = solved_flag;
						if (!terminal_node_flag) 
						  e->best_index = best_index;
						else
						  e->best_index = TERMINAL_INDEX;
					}
				}
				else {
				  if (solved_flag) {
					e->result = result;
					e->pseudo_result = result;
					e->solved_flag = solved_flag;
					if (!terminal_node_flag) 
					  e->best_index = best_index;
					else
					  e->best_index = TERMINAL_INDEX;
				  }
				}
				e->upperbound = min(element->upperbound, upperbound);
				e->subtree_size += increased_subtree_size;
				//If the same entry was written, do not use the entry
				written_by_another_flag = true;
				if (mapvar_flag) {
					if (!solved_flag && e->result != ELEM_INF && e->result + RBFAOO_EPS >= e->upperbound) {
						e->solved_flag = 1;
						e->dn = DNINF;
					}
					if (e->result == e->upperbound && e->result != ELEM_INF) {
						; assert(e->solved_flag);
						; assert(e->dn == DNINF);
					}
				}
				break;
			}
		}
		if (!written_by_another_flag) {
			element->next = m_table[index];
			m_table[index] = element;
		}
	}
	if (written_by_another_flag) {
		element->next = m_thread_local[thread_id]->m_free_list;
		m_thread_local[thread_id]->m_free_list = element;
	}
}

ParallelRbfaooCacheElement *ParallelRbfaooCacheTable::enter(size_t thread_id, const Zobrist &zob, int var, bool mapvar_flag, 
				     context_t &ctxt, int or_flag, double result) {
	size_t index = tableIndex(zob);

	m_thread_local[thread_id]->m_num_saved++;
	//KISHI: SAFE APPROACH BUT MIGHT NEED TO THINK ABOUT HOW TO MAKE THIS LOCK-FREE
	{
		size_t mutex_index = index % m_table_size;
		scoped_lock_t lock(m_mutex_table[mutex_index]);
		for (ParallelRbfaooCacheElement *element = m_table[index]; element != NULL;
		     element = element->next) {
			if (zob == element->zobkey && var == element->var
					&& element->or_node_flag == or_flag
					&& element->context == ctxt) {
				if (mapvar_flag) {
					element->result = max(element->result, result);
					element->pseudo_result += virtual_value;
				}
				//else 
				//	element->result = min(element->result, result);
				element->num_threads ++;
				//cerr << "Thread " << thread_id << ": enterA (" << zob.getC0() 
				//<< "," << zob.getC1() << "): nth=" << element->num_threads << "\n"; 
				return element;
			}
		}
	}

	if (m_thread_local[thread_id]->m_free_list == NULL)
		smallTreeGC(thread_id);

	ParallelRbfaooCacheElementPtr element = m_thread_local[thread_id]->m_free_list;
	m_thread_local[thread_id]->m_num_newly_saved++;
	m_thread_local[thread_id]->m_free_list = m_thread_local[thread_id]->m_free_list->next;
	element->zobkey = zob;
	element->subtree_size = 0;
	element->var = var;
	element->result = result;
	//element->pseudo_result = result;
	if (mapvar_flag) 
		element->pseudo_result = result + virtual_value;
	else
		element->pseudo_result = result;
	element->solved_flag = 0;
	element->context.swap(ctxt);
	element->best_index = 0;
	element->or_node_flag = or_flag;
	element->upperbound = ELEM_INF;
	element->num_threads = 1;
	element->dn = 1;
	//cerr << "Thread " << thread_id << ": enterB (" << zob.getC0() 
	//<< "," << zob.getC1() << "): nth=" << element->num_threads << "\n"; 
	bool written_by_another_flag = false;
	ParallelRbfaooCacheElement *r = element;
	{
		size_t mutex_index = index % m_table_size;
		scoped_lock_t lock(m_mutex_table[mutex_index]);
		//cerr << "Thread " << thread_id << ": acquire (" << zob.getC0() 
		//<< "," << zob.getC1() << "):mutex_index = " << mutex_index << "index=" << index 
		//<< ":tbl=" << m_table[index] << "\n"; 
		// Need to check again if the same entry is added when the lock was released. 
		for (ParallelRbfaooCacheElement *e = m_table[index]; e != NULL; e = e->next) {
			if (zob == e->zobkey && var == e->var
					&& e->or_node_flag == or_flag
					&& e->context == element->context) {
				if (mapvar_flag) 
					e->result = max(element->result, result);
				//else
				//	e->result = result;
				e->num_threads ++;
				//If the same entry was written, do not use the entry
				written_by_another_flag = true;
				//cerr << "Thread " << thread_id << ": enterC (" << zob.getC0() 
				//<< "," << zob.getC1() << "): nth=" << e->num_threads << "\n"; 
				r = e;
				break;
			}
		}
		if (!written_by_another_flag) {
			element->next = m_table[index];
			m_table[index] = element;
		}
		//cerr << "Thread " << thread_id << ": release (" << zob.getC0() 
		//<< "," << zob.getC1() << "):mutex_index = " << mutex_index << "index=" << index 
		//<< ":tbl=" << m_table[index] << "\n"; 
	}
	if (written_by_another_flag) {
	  //cerr << "Thread " << thread_id << " reaches here\n";
		element->next = m_thread_local[thread_id]->m_free_list;
		m_thread_local[thread_id]->m_free_list = element;
	}
	return r;
}

bool ParallelRbfaooCacheTable::setFirstThread(int thread_id) {
	if (m_first_thread >= 0) 
		return false;
	else {
		scoped_lock_t lock(m_first_thread_lock);
		if (m_first_thread >= 0)
			return false;
		m_first_thread = thread_id;
		return true;
	}
}

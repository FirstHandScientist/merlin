/*
 * RbfaooCacheTable.cpp
 *
 */

#include <assert.h>

#include <algorithm>
#include <iostream>
#include "rbfaoo_cache_table.h"

using namespace std;

RbfaooCacheTable::~RbfaooCacheTable() {
	printWriteStats();
	if (m_table) {
		delete[] m_table;
		delete[] m_entries;
		m_table = NULL;
		m_entries = NULL;
	}
}

void RbfaooCacheTable::printWriteStats() {
	cout << "# of times written in TT = " << m_num_saved << endl;
	cout << "# of nodes newly saved in TT = " << m_num_newly_saved << endl;
}

void RbfaooCacheTable::init(size_t size) {
	if (m_table) {
		delete[] m_table;
		delete[] m_entries;
		m_free_list = NULL;
	}

	// kishi: need to think about how much size is used for context
	m_table_size = size * 1024
			/ (sizeof(RbfaooCacheElement) + sizeof(RbfaooCacheElementPtr));
	m_entries = new RbfaooCacheElement[m_table_size];
	m_table = new RbfaooCacheElementPtr[m_table_size];
	cout << "Num of cache entries = " << m_table_size << "\n";

	fill(&m_table[0], &m_table[m_table_size],
			reinterpret_cast<RbfaooCacheElementPtr>(NULL));
	m_free_list = m_entries;
	for (unsigned int i = 0; i < m_table_size - 1; i++)
		m_entries[i].next = &m_entries[i + 1];
	m_entries[m_table_size - 1].next = NULL;
	m_gc_subtree_size = 0; // Kishi: I made up a number

}

void RbfaooCacheTable::smallTreeGC() {

	const double TABLE_OCCUPACY = 0.70;
	size_t num_stored = m_table_size;
	double occupacy;
	context_t empty_context;

	//cout << "Garbage collecting...\n";
	do {
		unsigned int smallest_subtree = 0x7fffffff;
		for (size_t i = 0; i < m_table_size; i++) {
			RbfaooCacheElementPtr prev, element;

			for (prev = NULL, element = m_table[i]; element != NULL;) {
				if (element->subtree_size <= m_gc_subtree_size) {
					RbfaooCacheElementPtr tmp;

					element->context.swap(empty_context);
					if (prev == NULL) {
						assert(element == m_table[i]);
						m_table[i] = m_table[i]->next;
					} else
						prev->next = element->next;
					num_stored--;
					tmp = element;
					element = element->next;
					tmp->next = m_free_list;
					m_free_list = tmp;
				} else {
					smallest_subtree = min(smallest_subtree,
							element->subtree_size);
					prev = element;
					element = element->next;
				}
			}
		}
		occupacy = static_cast<double>(num_stored) / m_table_size;
#ifdef DEBUG
		cout << "! Garbage collection: Subtree Size = " << std::setw(4) << m_gc_subtree_size
			<< " Occupacy = " << std::setw(4) << occupacy << " % !\n";
#endif
		m_gc_subtree_size = smallest_subtree;
	} while (occupacy >= TABLE_OCCUPACY);
}

RbfaooCacheElementPtr RbfaooCacheTable::read(const Zobrist &zob, int var,
		int or_flag, const context_t &ctxt) {

	size_t index = tableIndex(zob);

	for (RbfaooCacheElement *element = m_table[index]; element != NULL;
			element = element->next) {
		if (zob == element->zobkey && var == element->var
				&& element->or_node_flag == or_flag && element->context == ctxt)
			return element;
	}
	return NULL;
}

RbfaooCacheElementPtr RbfaooCacheTable::write(const Zobrist &zob, int var,
		context_t &ctxt, int or_flag, unsigned int best_index, double result,
		int dn, int solved_flag, unsigned int increased_subtree_size) {
	size_t index = tableIndex(zob);

	m_num_saved++;

	for (RbfaooCacheElement *element = m_table[index]; element != NULL;
			element = element->next) {
		if (zob == element->zobkey && var == element->var
				&& element->or_node_flag == element->or_node_flag
				&& element->context == ctxt) {
			element->best_index = best_index;
			element->result = result;
			element->dn = dn;
			element->solved_flag = solved_flag;
			element->subtree_size += increased_subtree_size;
			return element;
		}
	}

	if (m_free_list == NULL)
		smallTreeGC();

	RbfaooCacheElementPtr element = m_free_list;
	m_num_newly_saved++;
	m_free_list = m_free_list->next;
	element->zobkey = zob;
	element->subtree_size = increased_subtree_size;
	element->next = m_table[index];
	m_table[index] = element;
	element->var = var;
	element->result = result;
	element->dn = dn;
	element->solved_flag = solved_flag;
	element->context.swap(ctxt);
	element->best_index = best_index;
	element->or_node_flag = or_flag;

	return element;
}

// frees all cache table entries
void RbfaooCacheTable::clear() {

	if (m_table && m_entries) {
		fill(&m_table[0], &m_table[m_table_size],
				reinterpret_cast<RbfaooCacheElementPtr>(NULL));
		m_free_list = m_entries;
		for (unsigned int i = 0; i < m_table_size - 1; i++)
			m_entries[i].next = &m_entries[i + 1];
		m_entries[m_table_size - 1].next = NULL;
		m_gc_subtree_size = 0; // Kishi: I made up a number
	}
}



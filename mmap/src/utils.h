/*
 * utils.h
 *
 *  Created on: Mar 29, 2013
 *      Author: radu
 */

#ifndef IBM_ANYTIME_UTILS_H_
#define IBM_ANYTIME_UTILS_H_

#include "base.h"
#include <malloc.h>

using namespace std;

inline bool equals(double d1, double d2) {
	double eps1 = fabs(d1), eps2 = fabs(d2), eps;
	eps = (eps1 > eps2) ? eps1 : eps2;

	if (eps == 0.0)
		return true; //both numbers are 0.0
	//eps hold the minimum distance between the values
	//that will be considered as the numbers are equal
	//considering the magnitude of the numbers
	eps *= DOUBLE_PRECISION;
	return (fabs(d1 - d2) < eps);
}


int memoryusage();
inline int memoryusage() {

	struct mallinfo info;
	info = mallinfo();

	cout << info.arena << '\t' << info.uordblks << '\t' << info.fordblks << '\t'
			<< info.hblkhd << endl;

	return info.arena;

}

void myprint(std::string s);
void myerror(std::string s);
void err_txt(std::string s);

ostream& operator <<(ostream& os, const vector<int>& s);
ostream& operator <<(ostream& os, const vector<size_t>& s);

ostream& operator <<(ostream& os, const vector<signed short>& s);
ostream& operator <<(ostream& os, const vector<signed char>& s);

ostream& operator <<(ostream& os, const set<int>& s);
ostream& operator <<(ostream& os, const set<size_t>& s);

ostream& operator <<(ostream& os, const vector<double>& s);

string str_replace(string& s, const string& x, const string& y);

/*
 * increments the tuple value, up to each entry's limit. Returns false
 * iff no more tuples can be generated
 */
inline bool increaseTuple(size_t& idx, vector<int>& tuple,
		const vector<int>& limit) {
	assert(tuple.size() == limit.size());

	size_t i = tuple.size();
	while (i) {
		++tuple[i - 1];
		if (tuple[i - 1] == limit[i - 1]) {
			tuple[i - 1] = 0;
			--i;
		} else
			break;
	}
	++idx;
	if (i)
		return true;
	else
		return false;
}

/* same as above, but with int* as the tuple (used by mini bucket elimination) */
inline bool increaseTuple(size_t& idx, int* tuple,
		const vector<int>& limit) {

	size_t i = limit.size();
	while (i) {
		++(tuple[i - 1]);
		if (tuple[i - 1] == limit[i - 1]) {
			tuple[i - 1] = 0;
			--i;
		} else
			break;
	}
	++idx;
	if (i)
		return true;
	else
		return false;
}

#ifdef false
/* Convert an int/size_t into a std::string */
inline char* myitoa(size_t x) {
	size_t z = x;
	size_t size = 2;
	while (z/=10) ++size;
	char* s = new char[size];
	snprintf(s,size,"%d",x);
	return s;
//  return std::string(s);
}
inline char* myitoa(int x) {
	int z = x;
	size_t size = (x<0)?3:2;
	while (z/=10) ++size;
	char* s = new char[size];
	snprintf(s,size,"%d",x);
	return s;
//  return std::string(s);
}
#endif

#ifdef false
/* Convert an int/size_t into a std::string */
inline std::string myitoa(size_t x) {
	size_t z = x;
	size_t size = 2;
	while (z/=10) ++size;
	char s[size];
	snprintf(s,size,"%d",x);
	return std::string(s);
}
inline std::string myitoa(int x) {
	int z = x;
	size_t size = (x<0)?3:2;
	while (z/=10) ++size;
	char s[size];
	snprintf(s,size,"%d",x);
	return std::string(s);
}
#endif

/*
 * returns true iff the intersection of X and Y is empty
 */
inline bool intersectionEmpty(const set<int>& a, const set<int>& b) {
	set<int>::iterator ita = a.begin();
	set<int>::iterator itb = b.begin();
	while (ita != a.end() && itb != b.end()) {
		if (*ita == *itb) {
			return false;
		} else if (*ita < *itb) {
			++ita;
		} else { // *ita > *itb
			++itb;
		}
	}
	return true;
}

/*
 * computes the intersection of two set<int>
 */
inline set<int> intersection(const set<int>& a, const set<int>& b) {
	set<int> s;
	set<int>::iterator ita = a.begin();
	set<int>::iterator itb = b.begin();
	while (ita != a.end() && itb != b.end()) {
		if (*ita == *itb) {
			s.insert(*ita);
		} else if (*ita < *itb) {
			++ita;
		} else { // *ita > *itb
			++itb;
		}
	}
	return s;
}

/*
 * returns the size of (set A without set B)
 */
inline int setminusSize(const set<int>& a, const set<int>& b) {
	set<int>::iterator ita = a.begin();
	set<int>::iterator itb = b.begin();
	int s = 0;
	while (ita != a.end() && itb != b.end()) {
		if (*ita == *itb) {
			++ita;
			++itb;
		} else if (*ita < *itb) {
			++s;
			++ita;
		} else { // *ita > *itb
			++itb;
		}
	}
	while (ita != a.end()) {
		++s;
		++ita;
	}
	return s;
}

/* hex output of values (e.g., doubles) */
template<typename _T>
void print_hex(const _T* d) {
	const unsigned char* ar = (const unsigned char*) d;
	for (size_t i = 0; i < sizeof(_T); ++i)
		printf("%X", ar[i]);
}


#endif /* UTILS_H_ */

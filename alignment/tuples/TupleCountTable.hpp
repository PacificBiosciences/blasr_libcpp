#ifndef _BLASR_TUPLE_COUNT_TABLE_HPP_
#define _BLASR_TUPLE_COUNT_TABLE_HPP_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <assert.h>
#include "tuples/TupleMetrics.hpp"
using namespace std;

template<typename TSequence, typename TTuple>
class TupleCountTable {
public:
	int *countTable;
	TupleData countTableLength;
	int nTuples;
	TupleMetrics tm;
	bool deleteStructures;
	void InitCountTable(TupleMetrics &ptm);

	TupleCountTable();
   	~TupleCountTable();
	void Free();

	void IncrementCount(TTuple &tuple);
	void AddSequenceTupleCountsLR(TSequence &seq);
	void Write(ofstream &out);
	void Read(ifstream &in);
};

#include "tuples/TupleCountTableImpl.hpp"

#endif

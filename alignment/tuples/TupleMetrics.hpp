#ifndef _BLASR_TUPLE_METRICS_HPP_
#define _BLASR_TUPLE_METRICS_HPP_

#include "Types.h"

class TupleMetrics {
public:
    int tupleSize; // Number of bits in a tuple.
    TupleData tupleMask;

    TupleMetrics();

    void InitializeMask();

    void Initialize(int pTupleSize);
};


#endif //_BLASR_TUPLE_METRICS_HPP_

#ifndef _BLASR_BASE_TUPLE_HPP_
#define _BLASR_BASE_TUPLE_HPP_

#include "tuples/TupleMetrics.hpp"

class BaseTuple {
public:
    TupleData tuple;

    TupleData HashPowerOfFour(int nBits, TupleMetrics &tm);

    bool operator<(const BaseTuple &rhs) const; 

    bool operator==(const BaseTuple &rhs) const; 

    bool operator!=(const BaseTuple &rhs) const;

    BaseTuple & ShiftLeft(TupleMetrics &tm, int shift=1); 

    BaseTuple & ShiftRight(int shift=1); 

    BaseTuple & Append(TupleData val, TupleMetrics &tm, int nBits); 

    operator TupleData() const;
};

#endif

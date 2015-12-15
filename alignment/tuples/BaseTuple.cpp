#include <stdint.h>
#include "Types.h"
#include "TupleMask.h"
#include "TupleMetrics.hpp"
#include "BaseTuple.hpp"

TupleData BaseTuple::HashPowerOfFour(int nBits, TupleMetrics &tm) {
    //
    // When the hash can fit inside the entire tuple, just return the
    // tuple.
    //
    if (tm.tupleSize > nBits) {
        return tuple;
    }
    else {
        return ((tuple & TupleMask[nBits]) + (tuple % 1063)) % (1 << (nBits*2));
    }
}

bool BaseTuple::operator<(const BaseTuple &rhs) const {
    return tuple < rhs.tuple;
}

bool BaseTuple::operator==(const BaseTuple &rhs) const {
    return tuple == rhs.tuple;
}

bool BaseTuple::operator!= (const BaseTuple &rhs) const {
    return tuple != rhs.tuple;
}

BaseTuple & BaseTuple::ShiftLeft(TupleMetrics &tm, int shift) {
    tuple = tuple << shift;
    tuple = tuple & tm.tupleMask;
    return *this;
}

BaseTuple & BaseTuple::ShiftRight(int shift) {
    tuple = tuple >> shift;
    return *this;
}

BaseTuple & BaseTuple::Append(TupleData val, TupleMetrics &tm, int nBits) {
    tuple = tuple << nBits;
    tuple = tuple & tm.tupleMask;
    tuple = tuple + val;
    return *this;
}

BaseTuple::operator TupleData() const{
    return tuple;
}


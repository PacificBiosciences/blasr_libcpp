#include <stdint.h>
#include "TupleMetrics.hpp"
#include "TupleMask.h"

TupleMetrics::TupleMetrics()
    : tupleSize(0)
    , tupleMask(0)
{ }

void TupleMetrics::InitializeMask() {
    tupleMask = TupleMask[tupleSize];
}

void TupleMetrics::Initialize(int pTupleSize) {
    tupleSize = pTupleSize;
    InitializeMask();
}

#include "DNATuple.hpp"

DNATuple::DNATuple()
    : pos(0) 
{ }

DNATuple::DNATuple(const DNATuple & rhs)
    : BaseTuple(rhs)
    , pos(rhs.pos)
{ }

int DNATuple::FromStringRL(Nucleotide *strPtr, TupleMetrics &tm) {

    //
    // Tuples are created with the right-most character
    // in the most significant bit in the sequence.
    //
    DNASequence tmpSeq;
    tmpSeq.seq = strPtr;
    tmpSeq.length = tm.tupleSize;
    if (!OnlyACTG(tmpSeq))
        return 0;

    if (tm.tupleSize == 0)
        return 1;

    tuple = 0;
    Nucleotide *p;
    for (p = strPtr + tm.tupleSize - 1; p > strPtr; p--) {
        tuple += TwoBit[*p];
        tuple <<=2;
    }
    //
    // The tuple size is guaranteed to be at least 
    // 1, so it's safe to add the last value.
    // This cannot be in the previous loop since
    // the shift shouldn't happen.
    tuple += TwoBit[*p];
    return 1;
}


int DNATuple::MakeRC(DNATuple &dest, TupleMetrics &tm) {
    int i;
    TupleData tempTuple = tuple;
    dest.tuple = 0;
    TupleData mask = 0x3;
    if (tm.tupleSize == 0)
        return 0;
    for (i = 0; i < tm.tupleSize - 1; i++ ) {
        dest.tuple += (~tempTuple & mask);
        tempTuple >>= 2;
        dest.tuple <<= 2;
    }
    dest.tuple += (~tempTuple & mask);
    return 1;
}

std::string DNATuple::ToString(TupleMetrics &tm) {
    int i;
    std::string s;
    TupleData tempTuple = tuple;
    for (i = 0;i < tm.tupleSize;i++) {
        s.insert(s.begin(), TwoBitToAscii[tempTuple & 3]);
        tempTuple = tempTuple >> 2;
    }
    return s;
}

bool 
CompareByTuple::operator()(const DNATuple &lhs, const DNATuple &rhs) 
const {
    return lhs.tuple < rhs.tuple;
}

CountedDNATuple::CountedDNATuple(const CountedDNATuple & rhs) 
    : DNATuple(rhs)
    , count(rhs.count)
{ }

PositionDNATuple::PositionDNATuple()
    : DNATuple() 
{ 
    pos = static_cast<DNALength>(-1);
}

PositionDNATuple::PositionDNATuple(const PositionDNATuple & rhs)
    : DNATuple(rhs)
{ }

PositionDNATuple::PositionDNATuple(const PositionDNATuple &tupleP, const DNALength posP)
    : DNATuple(tupleP)
{
    pos   = posP;
}

int PositionDNATuple::operator==(const DNATuple &pTuple) const {
    return tuple == pTuple.tuple;
}

int 
OrderPositionDNATuplesByPosition::operator()(
    const PositionDNATuple &lhs, const PositionDNATuple &rhs) 
const {
    return lhs.pos < rhs.pos;
}

int
OrderPositionDNATuplesByTuple::operator()(
    const PositionDNATuple &lhs, const PositionDNATuple &rhs) 
const {
    return lhs.tuple < rhs.tuple;
}

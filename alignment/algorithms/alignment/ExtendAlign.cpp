#include "ExtendAlign.hpp"

RCToIndex::RCToIndex() {
    qStart = 0; tStart = 0;
    band = middleCol = nCols = 0;
}

int RCToIndex::operator()(int r, int c, int &index) {
    //
    // First do some error checking on the row and column to see if it
    // is within the band.
    //
    if (r < qStart) { return 0; }
    if (c < tStart) { return 0; }
    r -= qStart;
    c -= tStart;
    if (std::abs(r-c) > band) { return 0; } // outside band range.
    if (c < 0) { return 0; }

    if (middleCol - (r - c) >= nCols) { return  0; }
    index = (r*nCols) + (middleCol - (r - c));
    return 1;
}

int BaseIndex::QNotAtSeqBoundary(int q) {
    return q != queryAlignLength;
}
int BaseIndex::TNotAtSeqBoundary(int t) {
    return t != refAlignLength;
}	
int BaseIndex::QAlignLength() {
    return queryAlignLength;
}

int BaseIndex::TAlignLength() {
    return refAlignLength;
}


int ForwardIndex::QuerySeqPos(int q) {
    return queryPos + q;
}

int ForwardIndex::RefSeqPos(int t) {
    return refPos + t;
}

int ForwardIndex::GetQueryStartPos(int startQ, int endQ) {
    (void)(endQ);
    return queryPos + startQ + 1;
}

int ForwardIndex::GetRefStartPos(int startT, int endT) {
    (void)(endT);
    return refPos + startT + 1;
}

void ForwardIndex::OrderArrowVector(std::vector<Arrow> &mat) {
    reverse(mat.begin(), mat.end());
}


int ReverseIndex::QuerySeqPos(int q) {
    return queryPos - q;
}

int ReverseIndex::RefSeqPos(int t) {
    return refPos - t;
}

int ReverseIndex::GetQueryStartPos(int startQ, int endQ) {
    (void)(startQ);
    return queryPos - (endQ-1);
}

int ReverseIndex::GetRefStartPos(int startT, int endT) {
    (void)(startT);
    return refPos - (endT-1);
}

void ReverseIndex::OrderArrowVector(std::vector<Arrow> &mat) {
    (void)(mat);
}

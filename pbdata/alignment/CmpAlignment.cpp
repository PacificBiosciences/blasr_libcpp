#include "CmpAlignment.hpp"

std::map<std::string,int> CmpAlignmentBase::columnNameToIndex;

unsigned int* CmpAlignmentBase::GetAlignmentIndex() {
    return &alignmentIndex[0];
}
int CmpAlignmentBase::GetAlignmentIndexSize() {
    return alignmentIndex.size();
}
unsigned int CmpAlignmentBase::GetAlignedStrand() {
    return LookupColumnValue("AlignedStrand");
}

unsigned int CmpAlignmentBase::GetRCRefStrand() {
    return LookupColumnValue("RCRefStrand");
}

// synonym 
unsigned int CmpAlignmentBase::GetTStrand() {
    return GetRCRefStrand();
}

bool CmpAlignmentBase::GetX(int &xp) {
    if (alignmentIndex.size() > 0) {
        xp = alignmentIndex[columnNameToIndex["x"]];
        return true;
    }
    else {
        xp = -1;
        return false;
    }
}

unsigned int CmpAlignmentBase::GetAlignmentId() {
    return LookupColumnValue("AlnID");
}

unsigned int CmpAlignmentBase::GetX() {
    return LookupColumnValue("x");
}

unsigned int CmpAlignmentBase::GetY() {
    return LookupColumnValue("y");
}

unsigned int CmpAlignmentBase::GetMovieId() {
    return LookupColumnValue("MovieId");
}

unsigned int CmpAlignmentBase::GetAlnGroupId() {
    return LookupColumnValue("AlnGroupId");
}

unsigned int CmpAlignmentBase::GetReadGroupId() {
    return LookupColumnValue("ReadGroupId");
}


unsigned int CmpAlignmentBase::LookupColumnValue(const char * columnName) {
    if (columnNameToIndex.find(columnName) != columnNameToIndex.end()) {
        int columnIndex = columnNameToIndex[columnName];
        return alignmentIndex[columnIndex];
    }
    else {
        std::cout << "ERROR, For now cmp files must contain a column " << columnName << std::endl;
        std::cout << "size of columnNameToIndex: " << columnNameToIndex.size() << std::endl;
        assert(0);
    }
}		


void CmpAlignmentBase::InitializeColumnNameToIndex(std::vector<std::string> &columnNames) {
    for (size_t i = 0; i < columnNames.size(); i++ ){
        columnNameToIndex[columnNames[i]] = i;
    }
}


unsigned int CmpAlignmentBase::GetHoleNumber() {
    return LookupColumnValue("HoleNumber");
}

unsigned int CmpAlignmentBase::GetRefGroupId() {
    return LookupColumnValue("RefGroupId");
}

unsigned int CmpAlignmentBase::GetRefSeqId() {
    return LookupColumnValue("RefSeqId");
}

unsigned int CmpAlignmentBase::GetOffsetBegin() {
    return LookupColumnValue("offset_begin");
}

unsigned int CmpAlignmentBase::GetOffsetEnd() {
    return LookupColumnValue("offset_end");
}

unsigned int CmpAlignmentBase::GetQueryStart() {
    return LookupColumnValue("rStart");
}

unsigned int CmpAlignmentBase::GetQueryEnd() {
    return LookupColumnValue("rEnd");
}

unsigned int CmpAlignmentBase::GetRefStart() {
    return LookupColumnValue("tStart");
}

unsigned int CmpAlignmentBase::GetRefEnd() {
    return LookupColumnValue("tEnd");
}

unsigned int CmpAlignmentBase::GetNMatch() {
    return LookupColumnValue("nM");
}

unsigned int CmpAlignmentBase::GetNMismatch() {
    return LookupColumnValue("nMM");
}

unsigned int CmpAlignmentBase::GetNInsertions() {
    return LookupColumnValue("nIns");
}

unsigned int CmpAlignmentBase::GetNDeletions() {
    return LookupColumnValue("nDel");
}

unsigned int CmpAlignmentBase::GetMapQV() {
    return LookupColumnValue("MapQV");
}

unsigned int CmpAlignmentBase::GetSubreadId() {
    return LookupColumnValue("SubreadId");
}

unsigned int CmpAlignmentBase::GetStrobeNumber() {
    return LookupColumnValue("StrobeNumber");
}

unsigned int CmpAlignmentBase::GetSetNumber() {
    return LookupColumnValue("SetNumber");
}

CmpAlignmentBase::CmpAlignmentBase(PlatformId platformIdP) {
    platformId = platformIdP;
}

void CmpAlignmentBase::SetPlatformId(PlatformId platformIdP) {
    platformId = platformIdP;
}

// CmpAlignment

CmpAlignment::CmpAlignment(PlatformId pid) : CmpAlignmentBase(pid) {
}

void CmpAlignment::StoreAlignmentIndex(unsigned int *alignmentIndexPtr, 
    DSLength alignmentIndexLength) {

    alignmentIndex.clear();
    alignmentIndex.insert(alignmentIndex.begin(), &alignmentIndexPtr[0], 
        &alignmentIndexPtr[alignmentIndexLength]);
}

void CmpAlignment::StoreAlignmentArray(unsigned char* alignmentArrayPtr, 
    DSLength alignmentArrayLength) {
    alignmentArray.resize(alignmentArrayLength);
    for (DSLength a = 0; a < alignmentArrayLength; a++ ){
        alignmentArray[a] = alignmentArrayPtr[a];
    }
}

CmpAlignment &CmpAlignment::operator=(const CmpAlignment &rhs) {
    // deep copy the alignment index
    alignmentIndex.resize(rhs.alignmentIndex.size());
    copy(rhs.alignmentIndex.begin(), rhs.alignmentIndex.end(), alignmentIndex.begin());
    // deep copy the alignment array
    alignmentArray.resize(rhs.alignmentIndex.size());
    copy(rhs.alignmentArray.begin(), rhs.alignmentArray.end(), alignmentArray.begin());
    // copy fields
    Z = rhs.Z;
    index = rhs.index; readGroupId = rhs.readGroupId; movieId = rhs.movieId;
    refSeqId = rhs.refSeqId;
    expId = rhs.expId; runId = rhs.runId; panel = rhs.panel;
    x = rhs.x; y = rhs.y;
    rcRefStrand = rhs.rcRefStrand;
    holeNumber = rhs.holeNumber;
    offsetBegin = rhs.offsetBegin; offsetEnd = rhs.offsetEnd;
    setNumber = rhs.setNumber, strobeNumber = rhs.strobeNumber, mapQV = rhs.mapQV, nBackRead = rhs.nBackRead, nReadOverlap = rhs.nReadOverlap;
    subreadId = rhs.subreadId;
    nMatch = rhs.nMatch, nMismatch = rhs.nMismatch, nIns = rhs.nIns, nDel = rhs.nDel;
    return *this;
}

int CmpAlignment::operator<(const CmpAlignment &rhs) const {
    if (alignmentArray[1] == rhs.alignmentArray[1]) {
        if (alignmentArray[2] == rhs.alignmentArray[2]) {
            if (alignmentArray[10] == rhs.alignmentArray[10]) {
                return (alignmentArray[4] < rhs.alignmentArray[4]);
            }
            else {
                return alignmentArray[10] < rhs.alignmentArray[10];
            }
        }
        else {
            return alignmentArray[2] < rhs.alignmentArray[2];
        }
    }
    else {
        return alignmentArray[1] < rhs.alignmentArray[1];
    }
}

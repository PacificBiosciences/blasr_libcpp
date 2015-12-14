#ifndef DATASTRUCTURES_ALIGNMENT_CMP_ALIGNMENT_H_
#define DATASTRUCTURES_ALIGNMENT_CMP_ALIGNMENT_H_
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <assert.h>
#include <algorithm>

#include "Types.h"
#include "Enumerations.h"

class CmpAlignmentBase  {
public:
    //
    // For use in referencing alignment sets. TODO: subclass.
    //
    PlatformId platformId;
    int Z;
    unsigned int index, readGroupId, movieId, refSeqId; 
    unsigned int expId, runId, panel;
    unsigned int x, y;
    unsigned int rcRefStrand;
    unsigned int holeNumber;
    unsigned int offsetBegin, offsetEnd;
    unsigned int setNumber, strobeNumber, mapQV, nBackRead, nReadOverlap;
    unsigned int subreadId;
    unsigned int nMatch, nMismatch, nIns, nDel;
    std::vector<unsigned char> alignmentArray;
    std::vector<unsigned int> alignmentIndex;
    static std::map<std::string,int> columnNameToIndex;
    static bool initializedColumnNameToIndex;
    std::map<std::string, std::vector<UChar> > fields;

    unsigned int *GetAlignmentIndex(); 

    int GetAlignmentIndexSize(); 

    unsigned int GetAlignedStrand(); 

    unsigned int GetRCRefStrand(); 

    // synonym 
    unsigned int GetTStrand(); 

    bool GetX(int &xp); 

    unsigned int GetAlignmentId(); 

    unsigned int GetX(); 

    unsigned int GetY(); 

    unsigned int GetMovieId(); 

    unsigned int GetAlnGroupId(); 

    unsigned int GetReadGroupId(); 

    unsigned int LookupColumnValue(const char * columnName); 

    void InitializeColumnNameToIndex(std::vector<std::string> &columnNames); 

    unsigned int GetHoleNumber(); 

    unsigned int GetRefGroupId(); 

    unsigned int GetRefSeqId(); 

    unsigned int GetOffsetBegin(); 

    unsigned int GetOffsetEnd(); 

    unsigned int GetQueryStart();

    unsigned int GetQueryEnd(); 

    unsigned int GetRefStart(); 

    unsigned int GetRefEnd();

    unsigned int GetNMatch(); 

    unsigned int GetNMismatch();

    unsigned int GetNInsertions();

    unsigned int GetNDeletions(); 

    unsigned int GetMapQV(); 

    unsigned int GetSubreadId(); 

    unsigned int GetStrobeNumber(); 

    unsigned int GetSetNumber(); 

    CmpAlignmentBase(PlatformId platformIdP=Springfield); 

    void SetPlatformId(PlatformId platformIdP); 
};

class CmpAlignment : public CmpAlignmentBase {
public:
    int qStrand, tStrand;
    DNALength qStart, qLength;
    DNALength tStart, tLength;
    //
    // Default constructor just calls the base constructor to initialize platoformType
    CmpAlignment(PlatformId pid=Springfield);

    void StoreAlignmentIndex(unsigned int *alignmentIndexPtr, DSLength alignmentIndexLength);

    void StoreAlignmentArray(unsigned char* alignmentArrayPtr, DSLength alignmentArrayLength); 

    template<typename T_Field>
    void StoreField(std::string fieldName, T_Field* fieldValues, DSLength length); 

    CmpAlignment &operator=(const CmpAlignment &rhs); 

    int operator<(const CmpAlignment &rhs) const; 

};

#include "CmpAlignmentImpl.hpp"
#endif

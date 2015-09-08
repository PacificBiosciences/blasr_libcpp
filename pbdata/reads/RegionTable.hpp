// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Mark Chaisson

#ifndef _BLASR_REGION_TABLE_HPP_
#define _BLASR_REGION_TABLE_HPP_

#include <cassert>
#include <cstring>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <ostream>
#include "Types.h"
#include "Enumerations.h"
#include "PacBioDefs.h"


class HDFRegionTableReader;
class HDFRegionTableWriter;

class RegionTypeMap {
public:
    static std::string ToString(RegionType rt);

    static RegionType ToRegionType(const std::string & str);

private:
    static const std::map<RegionType, std::string> RegionTypeToString ;

    static const std::map<std::string, RegionType> StringToRegionType;
};

class RegionAnnotation {
friend class HDFRegionTableReader;
friend class HDFRegionTableWriter;
friend class HDFRegionsWriter;

public:
    static const int HOLENUMBERCOL = 0;
    static const int REGIONTYPEINDEXCOL = 1;
    static const int REGIONSTARTCOL = 2;
    static const int REGIONENDCOL = 3;
    static const int REGIONSCORECOL = 4;
    static const int NCOLS=5;

    int row[NCOLS];

public:
    inline RegionAnnotation(UInt holeNumber = 0, 
                            int typeIndex = 0,
                            int start = 0, int end = 0,
                            int score = -1);

    inline bool operator<(const RegionAnnotation &rhs) const;

    inline bool operator<(int holeNumber) const;

    inline RegionAnnotation& operator=(const RegionAnnotation &rhs); 

    inline int GetHoleNumber(void) const;

    inline RegionAnnotation & SetHoleNumber(int holeNumber); 

    inline int GetTypeIndex(void) const; 

    inline std::string GetTypeString(const std::vector<RegionType> & types) const; 

    inline RegionAnnotation & SetTypeIndex(int typeIndex); 

    inline int GetStart(void) const; 

    inline RegionAnnotation & SetStart(int start); 

    inline int GetEnd(void) const; 

    inline RegionAnnotation & SetEnd(int end); 

    inline int GetScore(void) const; 

    inline RegionAnnotation & SetScore(int score); 

public:
    friend std::ostream & operator << (std::ostream & os, const RegionAnnotation& ra);
};

class RegionTable {
public:
    std::vector<RegionAnnotation> table;
    std::vector<std::string> columnNames;
    std::vector<std::string> regionTypes;
    std::vector<std::string> regionDescriptions;
    std::vector<std::string> regionSources;
    std::vector<RegionType> regionTypeEnums;

    // Return default region types used in a region table 
    // Note that the ORDER of region types does matter.
    static std::vector<RegionType> DefaultRegionTypes(void);

    int LookupRegionsByHoleNumber(int holeNumber, int &low, int &high) const; 

    //
    // Define a bunch of accessor functions.
    //

    //
    // Different region tables have different ways of encoding regions.
    // This maps from the way they are encoded in the rgn table to a
    // standard encoding.
    //
    RegionType GetType(int regionIndex) const; 

    int GetStart(const int regionIndex) const; 

    void SetStart(int regionIndex, int start); 

    int GetEnd(const int regionIndex) const; 

    void SetEnd(int regionIndex, int end); 

    int GetHoleNumber(int regionIndex) const;

    void SetHoleNumber(int regionIndex, int holeNumber); 

    int GetScore(int regionIndex) const; 

    void SetScore(int regionIndex, int score); 

    void SortTableByHoleNumber(); 

    void Reset(); 

    void CreateDefaultAttributes(); 
};


inline
RegionAnnotation::RegionAnnotation(UInt holeNumber, 
        int typeIndex, int start, int end, int score) {
    SetHoleNumber(static_cast<int>(holeNumber));
    SetTypeIndex(typeIndex);
    SetStart(start);
    SetEnd(end);
    SetScore(score);
}

inline
bool RegionAnnotation::operator<(const RegionAnnotation &rhs) const
{ 
    if (GetHoleNumber() == rhs.GetHoleNumber()) 
        return GetStart() < rhs.GetStart();
    else
        return GetHoleNumber() < rhs.GetHoleNumber();
}

inline
bool RegionAnnotation::operator<(int holeNumber) const
{ return GetHoleNumber() < holeNumber; }


inline
RegionAnnotation& RegionAnnotation::operator=(const RegionAnnotation &rhs) {
    memcpy(row, rhs.row, sizeof(int)*NCOLS);
    return *this;
}

inline
int RegionAnnotation::GetHoleNumber(void) const {
    return row[HOLENUMBERCOL];
}

inline
RegionAnnotation & RegionAnnotation::SetHoleNumber(int holeNumber) {
    row[HOLENUMBERCOL] = holeNumber;
    return *this;
}

inline
int RegionAnnotation::GetTypeIndex(void) const {
    return row[REGIONTYPEINDEXCOL];
}

inline std::string RegionAnnotation::GetTypeString(const std::vector<RegionType> & typesTable) const {
    assert(GetTypeIndex() >= 0 and GetTypeIndex() < static_cast<int>(typesTable.size()));
    return RegionTypeMap::ToString(typesTable[GetTypeIndex()]);
}

inline
RegionAnnotation & RegionAnnotation::SetTypeIndex(int regionTypeIndex) {
    row[REGIONTYPEINDEXCOL] = regionTypeIndex;
    return *this;
}

inline
int RegionAnnotation::GetStart(void) const {
    return row[REGIONSTARTCOL];
}

inline
RegionAnnotation & RegionAnnotation::SetStart(int start) {
    row[REGIONSTARTCOL] = start;
    return *this;
}

inline
int RegionAnnotation::GetEnd(void) const {
    return row[REGIONENDCOL];
}

inline
RegionAnnotation & RegionAnnotation::SetEnd(int end) {
    row[REGIONENDCOL] = end;
    return *this;
}

inline
int RegionAnnotation::GetScore(void) const {
    return row[REGIONSCORECOL];
}

inline
RegionAnnotation & RegionAnnotation::SetScore(int score) {
    row[REGIONSCORECOL] = score;
    return *this;
}
#endif // _BLASR_REGION_TABLE_HPP_

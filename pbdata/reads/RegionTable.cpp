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


#include <algorithm>
#include <iostream>
#include <ostream>
#include "RegionTable.hpp"

using namespace std;

std::string RegionTypeMap::ToString(RegionType rt) {
    assert(RegionTypeToString.find(rt) != RegionTypeToString.end());
    return RegionTypeToString.find(rt)->second;
}

RegionType RegionTypeMap::ToRegionType(const std::string & str) {
    if (StringToRegionType.find(str) == StringToRegionType.end()) {
        std::cout << "Unsupported RegionType " << str << std::endl;
        assert(false);
    }
    return StringToRegionType.find(str)->second;
}

const std::map<RegionType, std::string> RegionTypeMap::RegionTypeToString = {
    {Adapter,  "Adapter"},
    {Insert,   "Insert"},
    {HQRegion, "HQRegion"},
    {BarCode,  "Barcode"}
};

const std::map<std::string, RegionType> RegionTypeMap::StringToRegionType = {
    {"Adapter",  Adapter},
    {"Insert",   Insert},
    {"HQRegion", HQRegion}, 
    {"Barcode",  BarCode},
};

std::ostream & operator << (std::ostream & os, const RegionAnnotation& ra) {
    os << "ZMW " << ra.GetHoleNumber() 
       << ", region type index " << ra.GetTypeIndex() 
       << " [" << ra.GetStart() 
       << ", " << ra.GetEnd()
       << "), " << ra.GetScore();
    return os;
}

int RegionTable::LookupRegionsByHoleNumber(int holeNumber, int &low, int &high) const {
    std::vector<RegionAnnotation>::const_iterator lowIt, highIt;
    lowIt  = std::lower_bound(table.begin(), table.end(), holeNumber);
    highIt = std::lower_bound(table.begin(), table.end(), holeNumber+1);
    low =  lowIt - table.begin();
    high = highIt - table.begin();
    return high-low;
}

//
// Define a bunch of accessor functions.
//

//
// Different region tables have different ways of encoding regions.
// This maps from the way they are encoded in the rgn table to a
// standard encoding.
//

RegionType RegionTable::GetType(int regionIndex) const {
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    return regionTypeEnums[table[regionIndex].GetTypeIndex()];
}

int RegionTable::GetStart(const int regionIndex) const {
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    return table[regionIndex].GetStart();
}

void RegionTable::SetStart(int regionIndex, int start) {
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    table[regionIndex].SetStart(start);
}

int RegionTable::GetEnd(const int regionIndex) const {
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    return table[regionIndex].GetEnd();
}

void RegionTable::SetEnd(int regionIndex, int end) {
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    table[regionIndex].SetEnd(end);
}

int RegionTable::GetHoleNumber(int regionIndex) const{
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    return table[regionIndex].GetHoleNumber();
}

void RegionTable::SetHoleNumber(int regionIndex, int holeNumber) {
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    table[regionIndex].SetHoleNumber(holeNumber);
}

int RegionTable::GetScore(int regionIndex) const{
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    return table[regionIndex].GetScore();//row[RegionAnnotationColumn::RegionScore];
}

void RegionTable::SetScore(int regionIndex, int score) {
    assert(regionIndex < table.size());
    assert(regionIndex >= 0);
    table[regionIndex].SetScore(score);//.row[RegionAnnotationColumn::RegionScore] = score;
}

void RegionTable::SortTableByHoleNumber() {
    std::stable_sort(table.begin(), table.end());
}

std::vector<RegionType> RegionTable::DefaultRegionTypes(void) {
    std::vector<RegionType> ret;
    for (std::string regionTypeString: PacBio::AttributeValues::Regions::regiontypes) {
        ret.push_back(RegionTypeMap::ToRegionType(regionTypeString));
    }
    return ret;
}

void RegionTable::Reset() {
    table.clear();
    columnNames.clear();
    regionTypes.clear();
    regionDescriptions.clear();
    regionSources.clear();
    regionTypeEnums.clear();
}

void RegionTable::CreateDefaultAttributes() {
    columnNames        = PacBio::AttributeValues::Regions::columnnames;
    regionTypes        = PacBio::AttributeValues::Regions::regiontypes;
    regionDescriptions = PacBio::AttributeValues::Regions::regiondescriptions;
    regionSources      = PacBio::AttributeValues::Regions::regionsources;
    regionTypeEnums    = DefaultRegionTypes();
}

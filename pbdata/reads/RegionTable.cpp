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
// Modified by: Yuan Li


#include <algorithm>
#include <iostream>
#include <ostream>
#include "RegionTable.hpp"

using namespace std;


RegionTable & RegionTable::Reset() {
    map_.clear();
    columnNames.clear();
    regionTypes.clear();
    regionDescriptions.clear();
    regionSources.clear();
    regionTypeEnums.clear();
    return *this;
}

std::vector<RegionType> RegionTable::RegionTypeEnums(void) const 
{ return regionTypeEnums; }

std::vector<std::string> RegionTable::RegionTypes(void) const
{ return regionTypes; }

std::vector<std::string> RegionTable::ColumnNames(void) const
{ return columnNames; }

std::vector<std::string> RegionTable::RegionDescriptions(void) const
{ return regionDescriptions; }

std::vector<std::string> RegionTable::RegionSources(void) const
{ return regionSources;}

RegionTable & RegionTable::ConstructTable(std::vector<RegionAnnotation> & table,
                                          const std::vector<std::string> & regionTypeStrs) {
    RegionTypes(regionTypeStrs); //< Set both regionTypes and regionTypeEnums. 

    // Must sort region annotations by HoleNumber, RegionTypeIndex, Start, End, and Score
    std::sort(table.begin(), table.end(), compare_region_annotation_by_type);

    // Construct map_<holeNumber, RegionAnnotations>
    if (table.size() > 0) {
        UInt pre_hn = table[0].GetHoleNumber();
        auto itBegin = table.begin();
        for (auto it = table.begin(); it != table.end(); it++) {
            if (it->GetHoleNumber() > pre_hn) {
                map_.insert(std::pair<UInt, RegionAnnotations>(pre_hn, 
                            RegionAnnotations(pre_hn, 
                                              std::vector<RegionAnnotation>(itBegin, it),
                                              regionTypeEnums)));
                pre_hn = it->GetHoleNumber();
                itBegin = it;
            }
        }

        map_.insert(std::pair<UInt, RegionAnnotations>(pre_hn, 
                    RegionAnnotations(pre_hn, 
                                      std::vector<RegionAnnotation>(itBegin, table.end()),
                                      regionTypeEnums)));
    }
}

std::vector<RegionType> RegionTable::DefaultRegionTypes(void) {
    std::vector<RegionType> ret;
    for (std::string regionTypeString: PacBio::AttributeValues::Regions::regiontypes) {
        ret.push_back(RegionTypeMap::ToRegionType(regionTypeString));
    }
    return ret;
}

RegionTable & RegionTable::RegionTypes(const std::vector<std::string> & regionTypeStrs) {
    regionTypes = regionTypeStrs;
    for (std::string regionTypeString: regionTypeStrs) {
        regionTypeEnums.push_back(RegionTypeMap::ToRegionType(regionTypeString));
    }
    return *this;
}

RegionTable & RegionTable::ColumnNames(const std::vector<std::string> & in)
{ columnNames = in; return *this; }

RegionTable & RegionTable::RegionDescriptions(const std::vector<std::string> & in)
{ regionDescriptions = in; return *this; }

RegionTable & RegionTable::RegionSources(const std::vector<std::string> & in)
{ regionSources = in; return *this; }

bool RegionTable::HasHoleNumber(const UInt holeNumber) const {
    return (map_.find(holeNumber) != map_.end());
}

RegionAnnotations RegionTable::operator [] (const UInt holeNumber) const {
    // Must check whether a zmw exists or not first.
    assert (HasHoleNumber(holeNumber) or 
            false == "Could not find zmw in region table.");
    return map_.find(holeNumber)->second;
}

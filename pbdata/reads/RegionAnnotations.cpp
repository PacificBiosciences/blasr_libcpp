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

// Author: Yuan Li


#include <algorithm>
#include "RegionAnnotations.hpp"


RegionAnnotations::RegionAnnotations(const UInt holeNumber, 
                                     const std::vector<RegionAnnotation> & annotations)
    : holeNumber_(holeNumber)
    , table_(annotations)
{
    for (auto annotation: annotations) {
        // Only allow RegionAnnotations of a single ZMW
        if (holeNumber_ != annotation.GetHoleNumber()) {
            assert(false == "RegionAnnotations must contain regions from a single ZMW");
        }
    }

    std::sort(table_.begin(), table_.end(), compare_region_annotation_by_type);
}

RegionAnnotations& RegionAnnotations::operator=(const RegionAnnotations & other) {
    holeNumber_ = other.HoleNumber();
    table_ = other.table_;
    return *this;
}

UInt RegionAnnotations::HoleNumber(void) const {
    return holeNumber_;
}

std::vector<RegionAnnotation> 
RegionAnnotations::RegionAnnotationsOfType(RegionType type, 
                                           const std::vector<std::string> & typeStrs) const {

    std::vector<RegionAnnotation> ret;
    int typeIndex = RegionTypeMap::ToIndex(type, typeStrs);
    if (typeIndex >= 0) {
        std::copy_if(table_.begin(), table_.end(), ret.begin(), 
                     [&](const RegionAnnotation & a)->bool{return a.GetTypeIndex()==typeIndex;});
        sort(ret.begin(), ret.end());
    }
    return ret;
}

std::vector<RegionAnnotation> RegionAnnotations::Adapters(const std::vector<std::string> & typeStrs) const {
    return RegionAnnotationsOfType(Adapter, typeStrs);
}

bool RegionAnnotations::HasHQRegion(const std::vector<std::string> & typeStrs) const {
    return (HQRegions(typeStrs).size() >= 1);
}

std::vector<RegionAnnotation>
RegionAnnotations::HQRegions(const std::vector<std::string> & typeStrs) const {
    return RegionAnnotationsOfType(HQRegion, typeStrs);
}

RegionAnnotation 
RegionAnnotations::TheHQRegion(const std::vector<std::string> & typeStrs) const {
    std::vector<RegionAnnotation> hqs_ = HQRegions(typeStrs);
    if (hqs_.size() == 0) 
        return RegionAnnotation(holeNumber_, RegionTypeMap::ToIndex(HQRegion, typeStrs), 0, 0, 0);
    else if (hqs_.size() == 1)
        return hqs_[0];
    else assert(false == "Zmw has more than one HQRegion."); 
}

DNALength RegionAnnotations::HQStart(const std::vector<std::string> & typeStrs) const {
    return TheHQRegion(typeStrs).GetStart();
}

DNALength RegionAnnotations::HQEnd(const std::vector<std::string> & typeStrs) const {
    return TheHQRegion(typeStrs).GetEnd();
}

std::vector<RegionAnnotation> 
RegionAnnotations::Inserts(const std::vector<std::string> & typeStrs) const {
    return RegionAnnotationsOfType(Insert, typeStrs);
}

std::vector<ReadInterval> 
RegionAnnotations::AdapterIntervals(const std::vector<std::string> & typeStrs) const {
    std::vector<ReadInterval> ret;
    for (auto adapter: Adapters(typeStrs)) {
        ret.push_back(ReadInterval(adapter));
    }
    return ret;
}

std::vector<ReadInterval> 
RegionAnnotations::SubreadIntervals(const std::vector<std::string> & typeStrs,
                                    const DNALength wholeLength,
                                    const bool byAdapter, 
                                    const bool byHQRegion) const {
    std::vector<RegionAnnotation> inserts;
    if (not byAdapter) {
        inserts = Inserts(typeStrs);
    } else {
        std::vector<DNALength> starts, ends; 
        starts.push_back(0);
        for(auto adapter: Adapters(typeStrs)) {
            starts.push_back(adapter.GetEnd());
            ends.push_back(adapter.GetStart());
        }
        ends.push_back(wholeLength);

        for (size_t i = 0; i < starts.size(); i++) {
            inserts.push_back(RegionAnnotation(holeNumber_,
                                               Insert, 
                                               starts[i],
                                               ends[i], 
                                               0));
        }
    }

    std::vector<ReadInterval> ret;
    for (auto insert: inserts) {
        if (byHQRegion and HasHQRegion(typeStrs)) {
            DNALength s = std::max(static_cast<UInt>(insert.GetStart()), HQStart(typeStrs));
            DNALength e = std::min(static_cast<UInt>(insert.GetEnd()), HQEnd(typeStrs));
            if (s < e) {
                ret.push_back(ReadInterval(s, e, insert.GetScore()));
            }
        } else {
           ret.push_back(ReadInterval(insert));
        }
    }
    return ret;
}
